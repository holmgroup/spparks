/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_moment.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "comm_lattice.h"
#include "comm_off_lattice.h"
#include "error.h"
#include "grain.h"
#include "math.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagMoment::DiagMoment(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass == App::LATTICE) latticeflag = 1;
  else if (app->appclass == App::OFF_LATTICE) latticeflag = 0;
  else error->all(FLERR,"Diag style incompatible with app style");
  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  fpdump = fopen(arg[iarg],"w");
	  if (!fpdump) error->one(FLERR,"Cannot open diag_style cluster output file");
	}
      } else error->all(FLERR,"Illegal diag_style cluster command");
    }
  }
}

DiagMoment::~DiagMoment() {
  if (me == 0)
    if (fpdump) fclose(fpdump);
}

/* ---------------------------------------------------------------------- */

void DiagMoment::init()
{
  if (latticeflag) applattice = (AppLattice *) app;
  else appofflattice = (AppOffLattice *) app;

  nlocal = applattice->nlocal;
  nghost = applattice->nghost;
  
  xyz = app->xyz;
  idsite = app->xyz;
  
  x_size = domain->boxxhi - domain->boxxlo;
  y_size = domain->boxyhi - domain->boxylo;
  z_size = domain->boxzhi - domain->boxzlo;
  
  // if domain is periodic, use minimum image convention
  // otherwise use plain distance
  if (domain->xperiodic == 1)
    x_dist = &DiagMoment::min_dist_x;
  else
    x_dist = &DiagMoment::dist;
  if (domain->yperiodic == 1)
    y_dist = &DiagMoment::min_dist_y;
  else
    y_dist = &DiagMoment::dist;
  if (domain->zperiodic == 1)
    z_dist = &DiagMoment::min_dist_z;
  else
    z_dist = &DiagMoment::dist;
  
}

/* ---------------------------------------------------------------------- */
/* Compute image moments of clusters. 
     Strong assumption that clusters have unique ivalue (ivalue == grain_id)  
     M_{ij} = \sum_x \sum_y x^i y^j I(x,y) 
*/
void DiagMoment::compute()
{
  int nlocal = app->nlocal;
  if (latticeflag) applattice->comm->all();
  else appofflattice->comm->all();

  /* Compute centroids (zeroth and first moments) in each processor domain,
       respecting PBC */
  /*
       - choose a reference point for each cluster
       - compute centroid in that reference frame with minimum image convention
       - translate centroid back into simulation reference frame
   */

  // get access to app-level data
  int* site = app->iarray[0];
  int* numneigh = app->numneigh;
  int** neighbor = app->neighbor;
  
  double x, y, z = 0;
  double dx, dy, dz = 0;
  Point3D reference;
  int grain_id = 0;
  
  // for each owned site i:
  for (int i = 0; i < nlocal; i++) {
    // get/set the reference site for that grain_id
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    grain_id = site[i];
    grain_iter = grains.find(grain_id);
    if (grain_iter == grains.end) {
      Grain new_grain = new Grain(grain_id, Point3D(x,y,z));
      dx = dy = dz = 0;
      grains[grain_id] = new_grain;
      grain_iter = grains.find(grain_id);
    }
    else {
      reference = grain_iter.second.reference;
      dx = (*x_dist)(x, reference.x);
      dy = (*y_dist)(y, reference.y);
      dz = (*z_dist)(z, reference.z);
    }

    grain_iter.second.update_centroid(Point3D(dx,dy,dz));
    grain_iter.second.volume += 1;

    // check neighboring site ids for grain boundaries
    for (int j = 0; j < numneigh[i], j++) {
      if (grain_id != neigbor[i][j]) 
	grain_iter.second.add_neigh(neighbor[i][j]);
    }
    
  }
  
  // move partial grain centroids back to simulation reference frame
  for (auto& grain : grains) {
    grain_id = grain.first;
    reference = grain.second.reference;
    Point3D delta = Point3D(-reference.x, -reference.y, -reference.z);
    grain.second.update_centroid(delta);
  }
  
  /* communicate partial moments to root process */
  /* pack grain data into buffer */

  int me_size,m,maxbuf;
  double* dbufclust;
  int tmp,nrecv;
  MPI_Status status;
  MPI_Request request;
  int nn;

  me_size = 0;
  for (auto& grain : grains) {
    me_size += 6+grain.second.nneigh;
  }
  if (me == 0) me_size = 0;

  MPI_Allreduce(&me_size,&maxbuf,1,MPI_INT,MPI_MAX,world);

  dbufclust = new double[maxbuf];

  if (me != 0) {
    m = 0;
    for (auto& grain : grains) {
      dbufclust[m++] = grain.second.ivalue;
      dbufclust[m++] = grain.second.volume;
      dbufclust[m++] = grain.second.centroid.x;
      dbufclust[m++] = grain.second.centroid.y;
      dbufclust[m++] = grain.second.centroid.z;
      dbufclust[m++] = grain.second.nneigh;
      for (int j = 0; j < grain.second.nneigh; j++) {
	dbufclust[m++] = grain.second.neighlist[j];
      }
    }
    
    if (me_size != m)
      error->one(FLERR,"Mismatch in counting for dbufclust");

  }

  // proc 0 pings each proc, receives it's data, adds it to list
  // all other procs wait for ping, send their data to proc 0

  if (me == 0) {
    for (int iproc = 1; iproc < nprocs; iproc++) {
      MPI_Irecv(dbufclust,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
      MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
      MPI_Wait(&request,&status);
      MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
      
      m = 0;
      while (m < nrecv) {
	iv = static_cast<int> (dbufclust[m++]);
	vol = dbufclust[m++];
	x = dbufclust[m++];
	y = dbufclust[m++];
	z = dbufclust[m++];
	nn = static_cast<int> (dbufclust[m++]);
	// add new grain, or merge respecting PBC
	merge_grain(iv,vol,x,y,z,nn,&dbufclust[m]);
	m+=nn;
	volsum+=vol;
      }
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(dbufclust,me_size,MPI_DOUBLE,0,0,world);
  }

  delete [] dbufclust;

  /* if higher-order moments wanted */
  /* Root process broadcasts cluster centroids */
  /* Compute partial higher moments relative to centroids, respecting PBC */
  /* communicate partial higher moments back to root process */
  /* root process combines partial higher moments */

  // root process saves to file
  if (me == 0) {
    if (fpdump) {
      fprintf(fpdump, '------------------------\n');
      fprintf(fpdump, 'ivalue volume Cent_x Cent_y Cent_z NumNeighs NeighList\n');
      for (auto& grain : grains) {
	fprintf(fpdump, "%d %f %f %f %f %d",
		grain.second.ivalue,
		grain.second.volume,
		grain.second.centroid.x,
		grain.second.centroid.y,
		grain.second.centroid.z,
		grain.second.nneigh);
	for (int ii = 0; ii < grain.second.nneigh; ii++)
	  fprintf(fpdump, " %d", grain.second.neighlist[ii]);
	fprintf(fpdump, "\n");
      }
    }
  }
}


/* ---------------------------------------------------------------------- */
void DiagMoment::merge_grain(int iv, double vol, double x, double y, double z, int nn, double* neighs) {
  // if grain does not exist, create a new grain.
  grain_iter = grains.find(iv);
  if (grain_iter == grains.end) {
    Grain new_grain = new Grain(iv, Point3D(x,y,z));
    new_grain.volume = vol;
    new_grain.centroid = Point3D(x,y,z);
    new_grain.nneigh = nn;
    for (int i = 0; i < nn; i++)
      new_grain.add_neigh(neighs[i]);
    grains[iv] = new_grain;
  }
  // else merge it
  else {
    double new_volume = grain_iter.second.volume + vol;
    double dx, dy, dz;
    // merge centroids
    dx = x_dist(x, grain_iter.second.centroid.x);
    grain_iter.second.centroid.x += (vol/new_volume)*dx;
    dy = y_dist(y, grain_iter.second.centroid.y);
    grain_iter.second.centroid.y += (vol/new_volume)*dy;
    dz = z_dist(z, grain_iter.second.centroid.z);
    grain_iter.second.centroid.z += (vol/new_volume)*dz;

    // merge neighbor list
    for (int i = 0; i < nn; i++)
      grain_iter.second.add_neigh(neighs[i]);
    // merge grain size
    grain_iter.second.volume = new_volume;
  }
}

/* ---------------------------------------------------------------------- */

void DiagMoment::stats(char *strtmp)
{
  sprintf(strtmp," %10g",energy);
}

/* ---------------------------------------------------------------------- */

void DiagMoment::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Energy");
}

double DiagMoment::dist(float p, float p_ref) {
  return p - p_ref;
}

double DiagMoment::min_dist_x(float p, float p_ref) {
  double dx = p - p_ref;
  dx = dx - x_size * floor(0.5 + (dx / x_size));
  return dx;
}

double DiagMoment::min_dist_y(float p, float p_ref) {
  double dy = p - p_ref;
  dy = dy - y_size * floor(0.5 + (dy / y_size));
  return dy;
}

double DiagMoment::min_dist_x(float p, float p_ref) {
  double dz = p - p_ref;
  dz = dy - z_size * floor(0.5 + (dz / z_size));
  return dz;
}
