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
#include "domain.h"
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
    iarg++;
  }
}

DiagMoment::~DiagMoment() {
  grains.clear();
  if (me == 0)
    if (fpdump) fclose(fpdump);
}

/* ---------------------------------------------------------------------- */

void DiagMoment::init()
{
  if (latticeflag) applattice = (AppLattice *) app;

  nlocal = applattice->nlocal;
  nghost = applattice->nghost;
  
  xyz = app->xyz;
  idsite = app->id;
  
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
  grains.clear();
  int nlocal = app->nlocal;
  applattice->comm->all();

  /* Compute centroids (zeroth and first moments) in each processor domain,
       respecting PBC */
  /*
       - choose a reference point for each cluster
       - compute centroid in that reference frame with minimum image convention
       - translate centroid back into simulation reference frame
   */

  // get access to app-level data
  int* site = app->iarray[0];
  int* numneigh = &(applattice->numneigh[0]);
  int** neighbor = &(applattice->neighbor[0]);
  
  double x, y, z = 0;
  double dx, dy, dz = 0;
  Point3D ref;
  int grain_id = 0;
  
  // for each owned site i:
  for (int i = 0; i < nlocal; i++) {
    // get/set the reference site for that grain_id
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    grain_id = site[i];

    // get the current grain, or insert a new grain
    grain_iter = grains.find(grain_id);
    if (grain_iter == grains.end()) {
      grain_iter = grains.insert(grain_iter, std::make_pair(grain_id,
							    Grain(grain_id, Point3D(x,y,z)) ));
    }
    Grain& grain = grain_iter->second;

    ref = grain.reference;
    dx = (this->*x_dist)(x, ref.x);
    dy = (this->*y_dist)(y, ref.y);
    dz = (this->*z_dist)(z, ref.z);
    grain.update_centroid(Point3D(dx,dy,dz));

    grain.volume += 1;

    // check neighboring site ids for grain boundaries
    int j_id = 0;
    for (int j = 0; j < numneigh[i]; j++) {
      j_id = site[neighbor[i][j]];
      if (grain_id != j_id) 
	grain.add_neigh(j_id);
    }
  }
  
  // move partial grain centroids back to simulation reference frame
  for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
    Grain& grain = grain_iter->second;
    ref = grain.reference;
    // first divide by area to get actual centroid
    grain.centroid.x = grain.centroid.x / grain.volume;
    grain.centroid.y = grain.centroid.y / grain.volume;
    grain.centroid.z = grain.centroid.z / grain.volume;
    Point3D delta = Point3D(ref.x, ref.y, ref.z);
    grain.update_centroid(delta);
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
  for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
    me_size += 6+grain_iter->second.nneigh;
  }
  if (me == 0) me_size = 0;

  MPI_Allreduce(&me_size,&maxbuf,1,MPI_INT,MPI_MAX,world);

  dbufclust = new double[maxbuf];

  if (me != 0) {
    m = 0;
    for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
      Grain& grain = grain_iter->second;
      dbufclust[m++] = grain.id;
      dbufclust[m++] = grain.volume;
      dbufclust[m++] = grain.centroid.x;
      dbufclust[m++] = grain.centroid.y;
      dbufclust[m++] = grain.centroid.z;
      dbufclust[m++] = grain.nneigh;
      for (int j = 0; j < grain.nneigh; j++) {
	dbufclust[m++] = grain.neighlist[j];
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
      double volsum = 0;
      double vol = 0;
      while (m < nrecv) {
	grain_id = static_cast<int> (dbufclust[m++]);
	vol = dbufclust[m++];
	x = dbufclust[m++];
	y = dbufclust[m++];
	z = dbufclust[m++];
	nn = static_cast<int> (dbufclust[m++]);
	// add new grain, or merge respecting PBC
	merge_grain(grain_id,vol,x,y,z,nn,&dbufclust[m]);
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
      fprintf(fpdump, "------------------------\n");
      fprintf(fpdump, "%d grains\n", grains.size());
      fprintf(fpdump, "id volume Cent_x Cent_y Cent_z NumNeighs NeighList\n");
      for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
	Grain& grain = grain_iter->second;
	grain.centroid.x -= x_size * floor(grain.centroid.x / x_size);
	grain.centroid.y -= y_size * floor(grain.centroid.y / y_size);
	grain.centroid.z -= z_size * floor(grain.centroid.z / z_size);
	fprintf(fpdump, "%d %f %f %f %f %d",
		grain.id,
		grain.volume,
		grain.centroid.x,
		grain.centroid.y,
		grain.centroid.z,
		grain.nneigh);
	for (int ii = 0; ii < grain.nneigh; ii++)
	  fprintf(fpdump, " %d", grain.neighlist[ii]);
	fprintf(fpdump, "\n");
      }
    }
  }
}


/* ---------------------------------------------------------------------- */
void DiagMoment::merge_grain(int iv, double vol, double x, double y, double z, int nn, double* neighs) {
  // if grain does not exist, create a new grain.
  grain_iter = grains.find(iv);
  if (grain_iter == grains.end()) {
    Grain new_grain = Grain(iv, Point3D(0,0,0));
    new_grain.volume = vol;
    new_grain.centroid = Point3D(x,y,z);
    for (int i = 0; i < nn; i++)
      new_grain.add_neigh(static_cast<int>(neighs[i]));
    grains.insert(std::make_pair(iv, new_grain));
  }
  // else merge it
  else {
    Grain& grain = grain_iter->second;
    double new_volume = grain.volume + vol;
    double dx, dy, dz;
    // merge centroids
    dx = (this->*x_dist)(x, grain.centroid.x);
    grain.centroid.x += (vol/new_volume)*dx;
    dy = (this->*y_dist)(y, grain.centroid.y);
    grain.centroid.y += (vol/new_volume)*dy;
    dz = (this->*z_dist)(z, grain.centroid.z);
    grain.centroid.z += (vol/new_volume)*dz;

    // merge neighbor list
    for (int i = 0; i < nn; i++)
      grain.add_neigh(static_cast<int>(neighs[i]));
    // merge grain size
    grain.volume = new_volume;
  }
}

/* ---------------------------------------------------------------------- */

void DiagMoment::stats(char *strtmp)
{
  sprintf(strtmp," %10d", grains.size());
}

/* ---------------------------------------------------------------------- */

void DiagMoment::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Ngrains");
}

inline double DiagMoment::dist(double p, double p_ref) {
  return p - p_ref;
}

inline double DiagMoment::min_dist(double p, double p_ref, double p_size) {
  double delta = dist(p,p_ref);
  delta = delta - (p_size * round(delta / p_size));
  return delta;
}


// inline double DiagMoment::min_dist(double p, double p_ref, double p_size) {
//   double delta = dist(p,p_ref);
//   if (abs(delta) > 0.5*p_size)
//     delta = delta - copysign(p_size,delta);
//   return delta;
// }

double DiagMoment::min_dist_x(double x, double x_ref) {
  return min_dist(x, x_ref, x_size);
}

double DiagMoment::min_dist_y(double y, double y_ref) {
  return min_dist(y, y_ref, y_size);
}

double DiagMoment::min_dist_z(double z, double z_ref) {
  return min_dist(z, z_ref, z_size);
}
