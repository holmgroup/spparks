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
#include "diag_grainsize.h"
#include "app.h"
#include "app_lattice.h"
#include "domain.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "comm_off_lattice.h"
#include "error.h"
#include "grain.h"
#include "math.h"
#include <iostream>

using namespace SPPARKS_NS;
using namespace H5;

/* lattice styles enum */
enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
       FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

#define DATA_PATH "grainsize"

/* ---------------------------------------------------------------------- */

DiagGrainsize::DiagGrainsize(SPPARKS *spk, int narg, char **arg) : 
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
	  std::cout << "creating hdf5 file " << arg[iarg] << std::endl;
	  create_hdf_file(arg[iarg]);
	  // TODO: error checking
	}
      } else error->all(FLERR,"Illegal diag_style cluster command");
    }
    iarg++;
  }
}

DiagGrainsize::~DiagGrainsize() {
  grains.clear();
  if (me == 0) {
    save_time(); // save std::vector for time
    prop.close();
    delete dataspace;
    delete dataset;
    output_file.close();
  }
}


void DiagGrainsize::save_time() {
  hsize_t time_dim = time.size();
  DataSpace timespace = DataSpace(1, &time_dim);
  DataSet timedata = output_file.createDataSet("time", PredType::NATIVE_DOUBLE, timespace);
  timedata.write(time.data(), PredType::NATIVE_DOUBLE);
  return;
}

/* ---------------------------------------------------------------------- */

void DiagGrainsize::init()
{
  if (latticeflag)
    applattice = (AppLattice *) app;

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
    x_dist = &DiagGrainsize::min_dist_x;
  else
    x_dist = &DiagGrainsize::dist;

  if (domain->yperiodic == 1)
    y_dist = &DiagGrainsize::min_dist_y;
  else
    y_dist = &DiagGrainsize::dist;

  if (domain->zperiodic == 1)
    z_dist = &DiagGrainsize::min_dist_z;
  else
    z_dist = &DiagGrainsize::dist;

  /* set up nearest neighbor lists for periodic lattices */
  n_nearest = 0;
  nearest_neighs = NULL;

  /* this initialization syntax only works in C++11 */
  if (domain->lattice->style == SC_26N) {
    // nearest neighbors in _periodic_ SC_26N lattice
    n_nearest = 6;
    nearest_neighs = new int[n_nearest] {4, 10, 12, 13, 15, 21};
  }
  else if (domain->lattice->style == SQ_8N) {
    // nearest neighbors in _periodic_ SQ_8N lattice
    n_nearest = 4;
    nearest_neighs = new int[n_nearest] {1, 3, 4, 6};
  }
  else
    error->all(FLERR,"diag_grainsize not implemented for this lattice.");

}

void DiagGrainsize::find_nspins() {
  // find max grain_id
  int *site = app->iarray[0];
  int max_id = 0;
  for (int i = 0; i < nlocal; i++) {
    max_id = std::max(max_id, site[i]);
  }

  MPI_Allreduce(&max_id,&nspins,1,MPI_INT,MPI_MAX,world);
  if (me == 0)
    std::cout << nspins << " grains." << std::endl;
}


void DiagGrainsize::create_hdf_file(char* filename) {
  init(); // needed to calculate nspins...
  if (me == 0) {
#ifdef SPPARKS_HDF5
    std::cout << "diag_grainsize filename: " << filename << std::endl;
    // output_file = H5Fcreate(filecurrent, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    output_file = H5File(filename, H5F_ACC_TRUNC);

    find_nspins();
    
    // shape: {0, app->nspins}; --> grow by one row on each iteration.
    dims[0] = 0;
    dims[1] = nspins+1;
    hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED}; // extensible dataset
    dataspace = new DataSpace(2, dims, maxdims);

    // enable hdf5 data chunking to get extensible dataset
    chunk_dims[0] = 2; chunk_dims[1] = 5;
    prop.setChunk(2, chunk_dims);

    dataset = new DataSet(output_file.createDataSet(DATA_PATH, PredType::NATIVE_INT,
						    *dataspace, prop));
#endif
  }
}


/* ---------------------------------------------------------------------- */
/* Compute image moments of clusters. 
     Strong assumption that clusters have unique ivalue (ivalue == grain_id)  
     M_{ij} = \sum_x \sum_y x^i y^j I(x,y) 
*/
void DiagGrainsize::compute()
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
  
  // compute partial sums for grain centroids
  // also build a list of grain neighbors
  for (int i = 0; i < nlocal; i++) {
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    grain_id = site[i];

    // get the current grain, or insert a new grain
    grain_iter = grains.find(grain_id);
    if (grain_iter == grains.end()) {
      Grain new_grain = Grain(grain_id, Point3D(x,y,z));
      grain_iter = grains.insert(grain_iter,
				 std::make_pair(grain_id, new_grain));
    }
    Grain& grain = grain_iter->second;

    // increment volume and centroid partial sum
    dx = (this->*x_dist)(x, grain.reference.x);
    dy = (this->*y_dist)(y, grain.reference.y);
    dz = (this->*z_dist)(z, grain.reference.z);
    grain.update_centroid(Point3D(dx,dy,dz));
    grain.volume += 1;

    // check neighboring site ids for grain boundaries
    // look at nearest neighbor sites, specific to lattice style
    int j_id = 0;    
    for (int j = 0; j < n_nearest; j++) {
      // nearest neighbors in _periodic_ SC_26N lattice
      j_id = site[neighbor[i][nearest_neighs[j]]];
      if (grain_id != j_id) 
	grain.add_neigh(j_id);
    }
  }
  
  // reduce centroid partial sums to centroids, move them to simulation reference frame
  for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
    Grain& grain = grain_iter->second;
    ref = grain.reference;
    // first divide by volume to get actual centroid
    grain.centroid.x = grain.centroid.x / grain.volume;
    grain.centroid.y = grain.centroid.y / grain.volume;
    grain.centroid.z = grain.centroid.z / grain.volume;
    grain.update_centroid(grain.reference);
  }
  
  /* communicate partial moments to root process */
  /* pack grain data into buffer */
  {
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
  }

  if (me == 0) {
    // move grain centroids back into the simulation box
    for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
      Grain& grain = grain_iter->second;
      grain.centroid.x -= x_size * floor(grain.centroid.x / x_size);
      grain.centroid.y -= y_size * floor(grain.centroid.y / y_size);
      grain.centroid.z -= z_size * floor(grain.centroid.z / z_size);
      grain.reference = grain.centroid;
    }
  }
  
  /* if higher-order moments wanted */
  /* Root process broadcasts cluster centroids */
  {
    int m,maxbuf;
    double* dbufcent;
    int tmp,nrecv;
    MPI_Status status;
    MPI_Request request;

    /* broadcast grain_id and centroid for each grain */
    maxbuf = 0;
    if (me == 0)
      maxbuf = 4 * grains.size();
    MPI_Bcast(&maxbuf, 1,MPI_INT,0,world);

    dbufcent = new double[maxbuf];

    if (me == 0) {
      /* root process packs buffer with grain ids and centroids */
      m = 0;
      for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
	Grain& grain = grain_iter->second;
	dbufcent[m++] = grain.id;
	dbufcent[m++] = grain.centroid.x;
	dbufcent[m++] = grain.centroid.y;
	dbufcent[m++] = grain.centroid.z;
	grain.reference = grain.centroid;
      }
    }
    
    /* root proc broadcasts centroid buffer */
    MPI_Bcast(dbufcent,maxbuf,MPI_DOUBLE,0,world);

    /* update centroids in place */
    if (me != 0) {
      for (m = 0; m < maxbuf; m+=4) {
	grain_id = static_cast<int>(dbufcent[m]);
	grain_iter = grains.find(grain_id);
	// skip grains not partially owned by this proc
	if (grain_iter == grains.end())
	  continue;
	Grain& grain = grain_iter->second;
	grain.centroid.x = dbufcent[m+1];
	grain.centroid.y = dbufcent[m+2];
	grain.centroid.z = dbufcent[m+3];
	grain.reference = grain.centroid;
      }
    }
    delete [] dbufcent;
  }
  
  /* Compute partial higher moments relative to centroids, respecting PBC */
  /* Note: volume and neighbor count/list not up to date */
  for (int i = 0; i < nlocal; i++) {
    grain_id = site[i];
    grain_iter = grains.find(grain_id);
    if (grain_iter == grains.end())
      error->all(FLERR,"inconsistent grain_id");
    Grain& grain = grain_iter->second;

    // increment volume and centroid partial sum
    x = (this->*x_dist)(xyz[i][0], grain.centroid.x);
    y = (this->*y_dist)(xyz[i][1], grain.centroid.y);
    z = (this->*z_dist)(xyz[i][2], grain.centroid.z);
    grain.update_central_moments(Point3D(x,y,z));
  }

  {
    /* communicate partial higher moments back to root process */
    int me_size,m,maxbuf;
    double* dbufmoment;
    int tmp,nrecv;
    MPI_Status status;
    MPI_Request request;

    me_size = 7 * grains.size();
    if (me == 0) me_size = 0;

    MPI_Allreduce(&me_size,&maxbuf,1,MPI_INT,MPI_MAX,world);

    dbufmoment = new double[maxbuf];

    if (me != 0) {
      m = 0;
      for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
	Grain& grain = grain_iter->second;
	dbufmoment[m++] = grain.id;
	for (int i = 0; i < NUM_MOMENTS; i++) {
	  dbufmoment[m++] = grain.second_moment[i];
	}
      }
    
      if (me_size != m)
	error->one(FLERR,"Mismatch in counting for dbufclust");
    
    }

    // proc 0 pings each proc, receives it's data, adds it to list
    // all other procs wait for ping, send their data to proc 0

    if (me == 0) {
      for (int iproc = 1; iproc < nprocs; iproc++) {
	MPI_Irecv(dbufmoment,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nrecv);
      
	m = 0;
	while (m < nrecv) {
	  /* root proc combines partial moments */
	  grain_id = static_cast<int> (dbufmoment[m++]);
	  merge_partial_moments(grain_id,&dbufmoment[m]);
	  m += NUM_MOMENTS;
	}
      }
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(dbufmoment,me_size,MPI_DOUBLE,0,0,world);
    }
    
    delete [] dbufmoment;
  }

  // root process saves to file
  if (me == 0) {
    time.push_back(applattice->time); // save time for later
    // fill the data buffer with grain sizes
    int *grainsize_buf = new int[nspins+1];
    for (int i = 0; i < nspins+1; i++)
      grainsize_buf[i] = 0;
    for (grain_iter = grains.begin(); grain_iter != grains.end(); grain_iter++) {
      Grain& grain = grain_iter->second;
      grainsize_buf[grain.id] = grain.volume;
    }

    // grow the hdf5 dataset
    hsize_t size[2];
    dims[0] = dims[0] + 1;
    size[0] = dims[0];
    size[1] = dims[1];
    dataset->extend(size);

    // select a hyperslab in the newly extended part of the dataset
    DataSpace *filespace = new DataSpace(dataset->getSpace());
    hsize_t offset[2];
    offset[0] = size[0] - 1;
    offset[1] = 0;
    hsize_t dimsext[2];
    dimsext[0] = 1;
    dimsext[1] = nspins+1;
    filespace->selectHyperslab(H5S_SELECT_SET, dimsext, offset);

    // define a memory space
    DataSpace *memspace = new DataSpace(2, dimsext, NULL);

    // actually write data
    dataset->write(grainsize_buf, PredType::NATIVE_INT, *memspace, *filespace);
    delete filespace;
    delete memspace;
  }
}


/* ---------------------------------------------------------------------- */
void DiagGrainsize::merge_grain(int iv, double vol, double x, double y, double z, int nn, double* neighs) {
  
  // 
  grain_iter = grains.find(iv);

  // merge incoming grain with pre-existing grain
  if (grain_iter != grains.end()) {
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

  // if incoming grain isn't already in the grains hashmap,
  // insert a new grain
  else if (grain_iter == grains.end()) {
    Grain new_grain = Grain(iv, Point3D(0,0,0));
    new_grain.volume = vol;
    new_grain.centroid = Point3D(x,y,z);
    for (int i = 0; i < nn; i++)
      new_grain.add_neigh(static_cast<int>(neighs[i]));
    grains.insert(std::make_pair(iv, new_grain));
  }
}


void DiagGrainsize::merge_partial_moments(int grain_id, double* moment) {
  grain_iter = grains.find(grain_id);
  
  if (grain_iter == grains.end()) {
    error->all(FLERR,"Merging moments: grain not found.");
  }
  Grain& grain = grain_iter->second;
  for (int i = 0; i < NUM_MOMENTS; i++) {
    grain.second_moment[i] += moment[i];
  }
}


/* ---------------------------------------------------------------------- */

void DiagGrainsize::stats(char *strtmp)
{
  sprintf(strtmp," %10d", grains.size());
}

/* ---------------------------------------------------------------------- */

void DiagGrainsize::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Ngrains");
}

inline double DiagGrainsize::dist(double p, double p_ref) {
  return p - p_ref;
}

inline double DiagGrainsize::min_dist(double p, double p_ref, double p_size) {
  double delta = dist(p,p_ref);
  delta = delta - (p_size * round(delta / p_size));
  return delta;
}

double DiagGrainsize::min_dist_x(double x, double x_ref) {
  return min_dist(x, x_ref, x_size);
}

double DiagGrainsize::min_dist_y(double y, double y_ref) {
  return min_dist(y, y_ref, y_size);
}

double DiagGrainsize::min_dist_z(double z, double z_ref) {
  return min_dist(z, z_ref, z_size);
}

