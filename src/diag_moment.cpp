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

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagMoment::DiagMoment(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass == App::LATTICE) latticeflag = 1;
  else if (app->appclass == App::OFF_LATTICE) latticeflag = 0;
  else error->all(FLERR,"Diag style incompatible with app style");
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
  
  boxxlo = domain->boxxlo;
  boxxhi = domain->boxxhi;
  boxylo = domain->boxylo;
  boxyhi = domain->boxyhi;
  boxzlo = domain->boxzlo;
  boxzhi = domain->boxzhi;
  
}

/* ---------------------------------------------------------------------- */
/* Compute image moments of clusters. 
     Strong assumption that clusters have unique ivalue (== grain_id)  
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
  int* site = app->iarray[0];
  double x, y, z = 0;
  int grain_id = 0;
  
  // for each owned site i:
  for (int i = 0; i < nlocal; i++) {
    // get/set the reference site for that grain_id
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    grain_id = site[i];
    
  }
  
  /* communicate partial moments to root process */
  
  /* Root process combines image moments for each cluster , respecting PBC */

  /* if higher-order moments wanted */
  /* Root process broadcasts cluster centroids */
  /* Compute partial higher moments relative to centroids, respecting PBC */
  /* communicate partial higher moments back to root process */
  /* root process combines partial higher moments */

  MPI_Allreduce(&etmp,&energy,1,MPI_DOUBLE,MPI_SUM,world);
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
