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

  energy = 0.0;
}

/* ---------------------------------------------------------------------- */

void DiagMoment::compute()
{
  int nlocal = app->nlocal;
  if (latticeflag) applattice->comm->all();
  else appofflattice->comm->all();

  /* Find a reference frame for each cluster */
  for (int i = 0; i < nlocal; i++) {
    ;
  }
  /* Compute partial image moment sums for each cluster */

  double etmp = 0.0;
  if (latticeflag)
    for (int i = 0; i < nlocal; i++) etmp += applattice->site_energy(i);
  else
    for (int i = 0; i < nlocal; i++) etmp += appofflattice->site_energy(i);

  /* communicate partial image moment sums to root process */
  
  /* Root process computes image moments for each cluster */
  

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
