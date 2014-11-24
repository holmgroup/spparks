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
    current_grain = grains.find(grain_id);
    if (current_grain == grains.end) {
      Grain new_grain = new Grain(grain_id, Point3D(x,y,z));
      dx = dy = dz = 0;
      grains[grain_id] = new_grain;
      current_grain = grains.find(grain_id);
    }
    else {
      reference = current_grain->reference;
      dx = (*x_dist)(point.x, reference.x);
      dy = (*y_dist)(point.y, reference.y);
      dz = (*z_dist)(point.z, reference.z);
    }

    current_grain.update_centroid(Point3D(dx,dy,dz));
    current_grain->volume += 1;
  }
  
  /* communicate partial moments to root process */
  // first move grain centroids back to simulation reference frame

  
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
