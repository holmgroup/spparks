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

#include "stdlib.h"
#include "stdio.h"
#include "grain.h"

using namespace SPPARKS_NS;

Grain::Grain(int id, double vol, int nn, int* neighs) {
  id = id;
  volume = vol;
  reference = Point3D(0,0,0);
  centroid = Point3D(0,0,0);
  for (int i = 0; i < NUM_MOMENTS; i++)
    second_moment[i] = 0;
  nneigh = nn;
  if (nneigh == 0) {
    neighlist = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = static_cast<int> (neighs[i]);
    }
  }
}

/* Create an 'empty' grain data structure */
Grain::Grain(int grain_id, Point3D ref) {
  id = grain_id;
  volume = 0;
  reference = ref;
  centroid = Point3D(0,0,0);
  for (int i = 0; i < NUM_MOMENTS; i++)
    second_moment[i] = 0;
  nneigh = 0;
  neighlist = NULL;
}

/* Copy constructor */
Grain::Grain(const Grain& g) {
  id = g.id;
  volume = g.volume;
  reference = g.reference;
  centroid = g.centroid;
  for (int i = 0; i < NUM_MOMENTS; i++) {
    second_moment[i] = g.second_moment[i];
  }
  nneigh = g.nneigh;
  if (nneigh == 0) {
    neighlist = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = g.neighlist[i];
    }
  }
}

// Define assigment operator.
// Because memory is allocated before
// the object is initialized, can not use constructor.

Grain& Grain::operator=(const Grain& g) {
  id = g.id;
  volume = g.volume;
  centroid = g.centroid;
  for (int i = 0; i < NUM_MOMENTS; i++) {
    second_moment[i] = g.second_moment[i];
  }
  nneigh = g.nneigh;
  if (nneigh == 0) {
    neighlist = NULL;
  } else {
    neighlist = (int*) malloc(nneigh*sizeof(int));
    for (int i = 0; i < nneigh; i++) {
      neighlist[i] = g.neighlist[i];
    }
  }
  return *this;
}
  
Grain::~Grain() {
  if (neighlist != NULL) free(neighlist);
}
  
void Grain::add_neigh(int id) {
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    if (id == neighlist[ineigh]) return;
  }
  nneigh++;
  neighlist = (int*) realloc(neighlist,nneigh*sizeof(int));
  if (neighlist == NULL)
    fprintf(stdout,"problem with realloc\n");
  // fprintf(stdout, "grain %d: nneigh = %d, neigh %d\n", ivalue, nneigh, id);
  neighlist[nneigh-1] = id;
}

void Grain::update_centroid(Point3D delta) {
  centroid.x += delta.x;
  centroid.y += delta.y;
  centroid.z += delta.z;
}

void Grain::update_central_moments(Point3D p) {
  /* second-order central moments: */
  
  second_moment[0] += p.x * p.x; // \mu_200
  second_moment[1] += p.y * p.y; // \mu_020
  second_moment[2] += p.z * p.z; // \mu_002
  second_moment[3] += p.x * p.y; // \mu_110
  second_moment[4] += p.x * p.z; // \mu_101
  second_moment[5] += p.y * p.z; // \mu_011
}

void Grain::print(FILE* fp) {
  fprintf(fp,"%d %g %d ",id,volume,nneigh);
  for (int ineigh = 0; ineigh < nneigh; ineigh++) {
    fprintf(fp,"%d ",neighlist[ineigh]);
  }
  fprintf(fp,"\n");
}
  
