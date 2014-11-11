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

#ifndef SPK_GRAIN_H
#define SPK_GRAIN_H

#include float3.h

namespace SPPARKS_NS {

class Grain {
 public:
  int global_id;
  int ivalue;
  double dvalue;
  double volume;
  Point3D reference;
  Point3D centroid;
  int nneigh;
  int* neighlist;

  Grain(int, int, double, double, int, double*);
  Grain& operator=(const Grain&);
  ~Grain();
  void add_neigh(int);
  void print(FILE*);
};

}

#endif
