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

#include "point3d.h"

#define NUM_MOMENTS 6

namespace SPPARKS_NS {

class Grain {
 public:
  int id;
  double volume;
  Point3D reference;
  Point3D centroid;
  int nneigh;
  int* neighlist;
  double second_moment[NUM_MOMENTS];

  Grain(int, double, int, int*);
  Grain(int, Point3D);
  Grain(const Grain&);
  Grain& operator=(const Grain&);
  ~Grain();
  void add_neigh(int);
  void update_centroid(Point3D);
  void update_central_moments(Point3D);
  void print(FILE*);
};

}

#endif
