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
#include "point3d.h"

using namespace SPPARKS_NS;

Point3D::Point3D() {
  x = 0;
  y = 0;
  z = 0;
}

Point3D::Point3D(float x_i, float y_i, float z_i) {
  x = x_i;
  y = y_i;
  z = z_i;
}
  
