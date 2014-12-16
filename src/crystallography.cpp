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
#include "math.h"
#include "crystallography.h"

using namespace SPPARKS_NS;

Crystallography::Crystallography() {
  quaternion_buf = NULL;
  
  /* default to uniform energy and mobility */
  this->energy_pt = &Crystallography::uniform_energy;
  this->mobility_pt = &Crystallography::uniform_mobility;

  /* default to cubic symmetry */
  load_symmetry_operator("cubic");
}

Crystallography::~Crystallography() {
  if (quaternion_buf != NULL) { 
    delete [] quaternion_buf;
    quaternion_buf = NULL;
  }
}

double Crystallography::energy(int i, int j) {
  return (*this.*energy_pt)(i,j);
}

double Crystallography::mobility(int i, int j) {
  return (*this.*mobility_pt)(i,j);
}

void Crystallography::set_options(int argc, char** argv) {
  int i = 0;
  while (i < argc) {
    if (strcmp(argv[i],"Energy") == 0) {
      i++;
      if (strcmp(argv[i],"ReadShockley") == 0) {
	i++;
	this->energy_pt = &Crystallography::read_shockley_energy;
      }
      else error->all(FLERR,"Crystallography energy function not implemented.");
    }
    else if (strcmp(argv[i],"Mobility") == 0) {
      i++;
      if (strcmp(argv[i],"HwangHumphreys") == 0) {
	i++;
	this->energy_pt = &Crystallography::hwang_humphreys_mobility;
      }
      else error->all(FLERR,"Crystallography mobility function not implemented.");
    }
    else error->all(FLERR,"illegal crystallography option.");
    i++
  }
}

double uniform_energy(int i, int j) {
  return 1.0;
}

double uniform_mobility(int i, int j) {
  return 1.0;
}

double read_shockley_energy(int i, int j) {
  double misori = get_misorientation_angle(i,j);
  // do read/shockley formula
  return 1.0;
}

double hwang_humphreys_mobility(int i, int j) {
  double misori = get_misorientation_angle(i,j);
  // do the hwang/humphreys formula
  return 1.0;
}

void Crystallography::copy_euler_angle_data(float *data, int num_orientations) {
  int data_size = 3 * num_orientations;
  float* euler_buf = new float[data_size];
  for (int i = 0; i < data_size; i++) {
    euler_buf[i] = data[i];
  }
}

void Crystallography::copy_quaternion_data(float *data, int num_orientations) {
  int data_size = 4 * num_orientations;
  quaternion_buf = new double[data_size];
  for (int i = 0; i < num_orientations; i+=4) {
    for (int j = 0; j < 4; j++) {
      quaternion_buf[i+j] = static_cast<double>(data[i+j]);
      QuatMap quat(&quaternion_buf[i], 4);
      quat.normalize();
    }
  }
}

