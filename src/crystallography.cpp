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
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "error.h"
#include "crystallography.h"

using namespace SPPARKS_NS;

Crystallography::Crystallography() {
  quaternion_buf = NULL;
  
  /* default to uniform energy and mobility */
  this->energy_pt = &Crystallography::uniform_energy;
  this->mobility_pt = &Crystallography::uniform_mobility;

  /* default to cubic symmetry */
  n_symm = 0;
  load_symmetry_operator("cubic");
}

Crystallography::~Crystallography() {
  if (quaternion_buf != NULL) { 
    delete [] quaternion_buf;
    quaternion_buf = NULL;
  }
  if (symmetry_buf != NULL) {
    delete [] symmetry_buf;
    symmetry_buf = NULL;
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
      else ; // error->all(FLERR,"Crystallography energy function not implemented.");
    }
    else if (strcmp(argv[i],"Mobility") == 0) {
      i++;
      if (strcmp(argv[i],"HwangHumphreys") == 0) {
	i++;
	this->energy_pt = &Crystallography::hwang_humphreys_mobility;
      }
      else ; // error->all(FLERR,"Crystallography mobility function not implemented.");
    }
    else ; // error->all(FLERR,"illegal crystallography option.");
    i++;
  }
}

double Crystallography::uniform_energy(int i, int j) {
  return 1.0;
}

double Crystallography::uniform_mobility(int i, int j) {
  return 1.0;
}

double Crystallography::read_shockley_energy(int i, int j) {
  double misori;
  misori = get_misorientation_angle(i,j);
  // do read/shockley formula
  return 1.0;
}

double Crystallography::hwang_humphreys_mobility(int i, int j) {
  double misori;
  misori = get_misorientation_angle(i, j);
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
      // Eigen::Map<Eigen::Quaternion<double> > q(&quaternion_buf[0]);
      QuatMap q(&quaternion_buf[i]);
      q.normalize();
    }
  }
}

void Crystallography::load_symmetry_operator(std::string symmetry_style) {
  /* Read Prof. Rollett's quaternion symmetry files into raw buffer */
  std::string blank;
  std::string filename = "quat.symm." + symmetry_style;
  std::ifstream infile("quat.symm.cubic");

  std::getline(infile, blank); // first line is a file header
  infile >> n_symm; // second line
  symmetry_buf = new double[4 * n_symm];
  
  for (int i = 0; i < n_symm; i++) {
    int j = 4 * i;
    infile >> symmetry_buf[j];   // x
    infile >> symmetry_buf[j+1]; // y
    infile >> symmetry_buf[j+2]; // z
    infile >> symmetry_buf[j+3]; // w
    QuatMap q(&symmetry_buf[i]);
    q.normalize();
  }
}

double Crystallography::get_misorientation_angle(int i, int j) {
  ;
} 
