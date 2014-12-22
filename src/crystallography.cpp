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

#include "signal.h"

using namespace SPPARKS_NS;

Crystallography::Crystallography() {

  quaternion_buf = NULL;
  
  /* default to uniform energy and mobility */
  this->energy_pt = &Crystallography::uniform_energy;
  this->mobility_pt = &Crystallography::uniform_mobility;

  /* default to cubic symmetry */
  n_symm = 0;
  e_theta_max = 15;
  m_theta_max = 15;

  // load_symmetry_operator(); // needs error checking code...
  initialize_cubic_symmetry();
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

void Crystallography::set_RS(double theta_max) {
  e_theta_max = theta_max;
  this->energy_pt = &Crystallography::read_shockley_energy;
}
void Crystallography::set_HH(double theta_max) {
  m_theta_max = theta_max;
  this->energy_pt = &Crystallography::hwang_humphreys_mobility;
}

double Crystallography::uniform_energy(int i, int j) {
  return 1.0;
}

double Crystallography::uniform_mobility(int i, int j) {
  return 1.0;
}

double Crystallography::read_shockley_energy(int i, int j) {
  double misori;
  if (i == j)
    return 0;
  misori = get_misorientation_angle(i,j);

  // high angle cutoff
  misori = std::min(misori, e_theta_max);
  // minimum allowed energy
  misori = std::max(misori, 0.1);
  
  misori = misori / e_theta_max;
  // do read/shockley formula
  return misori * (1 - log(misori));
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
  for (int i = 0; i < num_orientations; i++) {
    for (int j = 0; j < 4; j++) {
      quaternion_buf[(4*i)+j] = static_cast<double>(data[(4*i)+j]);
      // Eigen::Map<Eigen::Quaternion<double> > q(&quaternion_buf[0]);
      QuatMap q(&quaternion_buf[4*i]);
      q.normalize();
    }
  }
}

void Crystallography::load_symmetry_operator() {
  /* Read Prof. Rollett's quaternion symmetry files into raw buffer */
  fprintf(stdout,"loading symmetry operator\n");
  std::string blank;
  // std::string filename = "quat.symm." + symmetry_style;
  std::ifstream infile("quat.symm.cubic");
  // TODO: test the file stream. or get rid of reading from files...
  std::getline(infile, blank); // ignore first line
  std::cout << blank << std::endl;
  infile >> n_symm; // second line
  fprintf(stdout, "symmetry: %d\n", n_symm);
  symmetry_buf = new double[4 * n_symm];
  
  for (int i = 0; i < n_symm; i++) {
    int j = 4 * i;
    infile >> symmetry_buf[j];   // x
    infile >> symmetry_buf[j+1]; // y
    infile >> symmetry_buf[j+2]; // z
    infile >> symmetry_buf[j+3]; // w
    QuatMap q(&symmetry_buf[4*i]);
    q.normalize();
  }
  for (int i = 0; i < n_symm; i++) {
    int j = 4 * i;
    QuatMap q(&symmetry_buf[j]);
  }
}

double Crystallography::get_misorientation_angle(int i, int j) {
  // compute misorientation angle
  // axis is not necessarily in the FZ
  double wmin = 999999;
  AxAng a;
  Quat qmis, q, qmin;
  
  QuatMapConst q_i(&quaternion_buf[4*i]);
  QuatMapConst q_j(&quaternion_buf[4*j]);
  
  qmis = q_j * q_i.inverse();
  qmis.normalize();
  // loop over crystal symmetry operators for grain i
  for (int ii = 0; ii < n_symm; ii++) {
    QuatMapConst symm_i(&symmetry_buf[4*ii]);
    q = symm_i * qmis;
    q.normalize();
    a = AxAng(q);
    if (a.angle() > M_PI) {
      a.angle() = 2*M_PI - a.angle();
    }
    if (a.angle() < wmin) {
      qmin = q;
      wmin = a.angle();
    }
  }

  return wmin * static_cast<double>(180 / M_PI);
} 

double Crystallography::get_cubic_misorientation_angle(int i, int j) {
  Quat misori;
  double wmin = 999999;
  QuatMapConst q_i(&quaternion_buf[4*i]);
  QuatMapConst q_j(&quaternion_buf[4*j]);
  misori = q_i * q_j.conjugate();
  return wmin;
}

void Crystallography::initialize_cubic_symmetry() {
  n_symm = 24;
  symmetry_buf = new double[4 * n_symm];
  int i = 0;
  double* s = &symmetry_buf[0];
  initialize_quat(&s[4*i++], 0,0,0,1);
  initialize_quat(&s[4*i++], 1,0,0,0);
  initialize_quat(&s[4*i++], 0,1,0,0);
  initialize_quat(&s[4*i++], 0,0,1,0);
  
  initialize_quat(&s[4*i++],  sqrt(2),0,0,sqrt(2));
  initialize_quat(&s[4*i++],  0,sqrt(2),0,sqrt(2));
  initialize_quat(&s[4*i++],  0,0,sqrt(2),sqrt(2));
  initialize_quat(&s[4*i++], -sqrt(2),0,0,sqrt(2));
  
  initialize_quat(&s[4*i++],  0,-sqrt(2),0,sqrt(2));
  initialize_quat(&s[4*i++],  0,0,-sqrt(2),sqrt(2));
  initialize_quat(&s[4*i++],  sqrt(2),sqrt(2),0,0);
  initialize_quat(&s[4*i++], -sqrt(2),sqrt(2),0,1);
  
  initialize_quat(&s[4*i++], 0,sqrt(2),sqrt(2),0);
  initialize_quat(&s[4*i++], 0,-sqrt(2),sqrt(2),0);
  initialize_quat(&s[4*i++], sqrt(2),0,sqrt(2),0);
  initialize_quat(&s[4*i++], -sqrt(2),0,sqrt(2),0);

  initialize_quat(&s[4*i++],  0.5,  0.5,  0.5,  0.5);
  initialize_quat(&s[4*i++], -0.5, -0.5, -0.5,  0.5);
  initialize_quat(&s[4*i++],  0.5, -0.5,  0.5,  0.5);
  initialize_quat(&s[4*i++], -0.5,  0.5, -0.5,  0.5);

  initialize_quat(&s[4*i++], -0.5,  0.5,  0.5,  0.5);
  initialize_quat(&s[4*i++],  0.5, -0.5, -0.5,  0.5);
  initialize_quat(&s[4*i++], -0.5, -0.5,  0.5,  0.5);
  initialize_quat(&s[4*i++],  0.5,  0.5, -0.5, -0.5);

}

void Crystallography::initialize_quat(double *dbuf,
				      double x, double y, double z, double w) {
  QuatMap q(dbuf);
  q.x() = x;
  q.y() = y;
  q.z() = z;
  q.w() = w;
  q.normalize();
}
