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

#include <unordered_map>

using namespace SPPARKS_NS;

Crystallography::Crystallography() {

  quaternion_buf = NULL;
  
  /* default to uniform energy and mobility */
  this->energy_pt = &Crystallography::unity;
  this->mobility_pt = &Crystallography::unity;
  this->misorientation_pt = &Crystallography::compute_misorientation_angle;
   
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
  if (misorientation_buf != NULL) {
    delete [] misorientation_buf;
    misorientation_buf = NULL;
  }
}

/* ---- crystallographic energy/mobility functions ---- */
double Crystallography::unity(int i, int j) {
  return 1.0;
}

double Crystallography::read_shockley_energy(int i, int j) {
  double misori;
  if (i == j)
    return 0;
  misori = misorientation(i,j);

  // high angle cutoff
  misori = std::min(misori, e_theta_max);
  // minimum allowed energy
  misori = std::max(misori, 0.1);
  
  misori = misori / e_theta_max;
  // do read/shockley formula
  return std::max(misori * (1 - log(misori)), 0.4);
}

double Crystallography::hwang_humphreys_mobility(int i, int j) {
  double misori = 0;
  double mobility = 1;
  if (i == j)
    return 0;
  
  misori = misorientation(i, j);
  
  // high angle cutoff
  misori = std::min(misori, m_theta_max);
  // minimum allowed mobility
  misori = std::max(misori, 0.1);

  misori = misori / m_theta_max;
  
  // do the hwang/humphreys formula
  // HH_n and HH_d are the experimental parameters
  mobility = 1 - exp(-HH_n * pow(misori, HH_d));
  return mobility;
}

double Crystallography::binary_mobility(int i, int j) {
  double misori = 0;
  double mobility = 1;
  if (i == j)
    return 0;

  misori = misorientation(i,j);

  if (misori < m_theta_max)
    mobility = min_mobility;

  return mobility;  
}

/* ---------------------------------------------------- */


/* -- wrapper functions for function pointers -- */
double Crystallography::energy(int i, int j) {
  return (*this.*energy_pt)(i,j);
}

double Crystallography::mobility(int i, int j) {
  return (*this.*mobility_pt)(i,j);
}

double Crystallography::misorientation(int i, int j) {
  return (*this.*misorientation_pt)(i,j);
}

/* -------------------------------------------- */

/* ------ swap function pointers around ------- */
void Crystallography::use_read_shockley(double theta_max) {
  if (theta_max > 0)
    e_theta_max = theta_max;
  this->energy_pt = &Crystallography::read_shockley_energy;
}

void Crystallography::use_hwang_humphreys(double theta_max, double n, double d) {
  HH_n = 5;
  HH_d = 4;

  if (theta_max > 0)
    m_theta_max = theta_max;

  if (n > 0)
    HH_n = n;
  if (d > 0)
    HH_d = d;
  
  this->mobility_pt = &Crystallography::hwang_humphreys_mobility;
}

void Crystallography::use_binary_mobility(double theta_max, double m_min) {
  m_min = 0.001;
  
  if (theta_max > 0)
    m_theta_max = theta_max;

  if (min_mobility > 0)
    min_mobility = m_min;

  this->mobility_pt = &Crystallography::binary_mobility;
}

void Crystallography::setup_precomputed(char* style) {
  double (Crystallography::*compute)(int,int);
  if (strcmp(style,"misorientation") == 0) {
    fprintf(stdout,"allocating misorientation lookup table...\n");
    compute = &Crystallography::compute_misorientation_angle;
  }
  else if (strcmp(style,"energy") == 0) {
    fprintf(stdout,"allocating energy lookup table...\n");
    compute = &Crystallography::energy;
  }
  else if (strcmp(style,"mobility") == 0) {
    fprintf(stdout,"allocating mobility lookup table...\n");
    compute = &Crystallography::mobility;
  }

  // grain_ids are 1-indexed, and g
  int size = (n_orientations+1)*(n_orientations+1);
  double* buf = new double[size];
  
  int index = 0;
  for (int i = 0; i <= n_orientations; i++) {
    for (int j = 0; j <= n_orientations; j++) {
      index = i*(n_orientations+1) + j;
      if (i == j || i == 0 || j == 0)
	buf[index] = 0;
      else
	buf[index] = (*this.*compute)(i,j);
    }
  }
  if (strcmp(style,"misorientation") == 0) {
    misorientation_buf = buf;
    this->misorientation_pt = &Crystallography::precomputed_misorientation_angle;
  }
  else if (strcmp(style,"energy") == 0) {
    energy_buf = buf;
    this->energy_pt = &Crystallography::precomputed_energy;
  }
  else if (strcmp(style,"mobility") == 0) {
    mobility_buf = buf;
    this->mobility_pt = &Crystallography::precomputed_mobility;
  }
}

void Crystallography::setup_cached(char* style) {
  if (strcmp(style,"misorientation") == 0) {
    // set up cache data structure
    fprintf(stdout,"using hash table 'cache' for misorientations\n");
    this->misorientation_pt = &Crystallography::cached_misorientation_angle;
  }
}

/* -------------------------------------------- */

double Crystallography::precomputed_misorientation_angle(int i, int j) {
  int index = i * (n_orientations+1) + j;
  return misorientation_buf[index];
}

double Crystallography::cached_misorientation_angle(int i, int j) {
  Key_t key (std::min(i, j), std::max(i,j));
  cache_iter = mis_cache.find(key);
  if (cache_iter != mis_cache.end())
    return cache_iter->second;

  double misori = compute_misorientation_angle(i, j);
  mis_cache.insert(std::make_pair(key, misori));
  return misori;
}

double Crystallography::precomputed_energy(int i, int j) {
  int index = i * (n_orientations+1) + j;
  return energy_buf[index];
}

double Crystallography::cached_energy(int i, int j) {
  Key_t key (std::min(i, j), std::max(i,j));
  cache_iter = e_cache.find(key);
  if (cache_iter != e_cache.end())
    return cache_iter->second;

  double e = energy(i,j);
  e_cache.insert(std::make_pair(key, e));
  return e;
}

double Crystallography::precomputed_mobility(int i, int j) {
  int index = i * (n_orientations+1) + j;
  return mobility_buf[index];
}

double Crystallography::cached_mobility(int i, int j) {
  Key_t key (std::min(i, j), std::max(i,j));
  cache_iter = m_cache.find(key);
  if (cache_iter != m_cache.end())
    return cache_iter->second;

  double m = mobility(i,j);
  m_cache.insert(std::make_pair(key, m));
  return m;
}

/* --- misorientation calculations ---- */
double Crystallography::compute_misorientation_angle(int i, int j) {
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

double Crystallography::compute_cubic_misorientation_angle(int i, int j) {
  Quat misori;
  double wmin = 999999;
  QuatMapConst q_i(&quaternion_buf[4*i]);
  QuatMapConst q_j(&quaternion_buf[4*j]);
  misori = q_i * q_j.conjugate();
  return wmin;
}
/* -------------------------------------------------------- */

/* --- data loading stuff --- */
void Crystallography::copy_euler_angle_data(float *data, int num_orientations) {
  int data_size = 3 * num_orientations;
  float* euler_buf = new float[data_size];
  for (int i = 0; i < data_size; i++) {
    euler_buf[i] = data[i];
  }
}

void Crystallography::copy_quaternion_data(float *data, int n_buf) {
  n_orientations = n_buf - 1;
  int data_size = 4 * n_buf;
  quaternion_buf = new double[data_size];
  // grains are 1-indexed, but HDF5 is zero-indexed.
  // first quaternion in buffer is meaningless.
  for (int i = 0; i < 4; i++) {
    quaternion_buf[i] = 0;
  }
  for (int i = 1; i <= n_orientations; i++) {
    for (int j = 0; j < 4; j++) {
      quaternion_buf[(4*i)+j] = static_cast<double>(data[(4*i)+j]);
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
