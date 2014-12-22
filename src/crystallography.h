
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

#ifndef SPK_CRYSTALLOGRAPHY_H
#define SPK_CRYSTALLOGRAPHY_H

#include <Eigen/Geometry>

namespace SPPARKS_NS {

  typedef Eigen::Quaternion<double> Quat;
  typedef Eigen::Map<Quat> QuatMap; 
  typedef Eigen::Map<const Quat> QuatMapConst;
  typedef Eigen::AngleAxis<double> AxAng;
  
class Crystallography {
 public:
  double *quaternion_buf;
  int n_orientations;
  
  int n_symm;
  double *symmetry_buf;

  double e_theta_max;
  double m_theta_max;
  
  Crystallography();
  ~Crystallography();
  double energy(int,int);
  double mobility(int,int);
  double misorientation(int,int);
  void copy_euler_angle_data(float*, int);
  void copy_quaternion_data(float*, int);
  void use_read_shockley(double);
  void use_hwang_humphreys(double);
  void setup_misorientation(char*);
  
 protected:
  double* misorientation_buf;
  
  double (Crystallography::*energy_pt)(int,int);
  double (Crystallography::*mobility_pt)(int,int);
  double (Crystallography::*misorientation_pt)(int,int);
  
  double unity(int,int);

  double read_shockley_energy(int,int);
  double hwang_humphreys_mobility(int,int);

  double precomputed_misorientation_angle(int,int);
  double cached_misorientation_angle(int,int);
  double compute_misorientation_angle(int,int);
  double compute_cubic_misorientation_angle(int,int);
  
  Quat misori(Quat,Quat,int,double*);
  void load_symmetry_operator();
  void initialize_cubic_symmetry();
  void initialize_quat(double*,double,double,double,double);
  
  bool cubic_FZ_test(Quat);
  double* load_euler_orientations_as_quats(char*);
  void quat_from_Bunge(double,double,double,double*);
  void compute_misorientation_angles(std::string);
  double* load_symm_table(int &N_symm, std::string);
};

}

#endif
