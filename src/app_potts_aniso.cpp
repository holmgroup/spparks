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

#include "spktype.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "app_potts_aniso.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"

#include <map>
#include <algorithm>
#include <fstream>
#include "math_const.h"

#include <Eigen/Geometry>

using namespace SPPARKS_NS;
using Eigen::Quaternion;

/* ---------------------------------------------------------------------- */

AppPottsAniso::AppPottsAniso(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  if (narg != 2) error->all(FLERR,"Illegal app_style command");

  // default to isotropic potts model
  energy = &AppPottsAniso::uniform_energy;
  mobility = &AppPottsAniso::uniform_mobility;
  e_table = m_table = NULL;
  ori_table = misori_table = NULL;
}

/* ---------------------------------------------------------------------- */

AppPottsAniso::~AppPottsAniso()
{
  // sites and unique are deleted by the ~AppPotts virtual destructor
  // delete [] sites;
  // delete [] unique;
  delete [] e_table;
  delete [] m_table;
  delete [] ori_table;
  delete [] misori_table;
}


/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppPottsAniso::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"load_mobility") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal load_mobility command");
    // extract the filename for mobility lookup table
    char *m_filename = arg[0];
    // set fn pointer (*mobility) to the lookup function
    mobility = &AppPottsAniso::lookup_mobility;
    // load the mobility lookup table *m_table
    fprintf(stdout, "loading mobility lookup table...");
    m_table = load_table(m_filename);
    if (m_table == NULL) fprintf(stdout, "problem allocating m_table\n");
    fprintf(stdout, " finished.\n");
  }
  else if (strcmp(command, "load_energy") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal load_energy command");
    // extract the filename for energy lookup table
    char *e_filename = arg[0];
    // set fn pointer (*energy) to the lookup function
    energy = &AppPottsAniso::lookup_energy;
    // load the energy lookup table *e_table
    e_table = load_table(e_filename);
  }
  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsAniso::init_app()
{
  delete [] sites;
  delete [] unique;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (spin[i] < 1 || spin[i] > nspins+1) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute energy of site using the pairwise energy function
   pointed to by (*energy)
------------------------------------------------------------------------- */

double AppPottsAniso::site_energy(int i)
{
  int isite = spin[i];
  int jsite = 0;
  double eng = 0;
  double eng_inc = 0;
  for (int j = 0; j < numneigh[i]; j++) {
    jsite = spin[neighbor[i][j]];
    if (isite != jsite) {
      eng += (this->*energy)(isite, jsite);
    }
  }
  return eng;
}
// double AppPottsAniso::site_energy(int i)
// {
//   int isite = spin[i];
//   int eng = 0;
//   for (int j = 0; j < numneigh[i]; j++)
//     if (isite != spin[neighbor[i][j]]) eng++;
//   return (double) eng;
// }


/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   flip to random spin from 1 to nspins
------------------------------------------------------------------------- */

void AppPottsAniso::site_event_rejection(int i, RandomPark *random)
{
 
  int oldstate = spin[i];
  double einitial = site_energy(i);

  // event = random spin from 1 to nspins, including self

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  spin[i] = iran;
  double efinal = site_energy(i);
  double M = (this->*mobility)(spin[i], oldstate);
  // accept or reject via Boltzmann criterion
  // null bin extends to nspins

  if (efinal <= einitial && random->uniform() < M) {
    ;
  } else if (temperature == 0.0) {
    spin[i] = oldstate;
  } else if (random->uniform() > M * exp((einitial-efinal)*t_inverse)) {
    spin[i] = oldstate;
  }

  if (spin[i] != oldstate) naccept++;

  // set mask if site could not have changed
  // if site changed, unset mask of sites with affected propensity
  // OK to change mask of ghost sites since never used

  if (Lmask) {
    if (einitial < 0.5*numneigh[i]) mask[i] = 1;
    if (spin[i] != oldstate)
      for (int j = 0; j < numneigh[i]; j++)
	mask[neighbor[i][j]] = 0;
  }
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppPottsAniso::site_propensity(int i)
{
  // events = spin flips to neighboring site different than self
  // disallow wild flips = flips to value different than all neighs

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == spin[i] || value > nspins) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  // for each flip:
  // compute energy difference between initial and final state
  // if downhill or no energy change, propensity = 1
  // if uphill energy change, propensity = Boltzmann factor

  int oldstate = spin[i];
  double einitial = site_energy(i);
  double efinal;
  double prob = 0.0;
  double M = 0.0;

  for (m = 0; m < nevent; m++) {
    spin[i] = unique[m];
    efinal = site_energy(i);
    M = (this->*mobility)(spin[i], oldstate);
    if (efinal <= einitial) prob += M;
    else if (temperature > 0.0) prob += M * exp((einitial-efinal)*t_inverse);
  }

  spin[i] = oldstate;
  return prob;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppPottsAniso::site_event(int i, RandomPark *random)
{
  int j,m,value;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  // perform event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double efinal;

  int oldstate = spin[i];
  double einitial = site_energy(i);
  double prob = 0.0;
  int nevent = 0;
  double M = 0.0;

  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == oldstate || value > nspins) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;

    spin[i] = value;
    efinal = site_energy(i);
    M = (this->*mobility)(spin[i], oldstate);
    if (efinal <= einitial) prob += M;
    else if (temperature > 0.0) prob += M * exp((einitial-efinal)*t_inverse);
    if (prob >= threshhold) break;
  }

  // compute propensity changes for self and neighbor sites
  // ignore update of neighbor sites with isite < 0

  int nsites = 0;
  int isite = i2site[i];
  sites[nsites++] = isite;
  propensity[isite] = site_propensity(i);

  for (j = 0; j < numneigh[i]; j++) {
    m = neighbor[i][j];
    isite = i2site[m];
    if (isite < 0) continue;
    sites[nsites++] = isite;
    propensity[isite] = site_propensity(m);
  }

  solve->update(nsites,sites,propensity);
}

/* ----------------------------------------------------------------------
   return uniform grain boundary energy
------------------------------------------------------------------------- */

double AppPottsAniso::uniform_energy(int ispin, int jspin) {
  return 1;
}

/* ----------------------------------------------------------------------
   return uniform grain boundary mobility
------------------------------------------------------------------------- */

double AppPottsAniso::uniform_mobility(int ispin, int jspin) {
  return 1;
}

/* ----------------------------------------------------------------------
   Look up the energy of the (ispin, jspin) grain boundary
------------------------------------------------------------------------- */

double AppPottsAniso::lookup_energy(int ispin, int jspin) {
  int r,c;
  c = std::min(ispin, jspin);
  r = std::max(ispin, jspin);
  int address = ((r * (r+1)) / 2) + c;
  return e_table[address];
}

/* ----------------------------------------------------------------------
   Lookup the mobility of the (ispin, jspin) grain boundary
------------------------------------------------------------------------- */

double AppPottsAniso::lookup_mobility(int ispin, int jspin) {
  int r,c;
  c = std::min(ispin, jspin);
  r = std::max(ispin, jspin);
  int address = ((r * (r+1)) / 2) + c;
  return m_table[address];
}

/* ----------------------------------------------------------------------
   allocate sparse triangular lookup table and populate from file
   including the diagonal elements
   lookup table should contain entries for spin == 0
   this way special energies can be assigned to particle interfaces
   and it's easier to think about zero-based arrays
------------------------------------------------------------------------- */

double *AppPottsAniso::load_table(char *filename) {

  // fragile header-parsing code, updates nspins
  std::string line;
  std::ifstream infile(filename);
  if (!infile.is_open()) error->all(FLERR,"Could not open lookup table file");
  std::getline(infile, line); // ignore the first line
  infile >> nspins >> line; // second line of file header is "$nspins spins"

  int table_size = (nspins+1) + (nspins* (nspins+1)) / 2;
  double *table;
  table = new double[table_size];

  double value = -1;
  for (int r = 0; r < nspins+1; r++) {
    for (int c = 0; c <= r; c++) {
      int address = ((r * (r+1)) / 2) + c;
      infile >> value;
      if (value < 0.0 || value > 1.0) error->all(FLERR,"Improperly specified lookup table");
      table[address] = value;
    }
  }

  return table;
}

/* ----------------------------------------------------------------------
   Populate e_table with Read-Shockley energies
   based on misorientation angles in misori_table
------------------------------------------------------------------------- */

double* AppPottsAniso::read_shockley_table() {
  int table_size = (nspins+1) + (nspins* (nspins+1)) / 2;
  double* table = NULL;
  table = new double[table_size];
  memset(table, 0, table_size*sizeof(*table));

  for (int r = 0; r < nspins+1; r++) {
    for (int c = 0; c <= r; c++) {
      int address = ((r * (r+1)) / 2) + c;
      table[address] = read_shockley_energy(misori_table[address]);
    }
  }

  return table;
}

/* ----------------------------------------------------------------------
   Compute the Read-Shockley energy for a grain boundary with
   minimum misorientation angle misori_angle (radians)
------------------------------------------------------------------------- */
double AppPottsAniso::read_shockley_energy(double misori_angle) {
  double RS_energy = 0;
  return RS_energy;
}

/* ----------------------------------------------------------------------
   Compute a triangular table of misorientation angles in radians 
------------------------------------------------------------------------- */

void AppPottsAniso::compute_misorientation_angles(std::string symmetry_filename) {
  // load the symmetry operators from a file
  // containing (x,y,z,w) unit quaternions
  double* symm = NULL;
  int N_symm = 0;
  symm = load_symm_table(N_symm, symmetry_filename);

  int table_size = (nspins+1) + (nspins* (nspins+1)) / 2;
  misori_table = new double[table_size];
  memset(misori_table, 0, table_size*sizeof(*misori_table));

  Quaternion<double> qmis;
  double angle;
  // skip over entries involving spin = 0 -- start at 1
  for (int r = 1; r < nspins+1; r++) {
    for (int c = 1; c <= r; c++) {
      Eigen::Map<Quaternion<double> > ori_r(&ori_table[r]);
      Eigen::Map<Quaternion<double> > ori_c(&ori_table[c]);
      qmis = misori(ori_r, ori_c, N_symm, symm);
      angle = 2 * acos(qmis.w());
      int address = ((r * (r+1)) / 2) + c;
      misori_table[address] = angle;
    }
  }
  return;
}

/* ----------------------------------------------------------------------
   Compute the minimum misorientation using the unit quaternion method
   Returns the full misorientation in the fundamental zone.
------------------------------------------------------------------------- */

Quaternion<double> AppPottsAniso::misori(Quaternion<double> ori_a, Quaternion<double> ori_b, int N_symm, double *symm) {
  Quaternion<double> min_misori(0,0,0,0);
  Quaternion<double> m(0,0,0,0);
  Quaternion<double> temp;
  double max_w = 0.0;

  // calculate misorientation: g_b inv(g_a)
  m = ori_b * ori_a.inverse();
  m = m.normalized();
  // for the i-th crystal symmetry operator calculate O_i mis
  for (int i = 0; i < N_symm; i++) {
    Eigen::Map<Quaternion<double> > i_temp(&symm[i*4]); 
    // for the j-th crystal symmetry operator calculate O_i mis O_j
    for (int j = 0; j < N_symm; j++) {
      Eigen::Map<Quaternion<double> > j_temp(&symm[j*4]);
      temp = i_temp * m * j_temp;
      // (x,y,z,w) == (-x,-y,-z,-w)
      // looping over id_sign accounts for that
      for (int id_sign = 0; id_sign < 2; id_sign++) {
	temp.x() = -temp.x();
	temp.y() = -temp.y();
	temp.z() = -temp.z();
	temp.w() = -temp.w();
	// Apply switching symmetry
	for (int id_switch = 0; id_switch < 2; id_switch++) {
	  temp = temp.inverse();
	  if (cubic_FZ_test(temp)) {
	    // test the angle for minimum criterion
	    if (temp.w() > max_w) {
	      max_w = temp.w();
	      min_misori = temp;
	    }
	  }
	}
      }
    }
  }
  return min_misori;
}

/* ----------------------------------------------------------------------
   test the axis for cubic FZ membership
   this filters to the standard stereographic triangle for cubics
------------------------------------------------------------------------- */

bool AppPottsAniso::cubic_FZ_test(Quaternion<double> quat) {
  
  return ( 0 <= quat.x()
	   && quat.x() <= quat.y()
	   && quat.y() <= quat.z()
	   && quat.z() <= quat.w());
}

/* ----------------------------------------------------------------------
   read symmetry operator from file -- Professor Rollett's quat.symm.${type}
   order: "N_variants"
   body: "x y z w" quaternion components
------------------------------------------------------------------------- */

double *AppPottsAniso::load_symm_table(int &N_symm, std::string infile_name) {
  std::string header_line;
  std::fstream symmfile(infile_name.c_str());
  std::getline(symmfile, header_line);
  symmfile >> N_symm;

  double* symm;
  symm = new double[4*N_symm];
  std::memset(symm, 0, 4*N_symm*sizeof(*symm));

  double x,y,z,w = 0.0;
  for (int i = 0; i < N_symm; i++) {
    int offset = i*4;
    symmfile >> symm[offset] >> symm[offset+1] >> symm[offset+2] >> symm[offset+3]; 
    Eigen::Map<Quaternion<double> > quat(&symm[offset]);
    quat.normalize();
  }
	   
  return symm;
}


/* ----------------------------------------------------------------------
   load Bunge convention Euler angles (degrees) from file
   into a flat array of quaternions {x,y,z,w}
   fragile header-parsing code, updates nspins
------------------------------------------------------------------------- */

double *AppPottsAniso::load_euler_orientations_as_quats(char *filename) {
  std::string line;
  std::ifstream infile(filename);
  if (!infile.is_open()) error->all(FLERR,"Could not open lookup table file");
  std::getline(infile, line); // ignore the first line
  infile >> nspins >> line; // second line of file header is "$nspins spins"

  int table_size = 4 * (nspins+1);
  ori_table = new double[table_size];
  std::memset(ori_table, 0, table_size*sizeof(*ori_table));
  
  double phi_1, Phi, phi_2 = 0.0;
  // ori_table[0-3] is reserved for spin == 0
  for (int id_spin = 1; id_spin < nspins+1; id_spin++) {
    infile >> phi_1; infile >> Phi; infile >> phi_2;
    if (phi_1 < 0.0 || phi_1 > 360.0) error->all(FLERR,"Improperly specified Euler angle");
    if (Phi < 0.0 || Phi > 180.0) error->all(FLERR,"Improperly specified Euler angle");
    if (phi_2 < 0.0 || phi_2 > 360.0) error->all(FLERR,"Improperly specified Euler angle");

    int offset = 4 * id_spin;
    quat_from_Bunge(phi_1, Phi, phi_2, &(ori_table[offset]));
  }

  return ori_table;
}

/* ----------------------------------------------------------------------
quat_from_Bunge
accepts Bunge convention Euler angles in degrees
computes a quaternion representation 
------------------------------------------------------------------------- */

void AppPottsAniso::quat_from_Bunge(double phi_1, double Phi, double phi_2, double *xyzw) {
  const double rad = MathConst::MY_PI / 180.0;
  phi_1 = phi_1 * rad;
  Phi   = Phi   * rad;
  phi_2 = phi_2 * rad;
  
  double x, y, z, w = 0.0;
  x = sin(Phi/2) * cos((phi_1 - phi_2) / 2);
  y = sin(Phi/2) * sin((phi_1 - phi_2) / 2);
  z = cos(Phi/2) * sin((phi_1 + phi_2) / 2);
  w = cos(Phi/2) * cos((phi_1 + phi_2) / 2);

  xyzw[0] = x;
  xyzw[1] = y;
  xyzw[2] = z;
  xyzw[3] = w;

  return;
}
