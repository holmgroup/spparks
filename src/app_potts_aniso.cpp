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

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsAniso::AppPottsAniso(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  if (narg != 2) error->all(FLERR,"Illegal app_style command");

  // default to isotropic potts model
  energy = &AppPottsAniso::uniform_energy;
  mobility = &AppPottsAniso::uniform_mobility;
  e_table = m_table = NULL;
}

/* ---------------------------------------------------------------------- */

AppPottsAniso::~AppPottsAniso()
{
  // sites and unique are deleted by the ~AppPotts virtual destructor
  // delete [] sites;
  // delete [] unique;
  delete [] e_table;
  delete [] m_table;
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


double AppPottsAniso::uniform_energy(int ispin, int jspin) {
  return 1;
}

double AppPottsAniso::uniform_mobility(int ispin, int jspin) {
  return 1;
}

double AppPottsAniso::lookup_energy(int ispin, int jspin) {
  int r,c;
  c = std::min(ispin, jspin);
  r = std::max(ispin, jspin);
  int address = ((r * (r+1)) / 2) + c;
  return e_table[address];
}

double AppPottsAniso::lookup_mobility(int ispin, int jspin) {
  int r,c;
  c = std::min(ispin, jspin);
  r = std::max(ispin, jspin);
  int address = ((r * (r+1)) / 2) + c;
  return m_table[address];
}

double *AppPottsAniso::load_table(char *filename) {
  // allocate sparse triangular lookup table
  // including the diagonal elements
  // lookup table should contain entries for spin == 0
  // this way special energies can be assigned to particle interfaces
  // and it's easier to think about zero-based arrays
  int table_size = (nspins+1) + (nspins* (nspins+1)) / 2;
  double *table;
  table = new double[table_size];

  std::ifstream infile(filename);
  if (!infile.is_open()) error->all(FLERR,"Could not open lookup table file");
  double val = -1;
  for (int r = 0; r < nspins+1; r++) {
    for (int c = 0; c <= r; c++) {
      int address = ((r * (r+1)) / 2) + c;
      infile >> val;
      table[address] = val;
    }
  }

  return table;
}

