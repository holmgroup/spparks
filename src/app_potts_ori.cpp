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
#include "app_potts_ori.h"
#include "solve.h"
#include "random_mars.h"
#include "random_park.h"
#include "error.h"
#include "crystallography.h"

#include <map>
#include <algorithm>
#include <fstream>

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

AppPottsOri::AppPottsOri(SPPARKS *spk, int narg, char **arg) : 
  AppPotts(spk,narg,arg)
{
  if (narg != 2) error->all(FLERR,"Illegal app_style command");
  gb_props = new Crystallography();
}

/* ---------------------------------------------------------------------- */

AppPottsOri::~AppPottsOri()
{
  ;
}


/* ----------------------------------------------------------------------
   input script commands unique to this app
------------------------------------------------------------------------- */

void AppPottsOri::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"set_energy") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal energy command");
    if (strcmp(arg[0],"uniform") == 0 && narg == 1) {
      ;
    }
    else if (strcmp(arg[0],"read_shockley") == 0 && narg == 2) {
      double theta_max = static_cast<double>(atof(arg[1]));
      gb_props->use_read_shockley(theta_max);
      fprintf(stdout,"Using read shockley with high angle cutoff %s degrees\n", arg[1]);
    }
    else error->all(FLERR,"Illegal energy command");
  }
  else if (strcmp(command,"set_mobility") == 0 ) {
    if (narg < 1) error->all(FLERR,"Illegal mobility command");
    if (strcmp(arg[0],"uniform") == 0 && narg == 1) {
      // uniform mobility is the default behavior for SPPARKS::Crystallography
      ; 
    }
    else if (strcmp(arg[0],"hwang_humphreys") == 0) {
      // default values:
      double theta_max = 15;
      double n = 5;
      double d = 4;
      if (narg == 1) {
	fprintf(stdout,"Using Hwang/Humphreys with defaults: ");
      }
      else if (narg == 7) {
	fprintf(stdout,"Using Hwang/Humphreys with: ");
	int iarg = 1; // arg[0] was "hwang_humphreys"
	bool theta_set, n_set, d_set = false;

	while (iarg < narg) {
	  if (strcmp(arg[iarg],"theta_max") == 0){
	    iarg++;
	    theta_max = static_cast<double>(atof(arg[iarg++]));
	    fprintf(stdout,"theta_max = %f\n", theta_max);
	    theta_set = true;
	  }
	  else if (strcmp(arg[iarg],"n") == 0) {
	    iarg++;
	    n = static_cast<double>(atof(arg[iarg++]));
	    n_set = true;
	  }
	  else if (strcmp(arg[iarg],"d") == 0) {
	    iarg++;
	    d = static_cast<double>(atof(arg[iarg++]));
	    d_set = true;
	  }
	  else error->all(FLERR,"Error using Hwang/Humphreys: set theta_max, n, and d");
	}
	if (!(theta_set && n_set && d_set))
	  error->all(FLERR,"Error using Hwang/Humphreys: set each of  theta_max, n, and d");
      }
      else error->all(FLERR,"Illegal Hwang/Humphreys parameters");

      gb_props->use_hwang_humphreys(theta_max, n, d);
      fprintf(stdout,"theta_max = %f, ", theta_max);
      fprintf(stdout,"n = %f, ", n);
      fprintf(stdout,"d = %f\n", d);
    }
    else error->all(FLERR, "Illegal mobility command");
  }
  else if (strcmp(command,"precompute") == 0) {
    if (narg == 1) {
      if (strcmp(arg[0],"misorientation") == 0
	  || strcmp(arg[0],"energy") == 0
	  || strcmp(arg[0],"mobility") == 0) {
	gb_props->setup_precomputed(arg[0]);
      } else error->all(FLERR,"Illegal precompute option.");
    } else error->all(FLERR,"Illegal precompute command.");
  }
  else if (strcmp(command,"cache") == 0) {
    if (narg == 1) {
      if (strcmp(arg[0],"misorientation") == 0
	  || strcmp(arg[0],"energy") == 0
	  || strcmp(arg[0],"mobility") == 0)
	{
	  gb_props->setup_cached(arg[0]);
      } else error->all(FLERR,"Illegal cache option.");
    } else error->all(FLERR,"Illegal cache command.");
  }
  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppPottsOri::init_app()
{
  delete [] sites;
  delete [] unique;
  sites = new int[1 + maxneigh];
  unique = new int[1 + maxneigh];

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (spin[i] < 1 || spin[i] > nspins) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ----------------------------------------------------------------------
   compute energy of site using the pairwise energy function
   pointed to by (*energy)
------------------------------------------------------------------------- */

double AppPottsOri::site_energy(int i)
{
  int isite = spin[i];
  int jsite = 0;
  double eng = 0;
  for (int j = 0; j < numneigh[i]; j++) {
    jsite = spin[neighbor[i][j]];
    if (isite != jsite) {
      eng += gb_props->energy(isite, jsite);
    }
  }
  return eng;
}

/* ----------------------------------------------------------------------
   rKMC method
   perform a site event with null bin rejection
   flip to random spin from 1 to nspins
------------------------------------------------------------------------- */

void AppPottsOri::site_event_rejection(int i, RandomPark *random)
{
 
  int oldstate = spin[i];
  double einitial = site_energy(i);

  // event = random spin from 1 to nspins, including self

  int iran = (int) (nspins*random->uniform()) + 1;
  if (iran > nspins) iran = nspins;
  spin[i] = iran;
  double efinal = site_energy(i);
  double M = gb_props->mobility(spin[i], oldstate);
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

double AppPottsOri::site_propensity(int i)
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
  // if downhill or no energy change, propensity = Reduced mobility M
  // if uphill energy change, propensity = M * Boltzmann factor

  int oldstate = spin[i];
  double einitial = site_energy(i);
  double efinal;
  double prob = 0.0;
  double M = 0.0;

  for (m = 0; m < nevent; m++) {
    spin[i] = unique[m];
    efinal = site_energy(i);
    M = gb_props->mobility(spin[i], oldstate);
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

void AppPottsOri::site_event(int i, RandomPark *random)
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
    M = gb_props->mobility(spin[i], oldstate);
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
