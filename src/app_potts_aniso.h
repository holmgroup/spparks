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

#ifdef APP_CLASS
AppStyle(potts/aniso,AppPottsAniso)

#else

#include "app_potts.h"

namespace SPPARKS_NS {

class AppPottsAniso : public AppPotts {
 public:
  AppPottsAniso(class SPPARKS *, int, char **);
  ~AppPottsAniso() {}
  void input_app(char *, int, char **);
  void init_app();

  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);
  void push_new_site(int, int*, int, std::stack<int>*);

 private:
  double pfraction;
  int multi,nthresh;
  double* e_table,m_table;

  double (*energy)(int,int);
  double (*mobility)(int,int);
  double iso_energy(int,int);
  double iso_mobility(int,int);
  double lookup_energy(int,int);
  double lookup_mobility(int,int);

  void load_e_table(char* filename);
  void load_m_table(char* filename);
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

*/
