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
  ~AppPottsAniso();
  void input_app(char *, int, char **);
  void init_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *);
  double site_propensity(int);
  void site_event(int, class RandomPark *);

  double *e_table, *m_table, *ori_table, *symm_table;

 private:
  double (AppPottsAniso::*energy)(int,int);
  double (AppPottsAniso::*mobility)(int,int);
  double uniform_energy(int,int);
  double uniform_mobility(int,int);
  double lookup_energy(int,int);
  double lookup_mobility(int,int);

  double* load_table(char*);
  double* load_euler_orientations_as_quats(char*);
  void    quat_from_Bunge(double,double,double,double*);
  double* load_symm_table(char*);
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
