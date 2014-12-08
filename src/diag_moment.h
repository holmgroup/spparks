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

#ifdef DIAG_CLASS
DiagStyle(moment,DiagMoment)

#else

#ifndef SPK_DIAG_MOMENT_H
#define SPK_DIAG_MOMENT_H

#include "stdio.h"
#include "diag.h"
#include <map>
#include "grain.h"

namespace SPPARKS_NS {

class DiagMoment : public Diag {
 public:
  DiagMoment(class SPPARKS *, int, char **);
  virtual ~DiagMoment();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);
  double (DiagMoment::*x_dist)(double, double);
  double (DiagMoment::*y_dist)(double, double);
  double (DiagMoment::*z_dist)(double, double);

  double dist(double, double);
  double min_dist(double, double, double);
  double min_dist_x(double, double);
  double min_dist_y(double, double);
  double min_dist_z(double, double);
  void merge_grain(int,double,double,double,double,int,double*);
  
 protected:
  std::map<int, Grain> grains;
  std::map<int, Grain>::iterator grain_iter;
  FILE *fpdump;
  
 private:
  int latticeflag;
  class AppLattice *applattice;
  class CommLattice *comm;
  
  double x_size, y_size, z_size;
  double **xyz; // site coordinates
  tagint *idsite; // global site id

  int nlocal, nghost;
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag style incompatible with app style

The lattice styles of the diagnostic and the on-lattice application
must match.

*/
