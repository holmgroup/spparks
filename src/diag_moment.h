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

namespace SPPARKS_NS {

class DiagMoment : public Diag {
 public:
  DiagMoment(class SPPARKS *, int, char **);
  ~DiagMoment() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  int latticeflag;
  class AppLattice *applattice;
  class AppOffLattice *appofflattice;
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