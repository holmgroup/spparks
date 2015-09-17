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
DiagStyle(grainsize,DiagGrainsize)
#else

#ifndef SPK_DIAG_GRAINSIZE_H
#define SPK_DIAG_GRAINSIZE_H

#include "stdio.h"
#include "diag.h"
#include <map>
#include "grain.h"

#ifdef SPPARKS_HDF5
#include "H5Cpp.h"
#include "hdf5_hl.h"
#endif

namespace SPPARKS_NS {

class DiagGrainsize : public Diag {
 public:
  DiagGrainsize(class SPPARKS *, int, char **);
  virtual ~DiagGrainsize();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);
  double (DiagGrainsize::*x_dist)(double, double);
  double (DiagGrainsize::*y_dist)(double, double);
  double (DiagGrainsize::*z_dist)(double, double);

  double dist(double, double);
  double min_dist(double, double, double);
  double min_dist_x(double, double);
  double min_dist_y(double, double);
  double min_dist_z(double, double);
  void merge_grain(int,double,double,double,double,int,double*);
  void merge_partial_moments(int,double*);
  void create_hdf_file(char*);
  void find_nspins();
  
 protected:
  std::map<int, Grain> grains;
  std::map<int, Grain>::iterator grain_iter;
  
 private:
  int latticeflag;
  class AppLattice *applattice;
  class CommLattice *comm;
  
  double x_size, y_size, z_size;
  double **xyz; // site coordinates
  tagint *idsite; // global site id

  int nspins;
  int nlocal, nghost;
  int n_nearest;
  int* nearest_neighs;

#ifdef SPPARKS_HDF5
  hsize_t dims[2];
  hsize_t chunk_dims[2];
  H5::H5File output_file;
  H5::DataSpace *dataspace;
  H5::DSetCreatPropList prop;
  H5::DataSet *dataset;
  H5::DataSpace *filespace;
  H5::DataSpace *memspace;
  herr_t h5_status;
#endif
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Diag style incompatible with app style

The lattice styles of the diagnostic and the on-lattice application
must match.

*/
