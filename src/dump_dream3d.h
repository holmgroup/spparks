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

#ifdef DUMP_CLASS

DumpStyle(dream3d,DumpDream3d)

#else

#ifndef SPK_DUMP_DREAM3D_H
#define SPK_DUMP_DREAM3D_H

#include "math.h"
#include <string>
#include "dump_text.h"

#ifdef SPPARKS_HDF5
#include "H5Cpp.h"
#include "hdf5_hl.h"
#endif

namespace SPPARKS_NS {

class DumpDream3d : public DumpText{
 public:
  DumpDream3d(class SPPARKS *, int, char**);
  ~DumpDream3d();
  
 private:
  int me,nprocs;
  
  int filetype;
  std::string version_string;
  std::string dataset_name;
  std::string dataset_root;
  std::string voxels_path;
  std::string grain_ids_path;
  
  int io_style;
  int major_version;
  
#ifdef SPPARKS_HDF5
  hid_t output_file;
  herr_t h5_status;
#endif

  class AppLattice *applattice;

  void init_style();
  int modify_param(int, char **);
  void set_dataset_paths();
  void create_groups();
  void set_attrs();
  void write(double);
  void create_hdf5_file();
  void close_hdf5_file();
  void write_dream3d();
  void write_h5();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump image filename

UNDOCUMENTED

E: Cannot dump JPG file

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Dump image with quantity application does not support

UNDOCUMENTED

E: Invalid dump image theta value

UNDOCUMENTED

E: Dump image persp option is not yet supported

UNDOCUMENTED

W: Using dump image boundary with spheres

UNDOCUMENTED

E: Dump image boundary requires lattice app

UNDOCUMENTED

E: Dump image drange must be set

UNDOCUMENTED

E: Dump image crange must be set

UNDOCUMENTED

E: Dump image requires one snapshot per file

UNDOCUMENTED

E: Variable name for dump image theta does not exist

UNDOCUMENTED

E: Variable for dump image theta is invalid style

UNDOCUMENTED

E: Variable name for dump image phi does not exist

UNDOCUMENTED

E: Variable for dump image phi is invalid style

UNDOCUMENTED

E: Variable name for dump image center does not exist

UNDOCUMENTED

E: Variable for dump image center is invalid style

UNDOCUMENTED

E: Variable name for dump image zoom does not exist

UNDOCUMENTED

E: Variable for dump image zoom is invalid style

UNDOCUMENTED

E: Variable name for dump image persp does not exist

UNDOCUMENTED

E: Variable for dump image persp is invalid style

UNDOCUMENTED

E: Invalid dump image zoom value

UNDOCUMENTED

E: Invalid dump image persp value

UNDOCUMENTED

E: Invalid color in dump_modify command

UNDOCUMENTED

E: Dump_modify scolor requires integer attribute for dump image color

UNDOCUMENTED

E: Dump_modify sdiam requires integer attribute for dump image diameter

UNDOCUMENTED

*/
