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

#include "mpi.h"
#include "math.h"
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "dump_dream3d.h"
#include "image.h"
#include "app.h"
#include "app_lattice.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "variable.h"
#include "random_park.h"
#include "random_mars.h"
#include "math_const.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;
using namespace MathConst;

enum{DREAM3D,PLAIN_H5};
enum{PARALLEL,SERIAL};
enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};  // also in dump_text
enum{INT,DOUBLE,BIGINT};                              // also in dump_text


enum{SPHERE,SQUARE}; // to be removed
enum{JPG,PPM};

/* ---------------------------------------------------------------------- */

DumpDream3d::DumpDream3d(SPPARKS *spk, int narg, char **arg) :
  DumpText(spk, narg, arg)
{
  if (binary || multiproc) error->all(FLERR,"Invalid dump dream3d filename");

  // set filetype based on filename suffix

  int n = strlen(filename);
  if (strlen(filename) > 4 && strcmp(&filename[n-4],".dream3d") == 0)
    filetype = DREAM3D;
  else if (strlen(filename) > 5 && strcmp(&filename[n-5],".h5") == 0)
    filetype = PLAIN_H5;
  else filetype = DREAM3D;

#ifndef SPPARKS_HDF5
  error->all(FLERR,"Cannot dump dream3d file without HDF5 support");
#endif

  // if (size_one != 2) error->all(FLERR,"Illegal dump dream3d command");

  // set defaults for optional args
  write_style = SERIAL;

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"write_style") == 1) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal dump dream3d command");
      if (strcmp(arg[iarg+1],"parallel") == 0) write_style = PARALLEL;
      else if (strcmp(arg[iarg+1],"serial") == 0) write_style = SERIAL;
      iarg += 2;
    } else error->all(FLERR,"Illegal dump image command");
  }

  // error checks

  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Dump dream3d requires lattice app");
  applattice = (AppLattice *) app;

}

/* ---------------------------------------------------------------------- */

DumpDream3d::~DumpDream3d()
{
  // delete image;
}

/* ---------------------------------------------------------------------- */

void DumpDream3d::init_style()
{
  if (multifile == 0) error->all(FLERR,"Dump dream3d requires one snapshot per file");

  // keep this here for consistency, though it's probably superfluous to DumpDream3d
  DumpText::init_style();

  // check variables
  // no variables for DumpDream3d right now...
}

/* ---------------------------------------------------------------------- */

void DumpDream3d::write(double time) {
  // open new file
  create_hdf5_file();
  idump++;

  // nme = # of atoms this proc will contribute to dump
  // pack buf with x,y,z,color,diameter
  // set minmax color range if using color map
  // create my portion of image for my particles
  
  int nme = count();

  if (nme > maxbuf) {
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  pack();
  // if (scolor == DATTRIBUTE) image->color_minmax(nchoose,buf,size_one);

  // create image on each proc, then merge them

  // image->clear();
  // create_image();
  // image->merge();

  // write data to hdf5 file

  if (me == 0) {
    if (filetype == PLAIN_H5) write_h5();
    else write_dream3d();
#ifdef SPPARKS_HDF5
    h5_status = H5Fclose(output_file);
#endif
  }

}

/* Evil code replication, mostly to keep HDF5 code out of dump.h */

void DumpDream3d::create_hdf5_file() {

  /* replace '*' with current timestep */
  char *filecurrent;
  if (multifile == 0) filecurrent = filename;
  else {
    filecurrent = new char[strlen(filename) + 16];
    char *ptr = strchr(filename,'*');
    *ptr = '\0';
    if (padflag == 0) 
      sprintf(filecurrent,"%s%d%s",filename,idump,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,"%d");
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filename,idump,ptr+1);
    }
    *ptr = '*';
  }

  if (me == 0) {
    std::cout << filecurrent << std::endl;
#ifdef SPPARKS_HDF5
    output_file = H5Fcreate(filecurrent, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5LTset_attribute_string(output_file, "/", "FileVersion", "4.0");
    H5LTset_attribute_string(output_file, "/", "FileOrigin", "SPPARKS with HDF5");
#endif
  }
}

void DumpDream3d::write_dream3d() {
#ifdef SPPARKS_HDF5
  // Specify the dataset dimensions
  hsize_t dims[1] = {3};
  long extents[3] = {domain->boxxhi - domain->boxxlo,
		    domain->boxyhi - domain->boxylo,
		    domain->boxzhi - domain->boxzlo};
  h5_status = H5Gcreate(output_file, "/VoxelDataContainer", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5_status = H5LTmake_dataset_long(output_file, "/VoxelDataContainer/DIMENSIONS", 1, dims, extents);
 
  // Specify the origin
  float origin[3] = {0,0,0};
  h5_status = H5LTmake_dataset_float(output_file, "/VoxelDataContainer/ORIGIN", 1, dims, origin);

  // Specify spacing in microns (Analagous to DREAM.3D)
  // default to 0.5 micron spacing, as in DREAM.3D
  float spacing[3] = {0.5, 0.5, 0.5};
  h5_status = H5LTmake_dataset_float(output_file, "/VoxelDataContainer/SPACING", 1, dims, spacing);

#endif
}

void DumpDream3d::write_h5() {
  error->all(FLERR,"Plain hdf5 format not specified; use dream3d format for now.");
}


/* ---------------------------------------------------------------------- */

int DumpDream3d::modify_param(int narg, char **arg)
{
  int n = DumpText::modify_param(narg,arg);
  if (n) return n;
  
  return 0;
}

