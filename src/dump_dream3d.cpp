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
#include <iostream>

using namespace SPPARKS_NS;
using namespace MathConst;

enum{DREAM3D,PLAIN_H5};
enum{PARALLEL,SERIAL};
enum{DREAM3D4,DREAM3D6}; // DREAM3D file format major versions
enum{ID,SITE,X,Y,Z,ENERGY,PROPENSITY,IARRAY,DARRAY};  // also in dump_text
enum{INT,DOUBLE,BIGINT};                              // also in dump_text


enum{SPHERE,SQUARE}; // to be removed
enum{JPG,PPM};

/* ---------------------------------------------------------------------- */

DumpDream3d::DumpDream3d(SPPARKS *spk, int narg, char **arg) :
  DumpText(spk, narg, arg)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

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

  int def = 0;
  char **argcopy;
  
  if (narg == 4) {
    def = 1;
    narg = 9;
    argcopy = new char*[narg];
    argcopy[4] = (char *) "id";
    argcopy[5] = new char[6];
    strcpy(argcopy[5],"site");
    argcopy[6] = (char *) "x";
    argcopy[7] = (char *) "y";
    argcopy[8] = (char *) "z";
  } else argcopy = arg;
  
  size_one = 2;
  pack_choice = new FnPtrPack[size_one];
  ioptional = parse_fields(narg,argcopy);
  if (size_one != 2) error->all(FLERR,"Illegal dump dream3d command");

  // set defaults for optional args
  io_style = SERIAL;
  major_version = DREAM3D6;
  dataset_name = "SpparksVolume";

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"io_style") == 1) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal dump dream3d command");
      if (strcmp(arg[iarg+1],"parallel") == 0) io_style = PARALLEL;
      else if (strcmp(arg[iarg+1],"serial") == 0) io_style = SERIAL;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"version") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal dump dream3d command: missing version arg");
      if (strcmp(arg[iarg+1],"4.0") == 0) major_version = DREAM3D4;
      else if (strcmp(arg[iarg+1],"6.0") == 0) major_version = DREAM3D6;
      else error->all(FLERR,"only DREAM3D formats 4.0 and 6.0 supported");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"dataset") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal dump dream3d command: missing version arg");
      dataset_name = arg[iarg+1];
      iarg += 2;
    }
    else error->all(FLERR,"Illegal dump dream3d command");
  }

  set_dataset_paths();

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
  if (multifile == 0) error->all(FLERR,"Dump dream3d requires one snapshot per file until the DREAM3D format supports timeseries data");

  // keep this here for consistency, though it's probably superfluous to DumpDream3d
  // DumpText::init_style();

  // check variables
  // no variables for DumpDream3d right now...
}

/* ---------------------------------------------------------------------- */

void DumpDream3d::set_dataset_paths() {
  switch (major_version){
  case DREAM3D4 :
    version_string = "4.0";
    dataset_root = "/VoxelDataContainer";
    voxels_path = dataset_root + "/" + "CELL_DATA";
    grain_ids_path = voxels_path + "/" + "GrainIds";
    break;
  case DREAM3D6 :
    version_string = "6.0";
    dataset_root = "/DataContainers/" + dataset_name;
    voxels_path = dataset_root + "/" + "CellData";
    grain_ids_path = voxels_path + "/" + "FeatureIds";
    break;
  default:
    error->all(FLERR,"Only DREAM3D format versions 4 and 6 supported");
  }
}

void DumpDream3d::create_groups() {
#ifdef SPPARKS_HDF5
  if (major_version == DREAM3D6) {
    h5_status = H5Gcreate(output_file, "/DataContainers",
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  
  h5_status = H5Gcreate(output_file, dataset_root.c_str(),
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5_status = H5Gcreate(output_file, voxels_path.c_str(),
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
}

void DumpDream3d::set_attrs() {
#ifdef SPPARKS_HDF5
  // a hacked together way to set attributes DREAM3D needs
  long container_type[1] = {0};
  h5_status = H5LTset_attribute_long(output_file, dataset_root.c_str(),
				     "DataContainerType", container_type, 1);
  long matrix_type[1] = {3};
  h5_status = H5LTset_attribute_long(output_file, voxels_path.c_str(),
				     "AttributeMatrixType", matrix_type, 1);
  long dims[3] = {0,0,0};
  std::string dimensions_path = dataset_root + "/" + "DIMENSIONS";
  h5_status = H5LTread_dataset_long(output_file, dimensions_path.c_str(), dims);
  h5_status = H5LTset_attribute_long(output_file, voxels_path.c_str(),
				     "TupleDimensions", dims, 3);

  int data_version[1] = {2};
  h5_status = H5LTset_attribute_int(output_file, grain_ids_path.c_str(),
				     "DataArrayVersion", data_version, 1);
  std::string object_type = "DataArray<int32_t>";
  h5_status = H5LTset_attribute_string(output_file, grain_ids_path.c_str(),
				    "ObjectType", object_type.c_str());
  int num_cpts[1] = {1};
  h5_status = H5LTset_attribute_int(output_file, grain_ids_path.c_str(),
				    "NumComponents", num_cpts, 1);
  h5_status = H5LTset_attribute_long(output_file, grain_ids_path.c_str(),
				     "TupleDimensions", dims, 3);
  long cpt_dims[1] = {1};
  h5_status = H5LTset_attribute_long(output_file, grain_ids_path.c_str(),
				    "ComponentDimensions", cpt_dims, 1);
  char sbuf[100];
  sprintf(sbuf, "%s%d%s%d%s%d", "x=", dims[0], ",y=", dims[1], ",z=", dims[2]);
  h5_status = H5LTset_attribute_string(output_file, grain_ids_path.c_str(),
				    "Tuple Axis Dimensions", sbuf);
#endif
}

void DumpDream3d::write(double time) {
  // open new file
  create_hdf5_file();
  idump++;

  // write data to hdf5 file
  switch (filetype) {
  case PLAIN_H5:
    write_h5();
    break;
  case DREAM3D:
    write_dream3d();
    break;
  default:
    write_dream3d();
  }

  close_hdf5_file();
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
    std::cout << "dump_dream3d filename: " << filecurrent << std::endl;
#ifdef SPPARKS_HDF5
    output_file = H5Fcreate(filecurrent, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5LTset_attribute_string(output_file, "/", "FileVersion", version_string.c_str());
    H5LTset_attribute_string(output_file, "/", "FileOrigin", "SPPARKS with HDF5");
    std::string dream3d_build = "5.1.709.e69fafd";
    H5LTset_attribute_string(output_file, "/", "DREAM3D Version", dream3d_build.c_str());
    create_groups();
#endif
  }
}

void DumpDream3d::close_hdf5_file() {
  if (me == 0) {
#ifdef SPPARKS_HDF5
    h5_status = H5Fclose(output_file);
#endif
  }
}

// void write_data(n, double* buf) {
//   int i,j;

//   int m = 0;
//   for (i = 0; i < n; i++) {
//     for (j = 0; j < size_one; j++) {
//       if (vtype[j] == INT)
// 	fprintf(fp,vformat[j],static_cast<int> (buf[m]));
//       else if (vtype[j] == DOUBLE)
// 	fprintf(fp,vformat[j],buf[m]);
//       else if (vtype[j] == TAGINT) 
// 	fprintf(fp,vformat[j],static_cast<tagint> (buf[m]));
//       m++;
//     }
//   }
// }

void unpack_data(int* data, int nlines, double* buf) {
  int m = 0;
  for (int i = 0; i < nlines; i++) {
    int id = static_cast<int>(buf[m]) - 1; // site ids offset by one
    data[id] = static_cast<int>(buf[m+1]);
    m += 2;
  }
}

void DumpDream3d::write_dream3d() {
#ifdef SPPARKS_HDF5
  // Specify the dataset dimensions
  hsize_t dims[1] = {3};
  long extents[3] = {static_cast<long>(domain->boxxhi - domain->boxxlo),
		     static_cast<long>(domain->boxyhi - domain->boxylo),
		     static_cast<long>(domain->boxzhi - domain->boxzlo)};

  std::string dimensions_path = dataset_root + "/" + "DIMENSIONS";
  h5_status =  H5LTmake_dataset_long(output_file, dimensions_path.c_str(), 1, dims, extents);
  
  // Specify the origin
  float origin[3] = {0,0,0};
  std::string origin_path = dataset_root + "/" + "ORIGIN";
  h5_status = H5LTmake_dataset_float(output_file, origin_path.c_str(), 1, dims, origin);

  // Specify spacing in microns (Analagous to DREAM.3D)
  // default to 0.5 micron spacing, as in DREAM.3D
  float spacing[3] = {0.5, 0.5, 0.5};
  std::string spacing_path = dataset_root + "/" + "SPACING";
  h5_status = H5LTmake_dataset_float(output_file, spacing_path.c_str(), 1, dims, spacing);

  // write grain ids
  int* data = NULL;
  if (me == 0) {
    data = new int[app->nglobal];
    for (int i = 0; i < app->nglobal; i++)
      data[i] = 0;
  }
      // nme = # of atoms this proc will contribute to dump
  // pack buf with id, site
  
  int nme = count();

  if (nme > maxbuf) {
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  pack();

  // copy data
  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;
  
  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = nme;

      unpack_data(data,nlines,buf);
    }
    
  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,world);
  }
  
    // assuming SPPARKS and DREAM3D both use row-major indexing, GrainIds are indexed by global_id
  if (me == 0) {
    const int rank = 4;
    hsize_t shape[rank] = {extents[0], extents[1], extents[2], 1};
    h5_status = H5LTmake_dataset_int(output_file, grain_ids_path.c_str(), rank, shape, data);
  }
#endif
  set_attrs();
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

