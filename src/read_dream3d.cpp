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
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "read_dream3d.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "app_potts_ori.h"
#include "domain.h"
#include "create_sites.h"
#include "error.h"
#include "memory.h"
#include "crystallography.h"

#include "input.h"

#include <map>
#include <string>

using namespace SPPARKS_NS;

#define MAXLINE 256
#define CHUNK 1024
#define DELTA 4
#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

ReadDream3d::ReadDream3d(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  buffer = new int[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;
  load_orientations = false;
  dataset_name = "SyntheticVolume";

#ifdef SPPARKS_HDF5
  H5::Exception::dontPrint();
#endif
}

/* ---------------------------------------------------------------------- */

ReadDream3d::~ReadDream3d()
{
  delete [] buffer;
  memory->sfree(arg);
}

/* ---------------------------------------------------------------------- */

void ReadDream3d::command(int narg, char **arg)
{
  if (app == NULL) error->all(FLERR,"read_dream3d command before app_style set");
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"filename") == 0) {
      iarg++;
      if (iarg < narg) {
	if (me == 0) {
	  input_file = arg[iarg];
	}
      } else error->all(FLERR,"Illegal diag_style cluster command");
    }
    else if (strcmp(arg[iarg],"load_ori") == 0) {
      load_orientations = true;
    }
    else if (strcmp(arg[iarg],"dataset") == 0) {
      iarg++;
      dataset_name = arg[iarg++];
    }
    else {
      error->all(FLERR,"Illegal read_dream3d command");
    }
    iarg++;
  }
  
  if (domain->dimension == 2 && domain->zperiodic == 0)
    error->all(FLERR,"Cannot run 2d simulation with nonperiodic Z dimension");
  if (domain->dimension == 1 && 
      (domain->yperiodic == 0 || domain->zperiodic == 0))
    error->all(FLERR,
	       "Cannot run 1d simulation with nonperiodic Y or Z dimension");

  if (app->appclass == App::LATTICE) {
    applattice = (AppLattice *) app;
    latticeflag = 1;
  } else if (app->appclass == App::OFF_LATTICE) {
    error->all(FLERR, "Cannot read DREAM3D file for off-lattice apps.");
  }

  if (load_orientations) {
    if ((strcmp(app->style,"potts/ori") == 0)) {
	app_potts_ori = (AppPottsOri *) app;
      }
    else
      error->all(FLERR, "Must use app_style potts/ori to load orientations from dream3d");
  }
    
  // proc 0 opens the dream3d file
  if (me == 0) {
    if (screen) fprintf(screen,"Reading dream3d file ...\n");
#ifdef SPPARKS_HDF5
    // make sure the file exists first
    fprintf(stdout,"input file: %s\n", input_file);
    file_id = H5Fopen(input_file, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
      std::string msg = "Could not open HDF5 file: " + std::string(input_file);
      error->all(FLERR, msg.c_str());
    }

    // get the DREAM3D file format version and broadcast it
    // lazy way to get the version string...
    char* file_version = new char[16];
    h5_status = H5LTget_attribute_string(file_id,"/", "FileVersion", file_version);

    major_version = file_version[0];

    if (!(major_version == '4' || major_version == '6' || major_version == '7'))
      error->all(FLERR,"Only DREAM3D data format versions 4, 6, and 7 supported.");
#else
    error->all(FLERR,"Compile with HDF5 and set SPPARKS_HDF5 to read dream3d files.");
#endif
  }

  MPI_Bcast(&major_version, 1, MPI_CHAR, 0, world);
  set_dataset_paths();
  
  // extract dimensions of simulation volume
  get_dimensions();

  // create a simulation box
  if (!domain->box_exist) {
    domain->boxxlo = boxxlo;
    domain->boxylo = boxylo;
    domain->boxzlo = boxzlo;
    domain->boxxhi = boxxhi;
    domain->boxyhi = boxyhi;
    domain->boxzhi = boxzhi;

    domain->set_box();
    domain->box_exist = 1;
    if (domain->dimension == 1) domain->procs2domain_1d();
    if (domain->dimension == 2) domain->procs2domain_2d();
    if (domain->dimension == 3) domain->procs2domain_3d();
  }

  // Create sites
  if (domain->lattice == NULL)
    error->all(FLERR,"Cannot create sites with undefined lattice");
  // execute this input-script command to abstract it from the user
  input->one("create_sites box");

  // Read site values
  read_grain_ids();

  // Read orientation data
  // do this if specified and if app_style allows orientations
  if (load_orientations)
    read_average_quaternions();

  // close file
  std::cout << "Finished reading dream3d file" << std::endl;
#ifdef SPPARKS_HDF5
  if (me == 0) {
    h5_status = H5Fclose(file_id);
  }
#endif
}

void ReadDream3d::set_dataset_paths() {
  if (major_version == '4') {
    dimensions_path = "/VoxelDataContainer/DIMENSIONS";
    grain_ids_path = "/VoxelDataContainer/CELL_DATA/GrainIds";
    quaternions_path = "/VoxelDataContainer/FIELD_DATA/AvgQuats";
  }
  else if (major_version == '6') {
    std::string dataset_root = "/DataContainers/" + dataset_name + "/";
    dimensions_path = dataset_root + "DIMENSIONS";
    grain_ids_path = dataset_root + "CellData/FeatureIds";
    quaternions_path = dataset_root + "CellFeatureData/AvgQuats";
  }
  else if (major_version == '7') {
    std::string dataset_root = "/DataContainers/" + dataset_name + "/";
    dimensions_path = dataset_root + "/_SIMPL_GEOMETRY/DIMENSIONS";
    grain_ids_path = dataset_root + "CellData/FeatureIds";
    quaternions_path = dataset_root + "CellFeatureData/AvgQuats";
  }
  // TODO: check that dataset exists
  hbool_t check_object_valid = false;
  if (!(H5LTpath_valid(file_id, grain_ids_path.c_str(), check_object_valid)))
    error->all(FLERR,"DREAM3D dataset path invalid");    
}

void ReadDream3d::get_dimensions() {
  int dims_buf[3] = {0, 0, 0};
#ifdef SPPARKS_HDF5
  if (me == 0) 
    h5_status = H5LTread_dataset_int(file_id, dimensions_path.c_str(), dims_buf);
#endif

  MPI_Bcast(dims_buf, 3, MPI_INT, 0, world);
  
  std::cout << "dataset dimensions: " << dims_buf[0] << " x " << dims_buf[1] << " x " << dims_buf[2];
  std::cout << std::endl;

  boxxlo = boxylo = boxzlo = -0.5;
  boxxhi = boxyhi = boxzhi = 0.5;

  if (dims_buf[0] > 1) {
    boxxlo = 0;
    boxxhi = dims_buf[1];
  }
  if (dims_buf[1] > 1) {
    boxylo = 0;
    boxyhi = dims_buf[1];
  }
  if (dims_buf[2] > 1) {
    boxzlo = 0;
    boxzhi = dims_buf[2];
  }
  return;
}

void ReadDream3d::read_grain_ids() {

  int* data = NULL;
#ifdef SPPARKS_HDF5
  if (me == 0) {
    data = new int[app->nglobal];
    // assuming SPPARKS and DREAM3D both use row-major indexing, GrainIds are indexed by global_id
    h5_status = H5LTread_dataset_int(file_id, grain_ids_path.c_str(), data);
  }
#endif

  int i,m,nchunk;
  tagint site_id;
  char *next;;
  int *buf;

  // put my entire list of owned site IDs in a hashmap

  std::map<tagint,int>::iterator loc;
  std::map<tagint,int> hash;

  tagint *id = app->id;
  int nlocal = app->nlocal;

  for (i = 0; i < nlocal; i++)
    hash.insert(std::pair<tagint,int> (id[i],i));

  // broadcast one CHUNK of values at a time
  // store site's values if I own its ID

  int nvalues = app->ninteger + app->ndouble;
  if (nvalues != 1)
    error->all(FLERR, "read_dream3d implemented only for apps with 1 INT per site.");
  char **values = new char*[nvalues];
  for (i=0; i < nvalues; i++) {
    values[i] = new char[256];
  }

  tagint nread = 0;
  tagint nglobal = app->nglobal;

  while (nread < nglobal) {
    if (nglobal-nread > CHUNK) nchunk = CHUNK;
    else nchunk = nglobal-nread;
    if (me == 0) {
      m = 0;
      for (i = 0; i < nchunk; i++) {
	// hard-coded single integer per site...
	buffer[m] = data[nread+m];
	m++;
      }
      // m++;
    }
    MPI_Bcast(&m,1,MPI_INT,0,world);
    MPI_Bcast(buffer,m,MPI_INT,0,world);

    buf = buffer;

    for (int idx = 0; idx < nchunk; idx++) {
      // site id starts at 1
      site_id = nread + idx + 1;
      loc = hash.find(site_id);

      if (loc != hash.end()) {
	for (m = 0; m < nvalues; m++) {
	  sprintf(values[m], "%d", buf[idx]);
	  // std::cout << buf[i*nvalues+m] << " " << values[m] << std::endl;
	}
	if (latticeflag) applattice->add_values(loc->second,values);
	else appoff->add_values(loc->second,values);
      }
    }

    nread += nchunk;
  }
  for (i=0; i < nvalues; i++) delete [] values[i];
  delete [] values;

  bigint nbig = nglobal;
  if (me == 0) {
    if (screen)
      fprintf(screen,"  " BIGINT_FORMAT " values\n", nbig*nvalues);
    if (logfile)
      fprintf(logfile,"  " BIGINT_FORMAT " values\n", nbig*nvalues);
  }
  delete [] data;
  return;
}

void ReadDream3d::read_average_quaternions() {
  int num_grains = 0;
  int num_values = 0;
#ifdef SPPARKS_HDF5
  if (me == 0) {
    hsize_t dims[2];
    h5_status = H5LTget_dataset_info(file_id, quaternions_path.c_str(),
				     dims,NULL,NULL);
    num_grains = dims[0];
    num_values = dims[1];
  }
#endif
  
  MPI_Bcast(&num_grains,1,MPI_INT, 0, world);
  int data_size = 4 * num_grains;
  float* data = new float[data_size];
  
#ifdef SPPARKS_HDF5
  if (me == 0) {
    h5_status = H5LTread_dataset_float(file_id, quaternions_path.c_str(), data);
  }
#endif

  /* root proc broadcasts orientation data */
  MPI_Bcast(data, data_size, MPI_FLOAT, 0, world);

  // update nspins!
  app_potts_ori->update_nspins(num_grains-1);

  app_potts_ori->gb_props->copy_quaternion_data(data,num_grains);
}
