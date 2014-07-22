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
#include "domain.h"
#include "create_sites.h"
#include "error.h"
#include "memory.h"

#include "input.h"

#include <map>
#include <string>

using namespace SPPARKS_NS;

#define MAXLINE 256
#define CHUNK 1024
#define DELTA 4
#define EPSILON 1.0e-6

#define NSECTIONS 3       // change when add to header::section_keywords

/* ---------------------------------------------------------------------- */

ReadDream3d::ReadDream3d(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  buffer = new int[CHUNK*MAXLINE];
  narg = maxarg = 0;
  arg = NULL;
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

  if (narg != 1) error->all(FLERR,"Illegal read_dream3d command");

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

  // proc 0 opens the dream3d file
  if (me == 0) {
    if (screen) fprintf(screen,"Reading dream3d file ...\n");
    // make sure the file exists first
    file_id = H5Fopen(arg[0], H5F_ACC_RDONLY, H5P_DEFAULT);
  }
   
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
  extract_grain_ids();

  // Read orientation data

  // close file
  std::cout << "Finished reading dream3d file" << std::endl;
  if (me == 0) {
    h5_status = H5Fclose(file_id);
  }

}

/* ----------------------------------------------------------------------
   read free-format header of data file
   if flag = 0, only called by proc 0
   if flag = 1, called by all procs so bcast lines as read them
   1st line and blank lines are skipped
   non-blank lines are checked for header keywords and leading value is read
   header ends with EOF or non-blank line containing no header keyword
     if EOF, line is set to blank line
     else line has first keyword line for rest of file
------------------------------------------------------------------------- */

void ReadDream3d::header(){}

/* ----------------------------------------------------------------------
   read all sites
------------------------------------------------------------------------- */

void ReadDream3d::sites(){}

/* ----------------------------------------------------------------------
   read all per-site values
   to find atoms, must build atom map if not a molecular system 
------------------------------------------------------------------------- */

void ReadDream3d::values(){}

/* ----------------------------------------------------------------------
   proc 0 opens data file
   test if gzipped
------------------------------------------------------------------------- */

void ReadDream3d::open(char *file){}
void ReadDream3d::parse_keyword(int) {}
void ReadDream3d::parse_coeffs(int, char *) {}
int ReadDream3d::count_words(char *) {}

void ReadDream3d::get_dimensions() {
  int dims_buf[3] = {0, 0, 0};
  h5_status = H5LTread_dataset_int(file_id,"/VoxelDataContainer/DIMENSIONS", dims_buf);
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

void ReadDream3d::extract_grain_ids() {

  int* data = NULL;
  if (me == 0) {
    data = new int[app->nglobal];
    // assuming SPPARKS and DREAM3D both use row-major indexing, GrainIds are indexed by global_id
    h5_status = H5LTread_dataset_int(file_id,"/VoxelDataContainer/CELL_DATA/GrainIds",data);
  }

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
	  values[m] = new char[256];
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
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " values\n",
			nbig*nvalues);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " values\n",
			 nbig*nvalues);
  }
  delete [] data;
  return;
}
