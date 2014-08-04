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

#ifdef SPPARKS_HDF5
#include "H5Cpp.h"
#include "hdf5_hl.h"
#endif

using namespace SPPARKS_NS;
using namespace MathConst;

enum{DREAM3D,PLAIN_H5};
enum{PARALLEL,SERIAL};
enum{NUMERIC,IATTRIBUTE,DATTRIBUTE};
enum{STATIC,DYNAMIC};
enum{NO,YES};
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

  // site color,diameter settings

  // if (size_one != 2) error->all(FLERR,"Illegal dump dream3d command");

  if (vtype[0] == INT) scolor = IATTRIBUTE;
  else scolor = DATTRIBUTE;

  if (vtype[1] == INT) sdiam = IATTRIBUTE;
  else sdiam = DATTRIBUTE;

  // create Image class

  image = new Image(spk);

  // set defaults for optional args
  write_style = SERIAL;
  shape = SPHERE;
  boundflag = NO;
  crange = drange = NO;
  thetastr = phistr = NULL;
  cflag = STATIC;
  cx = cy = cz = 0.5;
  cxstr = cystr = czstr = NULL;

  if (domain->dimension == 3) {
    image->up[0] = 0.0; image->up[1] = 0.0; image->up[2] = 1.0;
  } else if (domain->dimension == 2) {
    image->up[0] = 0.0; image->up[1] = 1.0; image->up[2] = 0.0;
  } else if (domain->dimension == 1) {
    image->up[0] = 1.0; image->up[1] = 0.0; image->up[2] = 0.0;
  }

  upxstr = upystr = upzstr = NULL;
  zoomstr = NULL;
  perspstr = NULL;
  boxflag = YES;
  boxdiam = 0.02;
  axesflag = NO;

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


  // allocate image buffer now that image size is known
  image->buffers();

}

/* ---------------------------------------------------------------------- */

DumpDream3d::~DumpDream3d()
{
  delete image;

  memory->destroy(diamattribute);
  if (colorattribute) {
    for (int i = clo; i <= chi; i++)
      if (color_memflag[i-clo]) delete [] colorattribute[i-clo];
    memory->sfree(colorattribute);
  }
  memory->destroy(color_memflag);
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
  create_dream3d_file();
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
  if (scolor == DATTRIBUTE) image->color_minmax(nchoose,buf,size_one);

  // create image on each proc, then merge them

  image->clear();
  create_image();
  image->merge();

  // write image file

  if (me == 0) {
    if (filetype == JPG) image->write_JPG(fp);
    else image->write_PPM(fp);
    fclose(fp);
  }
}

void DumpDream3d::create_dream3d_file() {
  std::cout << "test\n";
  std::cout << filename;
}

/* ----------------------------------------------------------------------
   reset view parameters
   called once from constructor if view is STATIC
   called every snapshot from write() if view is DYNAMIC
------------------------------------------------------------------------- */

void DumpDream3d::box_center(){}

/* ----------------------------------------------------------------------
   create image for atoms on this proc
   every pixel has depth 
------------------------------------------------------------------------- */

void DumpDream3d::create_image()
{
  int i,j,m,n,ivalue;
  double diameter;
  double *color;

  // render my sites

  double **xyz = app->xyz;

  m = 0;
  for (i = 0; i < nchoose; i++) {
    j = clist[i];
    
    if (scolor == IATTRIBUTE) {
      ivalue = static_cast<int> (buf[m]);
      ivalue = MAX(ivalue,clo);
      ivalue = MIN(ivalue,chi);
      color = colorattribute[ivalue-clo];
    } else if (scolor == DATTRIBUTE) {
      color = image->value2color(buf[m]);
    }

    if (sdiam == NUMERIC) {
      diameter = sdiamvalue;
    } else if (sdiam == IATTRIBUTE) {
      ivalue = static_cast<int> (buf[m+1]);
      ivalue = MAX(ivalue,dlo);
      ivalue = MIN(ivalue,dhi);
      diameter = diamattribute[ivalue-dlo];
    } else if (sdiam == DATTRIBUTE) {
      diameter = buf[m+1];
    }

    if (shape == SPHERE) image->draw_sphere(xyz[j],color,diameter);
    else image->draw_cube(xyz[j],color,diameter);

    m += size_one;
  }

  // render my boundaries bewteen adjacent sites
  // loop over all chosen sites and all their neighbors
  // neighbor does not have to be chosen site
  // if 2 sites do not share adjacent face, do not draw boundary
  // if 2 sites have same value, do not draw boundary

  if (boundflag == YES) {
    int k,flag;
    double c1[3],c2[3],c3[4],c4[4];
    int dimension = domain->dimension;
    double dx = domain->lattice->xlattice;
    double dy = domain->lattice->ylattice;
    double dz = domain->lattice->zlattice;
    int *numneigh = applattice->numneigh;
    int **neighbor = applattice->neighbor;
    tagint *id = app->id;
    int **iarray = app->iarray;
    double **darray = app->darray;

    for (int ii = 0; ii < nchoose; ii++) {
      i = clist[ii];
      for (int jj = 0; jj < numneigh[i]; jj++) {
	j = neighbor[i][jj];
	if (boundvalue == ID) {
	  if (id[i] == id[j]) continue;
	} else if (boundvalue == IARRAY) {
	  if (iarray[boundindex][i] == iarray[boundindex][j]) continue;
	} else if (boundvalue == DARRAY) {
	  if (darray[boundindex][i] == darray[boundindex][j]) continue;
	} else if (boundvalue == X) {
	  if (xyz[i][0] == xyz[j][0]) continue;
	} else if (boundvalue == Y) {
	  if (xyz[i][1] == xyz[j][1]) continue;
	} else if (boundvalue == Z) {
	  if (xyz[i][2] == xyz[j][2]) continue;
	}

	flag = 0;
	if (xyz[i][0] != xyz[j][0]) flag++;
	if (xyz[i][1] != xyz[j][1]) flag++;
	if (xyz[i][2] != xyz[j][2]) flag++;
	if (flag >= 2) continue;

	if (xyz[i][0] != xyz[j][0]) {
	  if (fabs(xyz[i][0]-xyz[j][0]) < 0.5*domain->xprd)
	    c1[0] = c2[0] = c3[0] = c4[0] = 0.5*(xyz[i][0]+xyz[j][0]);
	  else 
	    c1[0] = c2[0] = c3[0] = c4[0] = 
	      0.5*(xyz[i][0]+xyz[j][0]-domain->xprd);
	  c1[1] = xyz[i][1] - 0.5*dy; c1[2] = xyz[i][2] - 0.5*dz;
	  c2[1] = xyz[i][1] - 0.5*dy; c2[2] = xyz[i][2] + 0.5*dz;
	  c3[1] = xyz[i][1] + 0.5*dy; c3[2] = xyz[i][2] + 0.5*dz;
	  c4[1] = xyz[i][1] + 0.5*dy; c4[2] = xyz[i][2] - 0.5*dz;
	} else if (xyz[i][1] != xyz[j][1]) {
	  if (fabs(xyz[i][1]-xyz[j][1]) < 0.5*domain->yprd)
	    c1[1] = c2[1] = c3[1] = c4[1] = 0.5*(xyz[i][1]+xyz[j][1]);
	  else 
	    c1[1] = c2[1] = c3[1] = c4[1] = 
	      0.5*(xyz[i][1]+xyz[j][1]-domain->yprd);
	  c1[0] = xyz[i][0] - 0.5*dx; c1[2] = xyz[i][2] - 0.5*dz;
	  c2[0] = xyz[i][0] - 0.5*dx; c2[2] = xyz[i][2] + 0.5*dz;
	  c3[0] = xyz[i][0] + 0.5*dx; c3[2] = xyz[i][2] + 0.5*dz;
	  c4[0] = xyz[i][0] + 0.5*dx; c4[2] = xyz[i][2] - 0.5*dz;
	} else {
	  if (fabs(xyz[i][2]-xyz[j][2]) < 0.5*domain->zprd)
	    c1[2] = c2[2] = c3[2] = c4[2] = 0.5*(xyz[i][2]+xyz[j][2]);
	  else
	    c1[2] = c2[2] = c3[2] = c4[2] = 
	      0.5*(xyz[i][2]+xyz[j][2]-domain->zprd);
	  c1[0] = xyz[i][0] - 0.5*dx; c1[1] = xyz[i][1] - 0.5*dy;
	  c2[0] = xyz[i][0] - 0.5*dx; c2[1] = xyz[i][1] + 0.5*dy;
	  c3[0] = xyz[i][0] + 0.5*dx; c3[1] = xyz[i][1] + 0.5*dy;
	  c4[0] = xyz[i][0] + 0.5*dx; c4[1] = xyz[i][1] - 0.5*dy;
	}
	
	image->draw_cylinder(c1,c2,boundcolor,bounddiam,3);
	image->draw_cylinder(c2,c3,boundcolor,bounddiam,3);
	image->draw_cylinder(c3,c4,boundcolor,bounddiam,3);
	image->draw_cylinder(c4,c1,boundcolor,bounddiam,3);
      }
    }
  }

  // render outline of simulation box

  if (boxflag) {
    double diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= boxdiam;

    double (*corners)[3];
    double corner[8][3];
    corner[0][0] = boxxlo; corner[0][1] = boxylo; corner[0][2] = boxzlo;
    corner[1][0] = boxxhi; corner[1][1] = boxylo; corner[1][2] = boxzlo;
    corner[2][0] = boxxlo; corner[2][1] = boxyhi; corner[2][2] = boxzlo;
    corner[3][0] = boxxhi; corner[3][1] = boxyhi; corner[3][2] = boxzlo;
    corner[4][0] = boxxlo; corner[4][1] = boxylo; corner[4][2] = boxzhi;
    corner[5][0] = boxxhi; corner[5][1] = boxylo; corner[5][2] = boxzhi;
    corner[6][0] = boxxlo; corner[6][1] = boxyhi; corner[6][2] = boxzhi;
    corner[7][0] = boxxhi; corner[7][1] = boxyhi; corner[7][2] = boxzhi;
    corners = corner;

    image->draw_box(corners,diameter);
  }

  // render XYZ axes in red/green/blue
  // offset by 10% of box size and scale by axeslen

  if (axesflag) {
    double diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= axesdiam;

    double (*corners)[3];
    double corner[8][3];
    corner[0][0] = boxxlo; corner[0][1] = boxylo; corner[0][2] = boxzlo;
    corner[1][0] = boxxhi; corner[1][1] = boxylo; corner[1][2] = boxzlo;
    corner[2][0] = boxxlo; corner[2][1] = boxyhi; corner[2][2] = boxzlo;
    corner[3][0] = boxxlo; corner[3][1] = boxylo; corner[3][2] = boxzhi;
    corners = corner;

    double offset = MAX(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) offset = MAX(offset,boxzhi-boxzlo);
    offset *= 0.1;
    corners[0][0] -= offset; corners[0][1] -= offset; corners[0][2] -= offset;
    corners[1][0] -= offset; corners[1][1] -= offset; corners[1][2] -= offset;
    corners[2][0] -= offset; corners[2][1] -= offset; corners[2][2] -= offset;
    corners[3][0] -= offset; corners[3][1] -= offset; corners[3][2] -= offset;

    corners[1][0] = corners[0][0] + axeslen*(corners[1][0]-corners[0][0]);
    corners[1][1] = corners[0][1] + axeslen*(corners[1][1]-corners[0][1]);
    corners[1][2] = corners[0][2] + axeslen*(corners[1][2]-corners[0][2]);
    corners[2][0] = corners[0][0] + axeslen*(corners[2][0]-corners[0][0]);
    corners[2][1] = corners[0][1] + axeslen*(corners[2][1]-corners[0][1]);
    corners[2][2] = corners[0][2] + axeslen*(corners[2][2]-corners[0][2]);
    corners[3][0] = corners[0][0] + axeslen*(corners[3][0]-corners[0][0]);
    corners[3][1] = corners[0][1] + axeslen*(corners[3][1]-corners[0][1]);
    corners[3][2] = corners[0][2] + axeslen*(corners[3][2]-corners[0][2]);

    image->draw_axes(corners,diameter);
  }
}

/* ---------------------------------------------------------------------- */

int DumpDream3d::modify_param(int narg, char **arg)
{
  int n = DumpText::modify_param(narg,arg);
  if (n) return n;
  
  if (strcmp(arg[0],"backcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    double *color = image->color2rgb(arg[1]);
    if (color == NULL) error->all(FLERR,"Invalid color in dump_modify command");
    image->background[0] = static_cast<int> (color[0]*255.0);
    image->background[1] = static_cast<int> (color[1]*255.0);
    image->background[2] = static_cast<int> (color[2]*255.0);
    return 2;
  }

  if (strcmp(arg[0],"boundcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    boundcolor = image->color2rgb(arg[1]);
    if (boundcolor == NULL) error->all(FLERR,
				       "Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"boxcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    image->boxcolor = image->color2rgb(arg[1]);
    if (image->boxcolor == NULL) 
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"color") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->addcolor(arg[1],atof(arg[2]),atof(arg[3]),atof(arg[4]));
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return 5;
  }

  if (strcmp(arg[0],"scolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    if (scolor != IATTRIBUTE)
      error->all(FLERR,"Dump_modify scolor requires integer attribute "
		 "for dump image color");

    int nlo,nhi;
    bounds(arg[1],clo,chi,nlo,nhi);

    // color arg = "random"
    // assign random RGB values to each attribute

    if (strcmp(arg[2],"random") == 0) {
      RandomPark *randomcolor = new RandomPark(ranmaster->uniform()); 
      for (int i = nlo; i <= nhi; i++) {
	double *rgb;
	if (color_memflag[i-clo] == 0) rgb = new double[3];
	else rgb = colorattribute[i-clo];
	rgb[0] = randomcolor->uniform();
	rgb[1] = randomcolor->uniform();
	rgb[2] = randomcolor->uniform();
	colorattribute[i-clo] = rgb;
	color_memflag[i-clo] = 1;
      }
      delete randomcolor;

    // color arg is a color name
    // ptrs = list of ncount colornames separated by '/'
    // assign each of ncount colors in round-robin fashion to attributes

    } else {
      int ncount = 1;
      char *nextptr;
      char *ptr = arg[2];
      while (nextptr = strchr(ptr,'/')) {
	ptr = nextptr + 1;
	ncount++;
      }
      char **ptrs = new char*[ncount+1];
      ncount = 0;
      ptrs[ncount++] = strtok(arg[2],"/");
      while (ptrs[ncount++] = strtok(NULL,"/"));
      ncount--;
      
      int m = 0;
      for (int i = nlo; i <= nhi; i++) {
	if (color_memflag[i-clo] == 1) delete [] colorattribute[i-clo];
	colorattribute[i-clo] = image->color2rgb(ptrs[m%ncount]);
	color_memflag[i-clo] = 0;
	if (colorattribute[i-clo] == NULL)
	  error->all(FLERR,"Invalid color in dump_modify command");
	m++;
      }

      delete [] ptrs;
    }

    return 3;
  }

  if (strcmp(arg[0],"sdiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    if (sdiam != IATTRIBUTE)
      error->all(FLERR,"Dump_modify sdiam requires integer attribute "
		 "for dump image diameter");

    int nlo,nhi;
    bounds(arg[1],dlo,dhi,nlo,nhi);

    double diam = atof(arg[2]);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) diamattribute[i-dlo] = diam;
    return 3;
  }

  if (strcmp(arg[0],"smap") == 0) {
    if (narg < 6) error->all(FLERR,"Illegal dump_modify command");
    if (strlen(arg[3]) != 2) error->all(FLERR,"Illegal dump_modify command");
    int factor = 2;
    if (arg[3][0] == 's') factor = 1;
    int nentry = atoi(arg[5]);
    if (nentry < 1) error->all(FLERR,"Illegal dump_modify command");
    int n = 6 + factor*nentry;
    if (narg < n) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->colormap(n-1,&arg[1]);
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return n;
  }
  
  return 0;
}

void DumpDream3d::bounds(char *str, int lo, int hi, int &nlo, int &nhi){}
