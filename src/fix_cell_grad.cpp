/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "string.h"
#include "domain.h"
#include "fix_cell_grad.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI,INTERIOR};         // same as Domain
enum{PKEEP,PINSERT,PDONE,PDISCARD,PENTRY,PEXIT,PSURF};   // several files
enum{INT,DOUBLE};                      // several files

// cell face quantities
// mass density, x-velocity, y-velocity, z-velocity

enum{RHO,U,V,W,UMAG,VMAG,WMAG,FACELASTSIZE};

// for outputs

enum{DX,DY,DZ,SUM};

#define MAXOUTPUTS 12

/* ---------------------------------------------------------------------- */

FixCellGrad::FixCellGrad(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal fix cell/grad command");

  // how frequently to update cell fluxes
  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix cell/grad command");
  T_interval = nevery * update->dt;

  // time average the quantities?

  if (strcmp(arg[3],"ave") != 0) 
    error->all(FLERR,"Expected 'ave' keyword\n");
  if (strcmp(arg[4],"running") == 0) aveflag = 1;
  else if (strcmp(arg[4],"one") == 0) aveflag = 0;
  else error->all(FLERR,"Unexpected input in fix cell/grad command");

  // number of quantities

  nvalues = atoi(arg[5]);

  // create direction and quantity

  memory->create(direction,nvalues,"fix/cell/grad:direction");
  memory->create(quantity,nvalues,"fix/cell/grad:quantity");

  // only need to track lower value (will handle both)
  int iarg = 6;
  size_per_grid_cols = 0; // number of face values
  while (iarg < narg) {
    if (strcmp(arg[iarg],"u_x") == 0) {
      direction[iarg-6] = XLO;
      quantity[iarg-6] = U;
    } else if (strcmp(arg[iarg],"u_y") == 0) {
      direction[iarg-6] = YLO;
      quantity[iarg-6] = U;
    } else if (strcmp(arg[iarg],"u_z") == 0) {
      direction[iarg-6] = ZLO;
      quantity[iarg-6] = U;
    } else if (strcmp(arg[iarg],"v_x") == 0) {
      direction[iarg-6] = XLO;
      quantity[iarg-6] = V;
    } else if (strcmp(arg[iarg],"v_y") == 0) {
      direction[iarg-6] = YLO;
      quantity[iarg-6] = V;
    } else if (strcmp(arg[iarg],"v_z") == 0) {
      direction[iarg-6] = ZLO;
      quantity[iarg-6] = V;
    } else if (strcmp(arg[iarg],"w_x") == 0) {
      direction[iarg-6] = XLO;
      quantity[iarg-6] = U;
    } else if (strcmp(arg[iarg],"w_y") == 0) {
      direction[iarg-6] = YLO;
      quantity[iarg-6] = V;
    } else if (strcmp(arg[iarg],"w_z") == 0) {
      direction[iarg-6] = ZLO;
      quantity[iarg-6] = W;
    } else error->all(FLERR,"Invalid fix cell/grad command");

    size_per_grid_cols++;
    iarg++;
  }

  if (nvalues != size_per_grid_cols) error->all(FLERR,"Should be same");

  // for outputting per-grid quantities
  per_grid_flag = 1;
  per_grid_freq = 1;
  nglocal = 0;
  array_grid = NULL;

  // count number of bulk and face values per-grid cell

  if (size_per_grid_cols == 0)
    error->one(FLERR,"No face values specified in fix cell grad");

  // store quantities in custom grid arrays
  // check if custom per-grid attribute exists

  cellbulkindex = grid->find_custom((char *) "cellbulk");
  cellfaceindex = grid->find_custom((char *) "cellface");

  if (cellfaceindex >= 0 || cellbulkindex >= 0)
    error->all(FLERR,"Fix cell gradient already exists");

  // cell bulk tracks density and number of bulk velocities (one extra)
  cellbulkindex = grid->add_custom((char *) "cellbulk", DOUBLE, nvalues+1);
  // each gradients needs 2 values for numer and denom
  // shouldn't need to do any syncronization between procs since this
  // .. quantity is measured during move where particles are migrated
  cellfaceindex = grid->add_custom((char *) "cellface", DOUBLE, 2*nvalues);
}

/* ---------------------------------------------------------------------- */

FixCellGrad::~FixCellGrad()
{
  if (copymode) return;

  memory->destroy(array_grid);
  memory->destroy(direction);
  memory->destroy(quantity);

  delete [] cellids;
  delete [] faceids;

  grid->remove_custom(cellbulkindex);
  grid->remove_custom(cellfaceindex);
}

/* ---------------------------------------------------------------------- */

int FixCellGrad::setmask()
{
  // averaging needs to be done during update::move()
  // gradient computations done at the end
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCellGrad::init()
{
  dim = domain->dimension;
  reallocate();
  nsample = 0;
}

/* ---------------------------------------------------------------------- */

void FixCellGrad::end_of_step()
{
  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;

  reallocate();

  int ivalue;
  double rho, vcell, igrad, volume;

  // grab the per-grid custom arrays

  double **face = grid->edarray[grid->ewhich[cellfaceindex]];
  double **cell = grid->edarray[grid->ewhich[cellbulkindex]];

  // cell geometry

  double *lo;
  double *hi;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!aveflag) {
      for (int i = 0; i < size_per_grid_cols; i++)
        array_grid[icell][i] = 0.0;
      nsample = 0;
    }

    lo = cells[icell].lo;
    hi = cells[icell].hi;
    volume = (hi[0]-lo[0])*(hi[1]-lo[1]);
    if (dim == 3) volume *= (hi[2]-lo[2]);

    rho = cell[icell][0]/(volume*T_interval);
    //printf("rho: %4.3e\n", rho);
    ivalue = 0;
    for (int i = 0; i < size_per_grid_cols; i++) {
      vcell = cell[icell][i+1]/(volume*T_interval);
      vcell /= rho;
      //printf("vcell: %4.3e\n", vcell);
      igrad = face[icell][ivalue] - face[icell][ivalue+1]*vcell;
      igrad /= (rho*volume*T_interval);
      array_grid[icell][i] = (array_grid[icell][i]*nsample+igrad)/(nsample+1);
      ivalue += 2;
    }

    // reset custom arrays for measuring over next interval

    ivalue = 0;
    cell[icell][0] = 0.0;
    for (int i = 0; i < size_per_grid_cols; i++) {
      cell[icell][i+1] = 0.0;
      face[icell][ivalue] = 0.0;
      face[icell][ivalue+1] = 0.0;
      ivalue += 2;
    }

  } // end cells

  //error->one(FLERR,"ck end of step");
}

/* ---------------------------------------------------------------------- */

void FixCellGrad::during_move(Particle::OnePart *ipart, int type, int icell, int outface, double dtremain)
{
  if (type == 0) face_flux_premove(ipart,icell);
  else if (type == 1) update_cell_bulk(ipart, icell, dtremain);
  else if (type == 2) face_flux_postmove(ipart, outface, icell);
  else error->all(FLERR,"Type not recognized");
}


/* ----------------------------------------------------------------------
   records crossing events in past
------------------------------------------------------------------------- */

void FixCellGrad::face_flux_premove(Particle::OnePart *ipart, int icell)
{
  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;

  // get particle attributes/properties

  int isp = ipart->ispecies;
  double pmass = particle->species[isp].mass;
  int pflag = ipart->flag;
  double dtremain = ipart->dtremain;

  // cell geometry

  double *lo = cells[icell].lo;
  double *hi = cells[icell].hi;

  // particle position/velocity

  double x[3],v[3];
  memcpy(x,ipart->x,3*sizeof(double));
  memcpy(v,ipart->v,3*sizeof(double));

  // check history if:
  // 1) particle has crossed from another cell (dtremain < update->dt)
  // ... or previously hit a boundary
  // 2) particle emitted from face (pflag == PINSERT)
  // 3) particle emitted from surf (pflag >= PSURF) (also includes boundaries)
  // second part of if statement checks for dtremain == update->dt so

  // note: pflag = 0 == PKEEP

  // note: if particle 'i' hits a surface, particle 'i' does not get 
  // ... PSURF flag. only particle 'j' (for reactive surfaces) get this flag

  int prev_cross_flag = 0;
  if (pflag == PINSERT || pflag >= PSURF) prev_cross_flag = 1;
  else if (dtremain > 0 && dtremain < update->dt) prev_cross_flag = 1;
  else {
    double xold[3];
    xold[0] = x[0] - v[0]*update->dt;
    xold[1] = x[1] - v[1]*update->dt;
    xold[2] = x[2] - v[2]*update->dt;

    if (xold[0] < lo[0] || xold[0] > hi[0]) prev_cross_flag = 1;
    if (xold[1] < lo[1] || xold[1] > hi[1]) prev_cross_flag = 1;
    if (dim == 3)
      if (xold[2] < lo[2] || xold[2] > hi[2]) prev_cross_flag = 1;
  }

  // Should never call
  if (x[0] < lo[0] || x[0] > hi[0]) error->one(FLERR,"x[0] outside cell");
  if (x[1] < lo[1] || x[1] > hi[1]) error->one(FLERR,"x[1] outside cell");
  if (dim == 3)
    if (x[2] < lo[2] || x[2] > hi[2]) error->one(FLERR,"x[2] outside cell");

  // if dtremain == 0, then it did NOT come from another cell or
  // ... was emitted from a face/surface
  // can ignore since no history to record
  if (!prev_cross_flag) return;

  // grab the per-grid custom arrays

  double **face = grid->edarray[grid->ewhich[cellfaceindex]];

  // find which face the particle is closest to
  // set guess as max cell dimension

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double mindist = MAX(boxhi[0]-boxlo[0],boxhi[1]-boxlo[1]);
  if (dim == 3) mindist = MAX(boxhi[2]-boxlo[2],mindist);

  // CAUTION: This won't account for particles equidistant from two
  //          cell edges or on cell corners (small/negligible error)

  int inface = -1;
  if (fabs(x[0]-lo[0]) < mindist && v[0] > 0.0) {
    inface = XLO;
    mindist = fabs(x[0]-lo[0]);
  }
  if (fabs(x[1]-lo[1]) < mindist && v[1] > 0.0) {
    inface = YLO;
    mindist = fabs(x[1]-lo[1]);
  }

  if (fabs(x[0]-hi[0]) < mindist && v[0] < 0.0) {
    inface = XHI;
    mindist = fabs(x[0]-hi[0]);
  }
  if (fabs(x[1]-hi[1]) < mindist && v[1] < 0.0) {
    inface = YHI;
    mindist = fabs(x[1]-hi[1]);
  }

  if (dim == 3) {
    if (fabs(x[2]-lo[2]) < mindist && v[2] > 0.0) {
      inface = ZLO;
      mindist = fabs(x[2]-lo[2]);
    }
    if (fabs(x[2]-hi[2]) < mindist && v[2] < 0.0) {
      inface = ZHI;
      mindist = fabs(x[2]-hi[2]);
    }
  }

  // error if no edge found (should never be called)

  if (inface < 0) error->one(FLERR,"Cannot find cell edge");

  // Heaviside function

  double theta = 1.0;
  if (inface == XLO || inface == YLO || inface == ZLO) theta = -1.0;

  // update cell face quantities of interest
  // note: face has Ngrid elements where each element contains Nface = 6 in 3D
  // ... (4 in 2D) cell faces. For each grid cell, the values are flattend
  // ... into a Nface x Nqoi where Nqoi is the number of quantities of interest

  int ivalue = 0;
  double denom, numer;
  for (int ival = 0; ival < nvalues; ival++) {
    if (quantity[ival] == U) numer = pmass*v[0];
    else if (quantity[ival] == V) numer = pmass*v[1];
    else if (quantity[ival] == W) numer = pmass*v[2];
    else error->all(FLERR,"Quantitiy not valid");

    if (direction[ival] == XLO || direction[ival]+1 == XHI)
      denom = fabs(v[0]);
    else if (direction[ival] == YLO || direction[ival]+1 == YHI)
      denom = fabs(v[1]);
    else if (direction[ival] == ZLO || direction[ival]+1 == ZHI)
      denom = fabs(v[2]);
    else error->all(FLERR,"Direction not valid");

    face[icell][ivalue] += (numer/denom)*theta;
    face[icell][ivalue+1] += (1.0/denom)*theta;
    ivalue += 2;
  }

  /*printf("x: %4.3e, %4.3e, %4.3e\n", x[0], x[1], x[2]);
  printf("v: %4.3e, %4.3e, %4.3e\n", v[0], v[1], v[2]);
  printf("xold: %4.3e, %4.3e, %4.3e\n",
    x[0]-v[0]*update->dt,
    x[1]-v[1]*update->dt,
    x[2]-v[2]*update->dt);
  printf("lo: %4.3e, %4.3e, %4.3e\n", lo[0], lo[1], lo[2]);
  printf("hi: %4.3e, %4.3e, %4.3e\n", hi[0], hi[1], hi[2]);
  printf("inface: %i\n", inface);
  printf("face_flux_premove\n");
  error->one(FLERR,"ck");*/

}

/* ----------------------------------------------------------------------
   update cell bulk based on relative time in cell
------------------------------------------------------------------------- */

void FixCellGrad::update_cell_bulk(Particle::OnePart *ipart, int icell, double dtremain)
{
  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;

  // get particle attributes/properties

  int isp = ipart->ispecies;
  double pmass = particle->species[isp].mass;

  // grab the per-grid custom arrays

  double **cell = grid->edarray[grid->ewhich[cellbulkindex]];

  // particle position/velocity

  double v[3];
  memcpy(v,ipart->v,3*sizeof(double));

  double dt_cell = dtremain/update->dt;
  cell[icell][0] += pmass*dt_cell;
  int ivalue = 1;
  for (int ival = 0; ival < nvalues; ival++) {
    if (quantity[ival] == U) cell[icell][ival+1] += pmass*v[0]*dt_cell;
    else if (quantity[ival] == V) cell[icell][ival+1] += pmass*v[1]*dt_cell;
    else if (quantity[ival] == W) cell[icell][ival+1] += pmass*v[2]*dt_cell;
    else error->all(FLERR,"Quantitiy not valid");

    /*if (icell == 0) {
      if (quantity[ival] == U)
        printf("ival: %i - U : %4.3e -> %4.3e\n", ival, pmass*v[0]*dt_cell, cell[icell][ival+1]);
      else if (quantity[ival] == V)
        printf("ival: %i - V : %4.3e -> %4.3e\n", ival, pmass*v[1]*dt_cell, cell[icell][ival+1]);
      else if (quantity[ival] == W)
        printf("ival: %i - W : %4.3e -> %4.3e\n", ival, pmass*v[2]*dt_cell, cell[icell][ival+1]);
    }*/

  } // END for cell values
}


/* ----------------------------------------------------------------------
   records current crossing events
------------------------------------------------------------------------- */

void FixCellGrad::face_flux_postmove(Particle::OnePart *ipart, int outface, int icell)
{
  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;

  // get particle mass

  int isp = ipart->ispecies;
  double pmass = particle->species[isp].mass;
  int pflag = ipart->flag;
  double dtremain = ipart->dtremain;

  // grab the per-grid custom arrays

  double **face = grid->edarray[grid->ewhich[cellfaceindex]];

  // particle props

  double v[3];
  memcpy(v,ipart->v,3*sizeof(double));

  // cell props

  double *lo = cells[icell].lo;
  double *hi = cells[icell].hi;

  if (outface == INTERIOR) error->all(FLERR,"Should never record outface");

  //printf("outface: %i\n", outface);

  // heaviside function
  double theta = 1.0;
  if (outface == XLO || outface == YLO || outface == ZLO) theta = -1.0;

  //printf("theta: %4.3e\n", theta);

  // update cell face quants
  // faces follow same enum order (range from 0 -> 2*DIM-1)
  // xlo -> 0, xhi -> 1, ylo -> 2, ...
  int ivalue = 0;
  double denom, numer;

  for (int ival = 0; ival < nvalues; ival++) {

    if (quantity[ival] == U) numer = pmass*v[0];
    else if (quantity[ival] == V) numer = pmass*v[1];
    else if (quantity[ival] == W) numer = pmass*v[2];
    else error->all(FLERR,"Quantitiy not valid");

    if (direction[ival] == XLO || direction[ival]+1 == XHI)
      denom = fabs(v[0]);
    else if (direction[ival] == YLO || direction[ival]+1 == YHI)
      denom = fabs(v[1]);
    else if (direction[ival] == ZLO || direction[ival]+1 == ZHI)
      denom = fabs(v[2]);
    else error->all(FLERR,"Direction not valid");

    face[icell][ivalue] += (numer/denom)*theta;
    face[icell][ivalue+1] += (1.0/denom)*theta;
    ivalue += 2;

    //if (quantity[ival] == U)
    //  printf("ival: %i - U : %4.3e\n", ival, v[0]);
    //else if (quantity[ival] == V)
    //  printf("ival: %i - V : %4.3e\n", ival, v[1]);
    //else if (quantity[ival] == W)
    //  printf("ival: %i - W : %4.3e\n", ival, v[2]);

    //if (direction[ival] == XLO || direction[ival]+1 == XHI)
    //  printf("ival: %i - X : %4.3e\n", ival, v[0]);
    //else if (direction[ival] == YLO || direction[ival]+1 == YHI)
    //  printf("ival: %i - Y : %4.3e\n", ival, v[1]);
    //else if (direction[ival] == ZLO || direction[ival]+1 == ZHI)
    //  printf("ival: %i - Z : %4.3e\n", ival, v[2]);
  }

  //printf("face_flux_postmove\n");
  //error->one(FLERR,"ck");
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void FixCellGrad::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(array_grid);
  nglocal = grid->nlocal;
  memory->create(array_grid,nglocal,size_per_grid_cols,"fix/cell/grad:array_grid");

  // initialize values
  for (int i = 0; i < nglocal; i++)
    for (int j = 0; j < size_per_grid_cols; j++)
      array_grid[i][j] = 0.0;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

//double FixCellGrad::memory_usage()
//{
//  double bytes = 0.0;
  //bytes += maxgrid*3 * sizeof(double);
//  return bytes;
//}
