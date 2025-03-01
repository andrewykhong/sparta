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

  // nface = number of faces in the grid cell

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
  error->one(FLERR,"ck end of step");

  Particle::OnePart *particles = particle->particles;
  Grid::ChildCell *cells = grid->cells;

  reallocate();

  int ivalue;
  double rho, vbulk, igrad, volume;

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

    rho = cell[icell][0];
    ivalue = 0;
    for (int i = 0; i < size_per_grid_cols; i++) {
      vbulk = cell[icell][i+1];
      igrad = face[icell][ivalue] - face[icell][ivalue+1]*vbulk;
      igrad *= rho/volume/T_interval;
      array_grid[icell][i] = (array_grid[icell][i]*nsample+igrad)/(nsample+1);
      ivalue += 2;
    }
  } // end cells

  return;
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

  // grab the per-grid custom arrays

  double **face = grid->edarray[grid->ewhich[cellfaceindex]];

  // particle position/velocity

  double x[3],v[3];
  memcpy(v,ipart->v,3*sizeof(double));

  // cell geometry

  double *lo = cells[icell].lo;
  double *hi = cells[icell].hi;

  // check history if:
  // 1) particle has crossed from another cell (dtremain < update->dt)
  // ... or previously hit a boundary
  // 2) particle emitted from face (pflag == PINSERT)
  // 3) particle emitted from surf (pflag >= PSURF) (also includes boundaries)
  // second part of if statement checks for dtremain == update->dt so

  // note: if particle 'i' hits a surface, particle 'i' does not get 
  // ... PSURF flag. only particle 'j' (for reactive surfaces) get this flag

  if (dtremain < update->dt || pflag == PINSERT || pflag >= PSURF) {

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
  } // END check previous crossing
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

  // heaviside function
  double theta = 1.0;
  if (outface == XLO || outface == YLO || outface == ZLO) theta = -1.0;

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
  }
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
