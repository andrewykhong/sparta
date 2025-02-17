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

enum{INT,DOUBLE};                      // several files

// cell face quantities
// mass density, x-velocity, y-velocity, z-velocity

enum{RHO,U,V,W,UMAG,VMAG,WMAG,FACELASTSIZE};

// cell bulk quantitis

enum{RHO_BULK,U_BULK,V_BULK,W_BULK,UVWSQ_BULK,CELLLASTSIZE};

/* ---------------------------------------------------------------------- */

FixCellGrad::FixCellGrad(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix cell/grad command");

  // how frequently to update cell fluxes
  nevery = atoi(arg[2]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix cell/grad command");
  T_interval = nevery * update->dt;

  // time average the quantities?
  int iarg = 3;
  aveflag = 0;
  if (strcmp(arg[iarg],"ave") == 0) {
    aveflag = 1;
    iarg++;
  }

  // specify pairs
  // first is the quantity (velocity, temperature, density, etc.)
  // second is the direction (x,y,z)
  
  // keep track which face values to find 
  faceids = new int[FACELASTSIZE];
  for (int i = 0; i < FACELASTSIZE; i++) faceids[i] = 0;

  // keep track which cell values to find 
  cellids = new int[CELLLASTSIZE];
  for (int i = 0; i < CELLLASTSIZE; i++) cellids[i] = 0;

  xface = yface = zface = 0.0; // which cell faces need to be computed
  while (iarg < narg) {
    if (strcmp(arg[iarg],"rho_x") == 0) {
      faceids[RHO] = 1;
      xface = 1;
    } else if (strcmp(arg[iarg],"rho_y") == 0) {
      faceids[RHO] = 1;
      yface = 1;
    } else if (strcmp(arg[iarg],"rho_z") == 0) {
      faceids[RHO] = 1;
      zface = 1;
    } else if (strcmp(arg[iarg],"u_x") == 0) {
      faceids[U] = 1;
      xface = 1;
      cellids[RHO_BULK] = 1;
      cellids[U_BULK] = 1;
    } else if (strcmp(arg[iarg],"u_y") == 0) {
      faceids[U] = 1;
      faceids[V] = 1;
      yface = 1;
      cellids[RHO_BULK] = 1;
      cellids[U_BULK] = 1;
    } else if (strcmp(arg[iarg],"u_z") == 0) {
      faceids[U] = 1;
      faceids[W] = 1;
      zface = 1;
      cellids[RHO_BULK] = 1;
      cellids[U_BULK] = 1;
    } else if (strcmp(arg[iarg],"v_x") == 0) {
      faceids[V] = 1;
      faceids[U] = 1;
      xface = 1;
      cellids[RHO_BULK] = 1;
      cellids[V_BULK] = 1;
    } else if (strcmp(arg[iarg],"v_y") == 0) {
      faceids[V] = 1;
      yface = 1;
      cellids[RHO_BULK] = 1;
      cellids[V_BULK] = 1;
    } else if (strcmp(arg[iarg],"v_z") == 0) {
      faceids[V] = 1;
      faceids[W] = 1;
      zface = 1;
      cellids[RHO_BULK] = 1;
      cellids[V_BULK] = 1;
    } else if (strcmp(arg[iarg],"w_x") == 0) {
      faceids[W] = 1;
      faceids[U] = 1;
      xface = 1;
      cellids[RHO_BULK] = 1;
      cellids[W_BULK] = 1;
    } else if (strcmp(arg[iarg],"w_y") == 0) {
      faceids[W] = 1;
      faceids[V] = 1;
      yface = 1;
      cellids[RHO_BULK] = 1;
      cellids[W_BULK] = 1;
    } else if (strcmp(arg[iarg],"w_z") == 0) {
      faceids[W] = 1;
      zface = 1;
      cellids[RHO_BULK] = 1;
      cellids[W_BULK] = 1;
    } else error->all(FLERR,"Invalid fix temp/rescale command");

    iarg += 2;
  }

  // count number of values to track on faces and in cells

  int nfacevalues = 0;
  for (int i = 0; i < FACELASTSIZE; i++)
    if (faceids[i]) nfacevalues++;

  int ncellvalues = 0;
  for (int i = 0; i < CELLLASTSIZE; i++)
    if (cellids[i]) ncellvalues++;

  if (nfacevalues == 0)
    error->one(FLERR,"No face values specified in fix cell grad");

  // nface = number of faces to sum over

  int nface = domain->dimension*2; // easiest to reference if all faces tracked
  //if (xface) nface+=2;
  //if (yface) nface+=2;
  //if (zface) nface+=2;
  //if (nface < 1 || nface > domain->dimension)
  //  error->one(FLERR,"Incorrect number of faces found in fix cell grad");

  // store quantities in custom grid arrays
  // check if custom per-particel attribute exists

  cellbulkindex = grid->find_custom((char *) "cellbulk");
  cellfaceindex = grid->find_custom((char *) "cellface");

  if (cellfaceindex >= 0 || cellbulkindex >= 0)
    error->all(FLERR,"Fix cell gradient already exists");

  cellbulkindex = grid->add_custom((char *) "cellbulk", DOUBLE, ncellvalues);
  cellfaceindex = grid->add_custom((char *) "cellface", DOUBLE, nface*nfacevalues);

  // per-cell array for aveflag = 1 case

  maxgrid = 0;
}

/* ---------------------------------------------------------------------- */

FixCellGrad::~FixCellGrad()
{
  if (copymode) return;

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
  return;
  //tprefactor = update->mvv2e / (3.0*update->boltz);
}

/* ---------------------------------------------------------------------- */

void FixCellGrad::end_of_step()
{
  if (update->ntimestep % nevery) return;

  // set current t_target

  //double delta = update->ntimestep - update->beginstep;
  //if (delta != 0.0) delta /= update->endstep - update->beginstep;
  //double t_target = tstart + delta * (tstop-tstart);

  // sort particles by grid cell if needed

  if (!particle->sorted) particle->sort();

  // 2 variants of thermostatting

  //if (!aveflag) end_of_step_no_average(t_target);
  //else end_of_step_average(t_target);
}

/* ---------------------------------------------------------------------- */
// ip - particle index
// pflag - type of particle
// icell - which grid cell 
// outface - which face particle is exitting to
// frac - fraction of dtremain particle is in cell icell
void FixCellGrad::mid_step(int ip, int pflag, int icell, int outface, double dtremain, double frac)
{
  // get particle mass
  int isp = particle[ip].species;
  double pmass = particle->species[isp].mass;

  // grab the perg-grid custom arrays
  double *cell = grid->edarray[particle->ewhich[cellbulkindex]];
  double *face = grid->edarray[particle->ewhich[cellfaceindex]];

  // update bulk properties in cell
  // consider residence time

  cell[icell][0] += mass*dtin;
  cell[icell][1] += mass*v[0]*dtin;
  cell[icell][2] += mass*v[1]*dtin;
  cell[icell][3] += mass*v[2]*dtin;

  // check face the particle WILL cross
  if (frac < 1.0) {

    // heaviside function
    double theta = 1.0;
    if (outface == XLO || outface == YLO || outface == ZLO) theta = -1.0;

    // update cell face quants
    // faces follow same enum order
    // xlo -> 0, xhi -> 1, ylo -> 2, ...
    int ivalue = 0;
    int i_index;
    for (int ival = 0; ival < FACELASTSIZE; ival++) {
      if (faceids[ival]) {
        i_index = ivalue+outface;
        if (ival == RHO) face[icell][i_index] += pmass;
        else if (ival == U) face[icell][i_index] += pmass*v[0];
        else if (ival == V) face[icell][i_index] += pmass*v[1];
        else if (ival == W) face[icell][i_index] += pmass*v[2];
        else if (ival == UMAG) face[icell][i_index] += pmass*fabs(v[0]);
        else if (ival == VMAG) face[icell][i_index] += pmass*fabs(v[1]);
        else if (ival == WMAG) face[icell][i_index] += pmass*fabs(v[2]);
        ivalue += (DIM*2);
      }
    } // END for face values

    // update bulk cell quants
    ivalue = 0;
    double dtcell = dtremain * frac;
    for (int ival = 0; ival < CELLLASTSIZE; ival++) {
      if (cellids[ival]) {
        if (ival == RHO_BULK) cell[icell][ivalue++] += pmass;
        else if (ival == U_BULK) cell[icell][ivalue++] += pmass*v[0];
        else if (ival == V_BULK) cell[icell][ivalue++] += pmass*v[1];
        else if (ival == W_BULK) cell[icell][ivalue++] += pmass*v[2];
        else if (ival == UVQSQ_BULK)
          cell[icell][ivalue++] +=
            pmass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
      }
    } // END for cell values

  } // END check future crossing

  // check history if particle HAS crossed a cell face
  // also considers any new particles created due to fix_emit
  if (particles[i].dtremain < dt || pflag == PINSERT) {

    // find which face the particle is closest to
    double mindist;
    int inface;

    inface = XLO;
    mindist = fabs(x[0] - lo[0]);
    if (fabs(x[0]-hi[0]) < mindist) inface = XHI;
    if (fabs(x[1]-lo[1]) < mindist) inface = YLO;
    if (fabs(x[1]-hi[1]) < mindist) inface = YHI;
    if (fabs(x[2]-lo[2]) < mindist) inface = ZLO;
    if (fabs(x[2]-hi[2]) < mindist) inface = ZHI;

    // heaviside function
    double theta = 1.0;
    if (inface == XLO || inface == YLO || inface == ZLO) theta = -1.0;

    // update cell face quants
    int ivalue = 0;
    int i_index;
    for (int ival = 0; ival < FACELASTSIZE; ival++) {
      if (faceids[ival]) {
        i_index = ivalue+inface;
        if (ival == RHO) face[icell][i_index] += pmass;
        else if (ival == U) face[icell][i_index] += pmass*v[0];
        else if (ival == V) face[icell][i_index] += pmass*v[1];
        else if (ival == W) face[icell][i_index] += pmass*v[2];
        else if (ival == UMAG) face[icell][i_index] += pmass*fabs(v[0]);
        else if (ival == VMAG) face[icell][i_index] += pmass*fabs(v[1]);
        else if (ival == WMAG) face[icell][i_index] += pmass*fabs(v[2]);
        ivalue += (DIM*2);
      }
    }

  } // END check previous crossing
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixCellGrad::memory_usage()
{
  double bytes = 0.0;
  //bytes += maxgrid*3 * sizeof(double);    // vcom
  return bytes;
}
