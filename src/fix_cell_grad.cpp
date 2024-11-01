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
#include "fix_cell_grad.h"
#include "update.h"
#include "grid.h"
#include "particle.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

/* ---------------------------------------------------------------------- */

FixCellGrad::FixCellGrad(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix cell/grad command");

  nevery = atoi(arg[2]);
  tstart = atof(arg[3]);
  tstop = atof(arg[4]);

  if (nevery <= 0) error->all(FLERR,"Illegal fix cell/grad command");
  if (tstart < 0.0 || tstop < 0.0)
    error->all(FLERR,"Illegal fix temp/rescale command");

  // optional keyword

  aveflag = 0;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix temp/rescale command");
      if (strcmp(arg[iarg+1],"yes") == 0) aveflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) aveflag = 0;
      else error->all(FLERR,"Invalid fix temp/rescale command");
      iarg += 2;
    } else error->all(FLERR,"Invalid fix temp/rescale command");
  }

  // per-cell array for aveflag = 1 case

  maxgrid = 0;
  vcom = NULL;
}

/* ---------------------------------------------------------------------- */

FixCellGrad::~FixCellGrad()
{
  if (copymode) return;

  //memory->destroy(vcom);
}

/* ---------------------------------------------------------------------- */

int FixCellGrad::setmask()
{
  int mask = 0;
  mask |= START_OF_STEP;
  mask |= MID_STEP;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCellGrad::init()
{
  //tprefactor = update->mvv2e / (3.0*update->boltz);
}

/* ---------------------------------------------------------------------- */

void FixCellGrad::start_of_step()
{

}

/* ---------------------------------------------------------------------- */

void FixCellGrad::mid_step()
{

}

/* ---------------------------------------------------------------------- */

void FixCellGrad::end_of_step()
{
  if (update->ntimestep % nevery) return;

  // set current t_target

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double t_target = tstart + delta * (tstop-tstart);

  // sort particles by grid cell if needed

  if (!particle->sorted) particle->sort();

  // 2 variants of thermostatting

  if (!aveflag) end_of_step_no_average(t_target);
  else end_of_step_average(t_target);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double FixCellGrad::memory_usage()
{
  double bytes = 0.0;
  bytes += maxgrid*3 * sizeof(double);    // vcom
  return bytes;
}
