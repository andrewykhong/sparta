/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
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
#include "fix_fnum_species.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "particle.h"
#include "mixture.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace SPARTA_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixFnumSpecies::FixFnumSpecies(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  if (particle->exist)
    error->all(FLERR,"Must call fix fnum/species before any particles exist");
  if (grid->cellweightflag)
    error->all(FLERR,"Fix fnum/species not currently supported for axi-symmetric");

  if (narg < 5) error->all(FLERR,"Not enough arguments for fix fnum/species command");
  
  // initialize array for storing weights of each species

  particle->specwtflag = 1;
  int nspecies = particle->nspecies;
  for (int i = 0; i < nspecies; i++)
    particle->species[i].specwt = 1.0;

  imix = particle->find_mixture(arg[2]);
  if (imix < 0)
    error->all(FLERR,"Fix fnum/species mixture ID does not exist");

  // relative to update->fnum or no?  

  int iarg = 3;
  int rflag = 0;
  if (strcmp(arg[iarg],"relative") == 0) {
    rflag = 1;
    iarg++;
  }

  // overwrite species dependent weights read in from species file
  // store raw fnum value (not relative to update->fnum)

  int isp;
  while (iarg < narg) {
    isp = particle->find_species(arg[iarg]);
    if (isp < 0) error->all(FLERR,"Undefined species specified in fix fnum/species");
    if (rflag) particle->species[isp].specwt = atof(arg[iarg+1]);
    else particle->species[isp].specwt = atof(arg[iarg+1])/update->fnum;
    iarg += 2;
  }

}

/* ---------------------------------------------------------------------- */

FixFnumSpecies::~FixFnumSpecies()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */
int FixFnumSpecies::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFnumSpecies::init()
{
  // comes after mixture init in particles so mixture will not reset
  // should not actually need to do this

  // reset weighted cummulatives
  int *fraction_flag = particle->mixture[imix]->fraction_flag;
  double *fraction_user = particle->mixture[imix]->fraction_user;
  double *fraction_wt = particle->mixture[imix]->fraction;
  double *cummulative_wt = particle->mixture[imix]->cummulative_wt;

  particle->mixture[imix]->init_fraction_wt(fraction_flag,fraction_user,fraction_wt,cummulative_wt);

  return;
}












