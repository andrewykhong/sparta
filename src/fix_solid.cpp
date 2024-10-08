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
#include "fix_solid.h"
#include "comm.h"
#include "update.h"
#include "grid.h"
#include "domain.h"
#include "particle.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "compute.h"
#include "fix.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files
enum{NOMOVE,EULER,LANGEVIN};                  // type of solid particle move
enum{GREEN,BURT,EMPIRICAL};            // type of solid particle force

// for compute_solid_grid

enum{SIZE,MASS,TEMP,FORCEX,FORCEY,FORCEZ,HEAT};
enum{PZERO,COMPUTE,FIX};

#define DELTADELETE 1024

/* ---------------------------------------------------------------------- */

FixSolid::FixSolid(SPARTA *sparta, int narg, char **arg) :
  Fix(sparta, narg, arg)
{
  // flag for update particle custom attributes

  flag_update_custom = 1;

  npmax = 0; // size of solid particle list
  nsample = 0; // number of samples
  dim = domain->dimension;
  id = NULL; // solid particle list
  dellist = NULL; // particles to delete

  // for storing per-grid solid-particle properties
  per_grid_flag = 1;
  per_grid_freq = 1;
  // average force (3), heat flux, and particle params (3)
  size_per_grid_cols = 7;
  nglocal = 0;
  array_grid = NULL;

  if (narg < 9) error->all(FLERR,"Not enough arguments for fix solid command");

  solid_species = particle->find_species(arg[2]);
  if (solid_species < 0) error->all(FLERR,"Fix solid drag species does not exist");

  // surface collision models

  alpha = atof(arg[3]);
  eps = atof(arg[4]);

  // initial solid particle temperature and specific heat

  Rp0 = atof(arg[5]);
  rhop0 = atof(arg[6]);
  if (dim == 2) mp0 = 3.141598*Rp0*Rp0*rhop0;
  else mp0 = 4./3.*3.14159*pow(Rp0,3.0)*rhop0;
  Tp0 = atof(arg[7]);
  in_csp = atof(arg[8]);

  // optional args

  ifix = -1;
  nevery = 1; // frequency to update particle vel + pos
  reduce_size_flag = 0; // reduce size of particle radii over time due to heat
  move_type = EULER;
  force_type = GREEN;
  pwhich = PZERO;
  uxp0 = uyp0 = uzp0 = 0.0;
  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"reduce") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Invalid fix solid command");
      reduce_size_flag = 1;

      // grab idsource from fix or compute (character between '[' and ']')
      int n = strlen(arg[iarg+1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg+1][2]);

      // check if there the fix/compute points to certain index
      int argindex = 0;
      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ablate command");
        argindex = atoi(ptr+1);
        *ptr = '\0';
      }

      n = strlen(suffix) + 1;
      char *idsource = new char[n];
      strcpy(idsource,suffix);
      delete [] suffix;

      // find compute and check if it is valid
      if(strncmp(arg[iarg+1],"c_",2 == 0)) {
        pwhich = COMPUTE;

        ifc = modify->find_compute(idsource);
        if (ifc < 0) // does compute exist?
          error->all(FLERR,"Compute or fix ID for fix solid does not exist");
        // does compute have per-grid values?
        if (modify->compute[ifc]->per_grid_flag == 0)
          error->all(FLERR,
                     "Fix solid compute does not calculate per-grid values");
        // if compute does not reference certain index, is it a 1xN vector?
        if (argindex == 0 &&
            modify->compute[ifc]->size_per_grid_cols != 0)
          error->all(FLERR,"Fix solid compute does not "
                     "calculate per-grid vector");
        // if compute index specified, is there at least one value per grid?
        if (argindex && modify->compute[ifc]->size_per_grid_cols == 0)
          error->all(FLERR,"Fix solid compute does not "
                     "calculate per-grid array");
        // if compute index specified, is index within range?
        if (argindex && argindex > modify->compute[ifc]->size_per_grid_cols)
          error->all(FLERR,"Fix solid compute array is accessed out-of-range");

      // find fix and check if it is valid
      } else if (strncmp(arg[iarg+1],"f_",2 == 0)) {
        pwhich = FIX;

        ifc = modify->find_fix(idsource);
        if (ifc < 0)
          error->all(FLERR,"Fix ID for fix solid does not exist");
        if (modify->fix[ifc]->per_grid_flag == 0)
          error->all(FLERR,"Fix solid fix does not calculate per-grid values");
        if (argindex == 0 && modify->fix[ifc]->size_per_grid_cols != 0)
          error->all(FLERR,
                     "Fix solid fix does not calculate per-grid vector");
        if (argindex && modify->fix[ifc]->size_per_grid_cols == 0)
          error->all(FLERR,
                     "Fix solid fix does not calculate per-grid array");
        if (argindex && argindex > modify->fix[ifc]->size_per_grid_cols)
          error->all(FLERR,"Fix solid fix array is accessed out-of-range");
        // are intervals for fix and fix solid compatible?
        if (nevery % modify->fix[ifc]->per_grid_freq)
          error->all(FLERR,
                     "Fix for fix solid not computed at compatible time");

      } else pwhich = PZERO;

      iarg += 3;
    } else if (strcmp(arg[iarg],"nevery") == 0) { 
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix solid command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Nevery must be greater than zero");

      // check if nevery compatible with pre-defined fixes
      if (pwhich == FIX)
        if (nevery % modify->fix[ifc]->per_grid_freq)
          error->all(FLERR,
                     "Fix solid nevery not compatible with freq of pre-defined fix");

      iarg += 2;
    } else if (strcmp(arg[iarg],"brownian") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix solid command");
      if (strcmp(arg[iarg+1],"yes") == 0) move_type = LANGEVIN;
      else if (strcmp(arg[iarg+1],"no") == 0) move_type = EULER;
      else error->all(FLERR,"Invalid fix solid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nomove") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid fix solid command");
      move_type = NOMOVE;
      uxp0 = atof(arg[iarg+1]);
      uyp0 = atof(arg[iarg+2]);
      uzp0 = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"green") == 0) {
      force_type = GREEN;
      iarg += 1;
    } else if (strcmp(arg[iarg],"empirical") == 0) {
      force_type = EMPIRICAL;
      error->all(FLERR,"Empirical-based drag forces not currently implemented");
      iarg += 1;
    } else if (strcmp(arg[iarg],"burt") == 0) { 
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix solid command");
      force_type = BURT;
      error->all(FLERR,"Burt-based drag forces not currently implemented");
      iarg += 1;
    } else error->all(FLERR,"Invalid fix temp/rescale command");
  }

  // check if custom particle parameters exist

  index_solid_params = particle->find_custom((char *) "solid_params");
  index_solid_force = particle->find_custom((char *) "solid_force");
  index_solid_bulk = particle->find_custom((char *) "solid_bulk");

  if (index_solid_params >= 0 || index_solid_force >= 0 || index_solid_bulk >= 0)
    error->all(FLERR,"Fix solid_drag custom attribute already exists");

  index_solid_params = particle->add_custom((char *) "solid_params",DOUBLE,5);
  index_solid_force = particle->add_custom((char *) "solid_force",DOUBLE,4);
  index_solid_bulk = particle->add_custom((char *) "solid_bulk",DOUBLE,5);

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

FixSolid::~FixSolid()
{
  if (copy || copymode) return;

  memory->destroy(id);
  memory->destroy(dellist);
  memory->destroy(array_grid);

  delete random;

  particle->remove_custom(index_solid_params);
  particle->remove_custom(index_solid_force);
  particle->remove_custom(index_solid_bulk);
}

/* ---------------------------------------------------------------------- */

int FixSolid::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSolid::init()
{
  // initialize arrays

  npmax = DELTAPART;
  memory->create(id,npmax,"soliddrag:id");

  // pre-compute prefactors for force and heat flux delta functions

  if (force_type == GREEN) {
    F1 = 1.0+4.0/9.0*(1.0-eps)*(1.0-alpha);
    F2 = (1.0-eps)*alpha*sqrt(MY_PI)/3.0;
    Q1 = (1.0-eps)*alpha;
  } else if (force_type == BURT) { 
    // alpha here is tau in Burt model
    F1 = 1.0;
    F2 = alpha*sqrt(2.0*MY_PI)/3.0;
    Q1 = alpha;
  }

  reallocate();
}

/* ----------------------------------------------------------------------
   set solid properties for new particles 
   ignore olds ones not to override existing velocities
------------------------------------------------------------------------- */

void FixSolid::update_custom(int index, double,
                                double, double,
                                double*)
{
  int ispecies = particle->particles[index].ispecies;
  if (ispecies != solid_species) return;

  double **solid_array = particle->edarray[particle->ewhich[index_solid_params]];
  double **solid_force = particle->edarray[particle->ewhich[index_solid_force]];
  double **solid_bulk  = particle->edarray[particle->ewhich[index_solid_bulk]];

  if (solid_array[index][0] != 0.0) return; // new particles will be zero

  solid_array[index][0] = 1.0; // existing particles are 1.0 (maybe not needed)
  solid_array[index][1] = Rp0; // radius
  solid_array[index][2] = mp0; // mass
  solid_array[index][3] = Tp0; // temperature
  solid_array[index][4] = in_csp; // specific heat

  solid_force[index][0] = 0.0; // Fx
  solid_force[index][1] = 0.0; // Fy
  solid_force[index][2] = 0.0; // Fz
  solid_force[index][3] = 0.0; // Q

  solid_bulk[index][0] = 0.0; // Ugx
  solid_bulk[index][1] = 0.0; // Ugy
  solid_bulk[index][2] = 0.0; // Ugz
  solid_bulk[index][3] = 0.0; // Tg
}

/* ---------------------------------------------------------------------- */

void FixSolid::end_of_step()
{
  if (!particle->sorted) particle->sort();
  reallocate(); // for outputting average force in each grid

  if (move_type == NOMOVE) reset_velocities(1);

  // force model type

  update_force();

  if (update->ntimestep % nevery) return;

  // update particles velocitis (and position if langevin)

  ndelete = 0;
  if (move_type == EULER) update_particle();
  //else if (move_type == LANGEVIN) move_langevin();
  else if (move_type == NOMOVE) reset_velocities(0); // zero out velocities

  // delete solid particles with no mass

  if (ndelete) particle->compress_reactions(ndelete,dellist);
}

/* ----------------------------------------------------------------------
   Momentum and heat fluxes based on Green's functions
   (Gallis, Torczynski, Rader - Phys. of Fluids - 2001)
---------------------------------------------------------------------- */

void FixSolid::update_force()
{
  // grab various particle and grid quantities

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // solid particle related vectors

  double **solid_array = particle->edarray[particle->ewhich[index_solid_params]];
  double **solid_force = particle->edarray[particle->ewhich[index_solid_force]];
  double **solid_bulk  = particle->edarray[particle->ewhich[index_solid_bulk]];

  int i,icell,ip,is; // dummy indices
  int ispecies, sid; // species of particle i; particle indices of solid particles
  int np, nsolid; // number of total simulators; number of solid simulators
  double Rp,mp,Tp,csp; // particle radius, mass, temperature, and specific heat
  double cx,cy,cz,cmag; // thermal velocity of gas and solid particles
  double *u,*up;  // velocity of gas and solid particles
  double um[3]; // drift velocity (zero if no charge or gravity) 
  double Fg[3],Qg; // force and heat flux as defined by Green function
  double totalmass; // for calculating temperature
  double mass,T; // gas particle mass, and temperature
  double csx,cp; // solid particle cross section and thermal velocity
  double frat; // ratio of gas to solid weight
  double prefactor; // prefactor

  for (icell = 0; icell < nglocal; icell++) {

    np = cinfo[icell].count;
    if (np <= 1) continue;

    // only need ids

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(id);
      memory->create(id,npmax,"soliddrag:id");
    }

    // get solid particle velocities

    nsolid = 0;
    um[0] = 0.0;
    um[1] = 0.0;
    um[2] = 0.0;
    totalmass = 0.0;

    // drift belocity defined as F / beta where
    // ... F is the nondrag force?

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      if (ispecies == solid_species) id[nsolid++] = ip;
      else {
        mass = species[ispecies].mass;
        u = particles[ip].v;
        totalmass += mass;
        um[0] += mass*u[0];
        um[0] += mass*u[1];
        um[0] += mass*u[2];
      }
      ip = next[ip];
    }

    um[0] /= totalmass;
    um[1] /= totalmass;
    um[2] /= totalmass;

    // calculate gas state

    T = 0;
    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      if (ispecies != solid_species) {
        mass = species[ispecies].mass;
        u = particles[ip].v;

        cx = u[0]-um[0];
        cy = u[1]-um[1];
        cz = u[2]-um[2];
        T += mass*(cx*cx+cy*cy+cz*cz);
      }
      ip = next[ip];
    }

    // removed half from numer. and denom.

    T /= (3.0*update->boltz*(np-nsolid));

    // calculate incident forces and heat flux

    for (is = 0; is < nsolid; is++) {
      sid = id[is];
      up = particles[sid].v;
      Rp = solid_array[sid][1];
      mp = solid_array[sid][2];
      Tp = solid_array[sid][3];
      csp = solid_array[sid][4];

      Fg[0] = Fg[1] = Fg[2] = Qg = 0.0;

      ip = cinfo[icell].first;
      while (ip >= 0) {
        ispecies = particles[ip].ispecies;
        if (ispecies != solid_species) {
          mass = species[ispecies].mass;

          // account for difference in species weight
          u = particles[ip].v;

          cx = u[0]-up[0];
          cy = u[1]-up[1];
          cz = u[2]-up[2];
          if (dim == 2) cz = 0.0;
          cmag = sqrt(cx*cx+cy*cy+cz*cz);

          // mean thermal speed of outgoing particles
          if (force_type == GREEN) {
            cp = sqrt(2.0*update->boltz*Tp/mass);

            Fg[0] += mass*cx*(F1*cmag+F2*cp);
            Fg[1] += mass*cy*(F1*cmag+F2*cp);
            if (dim == 3)
              Fg[2] += mass*cz*(F1*cmag+F2*cp);
            Qg    += mass*cmag*Q1*(0.5*cmag*cmag-cp*cp);
            //Qg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          } else if (force_type == BURT) {
            cp = sqrt(update->boltz*Tp*mass);

            Fg[0] += cx*(mass*cmag+F2*cp);
            Fg[1] += cy*(mass*cmag+F2*cp);
            if (dim == 3)
              Fg[2] += cz*(mass*cmag+F2*cp);
            Qg += cmag*Q1*(0.5*mass*cmag*cmag-2.0*update->boltz*Tp);
            //Qg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          }
        }
        ip = next[ip];
      } // end while

      csx = Rp*Rp*MY_PI;
      prefactor = csx * update->fnum / cinfo[icell].volume;
      Fg[0] *= prefactor;
      Fg[1] *= prefactor;
      Fg[2] *= prefactor;
      Qg    *= prefactor;

      // update per-particle forces

      solid_force[sid][0] = (solid_force[sid][0]*nsample+Fg[0])/(nsample+1.0);
      solid_force[sid][1] = (solid_force[sid][1]*nsample+Fg[1])/(nsample+1.0);
      solid_force[sid][2] = (solid_force[sid][2]*nsample+Fg[2])/(nsample+1.0);
      solid_force[sid][3] = (solid_force[sid][3]*nsample+Qg)/(nsample+1.0);

      solid_bulk[sid][0] = (solid_bulk[sid][0]*nsample+um[0])/(nsample+1.0);
      solid_bulk[sid][1] = (solid_bulk[sid][1]*nsample+um[1])/(nsample+1.0);
      solid_bulk[sid][2] = (solid_bulk[sid][2]*nsample+um[2])/(nsample+1.0);
      solid_bulk[sid][3] = (solid_bulk[sid][3]*nsample+T)/(nsample+1.0);

    } // end solid loop
  } // end cells

  // number of samples
  nsample++;

}

/* ----------------------------------------------------------------------
   Momentum and heat fluxes based on Burt Model
   (Burt, Boyd - AIAA - 2004)
---------------------------------------------------------------------- */

void FixSolid::force_burt)
{
  // grab various particle and grid quantities

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // solid particle related vectors

  double **solid_array = particle->edarray[particle->ewhich[index_solid_params]];
  double **solid_force = particle->edarray[particle->ewhich[index_solid_force]];
  double **solid_bulk  = particle->edarray[particle->ewhich[index_solid_bulk]];

  int i,icell,ip,is; // dummy indices
  int ispecies, sid; // species of particle i; particle indices of solid particles
  int np, nsolid; // number of total simulators; number of solid simulators
  double Rp,mp,Tp,csp; // particle radius, mass, temperature, and specific heat
  double cx,cy,cz,cmag; // thermal velocity of gas and solid particles
  double *u,*up;  // velocity of gas and solid particles
  double um[3]; // drift velocity (zero if no charge or gravity) 
  double Fg[3],Qg; // force and heat flux as defined by Green function
  double totalmass; // for calculating temperature
  double mass,T; // gas particle mass, and temperature
  double csx,cp; // solid particle cross section and thermal velocity
  double frat; // ratio of gas to solid weight
  double prefactor; // prefactor

  for (icell = 0; icell < nglocal; icell++) {

    np = cinfo[icell].count;
    if (np <= 1) continue;

    // only need ids

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(id);
      memory->create(id,npmax,"soliddrag:id");
    }

    // calculate mean gas velocity and record solid particle id's

    nsolid = 0;
    um[0] = 0.0;
    um[1] = 0.0;
    um[2] = 0.0;
    totalmass = 0.0;

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      if (ispecies == solid_species) id[nsolid++] = ip;
      else {
        mass = species[ispecies].mass;
        u = particles[ip].v;
        totalmass += mass;
        um[0] += mass*u[0];
        um[0] += mass*u[1];
        um[0] += mass*u[2];
      }
      ip = next[ip];
    }

    um[0] /= totalmass;
    um[1] /= totalmass;
    um[2] /= totalmass;

    // calculate gas state

    T = 0;
    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      if (ispecies != solid_species) {
        mass = species[ispecies].mass;
        u = particles[ip].v;

        cx = u[0]-um[0];
        cy = u[1]-um[1];
        cz = u[2]-um[2];
        T += mass*(cx*cx+cy*cy+cz*cz);
      }
      ip = next[ip];
    }

    // removed half from numer. and denom.

    T /= (3.0*update->boltz*(np-nsolid));

    // calculate incident forces and heat flux

    for (is = 0; is < nsolid; is++) {
      sid = id[is];
      up = particles[sid].v;
      Rp = solid_array[sid][1];
      mp = solid_array[sid][2];
      Tp = solid_array[sid][3];
      csp = solid_array[sid][4];

      Fg[0] = Fg[1] = Fg[2] = Qg = 0.0;

      ip = cinfo[icell].first;
      while (ip >= 0) {
        ispecies = particles[ip].ispecies;
        if (ispecies != solid_species) {
          mass = species[ispecies].mass;

          // account for difference in species weight
          u = particles[ip].v;

          cx = u[0]-up[0];
          cy = u[1]-up[1];
          cz = u[2]-up[2];
          if (dim == 2) cz = 0.0;
          cmag = sqrt(cx*cx+cy*cy+cz*cz);

          // mean thermal speed of outgoing particles
          cp = sqrt(2.0*MY_PI*update->boltz*Tp*mass);

          Fg[0] += mass*cx*(F1*cmag+F2*cp);
          Fg[1] += mass*cy*(F1*cmag+F2*cp);
          if (dim == 3)
            Fg[2] += mass*cz*(F1*cmag+F2*cp);
          Qg    += mass*cmag*Q1*(0.5*cmag*cmag-cp*cp);
        }
        ip = next[ip];
      } // end while

      csx = Rp*Rp*MY_PI;
      prefactor = csx * update->fnum / cinfo[icell].volume;
      Fg[0] *= prefactor;
      Fg[1] *= prefactor;
      Fg[2] *= prefactor;
      Qg    *= prefactor;

      // update per-particle forces

      solid_force[sid][0] = (solid_force[sid][0]*nsample+Fg[0])/(nsample+1.0);
      solid_force[sid][1] = (solid_force[sid][1]*nsample+Fg[1])/(nsample+1.0);
      solid_force[sid][2] = (solid_force[sid][2]*nsample+Fg[2])/(nsample+1.0);
      solid_force[sid][3] = (solid_force[sid][3]*nsample+Qg)/(nsample+1.0);

      solid_bulk[sid][0] = (solid_bulk[sid][0]*nsample+um[0])/(nsample+1.0);
      solid_bulk[sid][1] = (solid_bulk[sid][1]*nsample+um[1])/(nsample+1.0);
      solid_bulk[sid][2] = (solid_bulk[sid][2]*nsample+um[2])/(nsample+1.0);
      solid_bulk[sid][3] = (solid_bulk[sid][3]*nsample+T)/(nsample+1.0);

    } // end solid loop
  } // end cells

  // number of samples
  nsample++;

}


/* ----------------------------------------------------------------------
   Update velocities and positions using simple Euler scheme
   (no Brownian motion)
---------------------------------------------------------------------- */

void FixSolid::update_particle()
{
  // grab various particle and grid quantities

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // solid particle related vectors

  double **solid_array = particle->edarray[particle->ewhich[index_solid_params]];
  double **solid_force = particle->edarray[particle->ewhich[index_solid_force]];
  double **solid_bulk  = particle->edarray[particle->ewhich[index_solid_bulk]];

  // indices

  int i,icell,ip,ispecies,np;

  // particle related

  double v[3],cx,cy,cz,csq;
  double Rp,mp,Tp,csp;
  double Q,mp_loss,rho,new_vol;

  int nsolid;

  for (icell = 0; icell < nglocal; icell++) {

    // reset array_grid (for output)
    for (i = 0; i < size_per_grid_cols; i++)
      array_grid[icell][i] = 0.0;

    np = cinfo[icell].count;
    if (np <= 0) continue;

    nsolid = 0;

    // update only solid particles

    ip = cinfo[icell].first;
    while (ip >= 0) {
      if (particles[ip].ispecies == solid_species) {
        Rp = solid_array[ip][1];
        mp = solid_array[ip][2];
        Tp = solid_array[ip][3];
        csp = solid_array[ip][4];

        // velocities
        for (int d = 0; d < dim; d++)
          particles[ip].v[d] = particles[ip].v[d] + solid_force[ip][d]*update->dt;
        
        // temperature
        Tp = Tp + solid_force[ip][3]*update->dt/csp/mp;
        solid_array[ip][3] = Tp;

        // update radius and mass (TODO: add later)

        if (reduce_size_flag) {

          // grab pressure in cell
          double p; // pressure
          if (pwhich == PZERO) p = 0;
          //else if (pwhich == COMPUTE)

          // calcaulte saturatino temperautre (Fanale and Salvail 1984)
          double psat = (3.56e12) * exp(-(6147.667/Tp));

          // Hertz Knudsen Equation for flux of sulbiming molecules

          // correction coefficients essentially
          double beta = 1.0;
          double alpha = 0.83; // rec from Gadsden 1998 and Winkler 2012 
          double m_h2o = 1.0; // TODO: replace with the mass from species file
          // units = (1/A) * 1/s
          double flux = (psat - p)/sqrt(2.0*3.14159*m_h2o*update->boltz*Tp);

          // does it sublimate?
          if (flux > 0) {
            double area = 3.14159*4.0*Rp*Rp;
            double srho = 1.0; // TODO: replace with actual mass density
            double mloss = flux * area * srho * update->dt;

            // update mass
            mp -= mloss;

            // recalculate radius
            double vol = mp/srho;
            Rp = pow(vol*0.75/3.14159,1./3.);
          }

        }

        // update per-grid forces for outputting

        array_grid[icell][0] += Rp;
        array_grid[icell][1] += mp;
        array_grid[icell][2] += Tp;
        array_grid[icell][3] += solid_force[ip][0];
        array_grid[icell][4] += solid_force[ip][1];
        array_grid[icell][5] += solid_force[ip][2];
        array_grid[icell][6] += solid_force[ip][3];
        nsolid++;

      } // end check species

      // reset forces and heat fluxes
      for (i = 0; i < 4; i++) solid_force[ip][i] = 0.0;
      ip = next[ip];

    } // end particle

    // update per-grid forces for outputting
    if (nsolid > 0)
      for (i = 0; i < 7; i++)
        array_grid[icell][i] /= nsolid;
    else 
      for (i = 0; i < 7; i++)
        array_grid[icell][i] = 0.0;

  } // end cells

  // reset number of samples
  nsample = 0;
}

/* ----------------------------------------------------------------------
   Update velocities and positions due to Brownian motion

   Algorithm by Ermak / Buckholz  
   (Ermak and Buckholz - J. Comp. Phys. - 1979)
---------------------------------------------------------------------- */

void FixSolid::move_langevin()
{
  // grab various particle and grid quantities

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;

  // solid particle related vectors

  double **solid_array = particle->edarray[particle->ewhich[index_solid_params]];
  double **solid_force = particle->edarray[particle->ewhich[index_solid_force]];
  double **solid_bulk  = particle->edarray[particle->ewhich[index_solid_bulk]];

  // indices

  int i,icell,ip,ispecies,np;

  double v[3],cx,cy,cz,csq;
  double Rp,mp,Tp,csp;
  double Fx,Fy,Fz,Q,Ux,Uy,Uz,Tg;
  double beta, tau, delta;
  double cp, stdevB1, stdevB2, B1, B2;
  double rho, mp_loss, new_vol;

  for (icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 0) continue;

    // update only solid particles

    ip = cinfo[icell].first;
    while (ip >= 0) {
      if (particles[ip].ispecies == solid_species) {
        Rp = solid_array[ip][1];
        mp = solid_array[ip][2];
        Tp = solid_array[ip][3];
        csp = solid_array[ip][4];

        Fx = solid_force[ip][0];
        Fy = solid_force[ip][1];
        Fz = solid_force[ip][2];
        Q  = solid_force[ip][3];

        Ux = solid_bulk[ip][0];
        Uy = solid_bulk[ip][1];
        Uz = solid_bulk[ip][2];
        Tg = solid_bulk[ip][3];

        memcpy(v,particles[ip].v,3*sizeof(double));

        // peculiar velocities

        cx = v[0]-Ux;
        cy = v[1]-Uy;
        cz = v[2]-Uz;

        // find effective friction (drag) coefficient beta as:
        // abs(Fg . up) / (up . up)
        // then find stopping time: tau = mp / beta
 
        csq = (cx*cx+cy*cy+cz*cz);
        beta = fabs(Fx*cx+Fy*cy+Fz*cz)/csq;
        tau = mp / beta;

        // delta (~time step) = exp(-dt / tau)

        delta = exp(-update->dt/tau);

        // sample random normal variables

        cp = sqrt(2.0*update->boltz*Tg/mp);
        stdevB1 = sqrt(1.5*(1.0-delta*delta))*cp;
        stdevB2 = update->dt/tau - 2.0*(1.0-delta)/(1.0+delta);
        stdevB2 = sqrt(3.0)*cp*tau*sqrt(stdevB2);

        B1 = random_normal()*stdevB1;
        B2 = random_normal()*stdevB2;

        // velocities and positions
        for (int d = 0; d < dim; d++) {
          particles[ip].x[d] = particles[ip].x[d] + 
            (particles[ip].v[d] / beta)*(1.0-delta) + B2;
          particles[ip].v[d] = particles[ip].v[d]*delta + B1;
        }
        
        // temperature
        Tp = Tp + Q*update->dt/csp/mp;

        // update particle radii
        if (reduce_size_flag && solid_array[ip][3] > Tvap) {

          rho = mp/(4.0/3.0*MY_PI*Rp*Rp*Rp);
          solid_array[ip][3] = Tvap;
          Q = Q/(4.0*Rp*Rp*MY_PI);
          mp_loss = Q/(hvap-hsolid)*update->dt*nevery;

          // updated particle radius
          if (mp > mp_loss) {
            new_vol = (mp - mp_loss) / rho;
            solid_array[ip][1] = pow(0.75*new_vol/MY_PI,1./3.);
          } else {
            solid_array[ip][1] = 0.0;

            // delete particles with no mass
            if (ndelete == maxdelete) {
              maxdelete += DELTADELETE;
              memory->grow(dellist,maxdelete,"solid:dellist");
            }
            dellist[ndelete++] = ip;
          }
        }

      }

      // reset forces and heat fluxes
      for (i = 0; i < 8; i++) solid_force[ip][i] = 0.0;

      ip = next[ip];
    }

  } // end cells

  // reset number of samples
  nsample = 0;
}

/* ----------------------------------------------------------------------
   Update velocities and positions using simple Euler scheme
   (no Brownian motion)
---------------------------------------------------------------------- */

void FixSolid::reset_velocities(int reset_flag)
{
  // grab various particle and grid quantities

  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;
  int *next = particle->next;
  Grid::ChildInfo *cinfo = grid->cinfo;

  int icell,ip; // dummy indices
  int np; // number of particles
  double v[3]; // particle velocity

  for (icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 0) continue;

    // update only solid particles

    ip = cinfo[icell].first;
    while (ip >= 0) {
      if (particles[ip].ispecies == solid_species) {
        if (reset_flag) {
          particles[ip].v[0] = uxp0;
          particles[ip].v[1] = uyp0;
          particles[ip].v[2] = uzp0;
        } else {
          particles[ip].v[0] = 0.0;
          particles[ip].v[1] = 0.0;
          particles[ip].v[2] = 0.0;
        }
      }
      ip = next[ip];
    } // end while
  } // end cells
}

/* ----------------------------------------------------------------------
   Get solid particle propery
---------------------------------------------------------------------- */

double FixSolid::get_particle_property(int which, int ip)
{
  // grab various particle and grid quantities

  Particle::OnePart *particles = particle->particles;

  if (particles[ip].ispecies != solid_species) return -1;

  // solid particle related vectors

  double **solid_array = particle->edarray[particle->ewhich[index_solid_params]];
  double **solid_force = particle->edarray[particle->ewhich[index_solid_force]];

  if (which == SIZE)        return solid_array[ip][1];
  else if (which == MASS)   return solid_array[ip][2];
  else if (which == TEMP)   return solid_array[ip][3];
  else if (which == FORCEX) return solid_force[ip][0];
  else if (which == FORCEY) return solid_force[ip][1];
  else if (which == FORCEZ) return solid_force[ip][2];
  else if (which == HEAT)   return solid_force[ip][3];
  return -1;

}

/* ----------------------------------------------------------------------
   gaussian RN based on Box Mueller transform
------------------------------------------------------------------------- */

double FixSolid::random_normal()
{
  double U1 = random->uniform();
  double U2 = random->uniform();

  double randnorm = sqrt(-2.0*log(U1));
  if (random->uniform() > 0.5) randnorm *= cos(2.0*MY_PI*U2);
  else randnorm *= sin(2.0*MY_PI*U2);

  return randnorm;
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void FixSolid::reallocate()
{
  if (grid->nlocal == nglocal) return;

  memory->destroy(array_grid);
  nglocal = grid->nlocal;
  memory->create(array_grid,nglocal,size_per_grid_cols,"fix/solid:array_grid");

  // initialize values
  for (int i = 0; i < nglocal; i++)
    for (int j = 0; j < size_per_grid_cols; j++)
      array_grid[i][j] = 0.0;
}
