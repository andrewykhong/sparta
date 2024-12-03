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
enum{NOFORCE,GREEN,BURT,LOTH,SINGH};            // type of solid particle force

// for compute_solid_grid

enum{SIZE,MASS,TEMP,FORCEX,FORCEY,FORCEZ,HEAT};
enum{PZERO,AVERAGE};

#define DELTADELETE 1024
#define MAXLINE 1024

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

  if (narg < 7) error->all(FLERR,"Not enough arguments for fix solid command");

  // read in initial solid particle properties

  fp = fopen(arg[2],"r");
  if (!fp) error->one(FLERR,"Solid species file not found");
  read_solid();

  // initial solid particle temperature and specific heat

  Tp0 = atof(arg[3]);

  // compute for temperature and pressure

  // optional args

  ifix = -1;
  nevery = 1; // frequency to update particle vel + pos
  phase_flag = 0; // reduce size of particle radii over time due to heat
  move_type = EULER;
  force_type = NOFORCE;
  pwhich = PZERO;
  uxp0 = uyp0 = uzp0 = 0.0;

  int iarg = 4;
  while (iarg < narg) {
    // to model particle mass loss
    if (strcmp(arg[iarg],"reduce") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix solid command");
      phase_flag = 1;
      if (strcmp(arg[iarg+1],"zero") == 0) pwhich = PZERO;
      else if (strcmp(arg[iarg+1],"average") == 0) pwhich = AVERAGE;
      else error->all(FLERR,"Invalid pressure choice for fix solid command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nevery") == 0) { 
      if (iarg+2 > narg) error->all(FLERR,"Invalid fix solid command");
      nevery = atoi(arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Nevery must be greater than zero");
      iarg += 2;
    } else if (strcmp(arg[iarg],"brownian") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Invalid fix solid command");
      move_type = LANGEVIN;
      iarg += 1;
    } else if (strcmp(arg[iarg],"euler") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Invalid fix solid command");
      move_type = EULER;
      iarg += 1;
    } else if (strcmp(arg[iarg],"nomove") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid fix solid command");
      move_type = NOMOVE;
      uxp0 = atof(arg[iarg+1]);
      uyp0 = atof(arg[iarg+2]);
      uzp0 = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"green") == 0) {
      force_type = GREEN;
      iarg++;
    } else if (strcmp(arg[iarg],"burt") == 0) { 
      force_type = BURT;
      iarg++;
    } else error->all(FLERR,"Invalid fix solid command");
  }

  // check if custom particle parameters exist

  index_solid_params = particle->find_custom((char *) "solid_params");
  index_solid_force = particle->find_custom((char *) "solid_force");
  index_solid_bulk = particle->find_custom((char *) "solid_bulk");

  if (index_solid_params >= 0 || index_solid_force >= 0 || index_solid_bulk >= 0)
    error->all(FLERR,"Fix solid_drag custom attribute already exists");

  index_solid_params = particle->add_custom((char *) "solid_params",DOUBLE,5);
  index_solid_force = particle->add_custom((char *) "solid_force",DOUBLE,4);
  index_solid_bulk = particle->add_custom((char *) "solid_bulk",DOUBLE,6);

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

  fnum_rat = 1.0/particle->species[solid_species].specwt;

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
  solid_array[index][2] = Tp0; // temperature
  solid_array[index][3] = 1.0; // volume fraction which is solid

  solid_force[index][0] = 0.0; // Fx
  solid_force[index][1] = 0.0; // Fy
  solid_force[index][2] = 0.0; // Fz
  solid_force[index][3] = 0.0; // Q

  solid_bulk[index][0] = 0.0; // Ugx
  solid_bulk[index][1] = 0.0; // Ugy
  solid_bulk[index][2] = 0.0; // Ugz
  solid_bulk[index][3] = 0.0; // Tg
  solid_bulk[index][3] = 0.0; // Pg
}

/* ---------------------------------------------------------------------- */

void FixSolid::end_of_step()
{
  if (!particle->sorted) particle->sort();
  reallocate(); // for outputting average force in each grid

  if (move_type == NOMOVE) reset_velocities(1);

  // force model type

  if (force_type == GREEN || force_type == BURT) update_Fq_fm();
  else if (force_type == LOTH || force_type == SINGH) update_Fq_emp();

  if (update->ntimestep % nevery) return;

  // update particles velocitis (and position if langevin)

  ndelete = 0;
  if (move_type == EULER) update_particle();
  //else if (move_type == LANGEVIN) move_langevin();
  else if (move_type == NOMOVE) reset_velocities(0); // zero out velocities

  // delete solid particles with no mass

  if (ndelete) particle->compress_reactions(ndelete,dellist);
}

/* ---------------------------------------------------------------------- */

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
  double Rp,phi_s,phi_l,Tp,csp;
  double mp, mp_s, mp_l;
  double Tp_new,Rp_new;
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
        Tp = solid_array[ip][2];
        phi_s = solid_array[ip][3];
        phi_l = 1.0-phi_s;
        csp = phi_s*cp_solid + phi_l*cp_liquid;

        double Vp = 4./3.*3.14159*pow(Rp,3.0);
        mp = Vp*(phi_s*rho_solid + phi_l*rho_liquid);

        // velocities
        for (int d = 0; d < dim; d++) {
          particles[ip].v[d] = particles[ip].v[d] +
            solid_force[ip][d]*update->dt*nevery;
        }
        
        // joules
        double Ein = solid_force[ip][3];

        // only assume particle can sublimate

        if (phase_flag) {

          // First, determine saturation pressure
          // grab pressure in cell
          double p = 0; // pressure
          if (pwhich == AVERAGE) p = solid_bulk[ip][4];

          // use updated temp
          double T_degC = Tp - 273.15;

          // calculate saturation pressure of ice and water based on particle temp
          // ref: Huang 2018 
          // Simple Accurate Formula. for Calc. Saturation Vapor Pressure

          // if T_degC < 0, cannot be liquid
          // if T_degC > 0, cannot be solid
          double psat;
          if (T_degC > 0)
            psat = exp(34.494 - (4924.99)/(T_degC+273.1))/pow(T_degC+105,1.57);
          else
            psat = exp(43.494 - (6545.8)/(T_degC+278))/pow(T_degC+868,2.0);

          // Second determine mass lost (if any) according to Hertz-Knudsen
          double m_h2o = 2.988e-26; // [kg]
          double M_h2o = 0.018; // [kg/mol]
          double R_u = 8.314; // [J/(mol K)]
          // kg per second per area
          double flux = (psat - p)*sqrt(M_h2o/(2.0*3.14159*R_u*Tp));
          if (flux < 0) flux = 0.0;

          // determine mass lost and update radius based on
          // ... Hertz Knudsen Equation
          // ref: Kossacki and Leliwa-Kopystynski (2014) Icarus
          // if sublimation rate is "slow", then can use current surface area
          double area = 3.14159*4.0*Rp*Rp;
          double del_mp = flux * area * update->dt * nevery;
          mp -= del_mp;

          // new reduced particle size
          if (mp > 0) {
            Rp_new = pow( (mp/rho_solid)*0.75/3.14159, 1.0/3.0);

            double Eflux = del_mp*dHsub;
            double Enet = Ein - Eflux;

            // update particle temp based on net heat flux
            // based on sign of qnet, gas may or may not provide enough energy
            // .. to compensate energy lost due to phase change
            Tp_new = Tp + Enet/csp/(mp+del_mp)*update->dt*nevery; // use old mass

          // particle is gone
          } else Rp_new = Tp_new = mp = 0.0;
            

        } else {
          Tp_new =  Tp + Ein/csp/mp*update->dt*nevery;
          Rp_new = Rp;
        }

        // stpre new particle temp
        solid_array[ip][1] = Rp_new;
        solid_array[ip][2] = Tp_new;
        //solid_array[ip][3] = phi_s;

        if (Rp_new > 0.0) {

          // update per-grid forces for outputting
          array_grid[icell][0] += Rp_new;
          array_grid[icell][1] += mp;
          array_grid[icell][2] += Tp_new;
          array_grid[icell][3] += solid_force[ip][0];
          array_grid[icell][4] += solid_force[ip][1];
          array_grid[icell][5] += solid_force[ip][2];
          array_grid[icell][6] += solid_force[ip][3];
          nsolid++;
        }

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

/* --------------------------------------------------------------------- */

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
   read list of species defined in species file
   store info in filespecies and nfile
   only invoked by proc 0
------------------------------------------------------------------------- */

void FixSolid::read_solid()
{
  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have NWORDS

  int NWORDS = 8;
  char **words = new char*[NWORDS];
  char line[MAXLINE],copy[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    // check line length
    int nwords = wordcount(line,words);
    if (nwords != NWORDS)
      error->one(FLERR,"Incorrect line format in solid species file");

    char isp[16];
    strcpy(isp,words[0]);
    solid_species = particle->find_species(isp);
    if (solid_species < 0) error->all(FLERR,"Fix solid drag species does not exist");

    Rp0 = atof(words[1]);
    rho_solid = atof(words[2]);
    rho_liquid = atof(words[3]);
    if (dim == 2) mp0 = 3.141598*Rp0*Rp0*rho_solid;
    else mp0 = 4./3.*3.14159*pow(Rp0,3.0)*rho_solid;
    cp_solid = atof(words[4]);
    cp_liquid = atof(words[5]);
    alpha = atof(words[6]);
    eps = 1.0 - alpha;
    dHsub = atof(words[7]);
  }

  delete [] words;

  fclose(fp);
}

/* ----------------------------------------------------------------------
   count whitespace-delimited words in line
   line will be modified, since strtok() inserts NULLs
   if words is non-NULL, store ptr to each word
------------------------------------------------------------------------- */

int FixSolid::wordcount(char *line, char **words)
{
  int nwords = 0;
  char *word = strtok(line," \t\n");

  while (word) {
    if (words) words[nwords] = word;
    nwords++;
    word = strtok(NULL," \t\n");
  }

  return nwords;
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
