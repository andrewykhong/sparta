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
#include "math_extra.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{INT,DOUBLE};                      // several files
enum{EULER,LANGEVIN};                  // type of solid particle move
enum{NOFORCE,GREEN,LOTH,SINGH};            // type of solid particle force
enum{SPHERE,DISC,CYLINDER,CUSTOM};

// for compute_solid_grid

enum{SIZE,MASS,TEMP,FORCEX,FORCEY,FORCEZ,HEAT};
enum{PZERO,AVERAGE};

#define DELTADELETE 1024
#define MAXLINE 1024
#define SMALLANGLE 1E-10

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

  maxdelete = 0.0;
  dellist = NULL; // particles to delete

  // for storing per-grid solid-particle properties
  per_grid_flag = 1;
  per_grid_freq = 1;
  // average force (3), heat flux, and particle params (3)
  size_per_grid_cols = 7;
  nglocal = 0;
  array_grid = NULL;
  Sn = NULL;
  //Fq_grid = NULL;

  if (narg < 3) error->all(FLERR,"Not enough arguments for fix solid command");

  // read in initial solid particle properties

  fp = fopen(arg[2],"r");
  if (!fp) error->one(FLERR,"Solid species file not found");
  read_solid();

  // initial solid particle temperature and specific heat

  Tp0 = atof(arg[3]);

  // optional args

  ifix = -1;
  nevery = 1; // frequency to update particle vel + pos
  phase_flag = 0; // reduce size of particle radii over time due to heat
  move_type = EULER;
  force_type = NOFORCE;
  pwhich = PZERO;
  uxp0 = uyp0 = uzp0 = 0.0;
  conserve_flag = 0;
  reset_flag = 0;
  merge_flag = 0;

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
    } else if (strcmp(arg[iarg],"reset") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Invalid fix solid command");
      reset_flag = 1;
      uxp0 = atof(arg[iarg+1]);
      uyp0 = atof(arg[iarg+2]);
      uzp0 = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"green") == 0) {
      force_type = GREEN;
      iarg++;
    } else if (strcmp(arg[iarg],"loth") == 0) {
      force_type = LOTH;
      iarg++;
    } else if (strcmp(arg[iarg],"conserve") == 0) { 
      conserve_flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"surf") == 0) { 
      FILE *fsurf = fopen(arg[iarg+1],"r");
      if (!fsurf) error->one(FLERR,"Particle surface file not found");
      int type;
      if (strcmp(arg[iarg+2],"stl") == 0) read_surf_stl(fsurf);
      else if (strcmp(arg[iarg+2],"shape") == 0) read_surf_shape(fsurf);
      else error->all(FLERR,"Illegal option for custom surf");
      iarg += 3;
    } else error->all(FLERR,"Invalid fix solid command");
  }

  // check surface specificed if using custom
  if (shape == CUSTOM && !Sn) error->one(FLERR,"Fix solid: No surface provided");

  // check if custom particle parameters exist

  index_solid_params = particle->find_custom((char *) "solid_params");
  index_solid_force = particle->find_custom((char *) "solid_force");
  index_solid_bulk = particle->find_custom((char *) "solid_bulk");

  if (index_solid_params >= 0 || index_solid_force >= 0 || index_solid_bulk >= 0)
    error->all(FLERR,"Fix solid_drag custom attribute already exists");

  index_solid_params = particle->add_custom((char *) "solid_params",DOUBLE,8);
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
  memory->destroy(Sn);
  //memory->destroy(Fq_grid);

  //delete [] argindex;
  //delete [] value2index;
  //delete [] post_process;
  //for (int i = 0; i < 2; i++) delete [] ids[i];

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

  // molecules can interact with particles in three ways:
  // 1) specular (eps)
  // 2) diffuse (1-eps)(alpha)
  // 3) adiabatice (1-eps)(1-alpha)
  // pre-compute prefactors for force and heat flux delta functions

  if (force_type == GREEN) {
    c_diff = (1.0-eps)*alpha;
    c_spec = eps;
    c_adia = (1.0-eps)*(1.0-alpha);

    // normalize so sum is 1
    double sum = c_diff + c_spec + c_adia;
    c_diff /= sum;
    c_spec /= sum;
    c_adia /= sum;
  }

  fnum_rat = 1.0/particle->species[solid_species].specwt;
  ofnum_rat = particle->species[solid_species].specwt;

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
  solid_array[index][4] = Lp0; // length (for cylinder)
  solid_array[index][5] = theta0; // polar angle (to find normal)
  solid_array[index][6] = phi0; // azithumal angle (to find normal)

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

  // For debugging
  if (reset_flag) reset_velocities(1);

  // force model type

  if (force_type == GREEN) update_Fq_fm();
  //else if (force_type == LOTH) update_Fq_emp();

  if (update->ntimestep % nevery) return;

  // update particles velocitis (and position if langevin)

  ndelete = 0;
  if (move_type == EULER) update_particle();

  // For debugging
  if (reset_flag) reset_velocities(0);

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

  int i,icell,ip,ispecies,np,isp;

  // particle related

  double vold[3],cx,cy,cz,csq;
  double Rp,Lp,phi_s,phi_l,Tp,csp;
  double mp, mp_s, mp_l;
  double Tp_new,Rp_new,Lp_new,mp_new;
  double Q,mp_loss,rho,new_vol;

  double mass_g, Eth_g, Erot_g, U_g[3], total_rdof;
  double dj[3], dke, dT;

  int nsolid;

  // FOR DEBUG
  double dmom_part[3], dE_part;
  double dmom_gas[3], dE_gas;

  for (icell = 0; icell < nglocal; icell++) {

    // reset array_grid (for output)
    for (i = 0; i < size_per_grid_cols; i++)
      array_grid[icell][i] = 0.0;

    np = cinfo[icell].count;
    if (np <= 0) continue;
    nsolid = 0;

    // manually check conservation

    //dmom_gas[0] = dmom_gas[1] = dmom_gas[2] = dE_gas = 0.0;
    //dmom_part[0] = dmom_part[1] = dmom_part[2] = dE_part = 0.0;

    // update only solid particles

    mass_g = Eth_g = total_rdof = Erot_g = 0.0;
    dke = dT = 0.0;
    U_g[0] = U_g[1] = U_g[2] = 0.0;
    dj[0] = dj[1] = dj[2] = 0.0;

    ip = cinfo[icell].first;
    while (ip >= 0) {
      isp = particles[ip].ispecies;
      if (isp == solid_species) {
        Rp = solid_array[ip][1];
        Lp = solid_array[ip][4];
        Tp = solid_array[ip][2];
        phi_s = solid_array[ip][3];
        phi_l = 1.0-phi_s;
        csp = phi_s*cp_solid + phi_l*cp_liquid;

        double Vp;
        if (shape == SPHERE) Vp = 4./3.*MY_PI*pow(Rp,3.0);
        else if (shape == DISC) Vp = 2.0*MY_PI*Rp*Rp;
        else if (shape == CYLINDER) Vp = (MY_PI*Rp*Rp)*Lp;

        mp = Vp*(phi_s*rho_solid + phi_l*rho_liquid);

        // force on each real solid particle (do not need to update with fnum)
        double Fd[3];
        Fd[0] = solid_force[ip][0];
        Fd[1] = solid_force[ip][1];
        Fd[2] = solid_force[ip][2];

        // first solve for new temperature
        double Ein = solid_force[ip][3]; // joules

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
          double area;
          if (shape == SPHERE || shape == DISC) area = 4.0*MY_PI*Rp*Rp;
          else if (shape == DISC) area = 2.0*MY_PI*Rp*Rp;
          else if (shape == CYLINDER) area = 2.0*MY_PI*Rp*(Rp+Lp);

          double del_mp = flux * area * update->dt * nevery;
          mp_new = mp - del_mp;

          // energy loss by mass loss
          double Edm = del_mp*dHsub;
          Tp_new = Tp + (Ein*update->dt-Edm)/csp/mp_new*nevery;

          // new reduced particle size
          double Vscale = mp_new/mp; // assume constant density
          if (mp_new > 0) {
            if (shape == SPHERE) Rp_new = Rp*pow(Vscale,1.0/3.0);
            else if (shape == DISC) Rp_new = Rp*pow(Vscale,0.5);
            else if (shape == CYLINDER) {
              Rp_new = Rp*pow(Vscale,1.0/3.0);
              Lp_new = Lp*pow(Vscale,1.0/3.0);
            }
          } else Rp_new = Lp_new = Tp_new = mp_new = 0.0; // particle is gone  
        } else {
          // joules -> kg x m^2 / s^2
          // specific heat: J/(kg K) -> kg x m^2 / (s^2 x kg x K)
          // ... -> m^2 / (s^2 x K)
          // K = K + (kg x m^2/s^3) / (m^2 / (s^2 x K)) / kg * s
          // K = K + K x (1/s^3) / (1/(s^2)) * s
          // K = K + K
          Tp_new = Tp + Ein/csp/mp*update->dt*nevery;
          mp_new = mp;
          Rp_new = Rp;
          Lp_new = Lp;
        } // end phase change check

        /*if (fabs(Ein) > 0) {
          printf("Ein: %4.3e\n", Ein);
          printf("mp: %4.3e -> %4.3e\n", mp, mp_new);
          printf("Tp: %4.3e -> %4.3e\n", Tp, Tp_new);
          printf("Rp: %4.3e -> %4.3e\n", Rp, Rp_new);
          error->one(FLERR,"ck");
        }*/

        // stpre new particle temp
        solid_array[ip][1] = Rp_new;
        solid_array[ip][2] = Tp_new;
        solid_array[ip][4] = Lp_new;

        // record change in momentum and kinetic energy
        for (int d = 0; d < dim; d++) {
          // old
          dj[d] -= mp*particles[ip].v[d];
          dke -= 0.5*mp*(particles[ip].v[d]*particles[ip].v[d]);
          // update
          particles[ip].v[d] += Fd[d]/mp*update->dt*nevery;
          // new
          dj[d] += mp_new*particles[ip].v[d];
          dke += 0.5*mp_new*(particles[ip].v[d]*particles[ip].v[d]);
        }

        // record change in temperature
        dT += (csp*mp_new*Tp_new-csp*mp*Tp);

        //printf("mp: %4.3e; mp_new: %4.3e\n", mp, mp_new);
        //printf("Ein: %4.3e, Told: %4.3e; Tnew: %4.3e\n", Ein, Tp, Tp_new);
        //printf("dj: %4.3e, %4.3e, %4.3e; dke: %4.3e; dT: %4.3e\n",
        //  dj[0], dj[1], dj[2], dke, dT);
        //error->one(FLERR,"ck");

        // update per-grid forces for outputting
        if (Rp_new > 0.0) {
          array_grid[icell][0] += Rp_new;
          array_grid[icell][1] += mp_new;
          array_grid[icell][2] += Tp_new;
          array_grid[icell][3] += solid_force[ip][0];
          array_grid[icell][4] += solid_force[ip][1];
          array_grid[icell][5] += solid_force[ip][2];
          array_grid[icell][6] += solid_force[ip][3];
          nsolid++;
        // delete particle
        } else {
          if (ndelete == maxdelete) {
            maxdelete += DELTADELETE;
            memory->grow(dellist,maxdelete,"fix/solid:dellist");
          }
          dellist[ndelete++] = ip;
        }
      
      // record total gas mass and momentum
      } else if (conserve_flag) {
        double imass = species[isp].mass;
        mass_g += imass;
        total_rdof += species[isp].rotdof;
        Erot_g += particles[ip].erot;
        for (int d = 0; d < dim; d++) {
          U_g[d] += imass*particles[ip].v[d];
          Eth_g  += 0.5*imass*(particles[ip].v[d]*particles[ip].v[d]);
        }
      } // end species check


      // reset forces and heat fluxes
      for (i = 0; i < 4; i++) solid_force[ip][i] = 0.0;
      // reset pressures and temperatures
      for (i = 0; i < 5; i++) solid_bulk[ip][i] = 0.0;
      ip = next[ip];

    } // end particle

    // solid->gas to conserve momentum and energy
    if (conserve_flag) {
      //printf("conserve\n");
      double dU_g[3]; // new - old
      for (int d = 0; d < dim; d++) {
        dU_g[d] = -dj[d]*ofnum_rat/mass_g;
        U_g[d] /= mass_g;
      }
      dke *= ofnum_rat;
      dT *= ofnum_rat;

      double dE = -0.5*mass_g*(dU_g[0]*dU_g[0]+dU_g[1]*dU_g[1]+dU_g[2]*dU_g[2])-(dke+dT);

      Eth_g -= 0.5*mass_g*(U_g[0]*U_g[0]+U_g[1]*U_g[1]+U_g[2]*U_g[2]);
      double dEth  = dE * (3.0*(np-nsolid))/(3.0*(np-nsolid)+total_rdof);
      double dErot = dE - dEth;
      double phi = sqrt( (Eth_g+dE)/Eth_g );

      //if (dE<0) error->one(FLERR,"negative energy change");

      // update velocities and internal modes
      ip = cinfo[icell].first;
      while (ip >= 0) {
        isp = particles[ip].ispecies;
        if (isp != solid_species) {
          for (int d = 0; d < dim; d++) {
            particles[ip].v[d] =
              (particles[ip].v[d]-U_g[d])*phi + (U_g[d]+dU_g[d]);
          }
          //particles[ip].erot *= (dErot+Erot_g)/Erot_g;
        }
        ip = next[ip];
      }
    } // end conserve

    // update per-grid forces for outputting
    if (nsolid > 0) {
      for (i = 0; i < size_per_grid_cols; i++)
        array_grid[icell][i] /= nsolid;
    } else { 
      for (i = 0; i < size_per_grid_cols; i++)
        array_grid[icell][i] = 0.0;
    }

  } // end cells

  //printf("end - update\n");

  // reset number of samples
  nsample = 0;
}

/* --------------------------------------------------------------------- */

void FixSolid::reset_velocities(int reset)
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
        if (reset) {
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
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void FixSolid::read_surf_stl(FILE *fin)
{
  
  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have NWORDS

  int MAXWORDS = 5;
  char **words = new char*[MAXWORDS];
  char line[MAXLINE];

  int npts, ntris;
  int ipt = 0;
  int itri = 0;
  double **tmp_pts;
  int **tmp_tris;

  // flag for header
  int header = 0;

  // read file
  tmp_pts = NULL;
  tmp_tris = NULL;
  while (fgets(line,MAXLINE,fin)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    // check line length
    int nwords = wordcount(line,words);

    //printf("%s %s\n", words[0], words[1]);

    // read number of points
    if (header == 0) {
      if (strcmp(words[1],"points") != 0) continue;
        //error->all(FLERR,"Fix solid: First line expected to be number of points");
      npts = atoi(words[0]);

      memory->destroy(tmp_pts);
      memory->create(tmp_pts,npts,3,"fix/solid:tmp_pts");
      header++;
    // read number of triangles
    } else if (header == 1){
      if (strcmp(words[1],"triangles") != 0)
        error->all(FLERR,"Fix solid: Second line expected to be number of triangles");
      ntris = atoi(words[0]);

      memory->destroy(tmp_tris);
      memory->create(tmp_tris,ntris,3,"fix/solid:tmp_tris");
      header++;
    // skip section headers
    } else if (strcmp(words[0],"Points") == 0 || strcmp(words[0],"Triangles") == 0) {
      continue;
    } else if (ipt < npts) {
      int pt_id = atoi(words[0]) - 1; // ids start with 1
      tmp_pts[pt_id][0] = atof(words[1]);
      tmp_pts[pt_id][1] = atof(words[2]);
      tmp_pts[pt_id][2] = atof(words[3]);
      ipt++;
    } else if (itri < ntris) {
      // triangle ids not needed
      tmp_tris[itri][0] = atoi(words[1])-1;
      tmp_tris[itri][1] = atoi(words[2])-1;
      tmp_tris[itri][2] = atoi(words[3])-1;
      itri++;
    }
  }

  //printf("read file\n");

  // create surface list
  nsurfs = ntris;
  Sn = NULL;
  memory->destroy(Sn);
  memory->create(Sn,nsurfs,3,"fix/solid:Sn");

  int ind;
  double p1[3], p2[3], p3[3];
  double delta12[3], delta13[3], delta23[3], norm[3], new_norm[3];
  for (int i = 0; i < ntris; i++) {
    // grab vertices
    ind = tmp_tris[i][0];
    p1[0] = tmp_pts[ind][0];
    p1[1] = tmp_pts[ind][1];
    p1[2] = tmp_pts[ind][2];

    ind = tmp_tris[i][1];
    p2[0] = tmp_pts[ind][0];
    p2[1] = tmp_pts[ind][1];
    p2[2] = tmp_pts[ind][2];

    ind = tmp_tris[i][2];
    p3[0] = tmp_pts[ind][0];
    p3[1] = tmp_pts[ind][1];
    p3[2] = tmp_pts[ind][2];

    //printf("tri: %i: %i, %i, %i\n", i, tmp_tris[i][0], tmp_tris[i][1], tmp_tris[i][2]);
    //printf("p1: %i: %4.3e, %4.3e, %4.3e\n", i, p1[0], p1[1], p1[2]);
    //printf("p2: %i: %4.3e, %4.3e, %4.3e\n", i, p2[0], p2[1], p2[2]);
    //printf("p3: %i: %4.3e, %4.3e, %4.3e\n", i, p3[0], p3[1], p3[2]);
    //error->one(FLERR,"ck");

    MathExtra::sub3(p2,p1,delta12);
    MathExtra::sub3(p3,p1,delta13);
    MathExtra::sub3(p3,p2,delta23);
    MathExtra::cross3(delta12,delta13,norm);
    MathExtra::norm3(norm);
    if(abs(norm[0]) <= SMALLANGLE) norm[0] = 0.0;
    if(abs(norm[1]) <= SMALLANGLE) norm[1] = 0.0;
    if(abs(norm[2]) <= SMALLANGLE) norm[2] = 0.0;
    //printf("in norm: %4.3e, %4.3e, %4.3e\n", norm[0], norm[1], norm[2]);

    // convert to spherical coordinates
    double theta = atan(norm[1]/norm[0]);
    if (theta != theta) theta = 0.0; // x- and y- component are zero
    if (norm[0] < 0 && norm[1] >= 0) theta += MY_PI;
    else if (norm[0] < 0 && norm[1] < 0) theta -= MY_PI;
    else if (abs(norm[0]) == 0.0 && norm[1] > 0) theta = MY_PI*0.5;
    else if (abs(norm[0]) == 0.0 && norm[1] < 0) theta = -MY_PI*0.5;

    double phi = acos(norm[2]);

    // find area with Heron's formula
    double l1 = MathExtra::len3(delta12);
    double l2 = MathExtra::len3(delta13);
    double l3 = MathExtra::len3(delta23);
  
    double s = 0.5*(l1+l2+l3);
    double area = sqrt(s*(s-l1)*(s-l2)*(s-l3));

    // store
    Sn[i][0] = theta;
    Sn[i][1] = phi;
    Sn[i][2] = area;

    //printf("theta: %4.3e; phi: %4.3e; A: %4.3e\n", theta,phi,area);
  }

  /*double totalA = 0.0;
  for (int i = 0; i < nsurfs; i++) {
    printf("theta: %4.3e; phi: %4.3e; A: %4.3e\n", Sn[i][0],Sn[i][1],Sn[i][2]);
    totalA += Sn[i][2];
  }
  printf("totalA: %4.3e\n", totalA);*/
  //error->one(FLERR,"ck");

  // destroy

  memory->destroy(tmp_pts);
  memory->destroy(tmp_tris);
  delete [] words;

  fclose(fp);
}

/* ----------------------------------------------------------------------
   reallocate arrays if nglocal has changed
   called by init() and whenever grid changes
------------------------------------------------------------------------- */

void FixSolid::read_surf_shape(FILE *fin)
{
  
  // read file line by line
  // skip blank lines or comment lines starting with '#'
  // all other lines must have NWORDS

  int MAXWORDS = 5;
  char **words = new char*[MAXWORDS];
  char line[MAXLINE];

  int npts, ntris;
  int i = 0;
  double **tmp_pts;
  int **tmp_tris;

  // flag for header
  int header = 0;

  // read file
  tmp_pts = NULL;
  tmp_tris = NULL;
  while (fgets(line,MAXLINE,fin)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    // check line length
    int nwords = wordcount(line,words);

    //printf("%s %s\n", words[0], words[1]);

    // read number of points
    if (header == 0) {
      if (strcmp(words[1],"surfs") != 0) continue;
        //error->all(FLERR,"Fix solid: First line expected to be number of points");
      nsurfs = atoi(words[0]);

      Sn = NULL;
      memory->destroy(Sn);
      memory->create(Sn,nsurfs,3,"fix/solid:Sn");
      header++;
    } else {
      Sn[i][0] = atof(words[0]); // theta
      Sn[i][1] = atof(words[1]); // phi
      Sn[i][2] = atof(words[2]); // area
      i++;
    }
  }

  delete [] words;

  fclose(fp);
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

  int MAXWORDS = 13;
  char **words = new char*[MAXWORDS];
  char line[MAXLINE];

  while (fgets(line,MAXLINE,fp)) {
    int pre = strspn(line," \t\n\r");
    if (pre == strlen(line) || line[pre] == '#') continue;

    // check line length
    int nwords = wordcount(line,words);
    if (nwords != MAXWORDS)
      error->one(FLERR,"Incorrect line format in solid species file");

    char isp[16];
    strcpy(isp,words[0]);
    solid_species = particle->find_species(isp);
    if (solid_species < 0) error->all(FLERR,"Fix solid drag species does not exist");

    if (strcmp(words[1],"sphere") == 0) shape = SPHERE;
    else if (strcmp(words[1],"disc") == 0) shape = DISC;
    else if (strcmp(words[1],"cylinder") == 0) shape = CYLINDER;
    else if (strcmp(words[1],"custom") == 0) shape = CUSTOM;
    else error->all(FLERR,"Fix solid does not recognize shape");
    
    // geometry of particle (some may be ignored based on shape)
    Rp0     = atof(words[2]);
    Lp0     = atof(words[3]);
    theta0   = atof(words[4]);
    phi0     = atof(words[5]);

    if (theta0 < 0 || theta0 > 180.0)
      error->all(FLERR,"Fix solid: polar angle invalid");
    if (phi0 < 0 || phi0 > 360)
      error->all(FLERR,"Fix solid: azimuth angle invalid");

    // convert to radians
    theta0 *= MY_PI/180.0;
    phi0 *= MY_PI/180.0;

    // material properties
    int iarg = 6;
    rho_solid  = atof(words[iarg++]);
    rho_liquid = atof(words[iarg++]);
    cp_solid   = atof(words[iarg++]);
    cp_liquid  = atof(words[iarg++]);
    dHsub      = atof(words[iarg++]);
    alpha      = atof(words[iarg++]);
    eps        = atof(words[iarg++]);

    //printf("Rp: %4.3e; Lp: %4.3e; theta: %4.3e; phi: %4.3e\n", 
    //  Rp0, Lp0, theta0, phi0);
    //printf("rho: %4.3e; cp: %4.3e; alpha: %4.3e; eps: %4.3e; dH: %4.3e\n\n",
    //  rho_solid, cp_solid, alpha, eps, dHsub);
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
  //memory->destroy(Fq_grid);
  //memory->destroy(cell_Tp);
  nglocal = grid->nlocal;
  memory->create(array_grid,nglocal,size_per_grid_cols,"fix/solid:array_grid");
  //memory->create(Fq_grid,nglocal,4,"fix/solid:Fq_grid");
  //memory->create(cell_Tp,nglocal,2,"fix/solid:cell_Tp");

  // initialize values
  for (int i = 0; i < nglocal; i++) {
    for (int j = 0; j < size_per_grid_cols; j++) {
      array_grid[i][j] = 0.0;
      //Fq_grid[i][j] = 0.0;
    }
    //for (int j = 0; j < 2; j++)
    //  cell_Tp[i][j] = 0.0;
  }
}

