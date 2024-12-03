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

#ifdef FIX_CLASS

FixStyle(solid,FixSolid)

#else

#ifndef SPARTA_FIX_SOLID_H
#define SPARTA_FIX_SOLID_H

#include "stdio.h"
#include "fix.h"

namespace SPARTA_NS {

#define DELTAPART 128

class FixSolid : public Fix {
 public:
  FixSolid(class SPARTA *, int, char **);
  virtual ~FixSolid();
  int setmask();
  void init();
  virtual void update_custom(int, double, double, double, double *);
  void end_of_step();

  double get_particle_property(int, int);

  int solid_species;  // species for solid species
  double alpha,eps; // for defining surface-type (Lord) for Green

 protected:

  int dim; // dimension (2 or 3)
  int move_type; // type of move
  int force_type; // type of force
  int nsample; // number of samples
  int ifix; // fix for storing species-dependent
  int nglocal; // number of local grid cells

  FILE *fp; // file pointer for solid properties

  // index for custom per-particle solid propeties
  // custom array for solid params : radius, mass, specific heat, temperature
  // custom array for solid force  : Fx, Fy, Fz, heat flux

  int index_solid_params, index_solid_force, index_solid_bulk;
  int npmax; // track upper bound for number of particles
  double fnum_rat; // ratio of gas to solid particle Fnum

  // for now, all solid particles have same initial properties

  double Rp0, mp0; // initial particle radius and mass
  double rho_solid, rho_liquid; // density of ice and water
  double Tp0, cp_solid, cp_liquid; // initial temperature and specific heats
  double uxp0, uyp0, uzp0; // initial velocity of particle (for testing)
  double dHsub; // latent heat of sublimation

  // for approximating mass loss due to heating

  int phase_flag; // account for phase change
  int pwhich; // pressure
  int ifc; // id of fix or compute
  double hvap, hsolid; // specific enthalpy of vapor and solid
  double Tvap; // temperature of vapor phase

  int ndelete,maxdelete;      // # of particles removed by sublimation
  int *dellist;               // list of particle indices to delete

  // for evaluating force based on green's function

  int nspmax;
  double F1, F2; // prefactor for force defined by surface model
  double Q1; // prefactor for heat flux defined by surface model
  int *id; // solid particle id
  double *Tg, *Ug; // gas temperature, velocity

  void update_Fq_fm(); // free molecular force
  void update_Fq_emp(); // empirically derive force
  void update_particle();
  void move_langevin();
  void reset_velocities(int);
  double update_mass(double, double);

  // misc functions

  void reallocate();
  void read_solid();
  int wordcount(char *, char **);

  // random num

  class RanKnuth *random;
  double random_normal();

};

}

#endif
#endif
