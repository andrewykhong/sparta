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

 protected:

  int dim; // dimension (2 or 3)
  int move_type; // type of move
  int force_type; // type of force
  int nsample; // number of samples
  int ifix; // fix for storing species-dependent
  int nglocal; // number of local grid cells

  // index for custom per-particle solid propeties
  // custom array for solid params : radius, mass, specific heat, temperature
  // custom array for solid force  : Fx, Fy, Fz, heat flux

  int index_solid_params, index_solid_force, index_solid_bulk;
  int solid_species;  // species for solid species
  int npmax; // track upper bound for number of particles

  // for now, all solid particles have same initial properties

  double rhop0, Rp0, mp0; // initial particle density, radius, and mass
  double Tp0, in_csp; // initial temperature and specific heat
  double uxp0, uyp0, uzp0; // initial velocity of particle (for testing)

  // for approximating mass loss due to heating

  int reduce_size_flag;
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
  double alpha,eps; // for defining surface-type (Lord) for Green

  void force_burt();
  void force_green();
  void force_empirical();
  void update_particle();
  void move_langevin();
  void reset_velocities(int);
  double update_mass(double, double);

  void reallocate();

  // random num

  class RanKnuth *random;
  double random_normal();

};

}

#endif
#endif
