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

#ifdef SURF_COLLIDE_CLASS

SurfCollideStyle(mixed,SurfCollideMixed)

#else

#ifndef SPARTA_SURF_COLLIDE_MIXED_H
#define SPARTA_SURF_COLLIDE_MIXED_H

#include "surf_collide.h"
#include "surf.h"

namespace SPARTA_NS {

class SurfCollideMixed : public SurfCollide {
 public:
  SurfCollideMixed(class SPARTA *, int, char **);
  SurfCollideMixed(class SPARTA *sparta) : SurfCollide(sparta) {} // needed Kokkos
  virtual ~SurfCollideMixed();
  virtual void init();
  Particle::OnePart *collide(Particle::OnePart *&, double &,
                             int, double *, int, int &);
  void wrapper(Particle::OnePart *, double *, int *, double*);
  void flags_and_coeffs(int *, double *);

 protected:
  int *stype;                // surface collision type for each species
  // flags for diffuse, rotating, translating, rot-trans surface;
  int diffuse_flag, rotate_flag, translate_flag, tr_flag;

  double acc;                // surface accomodation coeff
  double vx,vy,vz;           // translational velocity of surface
  double wx,wy,wz;           // angular velocity of surface
  double px,py,pz;           // point to rotate surface around
  double dt;                 // time step
  double vwall;              // wall velocity for piston
  int tflag,rflag;           // flags for translation and rotation
  int trflag;                // 1 if either tflag or rflag is set
  int noslip_flag;           // 1 if no slip at wall
  int piston_flag;           // 1 if any species has a piston_flag assigned

  Surf::Line *lines;
  Surf::Tri *tris;

  double vstream[3];
  class RanKnuth *random;     // RNG for particle reflection

  void diffuse(Particle::OnePart *, double *);
  void specular(Particle::OnePart *, double *);
  void scatter_isotropic(Particle::OnePart *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Surf_collide diffuse rotation invalid for 2d

Specified rotation vector must be in z-direction.

E: Surf_collide diffuse variable name does not exist

Self-explanatory.

E: Surf_collide diffuse variable is invalid style

It must be an equal-style variable.

*/
