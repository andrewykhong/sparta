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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "surf_collide_mixed.h"
#include "surf.h"
#include "surf_react.h"
#include "input.h"
#include "variable.h"
#include "particle.h"
#include "domain.h"
#include "update.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "math_const.h"
#include "math_extra.h"
#include "error.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NUMERIC,CUSTOM,VARIABLE,VAREQUAL,VARSURF};   // surf_collide classes
enum{DIFFUSE,SPECULAR,ADIABATIC,CLL};
enum{NONE,DISCRETE,SMOOTH};

// straightforward to add other surface models:
// add surface tpye into loop and into above enum
/* ---------------------------------------------------------------------- */

SurfCollideMixed::SurfCollideMixed(SPARTA *sparta, int narg, char **arg) :
  SurfCollide(sparta, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal surf collide mixed command");

  int nspecies = particle->nspecies;
  memory->create(stype,nspecies,"surf_mixed:stype");
  for (int i = 0; i < nspecies; i++) stype[i] = -1;

  dflag = rflag = sflag = 0;
  tflag = rflag = 0;
  noslip_flag = 0;

  // find all species

  int iarg = 2;
  int lline;
  int isp;
  while (iarg < narg) {
    if(strcmp(arg[iarg],"diffuse") == 0) {
      dflag = 1;
      iarg++;
      while (1) {
        isp = particle->find_species(arg[iarg]);
        if (isp < 0) break;
        stype[isp] = DIFFUSE;
        iarg++;
      }
      parse_tsurf(arg[iarg++]);
      acc = atof(arg[iarg++]);
      if (acc < 0.0 || acc > 1.0)
        error->all(FLERR,"Illegal surf_collide mixed (diffuse) command");
    } else if (strcmp(arg[iarg],"specular") == 0) {
      sflag = 1;
      iarg++;
      while (1) {
        isp = particle->find_species(arg[iarg]);
        if (isp < 0) break;
        stype[isp] = SPECULAR;
        iarg++;
      }
      if (iarg < narg) { // need to check if any args left
        if (strcmp(arg[iarg],"noslip") == 0) {
          noslip_flag = 1;
          iarg ++;
        }
      }
    } else if (strcmp(arg[iarg],"adiabtic") == 0) {
      aflag = 1;
      iarg++;
      while (1) {
        isp = particle->find_species(arg[iarg]);
        if (isp < 0) break;
        stype[isp] = ADIABATIC;
        iarg++;
      }
    } else if (strcmp(arg[iarg],"translate") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal surf_collide mixed command");
      tflag = 1;
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+7 > narg)
        error->all(FLERR,"Illegal surf_collide mixed command");
      rflag = 1;
      px = atof(arg[iarg+1]);
      py = atof(arg[iarg+2]);
      pz = atof(arg[iarg+3]);
      wx = atof(arg[iarg+4]);
      wy = atof(arg[iarg+5]);
      wz = atof(arg[iarg+6]);

      if (domain->dimension == 2) {
        if (pz != 0.0)
          error->all(FLERR,"Surf_collide diffuse rotation invalid for 2d");
        if (!domain->axisymmetric && (wx != 0.0 || wy != 0.0))
          error->all(FLERR,"Surf_collide diffuse rotation invalid for 2d");
        if (domain->axisymmetric && (wy != 0.0 || wz != 0.0))
          error->all(FLERR,
                     "Surf_collide diffuse rotation invalid for 2d axisymmetric");
      }
      iarg += 7;
    } else error->all(FLERR,"Illegal surf_collide mixed command");
  }

  // check all species have a surface collision model assigned

  for (int i = 0; i < nspecies; i++)
    if (stype[i] < 0) error->all(FLERR,"No surface model assigned to a species");

  if (tflag && rflag) error->all(FLERR,"Illegal surf_collide mixed command");
  if (tflag || rflag) trflag = 1;
  else trflag = 0;

  // for updating custom for diffuse walls
  vstream[0] = vstream[1] = vstream[2] = 0.0;

  // initialize RNG

  random = new RanKnuth(update->ranmaster->uniform());
  double seed = update->ranmaster->uniform();
  random->reset(seed,comm->me,100);
}

/* ---------------------------------------------------------------------- */

SurfCollideMixed::~SurfCollideMixed()
{
  if (copy) return;

  memory->destroy(stype);
  delete random;
}

/* ---------------------------------------------------------------------- */

void SurfCollideMixed::init()
{
  SurfCollide::init();
  if (dflag) check_tsurf();
}

/* ----------------------------------------------------------------------
   particle collision with surface with optional chemistry
   ip = particle with current x = collision pt, current v = incident v
   isurf = index of surface element
   norm = surface normal unit vector
   isr = index of reaction model if >= 0, -1 for no chemistry
   ip = reset to NULL if destroyed by chemistry
   return jp = new particle if created by chemistry
   return reaction = index of reaction (1 to N) that took place, 0 = no reaction
   resets particle(s) to post-collision outward velocity

   can possible simplify this by directly calling collide from each type
------------------------------------------------------------------------- */

Particle::OnePart *SurfCollideMixed::
collide(Particle::OnePart *&ip, double &,
        int isurf, double *norm, int isr, int &reaction)
{
  nsingle++;

  // if surface chemistry defined, attempt reaction
  // reaction = 1 to N for which reaction took place, 0 for none
  // velreset = 1 if reaction reset post-collision velocity, else 0

  Particle::OnePart iorig;
  Particle::OnePart *jp = NULL;
  reaction = 0;
  int velreset = 0;

  if (isr >= 0) {
    if (modify->n_surf_react) memcpy(&iorig,ip,sizeof(Particle::OnePart));
    reaction = surf->sr[isr]->react(ip,isurf,norm,jp,velreset);
    if (reaction) surf->nreact_one++;
  }

  int isp;
  if (ip) {
    isp = ip->ispecies;
    // diffuse 
    if (stype[isp] == DIFFUSE) {
      // set temperature of isurf if VARSURF or CUSTOM
      if (persurf_temperature) {
        tsurf = t_persurf[isurf];
        if (tsurf <= 0.0) error->one(FLERR,"Surf_collide tsurf <= 0.0");
      }
      if (!velreset) diffuse(ip, norm);
      if (modify->n_update_custom) {
        int i = ip - particle->particles;
        modify->update_custom(i,tsurf,tsurf,tsurf,vstream);
      }
    // specular
    } else if (stype[isp] == SPECULAR && !velreset) {
      if (noslip_flag)  MathExtra::negate3(ip->v);
      else  MathExtra::reflect3(ip->v,norm);
    //adiabatic
    } else if (stype[isp] == ADIABATIC && !velreset) {
      scatter_isotropic(ip,norm);
    } else error->all(FLERR,"Could not find surface collision model");
  }

  if (jp) {
    isp = jp->ispecies;
    // diffuse 
    if (stype[isp] == DIFFUSE) {
      // set temperature of isurf if VARSURF or CUSTOM
      if (persurf_temperature) {
        tsurf = t_persurf[isurf];
        if (tsurf <= 0.0) error->one(FLERR,"Surf_collide tsurf <= 0.0");
      }
      if (!velreset) diffuse(jp, norm);
      if (modify->n_update_custom) {
        int i = jp - particle->particles;
        modify->update_custom(i,tsurf,tsurf,tsurf,vstream);
      }
    // specular
    } else if (stype[isp] == SPECULAR && !velreset) {
      if (noslip_flag)  MathExtra::negate3(jp->v);
      else  MathExtra::reflect3(jp->v,norm);
    //adiabatic
    } else if (stype[isp] == ADIABATIC && !velreset) {
      scatter_isotropic(jp,norm);
    } else error->all(FLERR,"Could not find surface collision model");
  }

  // call any fixes with a surf_react() method
  // they may reset j to -1, e.g. fix ambipolar
  //   in which case newly created j is deleted

  if (reaction && modify->n_surf_react) {
    int i = -1;
    if (ip) i = ip - particle->particles;
    int j = -1;
    if (jp) j = jp - particle->particles;
    modify->surf_react(&iorig,i,j);
    if (jp && j < 0) {
      jp = NULL;
      particle->nlocal--;
    }
  }

  return jp;
}

/* ----------------------------------------------------------------------
   taken from SurfCollideDiffuse (refer there)
------------------------------------------------------------------------- */

void SurfCollideMixed::diffuse(Particle::OnePart *p, double *norm)
{
  if (random->uniform() > acc) {
    MathExtra::reflect3(p->v,norm);
  } else {
    double tangent1[3],tangent2[3];
    Particle::Species *species = particle->species;
    int ispecies = p->ispecies;

    double vrm = sqrt(2.0*update->boltz * tsurf / species[ispecies].mass);
    double vperp = vrm * sqrt(-log(random->uniform()));

    double theta = MY_2PI * random->uniform();
    double vtangent = vrm * sqrt(-log(random->uniform()));
    double vtan1 = vtangent * sin(theta);
    double vtan2 = vtangent * cos(theta);

    double *v = p->v;
    double dot = MathExtra::dot3(v,norm);

    double beta_un,normalized_distbn_fn;

    tangent1[0] = v[0] - dot*norm[0];
    tangent1[1] = v[1] - dot*norm[1];
    tangent1[2] = v[2] - dot*norm[2];

    if (MathExtra::lensq3(tangent1) == 0.0) {
      tangent2[0] = random->uniform();
      tangent2[1] = random->uniform();
      tangent2[2] = random->uniform();
      MathExtra::cross3(norm,tangent2,tangent1);
    }

    MathExtra::norm3(tangent1);
    MathExtra::cross3(norm,tangent1,tangent2);

    if (trflag) {
      double vxdelta,vydelta,vzdelta;
      if (tflag) {
        vxdelta = vx; vydelta = vy; vzdelta = vz;
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];

        if (fabs(dot) > 0.001) {
          dot /= vrm;
          do {
            do {
              beta_un = (6.0*random->uniform() - 3.0);
            } while (beta_un + dot < 0.0);
            normalized_distbn_fn = 2.0 * (beta_un + dot) /
              (dot + sqrt(dot*dot + 2.0)) *
              exp(0.5 + (0.5*dot)*(dot-sqrt(dot*dot + 2.0)) -
                  beta_un*beta_un);
          } while (normalized_distbn_fn < random->uniform());
          vperp = beta_un*vrm;
        }

      } else {
        double *x = p->x;
        vxdelta = wy*(x[2]-pz) - wz*(x[1]-py);
        vydelta = wz*(x[0]-px) - wx*(x[2]-pz);
        vzdelta = wx*(x[1]-py) - wy*(x[0]-px);
        double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];
        vxdelta -= dot*norm[0];
        vydelta -= dot*norm[1];
        vzdelta -= dot*norm[2];
      }

      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0] + vxdelta;
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1] + vydelta;
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2] + vzdelta;

    } else {
      v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
      v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
      v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];
    }

    p->erot = particle->erot(ispecies,tsurf,random);
    p->evib = particle->evib(ispecies,tsurf,random);
  }
}

/* ----------------------------------------------------------------------
   taken from SurfCollideAdiabatic (refer there)
------------------------------------------------------------------------- */

void SurfCollideMixed::scatter_isotropic(Particle::OnePart *p, double *norm)
{
  double *v = p->v;
  double dot = MathExtra::dot3(v,norm);

  // tangent1/2 = surface tangential unit vectors

  double tangent1[3], tangent2[3];
  tangent1[0] = v[0] - dot*norm[0];
  tangent1[1] = v[1] - dot*norm[1];
  tangent1[2] = v[2] - dot*norm[2];

  // if mag(tangent1) == 0 mean normal collision, in that case choose
  // a random tangential vector
  if (MathExtra::lensq3(tangent1) == 0.0) {
    tangent2[0] = random->uniform();
    tangent2[1] = random->uniform();
    tangent2[2] = random->uniform();
    MathExtra::cross3(norm,tangent2,tangent1);
  }

  // normalize tangent1
  MathExtra::norm3(tangent1);
  // compute tangent2 as norm x tangent1, so that tangent1 and tangent2 are
  // orthogonal
  MathExtra::cross3(norm,tangent1,tangent2);

  // isotropic scattering
  // vmag = magnitude of incidient particle velocity vector
  // vperp = velocity component perpendicular to surface along norm (cy)
  // vtan1/2 = 2 remaining velocity components tangential to surface

  double vmag = MathExtra::len3(v);

  double theta = MY_2PI*random->uniform();
  double f_phi = random->uniform();
  double sqrt_f_phi = sqrt(f_phi);

  double vperp = vmag * sqrt(1.0 - f_phi);
  double vtan1 = vmag * sqrt_f_phi * sin(theta);
  double vtan2 = vmag * sqrt_f_phi * cos(theta);

  // update particle velocity
  v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0];
  v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1];
  v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2];

  // p->erot and p->evib stay identical
}

/* ----------------------------------------------------------------------
   combine wrappers from SurfCollideDiffuse and SurfCollideSpecular
------------------------------------------------------------------------- */

void SurfCollideMixed::wrapper(Particle::OnePart *p, double *norm,
                                 int *flags, double *coeffs)
{
  int isp = p->ispecies;
  if (stype[isp] == DIFFUSE) {
    if (coeffs) {
      tsurf = coeffs[0];
      acc = coeffs[1];
    }
    diffuse(p,norm);
  } else if (stype[isp] == SPECULAR) {
    if (flags) noslip_flag = flags[0];
    MathExtra::reflect3(p->v,norm);
  } else if (stype[isp] == ADIABATIC) {
    scatter_isotropic(p,norm);
  }
}

/* ----------------------------------------------------------------------
   return flags and coeffs for this SurfCollide instance to caller
------------------------------------------------------------------------- */

void SurfCollideMixed::flags_and_coeffs(int *flags, double *coeffs)
{
  if (tmode != NUMERIC)
    error->all(FLERR,"Surf_collide diffuse with non-numeric Tsurf "
               "does not support external caller");

  coeffs[0] = tsurf;
  coeffs[1] = acc;
  flags[0] = noslip_flag;
}
