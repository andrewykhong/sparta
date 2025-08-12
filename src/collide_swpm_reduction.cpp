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

#include <algorithm>

#include "domain.h"
#include "math.h"
#include "math_extra.h"
#include "math_eigen.h"
#include "math_eigen_impl.h"
#include "string.h"
#include "collide.h"
#include "particle.h"
#include "mixture.h"
#include "update.h"
#include "grid.h"
#include "comm.h"
#include "react.h"
#include "modify.h"
#include "fix.h"
#include "fix_ambipolar.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "memory.h"
#include "error.h"

using namespace SPARTA_NS;

enum{ENERGY,HEAT,STRESS};   // particle reduction choices
enum{BINARY,WEIGHT}; // grouping choices
enum{NO_BALANCE,NUMBER_BALANCE,WEIGHT_BALANCE};

#define DELTADELETE 1024
#define BIG 1.0e20
#define SMALL 1.0e-50

/* ----------------------------------------------------------------------
   Merge particles using energy scheme
------------------------------------------------------------------------- */
void Collide::reduce(int istart, int iend,
                     double rho, double *V,
                     double E, double Erot)
{
  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  int np = iend-istart;
  ip = MIN(np * random->uniform() + istart,iend-1);
  jp = MIN(np * random->uniform() + istart,iend-1);
  while (ip == jp) jp = MIN(np * random->uniform() + istart,iend-1);

  ipart = &particles[plist[ip]];
  jpart = &particles[plist[jp]];

  // find direction of velocity wrt CoM frame

  double theta = 2.0 * 3.14159 * random->uniform();
  double phi = acos(1.0 - 2.0 * random->uniform());
  double uvec[3];
  uvec[0] = sin(phi) * cos(theta);
  uvec[1] = sin(phi) * sin(theta);
  uvec[2] = cos(phi);

  // set reduced particle velocities
  // convert translational temperature to kinetic energy
  double sqE = sqrt(E);
  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + sqE*uvec[d];
    jpart->v[d] = V[d] - sqE*uvec[d];
  }

  // set reduced particle rotational energies

  ipart->erot = Erot/rho;
  jpart->erot = Erot/rho;

  // set reduced particle weights

  ipart->weight = rho*0.5/update->fnum;
  jpart->weight = rho*0.5/update->fnum;

  // if there is bad weight, only create one and clone

  if (E <= 0)
    error->one(FLERR,"Negative Energy");

  // delete other particles

  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[plist[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = plist[i];
    idelete++;
  }

  return;
}

/* ---------------------------------------------------------------------
   Merge particles using heat flux scheme
------------------------------------------------------------------------ */
void Collide::reduce(int istart, int iend,
                     double rho, double *V,
                     double E, double *q,
                     double Erot)
{
  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  int np = iend-istart;
  ip = np * random->uniform() + istart;
  jp = np * random->uniform() + istart;
  while (ip == jp) jp = np * random->uniform() + istart;

  ipart = &particles[plist[ip]];
  jpart = &particles[plist[jp]];

  // precompute vars from macro props
  double eps = sqrt(E);
  // has already been divided by rho from earlier
  double qmag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
  double qge = qmag / pow(eps,3.0) / rho;
  double itheta = qge + sqrt(1.0 + qge*qge);
  double alpha = eps*itheta;
  double beta = eps/itheta;

  // find direction of velocity wrt CoM frame

  double uvec[3];
  if (qmag < SMALL) { // revert back to energy scheme
    double theta = 2.0 * 3.14159 * random->uniform();
    double phi = acos(1.0 - 2.0 * random->uniform());
    uvec[0] = sin(phi) * cos(theta);
    uvec[1] = sin(phi) * sin(theta);
    uvec[2] = cos(phi);
    alpha = beta = eps;
  } else {
    for (int d = 0; d < 3; d++) uvec[d] = q[d]/qmag;
  }

  // set reduced particle velocities

  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + alpha*uvec[d];
    jpart->v[d] = V[d] - beta*uvec[d];
  }

  // set reduced particle weights

  double isw = rho/(1.0+itheta*itheta);
  double jsw = rho - isw;

  // set reduced particle rotational energies

  ipart->erot = Erot/rho;
  jpart->erot = Erot/rho;

  // set reduced particle weights

  ipart->weight = isw/update->fnum;
  jpart->weight = jsw/update->fnum;

  // if there is bad weight, only create one and clone

  if (isw != isw || isw <= 0.0 || jsw != jsw || jsw <= 0.0) {
    printf("rho: %g\n", rho);
    printf("V: %g %g %g\n", V[0], V[1], V[2]);
    printf("q: %g %g %g\n", q[0], q[1], q[2]);
    printf("eps: %g\n", eps);
    printf("qge: %g\n", qge);
    error->one(FLERR,"Bad weights");
  }

  // delete other particles

  for (int i = istart; i < iend; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[plist[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = plist[i];
    idelete++;
  }

  return;
}

/* ---------------------------------------------------------------------
   Merge particles using stress scheme
------------------------------------------------------------------------ */
// TODO : Check that stress tensor from group is conserved
void Collide::reduce(int istart, int iend,
                     double rho, double *V,
                     double *q, double pij[3][3],
                     double Erot)
{
  // find eigenpairs of stress tensor

  double ilambda[3], itheta[3][3];
  int ierror = MathEigen::jacobi3(pij,ilambda,itheta);

  // find number of non-zero eigenvalues
  // columns are eigenvectors

  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (fabs(ilambda[i]) > SMALL) {
      ilambda[nK] = ilambda[i];
      for (int d = 0; d < 3; d++) itheta[d][nK] = itheta[d][i];
      nK++;
    }
  }

  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  double qhat, c_theta, gam, igam;
  double isw, jsw;
  double uvec[3];

  for (int iK = 0; iK < nK; iK++) {

    // reduced particles chosen as first two

    ipart = &particles[plist[2*iK+istart]];
    jpart = &particles[plist[2*iK+1+istart]];

    // project heat flux to one of the eigenvectors
    qhat = itheta[0][iK]*q[0] + itheta[1][iK]*q[1] + itheta[2][iK]*q[2];
    if (qhat < 0)
      for (int d = 0; d < 3; d++) itheta[d][iK] *= -1.0;
    qhat = fabs(qhat);

    // factor to scale theta for setting final velocity
    c_theta = sqrt(nK*ilambda[iK]/rho);

    igam = sqrt(rho)*qhat/sqrt(nK)/pow(ilambda[iK],1.5);
    gam = igam + sqrt(1.0 + igam*igam);

    // set reduced particle velocities

    for (int d = 0; d < 3; d++) {
      ipart->v[d] = V[d] + gam*c_theta*itheta[d][iK];
      jpart->v[d] = V[d] - 1.0/gam*c_theta*itheta[d][iK];
    }

    // set reduced particle weights

    isw = rho/(nK*(1.0+gam*gam));
    jsw = rho/nK - isw;

    // if there is bad weight, only create one and clone

    if (isw != isw || isw <= 0.0 || jsw != jsw || jsw <= 0.0) {
      error->one(FLERR,"Bad Weights");
    }

    // set reduced particle rotational energies

    ipart->erot = Erot/rho/nK;
    jpart->erot = Erot/rho/nK;

    // set reduced weights

    ipart->weight = isw/update->fnum;
    jpart->weight = jsw/update->fnum;


  } // end nK
  
  // delete other particles
  for (int i = istart; i < iend; i++) {
    if (i < 2*nK + istart) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[plist[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = plist[i];
    idelete++;
  }

  return;
}
