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

  ipart->erot = Erot/(rho*0.5);
  jpart->erot = Erot/(rho*0.5);

  // set reduced particle weights

  ipart->weight = rho*0.5/update->fnum;
  jpart->weight = rho*0.5/update->fnum;

  // if there is bad weight, only create one and clone

  if (E < 0) error->one(FLERR,"Negative Energy");

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
  if (eps <= 1e-16) eps = 0.0;
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

  ipart->erot = Erot/isw;
  jpart->erot = Erot/jsw;

  // set reduced particle weights

  ipart->weight = isw/update->fnum;
  jpart->weight = jsw/update->fnum;

  // if there is bad weight, only create one and clone

  if (isw != isw || isw <= 0.0 || jsw != jsw || jsw <= 0.0) {
    printf("rho: %g\n", rho);
    printf("V: %g %g %g\n", V[0], V[1], V[2]);
    printf("q: %g %g %g\n", q[0], q[1], q[2]);
    printf("eps: %g: E: %g\n", eps, E);
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

  //MathEigen::jacobi3(pijn,ilambda,itheta);
  bool converged = jacobiEigen(pij,ilambda,itheta);
  //if (not converged) analytic(pij,ilambda,itheta);

  /*for (int i = 0; i < 3; i++) {
    if (ilambda[i] < 0.0) {
      printf("pij: %g %g %g\n", pij[0][0], pij[0][1], pij[0][2]);
      printf("   : %g %g %g\n", pij[1][0], pij[1][1], pij[1][2]);
      printf("   : %g %g %g\n", pij[2][0], pij[2][1], pij[2][2]);

      printf("%g - (%g %g %g)\n",
        ilambda[0], itheta[0][0], itheta[1][0], itheta[2][0]);
      printf("%g - (%g %g %g)\n",
        ilambda[1], itheta[0][1], itheta[1][1], itheta[2][1]);
      printf("%g - (%g %g %g)\n",
        ilambda[2], itheta[0][2], itheta[1][2], itheta[2][2]);
      error->one(FLERR,"Negative eigenvalue");
    }
  }*/

  // find number of non-zero eigenvalues
  // columns are eigenvectors
  // can at times find a "negative" eigenvalue
  // I believe these are due to those eigenvalues being very close to zero

  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (ilambda[i] > SMALL) {
      if (i != nK) {
        ilambda[nK] = ilambda[i];
        for (int d = 0; d < 3; d++) itheta[nK][d] = itheta[i][d];
      }
      nK++;
    }
  }
  
  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  double qhat, c_theta, gam, igam;
  double isw, jsw;
  double uvec[3];
  double ipair = 0;

  for (int iK = 0; iK < nK; iK++) {

    if (ilambda[iK] <= SMALL)
      error->one(FLERR,"Bad eigenvalue");

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
      printf("pij: %g %g %g\n", pij[0][0], pij[0][1], pij[0][2]);
      printf("   : %g %g %g\n", pij[1][0], pij[1][1], pij[1][2]);
      printf("   : %g %g %g\n", pij[2][0], pij[2][1], pij[2][2]);

      printf("gam: %g; rho: %g; qhat: %g\n",
        igam, rho, qhat);
      printf("%g - (%g %g %g)\n",
        ilambda[0], itheta[0][0], itheta[1][0], itheta[2][0]);
      printf("%g - (%g %g %g)\n",
        ilambda[1], itheta[0][1], itheta[1][1], itheta[2][1]);
      printf("%g - (%g %g %g)\n",
        ilambda[2], itheta[0][2], itheta[1][2], itheta[2][2]);
      error->one(FLERR,"Bad Weights");
    }

    // set reduced particle rotational energies

    ipart->erot = Erot/isw;
    jpart->erot = Erot/jsw;

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

/* ---------------------------------------------------------------------
   Jacobi rotation to find eigenpairs
------------------------------------------------------------------------ */
bool Collide::jacobiEigen(double const mat[3][3], double *eval, double evec[3][3], int max_iterations, double max_tolerance)
{
  // Initialize eigenvectors to the identity matrix
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      evec[i][j] = (i==j) ? 1.0 : 0.0;

  // copy the matrix
  double A[3][3] = {{mat[0][0], mat[0][1], mat[0][2]},
                    {mat[1][0], mat[1][1], mat[1][2]},
                    {mat[2][0], mat[2][1], mat[2][2]}};

  // normalize A and shift eigenvalues to help convergence

  double scale = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      scale = MAX(scale,fabs(A[i][j]));

  double lambda_shift = 1.5;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) A[i][j] = A[i][j]/scale;
    A[i][i] += lambda_shift;
  }

  bool converged = false;
  for (int iter = 0; iter < max_iterations; ++iter) {

    // Find the largest off-diagonal element
    double maxVal = 0.0;
    int p = 0, q = 1;
    for (int i = 0; i < 3; i++) {
      for (int j = i + 1; j < 3; j++) {
        if (fabs(A[i][j]) > maxVal) {
          maxVal = fabs(A[i][j]);
          p = i;
          q = j;
        }
      }
    }

    // Check for convergence
    if (maxVal < max_tolerance) {
      converged = true;
      break;
    }

    // Compute the Jacobi rotate
    jacobiRotate(A, evec, p, q);
  }

  // if not converged, use analytic
  if (not converged) return false;

  // Extract eigenvalues from the diagonal of A
  for (int i = 0; i < 3; i++) eval[i] = (A[i][i]-lambda_shift)*scale;
  return true;
}

/* ---------------------------------------------------------------------
   Jacobi rotation to find eigenpairs
------------------------------------------------------------------------ */
void Collide::jacobiRotate(double A[3][3], double V[3][3], int p, int q)
{
  if (A[p][q] == 0.0) return;

  double app = A[p][p];
  double aqq = A[q][q];
  double apq = A[p][q];

  double phi;
  if (fabs(aqq-app) < 1E-50) phi = (apq > 0) ? M_PI*0.25 : -M_PI*0.25;
  else phi = 0.5*atan2(2.0*apq,aqq-app);

  double c = cos(phi);
  double s = sin(phi);

  for (int j = 0; j < 3; j++) {
    double apj = A[p][j];
    double aqj = A[q][j];
    A[p][j] = c*apj-s*aqj;
    A[q][j] = s*apj+c*aqj;
  }

  for (int i = 0; i < 3; i++) {
    double aip = A[i][p];
    double aiq = A[i][q];
    A[i][p] = c*aip-s*aiq;
    A[i][q] = s*aip+c*aiq;
  }

  A[p][q] = A[q][p] = 0.0;

  for (int i = 0; i < 3; i++) {
    double vip = V[i][p];
    double viq = V[i][q];
    V[i][p] = c*vip-s*viq;
    V[i][q] = s*vip+c*viq;
  }
}

/* ---------------------------------------------------------------------
   Jacobi rotation to find eigenpairs
------------------------------------------------------------------------ */
/*void Collide::analytic(double const mat[3][3], double *eval, double evec[3][3])
{
  double m00 = mat[0][0], m01 = mat[0][1], m02 = mat[0][2];
  double m11 = mat[1][1], m12 = mat[1][2];
  double m22 = mat[2][2];

  // initialize eigenvectors with unit matrix
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      evec[i][j] = (i==j) ? 1.0 : 0.0;

  double p1 = m01*m01 + m02*m02 + m12*m12;
  if (p1 < 1E-50) {
    eval[0] = m00;
    eval[1] = m11;
    eval[2] = m22;
  }

  double q = (m00+m11+m22)/3.0;
  double p2 = (m00-q)*(m00-q)+(m11-q)*(m11-q)+(m22-q)*(m22-q))+2.0*p1;
  double p = sqrt(p2/6.0);

  double B[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      B[i][j] = (A[i][j] - ( (i==j) ? q : 0 ))/p;

  




}*/

/* --------------------------------------------------------------------- */

