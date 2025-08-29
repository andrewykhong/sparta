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
#define SMALL 1.0e-16

/* ----------------------------------------------------------------------
   Stochastic weighted algorithm
------------------------------------------------------------------------- */

void Collide::collisions_one_sw()
{
  int i,j,n,ip,np,newp;
  int nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*lpart,*mpart;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double isw;
  double cell_scale;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;

    volume = cinfo[icell].volume / cinfo[icell].weight;
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // setup particle list for this cell

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
      memory->create(pL,npmax,"collide:pL");
      memory->create(pLU,npmax,"collide:pLU");
    }

    // build particle list and find maximum particle weight
    // particle weights are relative to update->fnum

    ip = cinfo[icell].first;
    n = 0;
    sweight_max = 0.0;
    while (ip >= 0) {
      plist[n++] = ip;
      ipart = &particles[ip];
      isw = ipart->weight;
      sweight_max = MAX(sweight_max,isw);
      if (isw != isw) error->one(FLERR,"Particle has NaN weight");
      if (isw <= 0.0) error->one(FLERR,"Particle has nonpositive weight");

      ip = next[ip];
    }
    sweight_max *= update->fnum;

    // scale max in cell count with cell level
    if (adapt_flag) {
      int level = cells[icell].level;
      if (level == 1) cell_scale = 1.0;
      else if (domain->dimension == 2) cell_scale = pow(4,level-1);
      else cell_scale = pow(8,level-1);
    } else {
      cell_scale = 1.0;
    }

    // attempt = exact collision attempt count for all particles in cell
    // nattempt = rounded attempt with RN
    // if no attempts, continue to next grid cell

    if (np >= MAX(Ncgmin*0.5,Ncmin/cell_scale) && Ncmin > 0.0)
      pre_wtf = 0.0;
    else
      pre_wtf = 1.0;

    attempt = attempt_collision(icell,np,volume);
    nattempt = static_cast<int> (attempt);

    if (!nattempt) continue;
    nattempt_one += nattempt;

    for (int iattempt = 0; iattempt < nattempt; iattempt++) {

      i = np * random->uniform();
      j = np * random->uniform();
      while (i == j) j = np * random->uniform();
      ipart = &particles[plist[i]];
      jpart = &particles[plist[j]];

      if (!test_collision(icell,0,0,ipart,jpart)) continue;

      // split particles

      newp = split(ipart,jpart,kpart,lpart);

      // add new particles to particle list

      if (newp > 1) {
        if (np+2 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
          memory->grow(pL,npmax,"collide:pL");
          memory->grow(pLU,npmax,"collide:pLU");
        }
        plist[np++] = particle->nlocal-2;
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      } else if (newp > 0) {
        if (np+1 >= npmax) {
          npmax += DELTAPART;
          memory->grow(plist,npmax,"collide:plist");
          memory->grow(pL,npmax,"collide:pL");
          memory->grow(pLU,npmax,"collide:pLU");
        }
        plist[np++] = particle->nlocal-1;
        particles = particle->particles;
      }

      // since ipart and jpart have same weight, do not need
      // ... to account for weight during collision itself
      // also the splits are all handled beforehand

      mpart = NULL; // dummy particle
      setup_collision(ipart,jpart);
      perform_collision(ipart,jpart,mpart);
      ncollide_one++;

    } // end attempt loop
  } // loop for cells

  // remove tiny weighted particles

  //if (remove_min_flag) remove_tiny();

  return;
}

/* ----------------------------------------------------------------------
   Splits particles and generates two new particles (for SWPM)
------------------------------------------------------------------------- */

int Collide::split(Particle::OnePart *&ip, Particle::OnePart *&jp,
                   Particle::OnePart *&kp, Particle::OnePart *&lp)
{
  double xk[3],vk[3];
  double xl[3],vl[3];
  double erotk, erotl;
  int ks, ls;
  int kcell, lcell;

  // checks if particles properly deleted

  kp = NULL;
  lp = NULL;

  // weight transfer function is assumed to be
  // ... MIN(ip->sweight,jp->sweight)/(1 + pre_wtf * wtf)

  double isw = ip->weight;
  double jsw = jp->weight;
  double Gwtf, ksw, lsw;

  if (isw <= 0.0 || jsw <= 0.0)
    error->one(FLERR,"Zero or negative weight before split");

  // particle ip has larger weight

  if(isw >= jsw) {
    Gwtf = jsw/(1.0+pre_wtf*wtf);
    ksw  = isw-Gwtf;
    lsw  = jsw-Gwtf;

    ks = ip->ispecies;
    ls = jp->ispecies;

    kcell = ip->icell;
    lcell = jp->icell;

    memcpy(xk,ip->x,3*sizeof(double));
    memcpy(vk,ip->v,3*sizeof(double));
    memcpy(xl,jp->x,3*sizeof(double));
    memcpy(vl,jp->v,3*sizeof(double));

    erotk = ip->erot;
    erotl = jp->erot;

  // particle jp has larger weight

  } else {
    Gwtf = isw/(1.0+pre_wtf*wtf);
    ksw  = jsw-Gwtf;
    lsw  = isw-Gwtf;

    ks = jp->ispecies;
    ls = ip->ispecies;

    kcell = jp->icell;
    lcell = ip->icell;

    memcpy(xk,jp->x,3*sizeof(double));
    memcpy(vk,jp->v,3*sizeof(double));
    memcpy(xl,ip->x,3*sizeof(double));
    memcpy(vl,ip->v,3*sizeof(double));

    erotk = jp->erot;
    erotl = ip->erot;
  }

  // update weights

  ip->weight = Gwtf;
  jp->weight = Gwtf;

  // Gwtf should never be negative or zero

  if (Gwtf <= 0.0)
    error->one(FLERR,"Negative weight assigned after split");

  if (Gwtf > 0.0 && pre_wtf > 0.0)
    if (ksw <= 0.0 || lsw <= 0.0)
      error->one(FLERR,"Zero or negative weight after split");

  // number of new particles

  int newp = 0;

  // gk is always the bigger of the two

  if(ksw > 0) {
    int id = MAXSMALLINT*random->uniform();
    Particle::OnePart *particles = particle->particles;
    int reallocflag = particle->add_particle(id,ks,kcell,xk,vk,erotk,0.0);
    if (reallocflag) {
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
    }
    kp = &particle->particles[particle->nlocal-1];
    kp->weight = ksw;
    newp++;
  }

  if (kp) {
    if (kp->weight <= 0.0) {
      printf("ksw: %2.3e; kp->weight: %2.3e\n", ksw,kp->weight);
      error->one(FLERR,"New particle [k] has bad weight");
    }
  }

  // there should never be case where you add particle "l" if
  // ... you did not add particle "k"

  if(lsw > 0) {
    if(ksw <= 0) error->one(FLERR,"Bad addition to particle list");
    int id = MAXSMALLINT*random->uniform();
    Particle::OnePart *particles = particle->particles;
    int reallocflag = particle->add_particle(id,ls,lcell,xl,vl,erotl,0.0);
    if (reallocflag) {
      ip = particle->particles + (ip - particles);
      jp = particle->particles + (jp - particles);
      kp = particle->particles + (kp - particles);
    }
    lp = &particle->particles[particle->nlocal-1];
    lp->weight = lsw;
    newp++;
  }

  if (lp) {
    if (lp->weight <= 0.0) {
      printf("lsw: %2.3e; lp->weight: %2.3e\n", lsw,lp->weight);
      error->one(FLERR,"New particle [l] has bad weight");
    }
  }

  return newp;
}

/* ----------------------------------------------------------------------
   Reorder plist depending on grouping strategy used and prepare for
   grouping and particle reduction
------------------------------------------------------------------------- */

void Collide::group_reduce()
{
  int n,nold,np,ip;
  double isw;

  Particle::OnePart *ipart;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;

  double swmean, swvar, swstd;
  double d1, d2;
  double lLim, uLim;
  int npL, npLU;
  double cell_scale, Ncmax_scale;
  int total_iter;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    // scale max in cell count with cell level
    if (adapt_flag) {
      int level = cells[icell].level;
      if (level == 1) cell_scale = 1.0;
      else if (domain->dimension == 2) cell_scale = pow(4,level-1);
      else cell_scale = pow(8,level-1);
    } else {
      cell_scale = 1.0;
    }

    // upper bound to double
    Ncmax_scale = MAX(2.0*Ncgmin,Ncmax/cell_scale);

    if (np <= Ncmax_scale) continue;

    gbuf = 0;
    total_iter = 0;
    while (1) {

      // recreate particle list for reduction
      ip = cinfo[icell].first;
      n = 0;
      while (ip >= 0) {
        ipart = &particles[ip];
        isw = ipart->weight;
        if(isw > 0) plist[n++] = ip;
        ip = next[ip];
      }
      idelete = 0; // keep track of how many deleted this iteration
      group_bt(0,n);

      if (n-idelete < Ncmax_scale) break;

      // if no particles reduced, increase group size
      // if less than 20% particles reduced, increase group size
      // 20% is arbitrary

      if (idelete < 0.20*n) gbuf += 1;

      // only increase buffer up to twice max group size

      if (gbuf >= Ngmax) break;

      total_iter++;

      if (total_iter > 10) break;

    } // while loop for n > ncmax
  }// loop for cells

  return;
}

/* ----------------------------------------------------------------------
   Recursivley divides particles using the binary tree strategy
------------------------------------------------------------------------- */
void Collide::group_bt(int istart, int iend)
{
  Particle::OnePart *ipart;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  // ignore groups which have too few particles

  // loops don't include iend
  int np = iend - istart;
  if (np <= Ngmin) return;

  // compute stress tensor since it's needed for
  // .. further dividing and reduction
  // since in each group, there is only one species, don't need mass

  int count;
  double gsum, mV[3];
  double mVV[3][3]; // for covariance matrix
  gsum  = count = 0.0;
  for (int i = 0; i < 3; i++) {
    mV[i] = 0.0;
    for (int j = 0 ; j < 3; j++) mVV[i][j] = 0.0;
  }
    
  // find maximum particle weight

  int ispecies;
  double mass, iweight, imass;
  double Erot = 0.0;
  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    iweight = ipart->weight * update->fnum;

   	gsum += iweight;
    Erot += iweight*ipart->erot;
    for (int i = 0; i < 3; i++) {
      mV[i] += iweight*ipart->v[i];
      for (int j = i; j < 3; j++)
        mVV[i][j] += iweight*ipart->v[i]*ipart->v[j];
    }
  }
  mVV[1][0] = mVV[0][1];
  mVV[2][0] = mVV[0][2];
  mVV[2][1] = mVV[1][2];

  // mean velocity

	double V[3], M[3][3];
  for (int i = 0; i < 3; i++) {
    V[i] = mV[i]/gsum;
    for (int j = 0; j < 3; j++) M[i][j] = mVV[i][j]/gsum;
  }

  // stress tensor / heat flux

  double pij[3][3], q[3];
  for (int i = 0; i < 3; i++) {
    q[i] = 0.0;
    for (int j = 0; j < 3; j++) pij[i][j] = 0.0;
  }

  // manual commputation better
  double c[3], csq;
  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    iweight = ipart->weight * update->fnum;

    c[0] = ipart->v[0]-V[0];
    c[1] = ipart->v[1]-V[1];
    c[2] = ipart->v[2]-V[2];    
    csq = MathExtra::lensq3(c);

    for (int i = 0; i < 3; i++) {
      q[i] += 0.5*iweight*c[i]*csq;
      for (int j = i; j < 3; j++)
        pij[i][j] += iweight*c[i]*c[j];
    }
  }

  // fill in remaining indices of symmetric matrix
  pij[1][0] = pij[0][1];
  pij[2][0] = pij[0][2];
  pij[2][1] = pij[1][2];

  // first determine if group needs to be divided
  int cp = istart;
  // weights on each side
  double gL, gR;
  gL = gR = 0.0;
  if (np > Ngmax+gbuf) {
    // Compute covariance matrix
    double Rij[3][3];
    double scale = fabs(pij[0][0]+pij[1][1]+pij[2][2])/3.0;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        Rij[i][j] = pij[i][j]/scale;
        //Rij[i][j] = M[i][j] - V[i]*V[j];

    // Find eigenpairs

    // columns are eigenvectors
    // last negative one is to output in decreasing order
    double eval[3], evec[3][3];
    jacobiEigen(Rij,eval,evec);

    // Find largest eigenpair

    double maxeval;
    double maxevec[3]; // normal of splitting plane

    // ignore negative eigenvalues (should never occur)
    maxeval = 0;
    for (int i = 0; i < 3; i++) {
      if (eval[i] > maxeval) {
        maxeval = fabs(eval[i]);
        for (int j = 0; j < 3; j++) {
          maxevec[j] = evec[j][i];  
        }
      }
    }

    // Separate based on particle velocity
    // cp is not included in left group
    double center = V[0]*maxevec[0] + V[1]*maxevec[1] + V[2]*maxevec[2];
    for (int i = istart; i < iend; i++) {
      ipart = &particles[plist[i]];
      if (MathExtra::dot3(ipart->v,maxevec) < center) {
        if (i != cp)
          std::swap(plist[cp],plist[i]);
        if (balance_flag == NUMBER_BALANCE) gL += 1.0;
        else gL += ipart->weight;
        cp++;
      } else {
        if (balance_flag == NUMBER_BALANCE) gR += 1.0;
        else gR += ipart->weight;
      }
    }

    // balance groups by number per group
    // will look for which group is bigger
    // then will sort that part of plist according to distance
    // ... of each point to the center hyperplane
    // then shift the center point
    if (balance_flag == NUMBER_BALANCE) {
      double vmV[3];
      int dnpLR = static_cast <int> (gL-gR);
      // number of particles to "move"
      int npmove =
        static_cast<int>(std::floor(fabs(dnpLR)*0.5));

      // left group too big
      if (dnpLR > bst_thresh) {
        // compute distances
        int np_group = cp-istart;
        double vel_dist[np_group];
        double dist;
        for (int i = istart; i < cp; i++) {
          ipart = &particles[plist[i]];
          vmV[0] = ipart->v[0]-V[0];
          vmV[1] = ipart->v[1]-V[1];
          vmV[2] = ipart->v[2]-V[2];
          dist = fabs(MathExtra::dot3(vmV,maxevec));
          vel_dist[i-istart] = dist;
        }

        // selection sort from greatest to smallest distance
        int iswap;
        double max_dist;
        for (int i = istart; i < cp; i++) {
          iswap = i;
          max_dist = vel_dist[i-istart];
          for (int j = i+1; j < cp; j++) {
            if (vel_dist[j-istart] > max_dist) {
              iswap = j;
              max_dist = vel_dist[j-istart];
            }
          }
          if (iswap != i) {
            std::swap(plist[i],plist[iswap]);
            std::swap(vel_dist[i-istart],vel_dist[iswap-istart]);
          }
        }

        // shift center point to left
        cp -= npmove;

      // right group too big
      } else if (dnpLR < -bst_thresh) {
        int np_group = iend-cp;
        double vel_dist[np_group];
        double dist;
        for (int i = cp; i < iend; i++) {
          ipart = &particles[plist[i]];
          vmV[0] = ipart->v[0]-V[0];
          vmV[1] = ipart->v[1]-V[1];
          vmV[2] = ipart->v[2]-V[2];
          dist = fabs(MathExtra::dot3(vmV,maxevec));
          vel_dist[i-cp] = dist;
        }

        // selection sort from smallest to greatest distance
        int iswap;
        double min_dist;
        for (int i = cp; i < iend; i++) {
          iswap = i;
          min_dist = vel_dist[i-cp];
          for (int j = i+1; j < iend; j++) {
            if (vel_dist[j-cp] < min_dist) {
              iswap = j;
              min_dist = vel_dist[j-cp];
            }
          }
          if (iswap != i) {
            std::swap(plist[i],plist[iswap]);
            std::swap(vel_dist[i-cp],vel_dist[iswap-cp]);
          }
        }

        // shift center point to right
        cp += npmove;
      }
    // weight based balanced
    } else if (balance_flag == WEIGHT_BALANCE) {
      double vmV[3];
      // total weight to move
      double dgLR = (gL-gR)*0.5;

      // left group too big
      if (dgLR > bst_thresh) {
        // compute distances
        int np_group = cp-istart;
        double vel_dist[np_group];
        double dist;
        for (int i = istart; i < cp; i++) {
          ipart = &particles[plist[i]];
          vmV[0] = ipart->v[0]-V[0];
          vmV[1] = ipart->v[1]-V[1];
          vmV[2] = ipart->v[2]-V[2];
          dist = fabs(MathExtra::dot3(vmV,maxevec));
          vel_dist[i-istart] = dist;
        }

        // selection sort from greatest to smallest distance
        int iswap;
        double max_dist;
        for (int i = istart; i < cp; i++) {
          iswap = i;
          max_dist = vel_dist[i-istart];
          for (int j = i+1; j < cp; j++) {
            if (vel_dist[j-istart] > max_dist) {
              iswap = j;
              max_dist = vel_dist[j-istart];
            }
          }
          if (iswap != i) {
            std::swap(plist[i],plist[iswap]);
            std::swap(vel_dist[i-istart],vel_dist[iswap-istart]);
          }
        }

        // shift center point to left
        double moved_weight = 0.0;
        while (1) {
          ipart = &particles[plist[cp-1]];
          iweight = ipart->weight;
          // prevents one massive particle from moving
          if (moved_weight + iweight > fabs(dgLR)) break;
          moved_weight += iweight;
          cp -= 1;
        }

      // right group too big
      } else if (dgLR < -bst_thresh) {
        int np_group = iend-cp;
        double vel_dist[np_group];
        double dist;
        for (int i = cp; i < iend; i++) {
          ipart = &particles[plist[i]];
          vmV[0] = ipart->v[0]-V[0];
          vmV[1] = ipart->v[1]-V[1];
          vmV[2] = ipart->v[2]-V[2];
          dist = fabs(MathExtra::dot3(vmV,maxevec));
          vel_dist[i-cp] = dist;
        }

        // selection sort from smallest to greatest distance
        int iswap;
        double min_dist;
        for (int i = cp; i < iend; i++) {
          iswap = i;
          min_dist = vel_dist[i-cp];
          for (int j = i+1; j < iend; j++) {
            if (vel_dist[j-cp] < min_dist) {
              iswap = j;
              min_dist = vel_dist[j-cp];
            }
          }
          if (iswap != i) {
            std::swap(plist[i],plist[iswap]);
            std::swap(vel_dist[i-cp],vel_dist[iswap-cp]);
          }
        }

        // shift center point to left
        double moved_weight = 0.0;
        while (1) {
          ipart = &particles[plist[cp]];
          iweight = ipart->weight;
          // prevents one massive particle from moving
          if (moved_weight + iweight > fabs(dgLR)) break;
          moved_weight += iweight;
          cp += 1;
        }
      }
    }
  }

  // force merge
  if (gL == 0.0 or gR == 0.0) cp = istart;

  // if group is small enough or poor division, merge the particles
  if (cp <= istart) {
    /*
    printf("rho: %g; n: %i; Ngmin: %i\n", gsum, np, Ngmin);
    printf("KE: %g\n", mVV[0][0]+mVV[1][1]+mVV[2][2]);
    printf("V: %g %g %g\n", V[0], V[1], V[2]);
    printf("q: %g %g %g\n", q[0], q[1], q[2]);
    printf("pij: %g %g %g\n", pij[0][0], pij[0][1], pij[0][2]);
    printf("   : %g %g %g\n", pij[1][0], pij[1][1], pij[1][2]);
    printf("   : %g %g %g\n\n", pij[2][0], pij[2][1], pij[2][2]);
    */

    // reduce based on type
    if (reduction_type == ENERGY) {
      // specific kinetic energy in center of mass reference
      double E = (pij[0][0] + pij[1][1] + pij[2][2])/gsum;
      if (E <= 1e-16) E = 0.0;
      reduce(istart, iend, gsum, V, E, Erot);
    } else if (reduction_type == HEAT) {
      double E = (pij[0][0] + pij[1][1] + pij[2][2])/gsum;
      if (E <= 1e-16) E = 0.0;
      // divide heat flux by mass of species because only
      // ... total weight is passed.
      reduce(istart, iend, gsum, V, E, q, Erot);
    } else if (reduction_type == STRESS) {
      reduce(istart, iend, gsum, V, q, pij, Erot);
    }

    /*
    // check that the reduction is correct
    gsum = count = 0.0;
    for (int i = 0; i < 3; i++) {
      mV[i] = 0.0;
      for (int j = 0; j < 3; j++) mVV[i][j] = 0.0;
    }

    for (int p = istart; p < iend; p++) {
      ipart = &particles[plist[p]];
      ispecies = ipart->ispecies;
      if (ipart->weight < 0) continue;

      iweight = ipart->weight * update->fnum;
      count += 1;
     	gsum += iweight;
      for (int i = 0; i < 3; i++) {
        mV[i] += iweight*ipart->v[i];
        for (int j = 0; j < 3; j++)
          mVV[i][j] += iweight*ipart->v[i]*ipart->v[j];
      }
    }

    for (int i = 0; i < 3; i++) V[i] = mV[i]/gsum;

    for (int i = 0; i < 3; i++) {
      q[i] = 0.0;
      for (int j = 0; j < 3; j++) pij[i][j] = 0.0;
    }

    // manual commputation better
    for (int p = istart; p < iend; p++) {
      ipart = &particles[plist[p]];
      ispecies = ipart->ispecies;
      if (ipart->weight < 0) continue;

      iweight = ipart->weight * update->fnum;

      c[0] = ipart->v[0]-V[0];
      c[1] = ipart->v[1]-V[1];
      c[2] = ipart->v[2]-V[2];    
      csq = MathExtra::lensq3(c);

      for (int i = 0; i < 3; i++) {
        q[i] += 0.5*iweight*c[i]*csq;
        for (int j = 0; j < 3; j++)
          pij[i][j] += iweight*c[i]*c[j];
      }
    }

    printf("rho: %g; n: %i\n", gsum, count);
    printf("KE: %g\n", mVV[0][0]+mVV[1][1]+mVV[2][2]);
    printf("V: %g %g %g\n", V[0], V[1], V[2]);
    printf("q: %g %g %g\n", q[0], q[1], q[2]);
    printf("pij: %g %g %g\n", pij[0][0], pij[0][1], pij[0][2]);
    printf("   : %g %g %g\n", pij[1][0], pij[1][1], pij[1][2]);
    printf("   : %g %g %g\n\n\n", pij[2][0], pij[2][1], pij[2][2]);
    error->one(FLERR,"ck");
    */

  } else {
    group_bt(istart,cp);
    group_bt(cp,iend);
  }

  return;
}

