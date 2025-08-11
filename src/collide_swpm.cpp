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
  int level;
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
      if (isw <= 0.0) {
        printf("%g\n", isw);
        error->one(FLERR,"Particle has negative or zero weight");
      }

      ip = next[ip];
    }
    sweight_max *= update->fnum;

    // scale max in cell count with cell level
    cell_scale = 1.0;
    //level = cells[icell].level;
    //if (level == 1) cell_scale = 1.0;
    //else if (domain->dimension == 2) cell_scale = pow(4,level-1);
    //else cell_scale = pow(8,level-1);

    // attempt = exact collision attempt count for all particles in cell
    // nattempt = rounded attempt with RN
    // if no attempts, continue to next grid cell

    if (np >= Ncmin*cell_scale && Ncmin > 0.0) pre_wtf = 0.0;
    else pre_wtf = 1.0;

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

  if (remove_min_flag) remove_tiny();

  return;
}

/* ----------------------------------------------------------------------
   Stochastic weighted algorithm (mixtures)
------------------------------------------------------------------------- */

void Collide::collisions_group_sw()
{
  int i,j,n,ip,np,newp;
  int isp,ipair,igroup,jgroup,newgroup,ngmax;
  int ng,ii,jj;
  int *ni,*nj,*ilist,*jlist;
  int nattempt;
  double attempt,volume;
  Particle::OnePart *ipart,*jpart,*kpart,*lpart,*mpart;

  // loop over cells I own

  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::ChildCell *cells = grid->cells;

  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int *species2group = mixture->species2group;

  double isw;
  double sw_max_spec[ngroups];
  int level;
  double cell_scale;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;
    if (np <= 1) continue;
    ip = cinfo[icell].first;
    volume = cinfo[icell].volume / cinfo[icell].weight;
    if (volume == 0.0) error->one(FLERR,"Collision cell volume is zero");

    // reallocate plist and p2g if necessary

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(plist);
      memory->create(plist,npmax,"collide:plist");
      memory->destroy(p2g);
      memory->create(p2g,npmax,2,"collide:p2g");
    }

    for (i = 0; i < ngroups; i++) {
      ngroup[i] = 0;
      sw_max_spec[i] = 0.0;
    }
    n = 0;

    // create particle list

    while (ip >= 0) {
      isp = particles[ip].ispecies;

      // store max weight for each species
      ipart = &particles[ip];
      isw = ipart->weight;
      if (isw != isw) error->all(FLERR,"Particle has NaN weight");
      if (isw <= 0.0) error->all(FLERR,"Particle has negative or zero weight");
      sw_max_spec[isp] = MAX(sw_max_spec[isp],isw);

      igroup = species2group[isp];
      if (ngroup[igroup] == maxgroup[igroup]) {
        maxgroup[igroup] += DELTAPART;
        memory->grow(glist[igroup],maxgroup[igroup],"collide:glist");
      }
      ng = ngroup[igroup];
      glist[igroup][ng] = n;
      p2g[n][0] = igroup;
      p2g[n][1] = ng;
      plist[n] = ip;
      ngroup[igroup]++;
      n++;
      ip = next[ip];
    } // end ip while

    // calculate attempts

    npair = 0;
    for (igroup = 0; igroup < ngroups; igroup++) {
      for (jgroup = igroup; jgroup < ngroups; jgroup++) {
        sweight_max = MAX(sw_max_spec[igroup],sw_max_spec[jgroup]);
        attempt = attempt_collision(icell,igroup,jgroup,volume);
        nattempt = static_cast<int> (attempt);

        if (nattempt) {
          gpair[npair][0] = igroup;
          gpair[npair][1] = jgroup;
          gpair[npair][2] = nattempt;
          nattempt_one += nattempt;
          npair++;
        }
      } // end jgroup
    } // end igroup

    // iterate thru species pairs

    for (ipair = 0; ipair < npair; ipair++) {
      igroup = gpair[ipair][0];
      jgroup = gpair[ipair][1];
      nattempt = gpair[ipair][2];

      ni = &ngroup[igroup];
      nj = &ngroup[jgroup];
      ilist = glist[igroup];
      jlist = glist[jgroup];

      if (*ni == 0 || *nj == 0) continue;
      if (igroup == jgroup && *ni == 1) continue;

      sweight_max = MAX(sw_max_spec[igroup],sw_max_spec[jgroup]);
      for (int iattempt = 0; iattempt < nattempt; iattempt++) {

        i = *ni * random->uniform();
        j = *nj * random->uniform();
        if (igroup == jgroup)
          while (i == j) j = *nj * random->uniform();

        ipart = &particles[plist[ilist[i]]];
        jpart = &particles[plist[jlist[j]]];

        // test if collision actually occurs
        // continue to next collision if no reaction

        if (!test_collision(icell,igroup,jgroup,ipart,jpart)) continue;

        // split particles

        newp = split(ipart,jpart,kpart,lpart);

        // add new particles to particle list

        if (newp > 1) {
          if (np+2 >= npmax) {
            npmax += DELTAPART;
            memory->grow(plist,npmax,"collide:plist");
            memory->grow(p2g,npmax,2,"collide:p2g");
          }

          plist[np++] = particle->nlocal-2;
          newgroup = species2group[kpart->ispecies];
          addgroup(newgroup,np-1);
          ilist = glist[igroup];
          jlist = glist[jgroup];

          plist[np++] = particle->nlocal-1;
          newgroup = species2group[lpart->ispecies];
          addgroup(newgroup,np-1);
          ilist = glist[igroup];
          jlist = glist[jgroup];

          particles = particle->particles;
        } else if (newp > 0) {
          if (np+1 >= npmax) {
            npmax += DELTAPART;
            memory->grow(plist,npmax,"collide:plist");
            memory->grow(p2g,npmax,2,"collide:p2g");
          }

          plist[np++] = particle->nlocal-1;
          newgroup = species2group[kpart->ispecies];
          addgroup(newgroup,np-1);
          ilist = glist[igroup];
          jlist = glist[jgroup];

          particles = particle->particles;
        }

        // since ipart and jpart have same weight, do not need
        // ... to account for weight during collision itself
        // also the splits are all handled beforehand

        mpart = NULL; // dummy particle
        setup_collision(ipart,jpart);
        perform_collision(ipart,jpart,mpart);
        ncollide_one++;

        // test to exit attempt loop due to groups becoming too small

        if (*ni <= 1) {
          if (*ni == 0) break;
          if (igroup == jgroup) break;
        }
        if (*nj <= 1) {
          if (*nj == 0) break;
          if (igroup == jgroup) break;
        }

      } // end attempts
    } // end species pairs

  } // end cells

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
  int level;
  double cell_scale, Ncmax_scale;
  int total_iter;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    // scale max in cell count with cell level
    cell_scale = 1.0;
    //level = cells[icell].level;
    //if (level == 1) cell_scale = 1.0;
    //else if (domain->dimension == 2) cell_scale = pow(4,level-1);
    //else cell_scale = pow(8,level-1);

    // upper bound to 1.5 times the max group size
    Ncmax_scale = MAX(Ncgmin,Ncmax/cell_scale); // make a user input

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
  int np = iend-istart;
  if (np <= Ngmin) return;

  // compute stress tensor since it's needed for
  // .. further dividing and reduction

  double gsum, msum, mV[3], mVV[3][3], mVVV[3][3];
  double count, uV[3], uVV[3][3]; // for unweighted division
  gsum = msum = count = 0.0;
  for (int i = 0; i < 3; i++) {
    mV[i] = 0.0;
    uV[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      mVV[i][j] = 0.0;
      mVVV[i][j] = 0.0;
      uVV[i][j] = 0.0;
    }
  }

  // find maximum particle weight

  int ispecies;
  double mass, psw, pmsw, vp[3];
  double Erot = 0.0;
  double KE_tmp = 0.0;
  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = ipart->weight;
    count += 1;
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));
   	gsum += psw;

    msum += pmsw;
    Erot += psw*ipart->erot;
    for (int i = 0; i < 3; i++) {
      mV[i] += (pmsw*vp[i]);
      uV[i] += vp[i];
      KE_tmp += psw*vp[i]*vp[i];
      //for (int j = 0; j < 3; j++) {
      //  mVV[i][j] += (pmsw*vp[i]*vp[j]);
      //  mVVV[i][j] += (pmsw*vp[i]*vp[j]*vp[j]);
      //  uVV[i][j] += (vp[i]*vp[j]);
      //}
    }
  }

  // mean velocity

	double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // stress tensor

  double pij[3][3], q[3];
  for (int i = 0; i < 3; i++) {
    q[i] = 0.0;
    for (int j = 0; j < 3; j++) pij[i][j] = 0.0;
  }

  //pij[i][j] = mVV[i][j] - mV[i]*mV[j]/msum;

  // manual commputation better
  for (int p = istart; p < iend; p++) {
    ipart = &particles[plist[p]];
    ispecies = ipart->ispecies;
    mass = species[ispecies].mass;

    psw = ipart->weight;
    pmsw = psw * mass;
    memcpy(vp, ipart->v, 3*sizeof(double));

    double c[3];
    c[0] = vp[0]-V[0];
    c[1] = vp[1]-V[1];
    c[2] = vp[2]-V[2];    
    double csq = MathExtra::lensq3(c);

    for (int i = 0; i < 3; i++) {
      q[i] += 0.5*pmsw*c[i]*csq;
      for (int j = 0; j < 3; j++)
        pij[i][j] += pmsw*c[i]*c[j];
    }
  }

  // if group is small enough, merge the particles

  if (np <= Ngmax+gbuf) {
    // reduce based on type
    if (reduction_type == ENERGY) {
      // specific kinetic energy in center of mass reference
      double E = (pij[0][0] + pij[1][1] + pij[2][2])/msum;
      reduce(istart, iend, gsum, V, E, Erot);
    } else if (reduction_type == HEAT) {
      double E = (pij[0][0] + pij[1][1] + pij[2][2])/msum;
      // specific heat flux
      // q/(g*mass) where mass is particle mass in Rja paper
      ipart = &particles[plist[istart]];
      mass = species[ipart->ispecies].mass;
      for (int d = 0; d < 3; d++) q[d] /= mass;
      reduce(istart, iend, gsum, V, E, q, Erot);
    } else if (reduction_type == STRESS) {
      printf("rho: %g\n", gsum);
      printf("KE: %g\n", KE_tmp);
      printf("V: %g %g %g\n", mV[0]/mass, mV[1]/mass, mV[2]/mass);
      printf("q: %g %g %g\n", q[0], q[1], q[2]);
      printf("pij[0]: %g %g %g\n", pij[0][0], pij[0][1], pij[0][2]);
      printf("pij[0]: %g %g %g\n", pij[1][0], pij[1][1], pij[1][2]);
      printf("pij[0]: %g %g %g\n", pij[2][0], pij[2][1], pij[2][2]);

      // mass-normalized heat flux
      for (int d = 0; d < 3; d++) q[d] = q[d]*gsum/msum;
      // in Lama paper, difference of half factor and assumes equal mass
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          pij[i][j] = pij[i][j]/msum*gsum*2.0;
      reduce(istart, iend, gsum, V, q, pij, Erot);
    }
  // group still too large so divide further

  } else {

    // Compute covariance matrix
    // unweighted needs to be updated

    double Rij[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (weighted) Rij[i][j] = pij[i][j]/msum;
        else Rij[i][j] = uVV[i][j] - uV[i]*uV[j]/count;
      }
    }

    // Find eigenpairs

    double eval[3], evec[3][3];
    int ierror = MathEigen::jacobi3(Rij,eval,evec);

    // Find largest eigenpair

    double maxeval;
    double maxevec[3]; // normal of splitting plane

    maxeval = 0;
    for (int i = 0; i < 3; i++) {
      if (std::abs(eval[i]) > maxeval) {
        maxeval = std::abs(eval[i]);
        for (int j = 0; j < 3; j++) {
          maxevec[j] = evec[j][i];  
        }
      }
    }

    // Separate based on particle velocity

    if (!weighted)
      for (int i = 0; i < 3; i++) V[i] = uV[i]/count;

    double center = V[0]*maxevec[0] + V[1]*maxevec[1] + V[2]*maxevec[2];
    int cp = istart;
    double gL, gR;
    gL = gR = 0.0;
    for (int i = istart; i < iend; i++) {
      ipart = &particles[plist[i]];
      if (MathExtra::dot3(ipart->v,maxevec) < center) {
        std::swap(plist[cp++],plist[i]);
        gL += ipart->weight;
      } else {
        gR += ipart->weight;
      }
    }

    // may occur due to machine precision issues
    //if(npL < 1) error->all(FLERR,"No particles in left group");
    //if(npR < 1) error->all(FLERR,"No particles in right group");

    // balance groups by number per group
    // will look for which group is bigger
    // then will sort that part of plist according to distance
    // ... of each point to the center hyperplane
    // then shift the center point
    if (balance_swpm_flag == NUMBER_BALANCE) {
      double vmV[3];
      int dnpLR = 2*cp-iend-istart;
      // number of particles to move
      int npmove =
        static_cast<int>(std::floor(fabs(dnpLR)*0.5+random->uniform()));

      // left group too big
      if (dnpLR > bst_number_thresh) {
        int np_group = cp-istart;
        double vel_dist[np_group];
        double dist;
        for (int i = istart; i < cp; i++) {
          ipart = &particles[plist[i]];
          vmV[0] = ipart->v[0]-V[0];
          vmV[1] = ipart->v[1]-V[1];
          vmV[2] = ipart->v[2]-V[2];
          dist = fabs(MathExtra::dot3(vmV,maxevec));
          vel_dist[i] = dist;
        }

        // lambda function
        // sorts from greatest to least
        // doesn't work (indices not lining up) - this is optimal
        /*std::sort(plist+istart,plist+istart+cp,
                  [&vel_dist,istart](int a, int b) {
                    printf("%i/%i, %g,%g - %i\n",
                      a-istart,b-istart,
                      vel_dist[a-istart],
                      vel_dist[b-istart],
                      vel_dist[a-istart] > vel_dist[b-istart]);
                    return vel_dist[a-istart] > vel_dist[b-istart];
                  });*/

        // selection sort from greatest to smallest distance
        int iswap;
        double max_dist = 0.0;
        for (int i = istart; i < cp; i++) {
          for (int j = i; j < cp; j++) {
            if (j == i || vel_dist[j] > max_dist) {
              iswap = j;
              max_dist = vel_dist[j];
            }
          }
          if (iswap != i) {
            std::swap(plist[i],plist[iswap]);
            std::swap(vel_dist[i],vel_dist[iswap]);
          }
        }

        // shift center point to left
        cp -= npmove;

      // right group too big
      } else if (dnpLR < -bst_number_thresh) {
        int np_group = iend-cp;
        double vel_dist[np_group];
        double dist;
        for (int i = cp; i < iend; i++) {
          ipart = &particles[plist[i]];
          vmV[0] = ipart->v[0]-V[0];
          vmV[1] = ipart->v[1]-V[1];
          vmV[2] = ipart->v[2]-V[2];
          dist = fabs(MathExtra::dot3(vmV,maxevec));
          vel_dist[i] = dist;
        }

        // selection sort from smallest to greatest distance
        int iswap;
        double min_dist = 0.0;
        for (int i = cp; i < iend; i++) {
          for (int j = i; j < iend; j++) {
            if (j == i || vel_dist[j] < min_dist) {
              iswap = j;
              min_dist = vel_dist[j];
            }
          }
          if (iswap != i) {
            std::swap(plist[i],plist[iswap]);
            std::swap(vel_dist[i],vel_dist[iswap]);
          }
        }

        // shift center point to right
        cp += npmove;

      }
    } // do weight balance later

    group_bt(istart,cp);
    group_bt(cp,iend);
  }

  return;
}

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
  ip = np * random->uniform() + istart;
  jp = np * random->uniform() + istart;
  while (ip == jp) jp = np * random->uniform() + istart;

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

  ipart->weight = rho*0.5;
  jpart->weight = rho*0.5;

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
  double qge = qmag / pow(eps,3.0);
  double itheta = qge + sqrt(1.0 + qge*qge);
  double alpha = eps*itheta;
  double beta = eps/itheta;

  // find direction of velocity wrt CoM frame

  double uvec[3];
  if (itheta < SMALL) { // revert back to energy scheme
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

    if (ipart->v[d] != ipart->v[d]) {
      printf("alpha %g beta %g eps %g q %g %g %g\n", 
              alpha, beta, eps, q[0], q[1], q[2]);
      error->one(FLERR,"ck");
    }
  }



  // set reduced particle weights

  double isw = rho/(1.0+itheta*itheta);
  double jsw = rho - isw;

  // set reduced particle rotational energies

  ipart->erot = Erot/rho;
  jpart->erot = Erot/rho;

  // set reduced particle weights

  ipart->weight = isw;
  jpart->weight = jsw;

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
    if (fabs(ilambda[i]) >= SMALL && ilambda[i] > 0) {
      ilambda[nK] = ilambda[i];
      for (int d = 0; d < 3; d++) itheta[d][nK] = itheta[d][i];
      nK++;
    }
  }

  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  double qhat, c_theta, gam;
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

    gam = sqrt(rho)*qhat/sqrt(nK)/pow(ilambda[iK],1.5);
    gam = gam + sqrt(1.0 + gam*gam);

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

    ipart->weight = isw;
    jpart->weight = jsw;


  } // end nK
  

  // check if right

  double rho_tmp = 0.0;
  double V_tmp[3];
  V_tmp[0] = V_tmp[1] = V_tmp[2] = 0.0;
  double KE_tmp = 0.0;

  for (int iK = 0; iK < nK; iK++) {

    // reduced particles chosen as first two

    ipart = &particles[plist[2*iK+istart]];
    jpart = &particles[plist[2*iK+1+istart]];

    // set reduced particle velocities

    rho_tmp += ipart->weight + jpart->weight;
    for (int d = 0; d < 3; d++) {
      V_tmp[d] += ipart->v[d]*ipart->weight;
      V_tmp[d] += jpart->v[d]*jpart->weight;
    }
  }

  printf("%g = %g\n", rho, ipart->weight+jpart->weight);
  KE_tmp = ipart->weight*(
    ipart->v[0]*ipart->v[0] +
    ipart->v[1]*ipart->v[1] +
    ipart->v[2]*ipart->v[2]);
  KE_tmp += jpart->weight*(
    jpart->v[0]*jpart->v[0] +
    jpart->v[1]*jpart->v[1] +
    jpart->v[2]*jpart->v[2]);
  KE_tmp /= rho;

  V_tmp[0] = ipart->weight*ipart->v[0]+jpart->weight*jpart->v[0];
  V_tmp[1] = ipart->weight*ipart->v[1]+jpart->weight*jpart->v[1];
  V_tmp[2] = ipart->weight*ipart->v[2]+jpart->weight*jpart->v[2];

  V_tmp[0] /= rho;
  V_tmp[1] /= rho;
  V_tmp[2] /= rho;

  double q_tmp[3];
  q_tmp[0] = q_tmp[1] = q_tmp[2] = 0.0;

  double psw = ipart->weight;
  double c[3];
  c[0] = ipart->v[0]-V[0];
  c[1] = ipart->v[1]-V[1];
  c[2] = ipart->v[2]-V[2];    
  double csq = MathExtra::lensq3(c);

  for (int i = 0; i < 3; i++) q_tmp[i] += 0.5*psw*c[i]*csq;

  psw = jpart->weight;
  c[0] = jpart->v[0]-V[0];
  c[1] = jpart->v[1]-V[1];
  c[2] = jpart->v[2]-V[2];    
  csq = MathExtra::lensq3(c);

  for (int i = 0; i < 3; i++) q_tmp[i] += 0.5*psw*c[i]*csq;

  printf("Vin: %g %g %g\n", V[0], V[1], V[2]);
  printf("%g %g %g; %g\n", V_tmp[0], V_tmp[1], V_tmp[2], KE_tmp);
  printf("q: %g %g %g\n", q_tmp[0], q[1], q[2]);

  printf("pij[0]: %g %g %g\n", pij[0][0], pij[0][1], pij[0][2]);
  printf("pij[0]: %g %g %g\n", pij[1][0], pij[1][1], pij[1][2]);
  printf("pij[0]: %g %g %g\n", pij[2][0], pij[2][1], pij[2][2]);
  error->one(FLERR,"ck");

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

/* ----------------------------------------------------------------------
   Delete tiny weighted particles
------------------------------------------------------------------------- */
void Collide::remove_tiny()
{
  int np, ip, n;
  double isw, sw_max;
  Grid::ChildInfo *cinfo = grid->cinfo;
  int *next = particle->next;

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart;

  for (int icell = 0; icell < nglocal; icell++) {
    np = cinfo[icell].count;

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = ipart->weight;
      if (isw <= 0) {
        error->one(FLERR,"negatvie weight found after compression");
      }
      ip = next[ip];
    }


    /*ip = cinfo[icell].first;
    n = 0;
    sw_max = 0.0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = ipart->weight;
      if (isw > 0) sw_max = MAX(sw_max,isw);
      ip = next[ip];
      n++;
    }*/

    // delete tiny weights

    /*ip = cinfo[icell].first;
    while (ip >= 0) {
      if (isw < sw_max*min_weight) {
        if (ndelete == maxdelete) {
          maxdelete += DELTADELETE;
          memory->grow(dellist,maxdelete,"collide:dellist");
        }
        ipart = &particles[ip];
        ipart->weight = -1.0;
        dellist[ndelete++] = ip;
      }
      ip = next[ip];
    }*/
  } // loop for cells

  return;
}

