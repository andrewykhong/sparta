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

  //if (comm->me == 6) printf("c: %i -- begin collision loop \n", comm->me);
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
      if (isw <= 0.0) error->one(FLERR,"Particle has negative or zero weight");

      ip = next[ip];
    }
    sweight_max *= update->fnum;

    // scale max in cell count with cell level
    level = cells[icell].level;
    if (level == 1) cell_scale = 1.0;
    else if (domain->dimension == 2) cell_scale = pow(4,level-1);
    else cell_scale = pow(8,level-1);

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
    level = cells[icell].level;
    if (level == 1) cell_scale = 1.0;
    else if (domain->dimension == 2) cell_scale = pow(4,level-1);
    else cell_scale = pow(8,level-1);

    // upper bound to 1.5 times the max group size
    Ncmax_scale = MAX(1.5*Ngmax,Ncmax/cell_scale);

    if (np <= Ncmax_scale) continue;

    // create particle list

    ip = cinfo[icell].first;
    n = 0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = ipart->weight;
      if(isw > 0) plist[n++] = ip;
      ip = next[ip];
    }

    gbuf = 0;
    total_iter = 0;
    while (1) {
      nold = n;
      group_bt(plist,n);

      // recreate particle list after reduction

      ip = cinfo[icell].first;
      n = 0;
      while (ip >= 0) {
        ipart = &particles[ip];
        isw = ipart->weight;
        if(isw > 0) plist[n++] = ip;
        ip = next[ip];
      }

      if (n < Ncmax_scale) break;

      // if no particles reduced, increase group size
      // if less than 20% particles reduced, increase group size
      
      if (n == nold || n/nold < 0.20) gbuf += 1;

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
void Collide::group_bt(int *plist_leaf, int np)
{
  Particle::OnePart *ipart;
  Particle::OnePart *particles = particle->particles;
  Particle::Species *species = particle->species;

  // ignore groups which have too few particles

  if (np <= Ngmin) return;

  // compute stress tensor since it's needed for
  // .. further dividing and reduction

  double gsum, msum;
  double mV[3], mVV[3][3], mVVV[3][3];
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
  for (int p = 0; p < np; p++) {
    ipart = &particles[plist_leaf[p]];
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
      for (int j = 0; j < 3; j++) {
        mVV[i][j] += (pmsw*vp[i]*vp[j]);
        mVVV[i][j] += (pmsw*vp[i]*vp[j]*vp[j]);
        uVV[i][j] += (vp[i]*vp[j]);
      }
    }
  }

  // mean velocity

	double V[3];
  for (int i = 0; i < 3; i++) V[i] = mV[i]/msum;

  // stress tensor

  double pij[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      pij[i][j] = mVV[i][j] - mV[i]*mV[j]/msum;

  // if group is small enough, merge the particles

  if (np <= Ngmax+gbuf) {

    // temperature
    double T = (pij[0][0] + pij[1][1] + pij[2][2])/
      (3.0 * gsum * update->boltz);

    // heat flux
    double Vsq = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
    double h,h1,h2,q[3];
    int i1,i2;
    for (int i = 0; i < 3; i++) {
      if (i == 0) {
        i1 = 1;
        i2 = 2;
      } else if (i == 1) {
        i1 = 2;
        i2 = 0;
      } else {
        i1 = 0;
        i2 = 1;
      }

      h  = mVVV[i][i] - 3.0*mV[i]*mVV[i][i]/msum +
           2.0*mV[i]*mV[i]*mV[i]/msum/msum;
      h1 = mVVV[i][i1] - 2.0*mVV[i][i1]*mV[i1]/msum -
           mV[i]*mVV[i1][i1]/msum + 2.0*mV[i]*mV[i1]*mV[i1]/msum/msum;
      h2 = mVVV[i][i2] - 2.0*mVV[i][i2]*mV[i2]/msum -
           mV[i]*mVV[i2][i2]/msum + 2.0*mV[i]*mV[i2]*mV[i2]/msum/msum;
      // normalized later in RW reduction (no fnum needed)
      q[i] = (h + h1 + h2) * 0.5;
    }

    // negative temp from big difference in particle weights
    if (T < 0) {
      T = 0.0;
      q[0] = q[1] = q[2] = 0.0;
    }

    // scale values to be consistent with definitions in
    // .. stochastic numerics book

    T *= update->boltz/mass;
    for (int i = 0; i < 3; i++) {
      q[i] /= mass;
      for(int j = 0; j < 3; j++) pij[i][j] /= mass;
    }

    // reduce based on type
    if (reduction_type == ENERGY) {
      reduce(plist_leaf, np, gsum, V, T, Erot);
    } else if (reduction_type == HEAT) {
      reduce(plist_leaf, np, gsum, V, T, Erot, q);
    } else if (reduction_type == STRESS) {
      reduce(plist_leaf, np, gsum, V, T, Erot, q, pij);
    }

  // group still too large so divide further

  } else {

    // Compute covariance matrix

    double Rij[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        if (weighted) Rij[i][j] = pij[i][j]/gsum;
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
    int pid, pidL[np], pidR[np];
    int npL, npR;
    npL = npR = 0;
    for (int i = 0; i < np; i++) {
      pid = plist_leaf[i];
      ipart = &particles[pid];
      if (MathExtra::dot3(ipart->v,maxevec) < center)
        pidL[npL++] = pid;
      else
        pidR[npR++] = pid;
    }

    if(npL < 1) error->all(FLERR,"No particles in left group");
    if(npR < 1) error->all(FLERR,"No particles in right group");

    if(npL > Ngmin) group_bt(pidL,npL);
    if(npR > Ngmin) group_bt(pidR,npR);
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using energy scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  ip = np * random->uniform();
  jp = np * random->uniform();
  while (ip == jp) jp = np * random->uniform();

  ipart = &particles[pleaf[ip]];
  jpart = &particles[pleaf[jp]];

  // find direction of velocity wrt CoM frame

  double theta = 2.0 * 3.14159 * random->uniform();
  double phi = acos(1.0 - 2.0 * random->uniform());
  double uvec[3];
  uvec[0] = sin(phi) * cos(theta);
  uvec[1] = sin(phi) * sin(theta);
  uvec[2] = cos(phi);

  // set reduced particle velocities

  double sqT = sqrt(3.0*T);
  for (int d = 0; d < 3; d++) {
    ipart->v[d] = V[d] + sqT*uvec[d];
    jpart->v[d] = V[d] - sqT*uvec[d];
  }

  // set reduced particle rotational energies

  ipart->erot = Erot/rho;
  jpart->erot = Erot/rho;

  // set reduced particle weights

  ipart->weight = rho*0.5;
  jpart->weight = rho*0.5;

  // delete other particles

  for (int i = 0; i < np; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[pleaf[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using heat flux scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot, double *q)
{

  // reduced particles chosen as first two in group

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  int ip, jp;
  ip = np * random->uniform();
  jp = np * random->uniform();
  while (ip == jp) jp = np * random->uniform();

  ipart = &particles[pleaf[ip]];
  jpart = &particles[pleaf[jp]];

  // precompute

  double sqT = sqrt(3.0*T);
  double qmag = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
  double qge = qmag / (rho * pow(sqT,3.0));
  double itheta = qge + sqrt(1.0 + qge*qge);
  double alpha = sqT*itheta;
  double beta = sqT/itheta;

  // find direction of velocity wrt CoM frame

  double uvec[3];
  if (qmag < SMALL) {
    for (int d = 0; d < 3; d++) {
      double A = sqrt(-log(random->uniform()));
      double phi = 6.283185308 * random->uniform();
      if (random->uniform() < 0.5) uvec[d] = A * cos(phi);
      else uvec[d] = A * sin(phi);
    }
  } else 
    for (int d = 0; d < 3; d++) uvec[d] = q[d]/qmag;

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

  ipart->weight = isw;
  jpart->weight = jsw;

  if (isw != isw || isw <= 0.0 || jsw != jsw || jsw <= 0.0)
    error->one(FLERR,"bad weight");

  // delete other particles
  for (int i = 0; i < np; i++) {
    if (i == ip || i == jp) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[pleaf[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = pleaf[i];
  }

  return;
}

/* ----------------------------------------------------------------------
   Merge particles using stress scheme
------------------------------------------------------------------------- */
void Collide::reduce(int *pleaf, int np,
                     double rho, double *V, double T, double Erot,
                     double *q, double pij[3][3])
{

  // scale by fnum
  rho *= update->fnum;
  Erot *= update->fnum;
  for (int i = 0; i < 3; i++) {
    q[i] *= update->fnum;
    for (int j = 0; j < 3; j++)
      pij[i][j] *= update->fnum;
  }

  // find eigenpairs of stress tensor

  double eval[3], evec[3][3];
  int ierror = MathEigen::jacobi3(pij,eval,evec);

  // find number of non-zero eigenvalues

  int nK = 0;
  for (int i = 0; i < 3; i++) {
    if (fabs(eval[i]) >= SMALL && eval[i] > 0) {
      eval[nK] = eval[i];
      for (int d = 0; d < 3; d++) evec[nK][d] = evec[i][d];
      nK++;
    }
  }

  // iterate through nK pairs

  Particle::OnePart *particles = particle->particles;
  Particle::OnePart *ipart, *jpart;

  double qli, itheta;
  double isw, jsw;
  double uvec[3];

  for (int iK = 0; iK < nK; iK++) {

    // reduced particles chosen as first two

    ipart = &particles[pleaf[2*iK]];
    jpart = &particles[pleaf[2*iK+1]];

    qli = evec[0][iK]*q[0] + evec[1][iK]*q[1] + evec[2][iK]*q[2];
    if (qli < 0)
      for (int d = 0; d < 3; d++) evec[d][iK] *= -1.0;
    qli = fabs(qli);

    itheta = sqrt(rho) * qli / (sqrt(nK) * pow(eval[iK],1.5))
      + sqrt(1.0 + (rho*qli*qli)/(nK*pow(eval[iK],3.0)));

    // set reduced particle velocities

    for (int d = 0; d < 3; d++) {
      ipart->v[d] = V[d] + itheta*sqrt(nK*eval[iK]/rho)*evec[d][iK];
      jpart->v[d] = V[d] - 1.0/itheta*sqrt(nK*eval[iK]/rho)*evec[d][iK];
    }

    // set reduced particle weights

    isw = rho/(nK*(1.0+itheta*itheta));
    jsw = rho/nK - isw;

    // set reduced particle rotational energies

    ipart->erot = Erot/rho/nK;
    jpart->erot = Erot/rho/nK;

    // scale back weights

    ipart->weight = isw/update->fnum;
    jpart->weight = jsw/update->fnum;

  } // end nK
  
  // delete other particles
  for (int i = 0; i < np; i++) {
    if (i < 2*nK) continue;
    if (ndelete == maxdelete) {
      maxdelete += DELTADELETE;
      memory->grow(dellist,maxdelete,"collide:dellist");
    }
    ipart = &particles[pleaf[i]];
    ipart->weight = -1.0;
    dellist[ndelete++] = pleaf[i];
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
    n = 0;
    sw_max = 0.0;
    while (ip >= 0) {
      ipart = &particles[ip];
      isw = ipart->weight;
      if (isw > 0) sw_max = MAX(sw_max,isw);
      ip = next[ip];
      n++;
    }

    // delete tiny weights

    ip = cinfo[icell].first;
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
    }
  } // loop for cells

  return;
}

