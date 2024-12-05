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

using namespace SPARTA_NS;
using namespace MathConst;

enum{NOFORCE,GREEN,BURT,LOTH,SINGH};            // type of solid particle force

#define INVOKED_PER_GRID 16

/* ----------------------------------------------------------------------
   Momentum and heat fluxes based on Green's functions
   Assumes free molecular
---------------------------------------------------------------------- */

void FixSolid::update_Fq_fm()
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

  int i,j,k,icell,ip,is; // dummy indices
  int ispecies, sid; // species of particle i; particle indices of solid particles
  int np, nsolid; // number of total simulators; number of solid simulators
  double Rp,mp,Tp,csp; // particle radius, mass, temperature, and specific heat
  double cx,cy,cz,cmag; // thermal velocity of gas and solid particles
  double *u,*up;  // velocity of gas and solid particles
  double mv[3], mvsq, um[3]; // drift velocity (zero if no charge or gravity) 
  double Fg[3],Eg; // force and energy as defined by Green function
  double totalmass; // for calculating temperature
  double mass,T,p; // gas particle mass, gas temperature and pressure
  double csx,cp; // solid particle cross section and thermal velocity
  double frat; // ratio of gas to solid weight
  double prefactor; // prefactor

  double Fgtotal[3], Egtotal; // total force and energy change

  /*
  modify->clearstep_compute();

  for (i = 0; i < 2; i++) {
    Compute *compute = modify->compute[value2index[i]];
    if (!(compute->invoked_flag & INVOKED_PER_GRID)) {
      compute->compute_per_grid();
      compute->invoked_flag |= INVOKED_PER_GRID;
    }

    if(post_process[i])
      // index, nsample, etally, emap, vec, nstride
      compute->post_process_grid(argindex[i],1,NULL,NULL,NULL,2);

    if (argindex[i] == 0) {
      double *cvec = compute->vector_grid;
      for (j = 0; j < nglocal; j++) cell_Tp[j][i] = cvec[j];
    } else {
      double **carray = compute->array_grid;
      k = argindex[i] - 1;
      for (j = 0; j < nglocal; j++) cell_Tp[j][i] = carray[j][k];
    }

  }
  error->one(FLERR,"ck");
  */

  for (icell = 0; icell < nglocal; icell++) {

    np = cinfo[icell].count;
    if (np <= 1) continue;

    // only need ids

    if (np > npmax) {
      while (np > npmax) npmax += DELTAPART;
      memory->destroy(id);
      memory->create(id,npmax,"soliddrag:id");
    }

    // get solid particle velocities

    nsolid = 0;
    mv[0] = mv[1] = mv[2] = mvsq = 0.0;
    totalmass = 0.0;

    // drift belocity defined as F / beta where
    // ... F is the nondrag force?

    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      if (ispecies == solid_species) id[nsolid++] = ip;
      else {
        mass = species[ispecies].mass;
        u = particles[ip].v;
        totalmass += mass;
        mv[0] += mass*u[0];
        mv[1] += mass*u[1];
        mv[2] += mass*u[2];
        mvsq = mass*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
      }
      ip = next[ip];
    }

    // calculate macroscopic vars

    T = mvsq - (mv[0]*mv[0] + mv[1]*mv[1] + mv[2]*mv[2])/totalmass;
    p = T*update->fnum/3.0/cinfo[icell].volume*cinfo[icell].weight;
    T /= (3.0*update->boltz*(np-nsolid));
    um[0] = mv[0]/totalmass;
    um[1] = mv[1]/totalmass;
    um[2] = mv[2]/totalmass;

    // calculate incident forces and heat flux

    Fgtotal[0] = Fgtotal[1] = Fgtotal[2] = Egtotal = 0.0;
    for (is = 0; is < nsolid; is++) {
      sid = id[is];
      up = particles[sid].v;
      Rp = solid_array[sid][1];
      mp = solid_array[sid][2];
      Tp = solid_array[sid][3];
      csp = solid_array[sid][4];

      Fg[0] = Fg[1] = Fg[2] = Eg = 0.0;

      ip = cinfo[icell].first;
      while (ip >= 0) {
        ispecies = particles[ip].ispecies;
        if (ispecies != solid_species) {
          mass = species[ispecies].mass;

          // account for difference in species weight
          u = particles[ip].v;

          cx = u[0]-up[0];
          cy = u[1]-up[1];
          cz = u[2]-up[2];
          if (dim == 2) cz = 0.0;
          cmag = sqrt(cx*cx+cy*cy+cz*cz);

          if (force_type == GREEN) {
            // mean thermal speed of outgoing particles
            cp = sqrt(2.0*update->boltz*Tp/mass);

            Fg[0] += mass*cx*(F1*cmag+F2*cp);
            Fg[1] += mass*cy*(F1*cmag+F2*cp);
            if (dim == 3)
              Fg[2] += mass*cz*(F1*cmag+F2*cp);
            Eg    += mass*cmag*Q1*(0.5*cmag*cmag-cp*cp);
            //Eg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          } else if (force_type == BURT) {
            cp = sqrt(update->boltz*Tp*mass);

            Fg[0] += cx*(mass*cmag+F2*cp);
            Fg[1] += cy*(mass*cmag+F2*cp);
            if (dim == 3)
              Fg[2] += cz*(mass*cmag+F2*cp);
            Eg += cmag*Q1*(0.5*mass*cmag*cmag-2.0*update->boltz*Tp);
            //Eg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          }
        }
        ip = next[ip];
      } // end while

      csx = Rp*Rp*MY_PI;
      prefactor = csx * update->fnum / cinfo[icell].volume * fnum_rat;
      Fg[0] *= prefactor;
      Fg[1] *= prefactor;
      Fg[2] *= prefactor;
      Eg    *= prefactor;

      // keep track of total force and energy change in cell

      /*if (conserve_flag) {
        Fq_grid[icell][0] += Fg[0]*update->dt;
        Fq_grid[icell][1] += Fg[1]*update->dt;
        Fq_grid[icell][2] += Fg[2]*update->dt;
        Fq_grid[icell][3] += Eg*update->dt;
      }*/

      // update per-particle forces

      solid_force[sid][0] = (solid_force[sid][0]*nsample+Fg[0])/(nsample+1.0);
      solid_force[sid][1] = (solid_force[sid][1]*nsample+Fg[1])/(nsample+1.0);
      solid_force[sid][2] = (solid_force[sid][2]*nsample+Fg[2])/(nsample+1.0);
      solid_force[sid][3] = (solid_force[sid][3]*nsample+Eg)/(nsample+1.0);

      solid_bulk[sid][0] = (solid_bulk[sid][0]*nsample+um[0])/(nsample+1.0);
      solid_bulk[sid][1] = (solid_bulk[sid][1]*nsample+um[1])/(nsample+1.0);
      solid_bulk[sid][2] = (solid_bulk[sid][2]*nsample+um[2])/(nsample+1.0);
      solid_bulk[sid][3] = (solid_bulk[sid][3]*nsample+T)/(nsample+1.0);
      solid_bulk[sid][4] = (solid_bulk[sid][4]*nsample+p)/(nsample+1.0);

    } // end solid loop

  } // end cells

  // number of samples
  nsample++;

}

/* ----------------------------------------------------------------------
   Empirical drag and heat flux laws
---------------------------------------------------------------------- */

void FixSolid::update_Fq_emp()
{
  error->all(FLERR,"Empirical drag laws not included yet");
  return;
}
