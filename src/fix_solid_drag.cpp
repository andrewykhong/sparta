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

  int i,icell,ip,is; // dummy indices
  int ispecies, sid; // species of particle i; particle indices of solid particles
  int np, nsolid; // number of total simulators; number of solid simulators
  double Rp,mp,Tp,csp; // particle radius, mass, temperature, and specific heat
  double cx,cy,cz,cmag; // thermal velocity of gas and solid particles
  double *u,*up;  // velocity of gas and solid particles
  double um[3], usq[3]; // drift velocity (zero if no charge or gravity) 
  double Fg[3],Qg; // force and heat flux as defined by Green function
  double totalmass; // for calculating temperature
  double mass,T,p; // gas particle mass, gas temperature and pressure
  double csx,cp; // solid particle cross section and thermal velocity
  double frat; // ratio of gas to solid weight
  double prefactor; // prefactor

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
    um[0] = um[1] = um[2] = 0.0;
    usq[0] = usq[1] = usq[2] = 0.0;
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
        um[0] += mass*u[0];
        um[0] += mass*u[1];
        um[0] += mass*u[2];
        usq[0] += mass*u[0]*u[0];
        usq[1] += mass*u[1]*u[1];
        usq[2] += mass*u[2]*u[2];
      }
      ip = next[ip];
    }

    um[0] /= totalmass;
    um[1] /= totalmass;
    um[2] /= totalmass;

    // calculate gas state

    T = 0;
    ip = cinfo[icell].first;
    while (ip >= 0) {
      ispecies = particles[ip].ispecies;
      if (ispecies != solid_species) {
        mass = species[ispecies].mass;
        u = particles[ip].v;

        cx = u[0]-um[0];
        cy = u[1]-um[1];
        cz = u[2]-um[2];
        T += mass*(cx*cx+cy*cy+cz*cz);
      }
      ip = next[ip];
    }

    // removed half from numer. and denom.

    T /= (3.0*update->boltz*(np-nsolid));
    p = (np-nsolid)/cinfo[icell].volume*update->boltz*T;

    // calculate incident forces and heat flux

    for (is = 0; is < nsolid; is++) {
      sid = id[is];
      up = particles[sid].v;
      Rp = solid_array[sid][1];
      mp = solid_array[sid][2];
      Tp = solid_array[sid][3];
      csp = solid_array[sid][4];

      Fg[0] = Fg[1] = Fg[2] = Qg = 0.0;

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
            Qg    += mass*cmag*Q1*(0.5*cmag*cmag-cp*cp);
            //Qg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          } else if (force_type == BURT) {
            cp = sqrt(update->boltz*Tp*mass);

            Fg[0] += cx*(mass*cmag+F2*cp);
            Fg[1] += cy*(mass*cmag+F2*cp);
            if (dim == 3)
              Fg[2] += cz*(mass*cmag+F2*cp);
            Qg += cmag*Q1*(0.5*mass*cmag*cmag-2.0*update->boltz*Tp);
            //Qg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          }
        }
        ip = next[ip];
      } // end while

      csx = Rp*Rp*MY_PI;
      prefactor = csx * update->fnum / cinfo[icell].volume * fnum_rat;
      Fg[0] *= prefactor;
      Fg[1] *= prefactor;
      Fg[2] *= prefactor;
      Qg    *= prefactor;

      // update per-particle forces

      solid_force[sid][0] = (solid_force[sid][0]*nsample+Fg[0])/(nsample+1.0);
      solid_force[sid][1] = (solid_force[sid][1]*nsample+Fg[1])/(nsample+1.0);
      solid_force[sid][2] = (solid_force[sid][2]*nsample+Fg[2])/(nsample+1.0);
      solid_force[sid][3] = (solid_force[sid][3]*nsample+Qg)/(nsample+1.0);

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
