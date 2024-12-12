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
enum{SPHERE,DISC,CYLINDER};

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
  double Lp, theta, phi; // particle length and direction of norm
  double c[3],cmag; // thermal velocity of gas and solid particles
  double *u,*up;  // velocity of gas and solid particles
  double mv[3], mvsq, um[3]; // drift velocity (zero if no charge or gravity) 
  double Fg[3],Eg; // force and energy as defined by Green function
  double totalmass; // for calculating temperature
  double mass,T,p; // gas particle mass, gas temperature and pressure
  double csx,cp; // solid particle cross section and thermal velocity
  double frat; // ratio of gas to solid weight
  double prefactor; // prefactor
  double Fi[3], Fs[3], Fd[3], Fa[3]; // prefactors for forces
  double Qi, Qs, Qd, Qa; // prefactors for heat

  double Fgtotal[3], Egtotal; // total force and energy change

  for (icell = 0; icell < nglocal; icell++) {

    // DEBUG
    if (reset_flag) {
      array_grid[icell][0] = 0.0;
      array_grid[icell][1] = 0.0;
      array_grid[icell][2] = 0.0;
      array_grid[icell][3] = 0.0;
      array_grid[icell][4] = 0.0;
      array_grid[icell][5] = 0.0;
      array_grid[icell][6] = 0.0;
    }

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
        mvsq += mass*(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
      }

      ip = next[ip];
    }

    // calculate macroscopic vars

    T = mvsq - (mv[0]*mv[0] + mv[1]*mv[1] + mv[2]*mv[2])/totalmass;
    T /= (3.0*update->boltz*(np-nsolid));
    p = update->fnum*(np-nsolid)/cinfo[icell].volume
        *update->boltz*T*cinfo[icell].weight;

    if (T < 0 || p < 0) error->one(FLERR,"Negative temp or pressure");

    //printf("n - %i / %i\n", np, nsolid);
    //printf("nden_g: %4.3e\n", update->fnum*(np-nsolid)/cinfo[icell].volume*particle->species[0].specwt);
    //printf("nden_s: %4.3e\n", update->fnum*(nsolid)/cinfo[icell].volume*particle->species[1].specwt);
    //printf("spectwt: %4.3e, %4.3e\n", particle->species[0].specwt,particle->species[1].specwt);
    //printf("p - drag: %4.3e\n", p);
    //printf("T - drag: %4.3e\n", T);
    //error->one(FLERR,"Ck");

    um[0] = mv[0]/totalmass;
    um[1] = mv[1]/totalmass;
    um[2] = mv[2]/totalmass;

    // calculate incident forces and heat flux

    for (is = 0; is < nsolid; is++) {
      sid = id[is];
      up = particles[sid].v;
      Rp = solid_array[sid][1];
      Tp = solid_array[sid][2];
      Lp = solid_array[sid][4];
      theta = solid_array[sid][5];
      phi = solid_array[sid][6];

      Fg[0] = Fg[1] = Fg[2] = Eg = 0.0;

      ip = cinfo[icell].first;
      while (ip >= 0) {
        ispecies = particles[ip].ispecies;
        if (ispecies != solid_species) {
          mass = species[ispecies].mass;

          // account for difference in species weight
          u = particles[ip].v;

          c[0] = u[0]-up[0];
          c[1] = u[1]-up[1];
          c[2] = u[2]-up[2];
          if (dim == 2) c[2] = 0.0;
          cmag = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

          if (force_type == GREEN) {
            // mean thermal speed of outgoing particles
            cp = sqrt(2.0*update->boltz*Tp/mass);

            if (shape == SPHERE) {
              csx = Rp*Rp*MY_PI;
              for (int d = 0; d < 3; d++) {
                // incident
                Fi[d] = cmag*c[d];
                // specular
                Fs[d] = 0.0;
                // diffuse
                Fd[d] = sqrt(MY_PI)/3.0*cp*c[d];
                // adiabatic
                Fa[d] = 4.0/9.0*cmag*c[d];
              }
              // heat
              Qi = cmag*cmag*cmag;
              Qs = -cmag*cmag*cmag;
              Qd = -2.0*cp*cp*cmag;
              Qa = -cmag*cmag*cmag;

              for (int d = 0; d < dim; d++) 
                Fg[d] += mass*csx*(Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
              Eg    += 0.5*mass*csx*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);

            // two contributions from each face
            // only need to consider face where c . x is negative
            } else if (shape == DISC || shape == CYLINDER) {
              csx = Rp*Rp*MY_PI;
              // first facee
              double norm[3];
              norm[0] = cos(theta)*sin(phi);
              norm[1] = sin(theta)*sin(phi);
              norm[2] = cos(phi);
              if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
              if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
              if(abs(norm[2]) < 1e-16) norm[2] = 0.0;

              // cos(theta) should always be a positive value
              double cos_tht = -(norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2])/cmag;
              double orth[3]; // yhat
              orth[0] = cos(theta)*cos(phi);
              orth[1] = sin(theta)*sin(phi);
              orth[2] = -sin(phi);
              if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
              if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
              if(abs(orth[2]) < 1e-16) orth[2] = 0.0;
              double sin_tht = sqrt(1.0 - cos_tht*cos_tht);

              // use other face if negative
              if (cos_tht < 0) {
                norm[0] = -norm[0];
                norm[1] = -norm[1];
                norm[2] = -norm[2];

                orth[0] = -orth[0];
                orth[1] = -orth[1];
                orth[2] = -orth[2];
              }
              cos_tht = abs(cos_tht);

              for (int d = 0; d < 3; d++) {
                Fi[d] = cmag*cos_tht*c[d];
                Fs[d] = -cmag*cmag*cos_tht*(cos_tht*norm[d]+sin_tht*orth[d]);
                Fd[d] = -0.5*sqrt(MY_PI)*cp*cmag*cos_tht*norm[d];
                Fa[d] = -2.0/3.0*cmag*cmag*cos_tht*norm[d];
              }

              Qi = cmag*cmag*cmag*cos_tht;
              Qs = -cmag*cmag*cmag*cos_tht;
              Qd = -2.0*cp*cp*cmag*cos_tht;
              Qa = -cmag*cmag*cmag*cos_tht;

              for (int d = 0; d < dim; d++) 
                Fg[d] += mass*csx*(Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
              Eg    += 0.5*mass*csx*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);

            } 

            // also add the half cylinder portion as well
            if (shape == CYLINDER) {
              csx = 2.0*Rp*Lp;
              // normal aligned with axis (zhat in paper)
              double norm[3];
              norm[0] = cos(theta)*sin(phi);
              norm[1] = sin(theta)*sin(phi);
              norm[2] = cos(phi);
              if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
              if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
              if(abs(norm[2]) < 1e-16) norm[2] = 0.0;

              // cos(theta) should always be a positive value
              double cZ = norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2];
              double sq_c = sqrt(cmag*cmag-cZ*cZ);

              for (int d = 0; d < 3; d++) {
                Fi[d] = sq_c*c[d];
                Fs[d] = sq_c/3.0*(c[d]-4.0*cZ*norm[d]);
                Fd[d] = pow(MY_PI,1.5)/8.0*cp*(c[d]-cZ*norm[d]);;
                Fa[d] = MY_PI/6.0*cmag*(c[d]-cZ*norm[d]);;
              }

              Qi = cmag*cmag*sq_c;
              Qs = -cmag*cmag*sq_c;
              Qd = -2.0*cp*cp*sq_c;
              Qa = -cmag*cmag*sq_c;

              for (int d = 0; d < dim; d++) 
                Fg[d] += mass*csx*(Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
              Eg += 0.5*mass*csx*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);
            }

            //Fg[0] += mass*cx*(F1*cmag+F2*cp);
            //Fg[1] += mass*cy*(F1*cmag+F2*cp);
            //if (dim == 3)
            //  Fg[2] += mass*cz*(F1*cmag+F2*cp);
            //Eg    += mass*cmag*Q1*(0.5*cmag*cmag-cp*cp);
            //Eg += cmag*Q1*(erot - (0.5*nrot)*update->boltz*Tp);
          }
        }
        ip = next[ip];
      } // end while

      prefactor = update->fnum / cinfo[icell].volume;
      Fg[0] *= prefactor;
      Fg[1] *= prefactor;
      Fg[2] *= prefactor;
      Eg    *= prefactor;

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

      // DEBUG
      if (reset_flag) {
        array_grid[icell][3] += Fg[0]/nsolid;
        array_grid[icell][4] += Fg[1]/nsolid;
        array_grid[icell][5] += Fg[2]/nsolid;
        array_grid[icell][6] += Eg/nsolid;
      }
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
