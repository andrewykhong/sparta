#include "compute.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_solid.h"
#include "grid.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "particle.h"
#include "random_mars.h"
#include "random_knuth.h"
#include "stdlib.h"
#include "string.h"
#include "update.h"

using namespace SPARTA_NS;
using namespace MathConst;

enum{NOFORCE,GREEN,LOTH,SINGH};
enum{SPHERE,DISC,CYLINDER,CUSTOM};

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
  double erot;
  int nrot;

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

    for (int d = 0; d < dim; d++)  um[d] = mv[d]/totalmass;

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
          erot = particles[ip].erot;
          nrot = species[ispecies].rotdof;

          for (int d = 0; d < dim; d++) c[d] = u[d]-up[d];
          if (dim == 2) c[2] = 0.0;
          cmag = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

          if (force_type == GREEN) {
            // mean thermal speed of outgoing particles
            cp = sqrt(2.0*update->boltz*Tp/mass);
            // energy exchange due to rotational energy only pertinent for 
            // ... diffuse contribution
            Eg += c_diff*cmag*(erot - (0.5*nrot)*update->boltz*Tp);

            double F[3], Q;
            if (shape == SPHERE) {
              csx = Rp*Rp*MY_PI;
              Fsphere(c, cmag, cp, F, Q);
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += 0.5*mass*csx*Q;
            } else if (shape == DISC) {
              csx = Rp*Rp*MY_PI;
              Fdisc(c, cmag, cp, theta, phi, F, Q);
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += 0.5*mass*csx*Q;
            } else if (shape == CYLINDER) {
              // disc part first
              csx = Rp*Rp*MY_PI;
              Fdisc(c, cmag, cp, theta, phi, F, Q);
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += 0.5*mass*csx*Q;

              // then half cylinder
              csx = Rp*Lp;
              Fcyl(c, cmag, cp, theta, phi, F, Q);
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += 0.5*mass*csx*Q;
            } else if (shape == CUSTOM) {
              Fcustom(c, cmag, cp, F, Q);
              // cross section included in previous calc.
              for (int d = 0; d < dim; d++) Fg[d] += mass*F[d];
              Eg += 0.5*mass*Q;
            }
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
  double erot;
  int nrot;

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


  } // end cells

  // number of samples
  nsample++;
  return;
}

/* ----------------------------------------------------------------------
   Sphere Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fsphere(const double *c, const double cmag, const double cp, double *F, double Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

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
    F[d] = (Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
  Q = 0.5*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);
  return;
}

/* ----------------------------------------------------------------------
   Disc Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fdisc(const double *c, const double cmag, const double cp, const double theta, const double phi, double *F, double Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

  // first facee
  double norm[3];
  norm[0] = cos(theta)*sin(phi);
  norm[1] = sin(theta)*sin(phi);
  norm[2] = cos(phi);
  if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
  if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
  if(abs(norm[2]) < 1e-16) norm[2] = 0.0;

  // cos(theta) should always be a positive value
  double orth[3]; // yhat
  orth[0] = cos(theta)*cos(phi);
  orth[1] = sin(theta)*sin(phi);
  orth[2] = -sin(phi);
  if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
  if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
  if(abs(orth[2]) < 1e-16) orth[2] = 0.0;

  double orth1[3]; // zhat
  orth1[0] = -sin(theta);
  orth1[1] = cos(theta);
  orth1[2] = 0.0;
  if(abs(orth1[0]) < 1e-16) orth1[0] = 0.0;
  if(abs(orth1[1]) < 1e-16) orth1[1] = 0.0;
  if(abs(orth1[2]) < 1e-16) orth1[2] = 0.0;

  // need opposite sign for one associated with normal
  double alpha = -(norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2])/cmag;

  // use other face if negative
  if (alpha < 0) {
    norm[0] = -norm[0];
    norm[1] = -norm[1];
    norm[2] = -norm[2];

    // only need to rotate around one of the in-plane vectors
    orth[0] = -orth[0];
    orth[1] = -orth[1];
    orth[2] = -orth[2];
  }

  // direction cosines for in-plane vectors
  double beta  = (orth[0]*c[0]+orth[1]*c[1]+orth[2]*c[2])/cmag;
  double gamma = (orth1[0]*c[0]+orth1[1]*c[1]+orth1[2]*c[2])/cmag;
  alpha = abs(alpha);

  for (int d = 0; d < 3; d++) {
    Fi[d] = cmag*alpha*c[d];
    Fs[d] = -cmag*cmag*alpha*(alpha*norm[d]+beta*orth[d]+gamma*orth1[d]);
    Fd[d] = -0.5*sqrt(MY_PI)*cp*cmag*alpha*norm[d];
    Fa[d] = -2.0/3.0*cmag*cmag*alpha*norm[d];
  }

  Qi = cmag*cmag*cmag*alpha;
  Qs = -cmag*cmag*cmag*alpha;
  Qd = -2.0*cp*cp*cmag*alpha;
  Qa = -cmag*cmag*cmag*alpha;

  for (int d = 0; d < dim; d++) 
    F[d] = (Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
  Q = 0.5*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);
  return;
}

/* ----------------------------------------------------------------------
   Cylinder Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fcyl(const double *c, const double cmag, const double cp, const double theta, const double phi, double *F, double Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

  // first facee
  double norm[3];
  norm[0] = cos(theta)*sin(phi);
  norm[1] = sin(theta)*sin(phi);
  norm[2] = cos(phi);
  if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
  if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
  if(abs(norm[2]) < 1e-16) norm[2] = 0.0;

  // cos(theta) should always be a positive value
  double orth[3]; // yhat
  orth[0] = cos(theta)*cos(phi);
  orth[1] = sin(theta)*sin(phi);
  orth[2] = -sin(phi);
  if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
  if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
  if(abs(orth[2]) < 1e-16) orth[2] = 0.0;

  double orth1[3]; // zhat
  orth1[0] = -sin(theta);
  orth1[1] = cos(theta);
  orth1[2] = 0.0;
  if(abs(orth1[0]) < 1e-16) orth1[0] = 0.0;
  if(abs(orth1[1]) < 1e-16) orth1[1] = 0.0;
  if(abs(orth1[2]) < 1e-16) orth1[2] = 0.0;

  // need opposite sign for one associated with normal
  double alpha = -(norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2])/cmag;
  double beta  = (orth[0]*c[0]+orth[1]*c[1]+orth[2]*c[2])/cmag;
  double gamma = (orth1[0]*c[0]+orth1[1]*c[1]+orth1[2]*c[2])/cmag;

  // cos(theta) should always be a positive value
  double cZ = norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2];
  double sq_c = sqrt(cmag*cmag-cZ*cZ);

  for (int d = 0; d < 3; d++) {
    Fi[d] = 2.0*sq_c*c[d];
    Fs[d] = 2.0/3.0*sq_c*(c[d]-4.0*cZ*norm[d]);
    Fd[d] = pow(MY_PI,1.5)/4.0*cp*(c[d]-cZ*norm[d]);;
    Fa[d] = MY_PI/3.0*cmag*(c[d]-cZ*norm[d]);;
  }

  Qi = cmag*cmag*sq_c;
  Qs = -cmag*cmag*sq_c;
  Qd = -2.0*cp*cp*sq_c;
  Qa = -cmag*cmag*sq_c;

  for (int d = 0; d < dim; d++) 
    F[d] = (Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
  Q = 0.5*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);
  return;
}

/* ----------------------------------------------------------------------
   Cylinder Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fcustom(const double *c, const double cmag, const double cp, double *F, double Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

  double norm[3], orth[3], orth1[3];
  double theta, phi, dS;  
  double alpha, beta, gamma;

  double chat[3];
  for (int d = 0; d < 3; d++) {
    chat[d] = c[d]/cmag;
    Fi[d] = Fs[d] = Fd[d] = Fa[d] = 0.0;
  }

  // iterate thru all surfs
  for (int isurf = 0; isurf < nsurfs; isurf++) {
    theta = Sn[isurf][0];
    phi = Sn[isurf][1];
    dS = Sn[isurf][2];

    // grab surface element normal
    norm[0] = cos(theta)*sin(phi);
    norm[1] = sin(theta)*sin(phi);
    norm[2] = cos(phi);
    if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
    if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
    if(abs(norm[2]) < 1e-16) norm[2] = 0.0;

    orth[0] = cos(theta)*cos(phi);
    orth[1] = sin(theta)*sin(phi);
    orth[2] = -sin(phi);
    if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
    if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
    if(abs(orth[2]) < 1e-16) orth[2] = 0.0;

    orth1[0] = -sin(theta);
    orth1[1] = cos(theta);
    orth1[2] = 0.0;
    if(abs(orth1[0]) < 1e-16) orth1[0] = 0.0;
    if(abs(orth1[1]) < 1e-16) orth1[1] = 0.0;
    if(abs(orth1[2]) < 1e-16) orth1[2] = 0.0;

    // direction cosines
    alpha = MathExtra::dot3(norm,chat);
    beta = MathExtra::dot3(orth,chat);
    gamma = MathExtra::dot3(orth1,chat);    

    for (int d = 0; d < 3; d++) {
      Fi[d] -= cmag*c[d]*MIN(alpha,0.0)*dS;
      Fs[d] -= cmag*cmag*MIN(alpha,0.0)*dS*(alpha*norm[d]-beta*orth[d]-gamma*orth1[d]);
      Fd[d] += cmag*cp*MIN(alpha,0.0)*dS*norm[d];
      Fa[d] += cmag*cmag*MIN(alpha,0.0)*dS*norm[d];
    }

    Qi -= alpha*pow(cmag,3.0)*dS;
    Qs += alpha*pow(cmag,3.0)*dS;
    Qd += alpha*cmag*cp*cp*dS;
    Qa += alpha*cmag*cmag*cmag*dS;
  }

  // includes cross section
  for (int d = 0; d < dim; d++) 
    F[d] = (Fi[d]+c_spec*Fs[d]+sqrt(3.14159)*0.5*c_diff*Fd[d]+2.0/3.0*c_adia*Fa[d]);
  Q = 0.5*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);

  return;
}


