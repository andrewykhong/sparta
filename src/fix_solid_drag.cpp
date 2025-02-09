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
  double F[3], Q;
  double erot;
  int nrot;
  int ngas;

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

    nsolid = ngas = 0;
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
        ngas++;
      }

      ip = next[ip];
    }

    // calculate macroscopic vars

    T = mvsq - (mv[0]*mv[0] + mv[1]*mv[1] + mv[2]*mv[2])/totalmass;
    T /= (3.0*update->boltz*ngas);
    p = update->fnum*ngas/cinfo[icell].volume
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

            Q = 0.0;
            F[0] = F[1] = F[2] = 0.0;
            if (shape == SPHERE) {
              Fsphere(c, cmag, cp, F, Q);
              csx = Rp*Rp*MY_PI;
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += mass*csx*Q;
              Eg += c_diff*cmag*(erot - (0.5*nrot)*update->boltz*Tp)*csx;
            } else if (shape == DISC) {
              Fdisc(c, cmag, cp, theta, phi, F, Q);
              csx = Rp*Rp*MY_PI;
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += mass*csx*Q;
              Eg += c_diff*cmag*(erot - (0.5*nrot)*update->boltz*Tp)*csx;
            } else if (shape == CYLINDER) {
              // disc part first
              Fdisc(c, cmag, cp, theta, phi, F, Q);
              csx = Rp*Rp*MY_PI;
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += mass*csx*Q;
              Eg += c_diff*cmag*(erot - (0.5*nrot)*update->boltz*Tp)*csx;

              // then half cylinder
              Fcyl(c, cmag, cp, theta, phi, F, Q);
              csx = 2.0*Rp*Lp;
              for (int d = 0; d < dim; d++) Fg[d] += mass*csx*F[d];
              Eg += mass*csx*Q;
              Eg += c_diff*cmag*(erot - (0.5*nrot)*update->boltz*Tp)*csx;
            } else if (shape == CUSTOM) {
              Fcustom(c, cmag, cp, F, Q);
              // cross section included in previous calc.
              for (int d = 0; d < dim; d++) Fg[d] += mass*F[d];
              Eg += mass*Q;
              Eg += c_diff*cmag*(erot - (0.5*nrot)*update->boltz*Tp)*csx;
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

      // DEBUG: Testing implementation
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
   Sphere Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fsphere(const double *c, const double cmag, const double cp, double *F, double &Q)
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
   Cylinder Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fcyl(const double *c, const double cmag, const double cp, const double theta, const double phi, double *F, double &Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

  // aligned with cylinder axis
  double norm[3];
  norm[0] = cos(theta)*sin(phi);
  norm[1] = sin(theta)*sin(phi);
  norm[2] = cos(phi);
  if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
  if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
  if(abs(norm[2]) < 1e-16) norm[2] = 0.0;
  double mag = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
  norm[0] /= mag;
  norm[1] /= mag;
  norm[2] /= mag;

  double orth[3]; // xhat
  orth[0] = cos(theta)*cos(phi);
  orth[1] = sin(theta)*cos(phi);
  orth[2] = -sin(phi);
  if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
  if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
  if(abs(orth[2]) < 1e-16) orth[2] = 0.0;
  mag = sqrt(orth[0]*orth[0]+orth[1]*orth[1]+orth[2]*orth[2]);
  orth[0] /= mag;
  orth[1] /= mag;
  orth[2] /= mag;

  double orth1[3]; // yhat
  orth1[0] = -sin(theta);
  orth1[1] = cos(theta);
  orth1[2] = 0.0;
  if(abs(orth1[0]) < 1e-16) orth1[0] = 0.0;
  if(abs(orth1[1]) < 1e-16) orth1[1] = 0.0;
  if(abs(orth1[2]) < 1e-16) orth1[2] = 0.0;
  mag = sqrt(orth1[0]*orth1[0]+orth1[1]*orth1[1]+orth1[2]*orth1[2]);
  orth1[0] /= mag;
  orth1[1] /= mag;
  orth1[2] /= mag;

  // randomly rotate using rodrigues axis-angle
  double rot_theta = random->uniform()*2.0*MY_PI;
  double ctht = cos(rot_theta);
  double octht = 1.0-ctht;
  double stht = sin(rot_theta);

  double qorth = norm[0]*orth[0]+norm[1]*orth[1]+norm[2]*orth[2];
  double qorth1 = norm[0]*orth1[0]+norm[1]*orth1[1]+norm[2]*orth1[2];

  double csx3[3];

  MathExtra::cross3(norm,orth,csx3);
  orth[0] = ctht*orth[0]+stht*csx3[0]+qorth*norm[0]*octht;
  orth[1] = ctht*orth[1]+stht*csx3[1]+qorth*norm[1]*octht;
  orth[2] = ctht*orth[2]+stht*csx3[2]+qorth*norm[2]*octht;

  MathExtra::cross3(norm,orth1,csx3);
  orth1[0] = ctht*orth1[0]+stht*csx3[0]+qorth1*norm[0]*octht;
  orth1[1] = ctht*orth1[1]+stht*csx3[1]+qorth1*norm[1]*octht;
  orth1[2] = ctht*orth1[2]+stht*csx3[2]+qorth1*norm[2]*octht;

  // cos(theta) should always be a positive value
  double cnorm = norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2];
  double corth = (c[0]*orth[0]+c[1]*orth[1]+c[2]*orth[2])+
                 (c[0]*orth1[0]+c[1]*orth1[1]+c[2]*orth1[2]);
  double sqccz = sqrt(cmag*cmag-cnorm*cnorm);

  for (int d = 0; d < 3; d++) {
    //Fi[d] = corth*c[d]; // ~ |c| (c.X) c
    Fi[d] = sqccz*c[d];
    //Fs[d] = 1.0/3.0*corth*(c[d]-4.0*cnorm*norm[d]); // ~ |c| (c.X) c
    Fs[d] = 1.0/3.0*sqccz*(c[d]-4.0*cnorm*norm[d]);
    Fd[d] = pow(MY_PI,1.5)/8.0*cp*(c[d]-cnorm*norm[d]);
    Fa[d] = MY_PI/6.0*cmag*(c[d]-cnorm*norm[d]);
  }

  Qi = cmag*cmag*sqccz;
  Qs = -cmag*cmag*sqccz;
  Qd = -2.0*cp*cp*sqccz;
  Qa = -cmag*cmag*sqccz;

  for (int d = 0; d < dim; d++) 
    F[d] = (Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
  Q = 0.5*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);
  return;
}

/* ----------------------------------------------------------------------
   Disc Green's Function
---------------------------------------------------------------------- */

void FixSolid::Fdisc(const double *c, const double cmag, const double cp, const double theta, const double phi, double *F, double &Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

  double norm[3]; // disc normal
  norm[0] = cos(theta)*sin(phi);
  norm[1] = sin(theta)*sin(phi);
  norm[2] = cos(phi);
  if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
  if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
  if(abs(norm[2]) < 1e-16) norm[2] = 0.0;
  double mag = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
  norm[0] /= mag;
  norm[1] /= mag;
  norm[2] /= mag;

  double orth[3]; // xhat
  orth[0] = cos(theta)*cos(phi);
  orth[1] = sin(theta)*cos(phi);
  orth[2] = -sin(phi);
  if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
  if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
  if(abs(orth[2]) < 1e-16) orth[2] = 0.0;
  mag = sqrt(orth[0]*orth[0]+orth[1]*orth[1]+orth[2]*orth[2]);
  orth[0] /= mag;
  orth[1] /= mag;
  orth[2] /= mag;

  double orth1[3]; // yhat
  orth1[0] = -sin(theta);
  orth1[1] = cos(theta);
  orth1[2] = 0.0;
  if(abs(orth1[0]) < 1e-16) orth1[0] = 0.0;
  if(abs(orth1[1]) < 1e-16) orth1[1] = 0.0;
  if(abs(orth1[2]) < 1e-16) orth1[2] = 0.0;
  mag = sqrt(orth1[0]*orth1[0]+orth1[1]*orth1[1]+orth1[2]*orth1[2]);
  orth1[0] /= mag;
  orth1[1] /= mag;
  orth1[2] /= mag;

  // randomly rotate using rodrigues axis-angle
  double rot_theta = random->uniform()*2.0*MY_PI;
  double ctht = cos(rot_theta);
  double octht = 1.0-ctht;
  double stht = sin(rot_theta);

  double qorth = norm[0]*orth[0]+norm[1]*orth[1]+norm[2]*orth[2];
  double qorth1 = norm[0]*orth1[0]+norm[1]*orth1[1]+norm[2]*orth1[2];

  double csx3[3];

  MathExtra::cross3(norm,orth,csx3);
  orth[0] = ctht*orth[0]+stht*csx3[0]+qorth*norm[0]*octht;
  orth[1] = ctht*orth[1]+stht*csx3[1]+qorth*norm[1]*octht;
  orth[2] = ctht*orth[2]+stht*csx3[2]+qorth*norm[2]*octht;

  MathExtra::cross3(norm,orth1,csx3);
  orth1[0] = ctht*orth1[0]+stht*csx3[0]+qorth1*norm[0]*octht;
  orth1[1] = ctht*orth1[1]+stht*csx3[1]+qorth1*norm[1]*octht;
  orth1[2] = ctht*orth1[2]+stht*csx3[2]+qorth1*norm[2]*octht;

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

void FixSolid::Fcustom(const double *c, const double cmag, const double cp, double *F, double &Q)
{
  double Fi[3], Fs[3], Fd[3], Fa[3];
  double Qi, Qs, Qd, Qa;

  double norm[3];
  double theta, phi, dS;

  // iterate over all surfs
  Fi[0] = Fi[1] = Fi[2] = 0.0;
  Fs[0] = Fs[1] = Fs[2] = 0.0;
  Fd[0] = Fd[1] = Fd[2] = 0.0;
  Fa[0] = Fa[1] = Fa[2] = 0.0;
  Qi = Qs = Qd = Qa = 0.0;

  for (int isurf = 0; isurf < nsurfs; isurf++)
  {
    theta = Sn[isurf][0];
    phi = Sn[isurf][1];
    dS = Sn[isurf][2];

    norm[0] = cos(theta)*sin(phi);
    norm[1] = sin(theta)*sin(phi);
    norm[2] = cos(phi);
    if(abs(norm[0]) < 1e-16) norm[0] = 0.0;
    if(abs(norm[1]) < 1e-16) norm[1] = 0.0;
    if(abs(norm[2]) < 1e-16) norm[2] = 0.0;
    double mag = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);
    norm[0] /= mag;
    norm[1] /= mag;
    norm[2] /= mag;

    double orth[3]; // xhat
    orth[0] = cos(theta)*cos(phi);
    orth[1] = sin(theta)*cos(phi);
    orth[2] = -sin(phi);
    if(abs(orth[0]) < 1e-16) orth[0] = 0.0;
    if(abs(orth[1]) < 1e-16) orth[1] = 0.0;
    if(abs(orth[2]) < 1e-16) orth[2] = 0.0;
    mag = sqrt(orth[0]*orth[0]+orth[1]*orth[1]+orth[2]*orth[2]);
    orth[0] /= mag;
    orth[1] /= mag;
    orth[2] /= mag;

    double orth1[3]; // yhat
    orth1[0] = -sin(theta);
    orth1[1] = cos(theta);
    orth1[2] = 0.0;
    if(abs(orth1[0]) < 1e-16) orth1[0] = 0.0;
    if(abs(orth1[1]) < 1e-16) orth1[1] = 0.0;
    if(abs(orth1[2]) < 1e-16) orth1[2] = 0.0;
    mag = sqrt(orth1[0]*orth1[0]+orth1[1]*orth1[1]+orth1[2]*orth1[2]);
    orth1[0] /= mag;
    orth1[1] /= mag;
    orth1[2] /= mag;

    // randomly rotate using rodrigues axis-angle
    double rot_theta = random->uniform()*2.0*MY_PI;
    double ctht = cos(rot_theta);
    double octht = 1.0-ctht;
    double stht = sin(rot_theta);

    double qorth = norm[0]*orth[0]+norm[1]*orth[1]+norm[2]*orth[2];
    double qorth1 = norm[0]*orth1[0]+norm[1]*orth1[1]+norm[2]*orth1[2];

    double csx3[3];

    MathExtra::cross3(norm,orth,csx3);
    orth[0] = ctht*orth[0]+stht*csx3[0]+qorth*norm[0]*octht;
    orth[1] = ctht*orth[1]+stht*csx3[1]+qorth*norm[1]*octht;
    orth[2] = ctht*orth[2]+stht*csx3[2]+qorth*norm[2]*octht;

    MathExtra::cross3(norm,orth1,csx3);
    orth1[0] = ctht*orth1[0]+stht*csx3[0]+qorth1*norm[0]*octht;
    orth1[1] = ctht*orth1[1]+stht*csx3[1]+qorth1*norm[1]*octht;
    orth1[2] = ctht*orth1[2]+stht*csx3[2]+qorth1*norm[2]*octht;

    // need opposite sign for one associated with normal
    double alpha = -(norm[0]*c[0]+norm[1]*c[1]+norm[2]*c[2])/cmag;
    // direction cosines for in-plane vectors
    double beta  = (orth[0]*c[0]+orth[1]*c[1]+orth[2]*c[2])/cmag;
    double gamma = (orth1[0]*c[0]+orth1[1]*c[1]+orth1[2]*c[2])/cmag;

    // use other face if negative
    if (alpha > 0) {
      for (int d = 0; d < 3; d++) {
        Fi[d] += cmag*alpha*c[d]*dS;
        Fs[d] += cmag*cmag*alpha*(alpha*norm[d]+beta*orth[d]+gamma*orth1[d])*dS;
        Fd[d] -= 0.5*sqrt(MY_PI)*cp*cmag*alpha*norm[d]*dS;
        Fa[d] -= 2.0/3.0*cmag*cmag*alpha*norm[d]*dS;
      }

      Qi += cmag*cmag*cmag*alpha*dS;
      Qs *= cmag*cmag*cmag*alpha*dS;
      Qd -= 2.0*cp*cp*cmag*alpha*dS;
      Qa -= cmag*cmag*cmag*alpha*dS;
    }
  }

  for (int d = 0; d < dim; d++) 
    F[d] = (Fi[d]+c_spec*Fs[d]+c_diff*Fd[d]+c_adia*Fa[d]);
  Q = 0.5*(Qi+c_spec*Qs+c_diff*Qd+c_adia*Qa);

  return;
}


