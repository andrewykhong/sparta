/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com
   Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(create_isurf,CreateISurf)

#else

#ifndef SPARTA_CREATE_ISURF_H
#define SPARTA_CREATE_ISURF_H

#include "stdio.h"
#include "pointers.h"
#include "surf.h"

namespace SPARTA_NS {

class CreateISurf : protected Pointers {
 public:
  CreateISurf(class SPARTA *);
  virtual ~CreateISurf();
  virtual void command(int, char **);

 protected:
  int me,nprocs;

	// For generating implicit surfaces
	int ggroup;							  // group id for grid cells
	double thresh;            // lower threshold for corner values
  double corner[3];         // corners of grid group
  double xyzsize[3];        // size of lowest level cell (must be uniform grid)
  int nxyz[3], Nxyz;        // dimensions of grid
  double *cvalues;          // array of corner point values
  double *mvalues;          // minimum intersection value
  int *svalues;             // marks corners as in or out
  double **ivalues;         // point of intersection between corner points

  double **icvalues;        // corner values for Fix Ablate
  int *tvalues;             // vector of per grid cell surf types

  int aveFlag;              // flag for how corners in unknown cells are set
  double mind;              // minimum cell length
  double cin, cout;         // in and out corner values
  class FixAblate *ablate;  // ablate fix

	// functions to set corner values

	void set_corners();

  // sets corner values whose cells have surfaces
  void surface_edge2d();
  void surface_edge3d();
  // marks corners which have no surfaces
  void set_inout();
  // find remaining corners
  void cleanup();
  // detects intersection between surfaces and cell edges
  int corner_hit2d(double*, double*, Surf::Line*, double&, int&);
  int corner_hit3d(double*, double*, Surf::Tri*, double&, int&);
  // remove old surfaces
  void remove_old();
  // misc functions
  void corner2cell();
  int get_cxyz(int *, double *);
  int get_cell(int, int, int);
  int get_corner(int, int, int);
  int get_corner(double, double, double);
  double param2in(double, double);
  double param2out(double, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPARTA to see the offending line.

E: Could not find surf_modify surf-ID

Self-explanatory.

E: Could not find surf_modify sc-ID

Self-explanatory.

E: %d surface elements not assigned to a collision model

All surface elements must be assigned to a surface collision model via
the surf_modify command before a simulation is perforemd.

E: Reuse of surf_collide ID

A surface collision model ID cannot be used more than once.

E: Invalid surf_collide style

Self-explanatory.

*/
