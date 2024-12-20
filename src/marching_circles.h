/* ----------------------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.github.io
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPARTA directory.
------------------------------------------------------------------------- */

#ifndef SPARTA_MARCHING_CIRCLES_H
#define SPARTA_MARCHING_CIRCLES_H

#include "pointers.h"

namespace SPARTA_NS {

class MarchingCircles : protected Pointers {
 public:
  MarchingCircles(class SPARTA *, int);
  ~MarchingCircles() {}
  void invoke(double **, double ***, int *);
  double mindist;

 private:
  int ggroup;
  double thresh;

  double extrapolate(double, double, double, double);
};

}

#endif

/* ERROR/WARNING messages:

*/
