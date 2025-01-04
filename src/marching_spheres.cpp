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

#include "math.h"
#include "math_extra.h"
#include "string.h"
#include "marching_spheres.h"
#include "grid.h"
#include "surf.h"
#include "irregular.h"
#include "lookup_table.h"
#include "geometry.h"
#include "my_page.h"
#include "memory.h"
#include "error.h"

// DEBUG
#include "update.h"

using namespace SPARTA_NS;

// prototype for non-class function

int compare_indices(const void *, const void *);

enum{UNKNOWN,OUTSIDE,INSIDE,OVERLAP};           // several files
enum{NCHILD,NPARENT,NUNKNOWN,NPBCHILD,NPBPARENT,NPBUNKNOWN,NBOUND};  // Grid

#define DELTA 128
#define BIG 1.0e20
#define EPSILON 1.0e-16

/* ---------------------------------------------------------------------- */

MarchingSpheres::MarchingSpheres(SPARTA *sparta, int ggroup_caller) :
  Pointers(sparta)
{
  MPI_Comm_rank(world,&me);

  ggroup = ggroup_caller;
}

/* ----------------------------------------------------------------------
   create 2d implicit surfs from grid point values
   follows https://en.wikipedia.org/wiki/Marching_squares
   see 2 sections: Basic algorithm and Disambiguation of saddle points
     treating open circles as flow volume, solid circles as material
     NOTE: Wiki page numbers points counter-clockwise
           SPARTA numbers them in x, then in y
           so bit2 and bit3 are swapped below
           this gives case #s here consistent with Wiki page
   process each grid cell independently
   4 corner points open/solid -> 2^4 = 16 cases
   cases infer 0,1,2 line segments in each grid cell
   order 2 points in each line segment to give normal into flow volume
   treat two saddle point cases (my 9,6) (Wiki 5,10)
     based on ave value at cell center
------------------------------------------------------------------------- */

void MarchingSpheres::invoke(double **cvalues, double ***mvalues, int *svalues, int **mcflags)
{
  int i,j,ipt,isurf,nsurf,icase,which;
  surfint surfID;
  surfint *ptr;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  MyPage<surfint> *csurfs = grid->csurfs;
  int nglocal = grid->nlocal;
  int groupbit = grid->bitmask[ggroup];
  double eps = 0.001;

  bigint maxsurfID = 0;
  if (sizeof(surfint) == 4) maxsurfID = MAXSMALLINT;
  if (sizeof(surfint) == 8) maxsurfID = MAXBIGINT;

  for (int icell = 0; icell < nglocal; icell++) {
    if (!(cinfo[icell].mask & groupbit)) continue;
    if (cells[icell].nsplit <= 0) continue;
    lo = cells[icell].lo;
    hi = cells[icell].hi;


    if (cvalues) {
      v000 = cvalues[icell][0];
      v001 = cvalues[icell][1];
      v010 = cvalues[icell][2];
      v011 = cvalues[icell][3];
      v100 = cvalues[icell][4];
      v101 = cvalues[icell][5];
      v110 = cvalues[icell][6];
      v111 = cvalues[icell][7];

      // vertex locations
      i0  = extrapolate(v000,v001,lo[0],hi[0]);
      i1  = extrapolate(v001,v011,lo[1],hi[1]);
      i2  = extrapolate(v010,v011,lo[0],hi[0]);
      i3  = extrapolate(v000,v010,lo[1],hi[1]);

      i4  = extrapolate(v100,v101,lo[0],hi[0]);
      i5  = extrapolate(v101,v111,lo[1],hi[1]);
      i6  = extrapolate(v110,v111,lo[0],hi[0]);
      i7  = extrapolate(v100,v110,lo[1],hi[1]);

      i8  = extrapolate(v000,v100,lo[2],hi[2]);
      i9  = extrapolate(v001,v101,lo[2],hi[2]);
      i10 = extrapolate(v011,v111,lo[2],hi[2]);
      i11 = extrapolate(v010,v110,lo[2],hi[2]);

    } else {
      // temporarily store all inner values

      for (i = 0; i < 8; i++)
        for (j = 0; j < 6; j++)
          inval[i][j] = mvalues[icell][i][j];

      // use averages for now

      v000 = v001 = v010 = v011 = 0.0;
      v100 = v101 = v110 = v111 = 0.0;

      // ordering in cvalues different from loop up table
      // manually change for consistency

      for (j = 0; j < 6; j++) {
        v000 += inval[0][j];
        v001 += inval[1][j];
        v010 += inval[2][j];
        v011 += inval[3][j];
        v100 += inval[4][j];
        v101 += inval[5][j];
        v110 += inval[6][j];
        v111 += inval[7][j];
      }

      v000 /= 6.0;
      v001 /= 6.0;
      v010 /= 6.0;
      v011 /= 6.0;
      v100 /= 6.0;
      v101 /= 6.0;
      v110 /= 6.0;
      v111 /= 6.0;

      i0  = extrapolate(inval[0][1],inval[1][0],lo[0],hi[0]);
      i1  = extrapolate(inval[1][3],inval[3][2],lo[1],hi[1]);
      i2  = extrapolate(inval[2][1],inval[3][0],lo[0],hi[0]);
      i3  = extrapolate(inval[0][3],inval[2][2],lo[1],hi[1]);

      i4  = extrapolate(inval[4][1],inval[5][0],lo[0],hi[0]);
      i5  = extrapolate(inval[5][3],inval[7][2],lo[1],hi[1]);
      i6  = extrapolate(inval[6][1],inval[7][0],lo[0],hi[0]);
      i7  = extrapolate(inval[4][3],inval[6][2],lo[1],hi[1]);

      i8  = extrapolate(inval[0][5],inval[4][4],lo[2],hi[2]);
      i9  = extrapolate(inval[1][5],inval[5][4],lo[2],hi[2]);
      i10 = extrapolate(inval[3][5],inval[7][4],lo[2],hi[2]);
      i11 = extrapolate(inval[2][5],inval[6][4],lo[2],hi[2]);

    }

    // for ambiguity tests
    v000iso = v000/255.0;
    v001iso = v001/255.0;
    v010iso = v010/255.0;
    v011iso = v011/255.0;
    v100iso = v100/255.0;
    v101iso = v101/255.0;
    v110iso = v110/255.0;
    v111iso = v111/255.0;

    // make bits 2, 3, 6 and 7 consistent with Lewiner paper (see NOTE above)

    bit0 = v000 <= thresh ? 0 : 1;
    bit1 = v001 <= thresh ? 0 : 1;
    bit2 = v011 <= thresh ? 0 : 1;
    bit3 = v010 <= thresh ? 0 : 1;
    bit4 = v100 <= thresh ? 0 : 1;
    bit5 = v101 <= thresh ? 0 : 1;
    bit6 = v111 <= thresh ? 0 : 1;
    bit7 = v110 <= thresh ? 0 : 1;

    which = (bit7 << 7) + (bit6 << 6) + (bit5 << 5) + (bit4 << 4) +
      (bit3 << 3) + (bit2 << 2) + (bit1 << 1) + bit0;

    // icase = case of the active cube in [0..15]

    icase = cases[which][0];
    config = cases[which][1];
    subconfig = 0;

    switch (icase) {
    case  0:
      nsurf = 0;
      break;

    case  1:
      nsurf = add_triangle(tiling1[config], 1);
      break;

    case  2:
      nsurf = add_triangle(tiling2[config], 2);
      break;

    case  3:
      if (test_face(test3[config]))
        nsurf = add_triangle(tiling3_2[config], 4); // 3.2
      else
        nsurf = add_triangle(tiling3_1[config], 2); // 3.1
      break;

    case  4:
      if (test_interior())
        nsurf = add_triangle(tiling4_1[config], 2); // 4.1.1
      else
        nsurf = add_triangle(tiling4_2[config], 6); // 4.1.2
      break;

    case  5:
      nsurf = add_triangle(tiling5[config], 3);
      break;

    case  6:
      if (test_face(test6[config][0]))
        nsurf = add_triangle(tiling6_2[config], 5); // 6.2
      else {
        if (test_interior())
          nsurf = add_triangle(tiling6_1_1[config], 3); // 6.1.1
        else {
          nsurf = add_triangle(tiling6_1_2[config], 9); // 6.1.2
        }
      }
      break;

    case  7:
      if (test_face(test7[config][0])) subconfig +=  1;
      if (test_face(test7[config][1])) subconfig +=  2;
      if (test_face(test7[config][2])) subconfig +=  4;
      switch (subconfig) {
      case 0:
        nsurf = add_triangle(tiling7_1[config], 3); break;
      case 1:
        nsurf = add_triangle(tiling7_2[config][0], 5); break;
      case 2:
        nsurf = add_triangle(tiling7_2[config][1], 5); break;
      case 3:
        nsurf = add_triangle(tiling7_3[config][0], 9); break;
      case 4:
        nsurf = add_triangle(tiling7_2[config][2], 5); break;
      case 5:
        nsurf = add_triangle(tiling7_3[config][1], 9); break;
      case 6:
        nsurf = add_triangle(tiling7_3[config][2], 9); break;
      case 7:
        if (test_interior())
          nsurf = add_triangle(tiling7_4_2[config], 9);
        else
          nsurf = add_triangle(tiling7_4_1[config], 5);
        break;
      };
      break;

    case  8:
      nsurf = add_triangle(tiling8[config], 2);
      break;

    case  9:
      nsurf = add_triangle(tiling9[config], 4);
      break;

    case 10:
      if (test_face(test10[config][0])) {
        if (test_face(test10[config][1]))
          nsurf = add_triangle(tiling10_1_1_[config], 4); // 10.1.1
        else {
          nsurf = add_triangle(tiling10_2[config], 8); // 10.2
        }
      } else {
        if (test_face(test10[config][1])) {
          nsurf = add_triangle(tiling10_2_[config], 8); // 10.2
        } else {
          if (test_interior())
            nsurf = add_triangle(tiling10_1_1[config], 4); // 10.1.1
          else
            nsurf = add_triangle(tiling10_1_2[config], 8); // 10.1.2
        }
      }
      break;

    case 11:
      nsurf = add_triangle(tiling11[config], 4);
      break;

    case 12:
      if (test_face(test12[config][0])) {
        if (test_face(test12[config][1]))
          nsurf = add_triangle(tiling12_1_1_[config], 4); // 12.1.1
        else {
          nsurf = add_triangle(tiling12_2[config], 8); // 12.2
        }
      } else {
        if (test_face(test12[config][1])) {
          nsurf = add_triangle(tiling12_2_[config], 8); // 12.2
        } else {
          if (test_interior())
            nsurf = add_triangle(tiling12_1_1[config], 4); // 12.1.1
          else
            nsurf = add_triangle(tiling12_1_2[config], 8); // 12.1.2
        }
      }
      break;

    case 13:
      if (test_face(test13[config][0])) subconfig +=  1;
      if (test_face(test13[config][1])) subconfig +=  2;
      if (test_face(test13[config][2])) subconfig +=  4;
      if (test_face(test13[config][3])) subconfig +=  8;
      if (test_face(test13[config][4])) subconfig += 16;
      if (test_face(test13[config][5])) subconfig += 32;

      switch (subconfig13[subconfig]) {
      case 0:/* 13.1 */
        nsurf = add_triangle(tiling13_1[config], 4); break;

      case 1:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][0], 6); break;
      case 2:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][1], 6); break;
      case 3:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][2], 6); break;
      case 4:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][3], 6); break;
      case 5:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][4], 6); break;
      case 6:/* 13.2 */
        nsurf = add_triangle(tiling13_2[config][5], 6); break;

      case 7:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][0], 10); break;
      case 8:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][1], 10); break;
      case 9:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][2], 10); break;
      case 10:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][3], 10); break;
      case 11:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][4], 10); break;
      case 12:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][5], 10); break;
      case 13:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][6], 10); break;
      case 14:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][7], 10); break;
      case 15:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][8], 10); break;
      case 16:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][9], 10); break;
      case 17:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][10], 10); break;
      case 18:/* 13.3 */
        nsurf = add_triangle(tiling13_3[config][11], 10); break;

      case 19:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][0], 12); break;
      case 20:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][1], 12); break;
      case 21:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][2], 12); break;
      case 22:/* 13.4 */
        nsurf = add_triangle(tiling13_4[config][3], 12); break;

      case 23:/* 13.5 */
        subconfig = 0;
        if (test_interior())
          nsurf = add_triangle(tiling13_5_1[config][0], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][0], 10);
        break;

      case 24:/* 13.5 */
        subconfig = 1;
        if (test_interior())
          nsurf = add_triangle(tiling13_5_1[config][1], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][1], 10);
        break;

      case 25:/* 13.5 */
        subconfig = 2;
        if (test_interior())
          nsurf = add_triangle(tiling13_5_1[config][2], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][2], 10);
        break;

      case 26:/* 13.5 */
        subconfig = 3;
        if (test_interior())
          nsurf = add_triangle(tiling13_5_1[config][3], 6);
        else
          nsurf = add_triangle(tiling13_5_2[config][3], 10);
        break;

      case 27:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][0], 10); break;
      case 28:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][1], 10); break;
      case 29:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][2], 10); break;
      case 30:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][3], 10); break;
      case 31:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][4], 10); break;
      case 32:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][5], 10); break;
      case 33:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][6], 10); break;
      case 34:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][7], 10); break;
      case 35:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][8], 10); break;
      case 36:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][9], 10); break;
      case 37:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][10], 10); break;
      case 38:/* 13.3 */
        nsurf = add_triangle(tiling13_3_[config][11], 10); break;

      case 39:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][0], 6); break;
      case 40:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][1], 6); break;
      case 41:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][2], 6); break;
      case 42:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][3], 6); break;
      case 43:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][4], 6); break;
      case 44:/* 13.2 */
        nsurf = add_triangle(tiling13_2_[config][5], 6); break;

      case 45:/* 13.1 */
        nsurf = add_triangle(tiling13_1_[config], 4); break;

      default:
        error->one(FLERR,"Marching cubes - impossible case 13");
      }
      break;

    case 14:
      nsurf = add_triangle(tiling14[config], 4);
      break;
    };

    // store 4 MC labels for FixAblate caller

    mcflags[icell][0] = icase;
    mcflags[icell][1] = config;
    mcflags[icell][2] = subconfig;
    mcflags[icell][3] = nsurf;

    // populate Grid and Surf data structs
    // points will be duplicated, not unique
    // surf ID = cell ID for all surfs in cell
    // check if uint cell ID overflows int surf ID

    if (nsurf) {
      if (cells[icell].id > maxsurfID)
        error->one(FLERR,"Grid cell ID overflows implicit surf ID");
      surfID = cells[icell].id;
    }

    ptr = csurfs->get(nsurf);

    ipt = 0;
    for (i = 0; i < nsurf; i++) {
      if (svalues) surf->add_tri(surfID,svalues[icell],
                                 pt[ipt+2],pt[ipt+1],pt[ipt]);
      else surf->add_tri(surfID,1,pt[ipt+2],pt[ipt+1],pt[ipt]);
      ipt += 3;
      isurf = surf->nlocal - 1;
      ptr[i] = isurf;
    }

    cells[icell].nsurf = nsurf;
    if (nsurf) {
      cells[icell].csurfs = ptr;
      cinfo[icell].type = OVERLAP;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate function used by both marching squares and cubes
   lo/hi = coordinates of end points of edge of square
   v0/v1 = values at lo/hi end points
   value = interpolated coordinate for thresh value
------------------------------------------------------------------------- */

double MarchingSpheres::extrapolate(double v0, double v1, double lo, double hi)
{
  double cmax = 255.0;
  if (v0 < 0 || v1 < 0) error->one(FLERR,"negative val");
  if (v0 > cmax || v1 > cmax) error->one(FLERR,"big val");

  // both inside or both outside
  if (v0 > 0 && v1 > 0) return 0;
  else if (v0 == 0 && v1 == 0) return 0; 

  // extrapolate from inside
  double value;
  if (v0 > v1) value = lo + (hi-lo)*(v0/cmax);
  else value = lo + (hi-lo)*(1.0-v0/cmax);

  if (value > hi || value < lo) error->one(FLERR,"Vertex off edge");

  // buffer to avoid degenerate triangles
  double ibuffer = (hi-lo)*mindist;
  value = MAX(value,lo+ibuffer);
  value = MIN(value,hi-ibuffer);

  return value;
}

/* ----------------------------------------------------------------------
   clean up issues that marching cubes occasionally generates
     that cause problems for SPARTA
   what MC does:
     may generate 0 or 2 triangles on the face of a cell
     the cell sharing the face may also generate 0 or 2 triangles
     the normals for the 2 triangles may be into or out of the owning cell
   what SPARTA needs:
     let cell1 and cell2 be two cells that share a face
     if cell1 has 2 tris on face and cell2 has none:
       if norm is into cell1: keep them in cell1
       if norm is into cell2: assign both tris to cell2
     if both cell1 and cell2 have 2 tris on face: delete all 4 tris
   algorithm to do this:
     loop over all my cells with implicit tris:
       count how many surfs on each face
     loop over all my cells with implicit tris:
       loop over faces with 2 tris:
         if I own adjoining cell:
           check its tally on shared face
           reassign and/or delete triangles as necessary
         if I do not own adjoining cell:
           add 2 tris to send list for this proc
     irregular comm of send list to nearby procs (share faces of my cells)
     each proc loops over its recv list:
       if my cell face has 2 tris: delete them
       if my cell face has 0 tris: skip or add 2 tris depending on norm
 ------------------------------------------------------------------------- */

void MarchingSpheres::cleanup()
{
  int i,j,k,m,icell,iface,nsurf,idim,nflag,inwardnorm;
  int ntri_other,othercell,otherface,otherproc,otherlocal,othernsurf;
  surfint *oldcsurfs;
  surfint *ptr;
  double *lo,*hi;
  double *norm;

  Surf::Tri *tris = surf->tris;
  Grid::ChildCell *cells = grid->cells;
  MyPage<surfint> *csurfs = grid->csurfs;
  int nglocal = grid->nlocal;

  // DEBUG

  //int nstotal;
  //MPI_Allreduce(&surf->nlocal,&nstotal,1,MPI_INT,MPI_SUM,world);
  //if (me == 0) printf("TOTAL TRI before count: %d\n",nstotal);

  // END of DEBUG

  // count # of tris on each face of every cell I own

  int **nfacetri;
  int ***facetris;
  memory->create(nfacetri,nglocal,6,"readisurf:nfacetri");
  memory->create(facetris,nglocal,6,2,"readisurf:facetris");

  for (icell = 0; icell < nglocal; icell++) {
    nfacetri[icell][0] = nfacetri[icell][1] = nfacetri[icell][2] =
      nfacetri[icell][3] = nfacetri[icell][4] = nfacetri[icell][5] = 0;

    if (cells[icell].nsplit <= 0) continue;
    nsurf = cells[icell].nsurf;
    if (nsurf == 0) continue;

    lo = cells[icell].lo;
    hi = cells[icell].hi;

    for (j = 0; j < nsurf; j++) {
      m = cells[icell].csurfs[j];
      iface = Geometry::tri_on_hex_face(tris[m].p1,tris[m].p2,tris[m].p3,lo,hi);
      if (iface < 0) continue;
      if (nfacetri[icell][iface] < 2)
        facetris[icell][iface][nfacetri[icell][iface]] = m;
      nfacetri[icell][iface]++;
    }
  }

  // check that every face has 0 or 2 tris

  int flag = 0;
  for (icell = 0; icell < nglocal; icell++)
    for (iface = 0; iface < 6; iface++)
      if (nfacetri[icell][iface] != 0 && nfacetri[icell][iface] != 2)
        flag++;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall)
    error->all(FLERR,"Some cell faces do not have zero or 2 triangles");

  // loop over all cell faces
  // check tri count for that face for both adjoining cells

  int *proclist = NULL;
  SendDatum *bufsend = NULL;
  int nsend = 0;
  int maxsend = 0;

  int *dellist = NULL;
  int ndelete = 0;
  int maxdelete = 0;

  // DEBUG
  //int ntotal = 0;
  //int nadd = 0;
  //int ndel = 0;

  for (icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;
    nsurf = cells[icell].nsurf;
    if (nsurf == 0) continue;

    for (iface = 0; iface < 6; iface++) {
      if (nfacetri[icell][iface] != 2) continue;
      //ntotal += 2;

      // other cell/face/proc = info for matching face in adjacent cell

      nflag = grid->neigh_decode(cells[icell].nmask,iface);
      if (nflag != NCHILD && nflag != NPBCHILD)
        error->one(FLERR,"Invalid neighbor cell in cleanup_MC()");

      norm = tris[facetris[icell][iface][0]].norm;
      idim = iface/2;
      if (iface % 2 && norm[idim] < 0.0) inwardnorm = 1;
      else if (iface % 2 == 0 && norm[idim] > 0.0) inwardnorm = 1;
      else inwardnorm = 0;
      if (iface % 2) otherface = iface-1;
      else otherface = iface+1;
      othercell = (int) cells[icell].neigh[iface];
      otherproc = cells[othercell].proc;
      otherlocal = cells[othercell].ilocal;

      // if I own the adjacent cell, make decision about shared tris
      // if both cells have 2 tris on face, delete all of them
      // otherwise cell that matches inward normal is assigned the 2 tris

      if (otherproc == me) {
        ntri_other = nfacetri[othercell][otherface];

        // icell keeps the 2 tris

        if (ntri_other == 0 && inwardnorm) continue;

        // add 2 tris to othercell
        // reset tri IDs to new owning cell

        if (ntri_other == 0) {
          othernsurf = cells[othercell].nsurf;
          oldcsurfs = cells[othercell].csurfs;
          ptr = csurfs->get(othernsurf+2);
          for (k = 0; k < othernsurf; k++)
            ptr[k] = oldcsurfs[k];
          ptr[othernsurf] = facetris[icell][iface][0];
          ptr[othernsurf+1] = facetris[icell][iface][1];
          cells[othercell].nsurf += 2;
          cells[othercell].csurfs = ptr;
          tris[facetris[icell][iface][0]].id = cells[othercell].id;
          tris[facetris[icell][iface][1]].id = cells[othercell].id;
          //printf("MC add1 %d %d\n",cells[icell].id,cells[othercell].id);
          //nadd += 2;
        }

        // delete 2 tris from othercell
        // set nfacetri[othercell] = 0, so won't delete again when it is icell

        if (ntri_other == 2) {
          nfacetri[othercell][otherface] = 0;
          othernsurf = cells[othercell].nsurf;
          ptr = cells[othercell].csurfs;
          m = facetris[othercell][otherface][0];
          for (k = 0; k < othernsurf; k++)
            if (ptr[k] == m) break;
          if (k == othernsurf)
            error->one(FLERR,"Could not find surf in cleanup_MC");
          cells[othercell].csurfs[k] = cells[othercell].csurfs[othernsurf-1];
          othernsurf--;
          m = facetris[othercell][otherface][1];
          for (k = 0; k < othernsurf; k++)
            if (ptr[k] == m) break;
          if (k == othernsurf)
            error->one(FLERR,"Could not find surf in cleanup_MC");
          cells[othercell].csurfs[k] = cells[othercell].csurfs[othernsurf-1];
          othernsurf--;
          cells[othercell].nsurf -= 2;
          //printf("MC del1 %d %d\n",cells[icell].id,cells[othercell].id);
          //ndel += 2;
        }

        // delete 2 tris from icell

        ptr = cells[icell].csurfs;
        m = facetris[icell][iface][0];
        for (k = 0; k < nsurf; k++)
          if (ptr[k] == m) break;
        if (k == nsurf) error->one(FLERR,"Could not find surf in cleanup_MC");
        cells[icell].csurfs[k] = cells[icell].csurfs[nsurf-1];
        nsurf--;
        m = facetris[icell][iface][1];
        for (k = 0; k < nsurf; k++)
          if (ptr[k] == m) break;
        if (k == nsurf) error->one(FLERR,"Could not find surf in cleanup_MC");
        cells[icell].csurfs[k] = cells[icell].csurfs[nsurf-1];
        nsurf--;
        cells[icell].nsurf -= 2;
        //printf("MC dele %d %d\n",cells[icell].id,cells[othercell].id);
        //ndel += 2;

        // add 4 tris to delete list if both cells deleted them

        if (ntri_other == 2) {
          if (ndelete+4 > maxdelete) {
            maxdelete += DELTA;
            memory->grow(dellist,maxdelete,"readisurf:dellist");
          }
          dellist[ndelete++] = facetris[icell][iface][0];
          dellist[ndelete++] = facetris[icell][iface][1];
          dellist[ndelete++] = facetris[othercell][otherface][0];
          dellist[ndelete++] = facetris[othercell][otherface][1];
        }

      // cell face is shared with another proc
      // send it the cell/face indices and the 2 tris,
      //   in case they need to be assigned to the other cell based on norm

      } else {
        if (nsend == maxsend) {
          maxsend += DELTA;
          proclist = (int *)
            memory->srealloc(proclist,maxsend*sizeof(int),
                             "readisurf:proclist");
          bufsend = (SendDatum *)
            memory->srealloc(bufsend,maxsend*sizeof(SendDatum),
                             "readisurf:bufsend");
        }
        proclist[nsend] = otherproc;
        bufsend[nsend].sendcell = icell;
        bufsend[nsend].sendface = iface;
        bufsend[nsend].othercell = otherlocal;
        bufsend[nsend].otherface = otherface;
        bufsend[nsend].inwardnorm = inwardnorm;
        memcpy(&bufsend[nsend].tri1,&tris[facetris[icell][iface][0]],
               sizeof(Surf::Tri));
        memcpy(&bufsend[nsend].tri2,&tris[facetris[icell][iface][1]],
               sizeof(Surf::Tri));
        nsend++;

        // if not inwardnorm, delete 2 tris from this cell
        // also add them to delete list

        if (!inwardnorm) {
          ptr = cells[icell].csurfs;
          m = facetris[icell][iface][0];
          for (k = 0; k < nsurf; k++)
            if (ptr[k] == m) break;
          if (k == nsurf) error->one(FLERR,"Could not find surf in cleanup_MC");
          cells[icell].csurfs[k] = cells[icell].csurfs[nsurf-1];
          nsurf--;
          m = facetris[icell][iface][1];
          for (k = 0; k < nsurf; k++)
            if (ptr[k] == m) break;
          if (k == nsurf) error->one(FLERR,"Could not find surf in cleanup_MC");
          cells[icell].csurfs[k] = cells[icell].csurfs[nsurf-1];
          nsurf--;
          cells[icell].nsurf -= 2;
          //ndel += 2;

          if (ndelete+2 > maxdelete) {
            maxdelete += DELTA;
            memory->grow(dellist,maxdelete,"readisurf:dellist");
          }
          dellist[ndelete++] = facetris[icell][iface][0];
          dellist[ndelete++] = facetris[icell][iface][1];
        }
      }
    }
  }

  // perform irregular communication of list of cell faces and tri pairs

  Irregular *irregular = new Irregular(sparta);
  int nrecv = irregular->create_data_uniform(nsend,proclist,1);

  SendDatum *bufrecv = (SendDatum *)
    memory->smalloc(nrecv*sizeof(SendDatum),"readisurf:bufrecv");

  irregular->exchange_uniform((char *) bufsend,sizeof(SendDatum),
                              (char *) bufrecv);
  delete irregular;
  memory->sfree(proclist);
  memory->sfree(bufsend);

  // loop over list of received face/tri info
  // if my matching face has 2 tris, delete them
  // if my matching face has 0 tris, skip or add 2 tris depending on norm

  for (i = 0; i < nrecv; i++) {
    icell = bufrecv[i].othercell;
    iface = bufrecv[i].otherface;

    // my icell is not affected, sender cell keeps its 2 tris

    if (nfacetri[icell][iface] == 0 && bufrecv[i].inwardnorm) continue;

    // add 2 tris to icell and this processor's Surf::tris list
    // set tri IDs to new owning cell, must be done after memcpy()
    // NOTE: what about tri types?

    if (nfacetri[icell][iface] == 0) {
      int nslocal = surf->nlocal;
      surf->add_tri(cells[icell].id,1,
                    bufrecv[i].tri1.p1,bufrecv[i].tri1.p2,bufrecv[i].tri1.p3);
      memcpy(&surf->tris[nslocal],&bufrecv[i].tri1,sizeof(Surf::Tri));
      surf->tris[nslocal].id = cells[icell].id;
      surf->add_tri(cells[icell].id,1,
                    bufrecv[i].tri2.p1,bufrecv[i].tri2.p2,bufrecv[i].tri2.p3);
      memcpy(&surf->tris[nslocal+1],&bufrecv[i].tri2,sizeof(Surf::Tri));
      surf->tris[nslocal+1].id = cells[icell].id;

      nsurf = cells[icell].nsurf;
      oldcsurfs = cells[icell].csurfs;
      ptr = csurfs->get(nsurf+2);
      for (k = 0; k < nsurf; k++)
        ptr[k] = oldcsurfs[k];
      ptr[nsurf] = nslocal;
      ptr[nsurf+1] = nslocal+1;
      cells[icell].nsurf += 2;
      cells[icell].csurfs = ptr;
      //nadd += 2;
    }

    // both cells have 2 tris on common face
    // need to delete my 2 tris from icell
    // sender will get similar message from me and delete
    // inwardnorm check to see if I already deleted when sent a message,
    // else delete now and add 2 tris to delete list

    if (nfacetri[icell][iface] == 2) {
      norm = tris[facetris[icell][iface][0]].norm;
      idim = iface/2;
      if (iface % 2 && norm[idim] < 0.0) inwardnorm = 1;
      else if (iface % 2 == 0 && norm[idim] > 0.0) inwardnorm = 1;
      else inwardnorm = 0;
      if (!inwardnorm) continue;

      nsurf = cells[icell].nsurf;
      ptr = cells[icell].csurfs;
      m = facetris[icell][iface][0];
      for (k = 0; k < nsurf; k++)
        if (ptr[k] == m) break;
      if (k == nsurf) error->one(FLERR,"Could not find surf in cleanup_MC");
      cells[icell].csurfs[k] = cells[icell].csurfs[nsurf-1];
      nsurf--;
      m = facetris[icell][iface][1];
      for (k = 0; k < nsurf; k++)
        if (ptr[k] == m) break;
      if (k == nsurf) error->one(FLERR,"Could not find surf in cleanup_MC");
      cells[icell].csurfs[k] = cells[icell].csurfs[nsurf-1];
      nsurf--;
      cells[icell].nsurf -= 2;
      //ndel += 2;

      if (ndelete+2 > maxdelete) {
        maxdelete += DELTA;
        memory->grow(dellist,maxdelete,"readisurf:dellist");
      }
      dellist[ndelete++] = facetris[icell][iface][0];
      dellist[ndelete++] = facetris[icell][iface][1];
    }
  }

  memory->sfree(bufrecv);
  memory->destroy(nfacetri);
  memory->destroy(facetris);

  // compress Surf::tris list to remove deleted tris
  // must sort dellist, so as to compress tris in DESCENDING index order
  // descending, not ascending, so that a surf is not moved from end-of-list
  //   that is flagged for later deletion
  // must repoint one location in cells->csurfs to moved surf
  //   requires grid hash to find owning cell of moved surf
  // note that ghost surfs exist at this point, but caller will clear them

  if (!grid->hashfilled) grid->rehash();

  qsort(dellist,ndelete,sizeof(int),compare_indices);

  tris = surf->tris;
  int nslocal = surf->nlocal;
  for (i = 0; i < ndelete; i++) {
    m = dellist[i];
    if (m != nslocal-1) memcpy(&tris[m],&tris[nslocal-1],sizeof(Surf::Tri));
    nslocal--;

    icell = (*grid->hash)[tris[m].id];
    nsurf = cells[icell].nsurf;
    ptr = cells[icell].csurfs;
    for (k = 0; k < nsurf; k++)
      if (ptr[k] == nslocal) {
        ptr[k] = m;
        break;
      }
    if (k == nsurf) error->one(FLERR,"Did not find moved tri in cleanup_MC()");
  }

  surf->nlocal = nslocal;
  memory->destroy(dellist);
}

/* ----------------------------------------------------------------------
   adding triangles
------------------------------------------------------------------------- */

int MarchingSpheres::add_triangle(int *trig, int n)
{
  for(int t = 0; t < 3*n; t++) {
    switch (trig[t]) {
    case 0:
      pt[t][0] = i0;
      pt[t][1] = lo[1];
      pt[t][2] = lo[2];
      break;
    case 1:
      pt[t][0] = hi[0];
      pt[t][1] = i1;
      pt[t][2] = lo[2];
      break;
    case 2:
      pt[t][0] = i2;
      pt[t][1] = hi[1];
      pt[t][2] = lo[2];
      break;
    case 3:
      pt[t][0] = lo[0];
      pt[t][1] = i3;
      pt[t][2] = lo[2];
      break;
    case 4:
      pt[t][0] = i4;
      pt[t][1] = lo[1];
      pt[t][2] = hi[2];
      break;
    case 5:
      pt[t][0] = hi[0];
      pt[t][1] = i5;
      pt[t][2] = hi[2];
      break;
    case 6:
      pt[t][0] = i6;
      pt[t][1] = hi[1];
      pt[t][2] = hi[2];
      break;
    case 7:
      pt[t][0] = lo[0];
      pt[t][1] = i7;
      pt[t][2] = hi[2];
      break;
    case 8:
      pt[t][0] = lo[0];
      pt[t][1] = lo[1];
      pt[t][2] = i8;
      break;
    case 9:
      pt[t][0] = hi[0];
      pt[t][1] = lo[1];
      pt[t][2] = i9;
      break;
    case 10:
      pt[t][0] = hi[0];
      pt[t][1] = hi[1];
      pt[t][2] = i10;
      break;
    case 11:
      pt[t][0] = lo[0];
      pt[t][1] = hi[1];
      pt[t][2] = i11;
      break;
    case 12: {
      int u = 0;
      pt[t][0] = pt[t][1] = pt[t][2] = 0.0;
      if (bit0 ^ bit1) {
        ++u;
        pt[t][0] += i0;
        pt[t][1] += lo[1];
        pt[t][2] += lo[2];
      }
      if (bit1 ^ bit2) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += i1;
        pt[t][2] += lo[2];
      }
      if (bit2 ^ bit3) {
        ++u;
        pt[t][0] += i2;
        pt[t][1] += hi[1];
        pt[t][2] += lo[2];
      }
      if (bit3 ^ bit0) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += i3;
        pt[t][2] += lo[2];
      }
      if (bit4 ^ bit5) {
        ++u;
        pt[t][0] += i4;
        pt[t][1] += lo[1];
        pt[t][2] += hi[2];
      }
      if (bit5 ^ bit6) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += i5;
        pt[t][2] += hi[2];
      }
      if (bit6 ^ bit7) {
        ++u;
        pt[t][0] += i6;
        pt[t][1] += hi[1];
        pt[t][2] += hi[2];
      }
      if (bit7 ^ bit4) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += i7;
        pt[t][2] += hi[2];
      }
      if (bit0 ^ bit4) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += lo[1];
        pt[t][2] += i8;
      }
      if (bit1 ^ bit5) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += lo[1];
        pt[t][2] += i9;
      }
      if (bit2 ^ bit6) {
        ++u;
        pt[t][0] += hi[0];
        pt[t][1] += hi[1];
        pt[t][2] += i10;
      }
      if (bit3 ^ bit7) {
        ++u;
        pt[t][0] += lo[0];
        pt[t][1] += hi[1];
        pt[t][2] += i11;
      }

      pt[t][0] /= static_cast<double> (u);
      pt[t][1] /= static_cast<double> (u);
      pt[t][2] /= static_cast<double> (u);
      break;
    }

    default:
      break;
    }
  }

  return n;
}

/* ----------------------------------------------------------------------
   test a face
   only need to check if sum of radii is greater than diag
------------------------------------------------------------------------- */

bool MarchingSpheres::test_face(int face)
{
  double A,B,C,D;
  double diag;

  switch (face) {
  // lower x-z
  case -1:
  case 1:
    A = v000iso;
    B = v100iso;
    C = v101iso;
    D = v001iso;
    break;
  // upper y-z
  case -2:
  case 2:
    A = v001iso;
    B = v101iso;
    C = v111iso;
    D = v011iso;
    break;
  // uppper x-z
  case -3:
  case 3:
    A = v011iso;
    B = v111iso;
    C = v110iso;
    D = v010iso;
    break;
  // lower y-z
  case -4:
  case 4:
    A = v010iso;
    B = v110iso;
    C = v100iso;
    D = v000iso;
    break;
  // lower x-y
  case -5:
  case 5:
    A = v000iso;
    B = v010iso;
    C = v011iso;
    D = v001iso;
    break;
  // upper x-y
  case -6:
  case 6:
    A = v100iso;
    B = v110iso;
    C = v111iso;
    D = v101iso;
    break;

  default:
    A = B = C = D = 0.0;
    error->one(FLERR,"Invalid face code");
  };

  if (A+C >= sqrt(2) || B+D >= sqrt(2)) return 1;
  else return 0;
}

/* ----------------------------------------------------------------------
   test the interior of a cube
   only need to check diagonls that cross interior point
------------------------------------------------------------------------- */

bool MarchingSpheres::test_interior()
{
  // check all large diagonals
  if (v000iso+v111iso > sqrt(3)) return 1;
  if (v001iso+v110iso > sqrt(3)) return 1;
  if (v101iso+v010iso > sqrt(3)) return 1;
  if (v100iso+v011iso > sqrt(3)) return 1;
  return 0;
}

