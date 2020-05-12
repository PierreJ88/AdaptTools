#ifndef _MSHMET_H
#define _MSHMET_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>

#include "chrono.h"
#include "libmesh5.h"
#include "ms_calls.h"


#define MS_VER   "4.0a"
#define MS_REL   "Mar. 12, 2016"
#define MS_CPY   "(C) Copyright 2007- , ICS-SU"

#define MS_EPS    1.e-6
#define MS_EPSD   1.e-200
#define LONMAX   4096

#define CTE2D    2.0 / 9.0
#define CTE3D    9.0 / 32.0

#define CTE_GEOPHY_FREQ  10.0
#define CTE_GEOPHY_ORDER 1
#define CTE_GEOPHY_PPW 20

#define CTE_GEOPHY_RMIN  1.0
#define CTE_GEOPHY_HMIN  50.0
#define CTE_GEOPHY_HMAX 200.0

#define MS_MAX(a,b)   ( ((a) < (b)) ? (b) : (a) )
#define MS_MIN(a,b)   ( ((a) < (b)) ? (a) : (b) )
#define MS_MIN3(a,b,c) ( (a) < (b) ? ((a)<(c) ? (a) : (c)) : ((b)<(c) ? (b) : (c)) )
#define MS_MAX3(a,b,c) ( (a) > (b) ? ((a)>(c) ? (a) : (c)) : ((b)>(c) ? (b) : (c)) )


typedef struct {
  double         c[3];
  int            s;
} Point;
typedef Point * pPoint;

typedef struct {
  int     v[3],adj[3],mark;
} Tria;
typedef Tria * pTria;

typedef struct {
  int     v[4],adj[4],mark;
} Tetra;
typedef Tetra * pTetra;

typedef struct {
  double   hmin,hmax,hgrad,err,rmin,freq;
  int      np,nt,ne,dim,ver,order,ppw;
  char     verb,iso,ls,nrm,grad;
  mytime   ctim[TIMEMAX];
} Info;

typedef struct {
  int      mark;
  char    *name;
  pPoint   point;
  pTria    tria;
  pTetra   tetra;
} Mesh;
typedef Mesh * pMesh;

typedef struct {
  double     *u,*g,*h,*m,umax;
  char       *name,*namepar;
} Sol;
typedef Sol * pSol;

struct _MSst {
  Mesh  mesh;
  Sol   sol;
  Info  info;
};

/* prototypes */
int  loadMesh(MSst *msst);
int  loadSol(MSst *msst);
int  saveMet(MSst *msst);
int  hashel_2d(MSst *mist);
int  hashel_3d(MSst *mist);
int  mshmet1_2d(MSst *mist);
int  mshmet1_3d(MSst *mist);
int  nrmgrad(MSst *msst);
int  defgrad2met(MSst *msst);


#endif
