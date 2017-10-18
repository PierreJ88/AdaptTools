#include "mshmet.h"
#include "ms_calls.h"


MSst *MS_init(int dim,int ver) {
  MSst   *msst;
  
  /* default values */
  msst = (MSst *)calloc(1,sizeof(MSst));
  memset(&msst->mesh,0,sizeof(Mesh));
  memset(&msst->sol,0,sizeof(Sol));
  
  msst->info.verb = '1';
  msst->info.dim  = dim;
  msst->info.ver  = ver;

  /* init timer */
  tminit(msst->info.ctim,TIMEMAX);
  chrono(ON,&msst->info.ctim[0]);

  return(msst);
}


int MS_stop(MSst *msst) {
  char   stim[32];
  
  free(msst->mesh.point);
  if ( msst->info.nt > 0 )  free(msst->mesh.tria);
  if ( msst->info.ne > 0 )  free(msst->mesh.tetra);
  free(msst->sol.u);
  free(msst->sol.m);
  
  chrono(OFF,&msst->info.ctim[0]);
  if ( msst->info.verb != '0' ) {
	  printim(msst->info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s sec.\n",stim);
  }
  
  return(1);
}


int MS_mshmet(MSst *msst) {
  int   ier;
  
  if ( msst->info.dim == 2 )
    ier = mshmet1_2d(msst);
  else
    ier = mshmet1_3d(msst);

  return(ier);
}

