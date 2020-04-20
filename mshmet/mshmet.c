/*
 * main program file for mshmet
 * (C) Copyright 1997 - , ICS-SU
 *
 * This file is part of mshmet.
 *
 * mshmet is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * mshmet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with AdaptTools.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "mshmet.h"
#include "ms_calls.h"


static void excfun(int sigid) {
  fprintf(stdout,"\n # unexpected error:");  fflush(stdout);
  switch(sigid) {
    case SIGABRT:
      fprintf(stdout," abnormal stop\n");  break;
    case SIGBUS:
      fprintf(stdout," code error...\n");  break;
    case SIGFPE:
      fprintf(stdout," floating-point exception\n"); break;
    case SIGILL:
      fprintf(stdout," illegal instruction\n"); break;
    case SIGSEGV:
      fprintf(stdout," segmentation fault.\n");  break;
    case SIGTERM:
    case SIGINT:
      fprintf(stdout," programm killed.\n");  break;
  }
  fprintf(stdout," # no data file saved.\n");
  exit(1);
}


static void usage(char *prog) {
  fprintf(stdout,"\n usage: %s [+/-v | -h | -i | -l] source[.mesh] [-e err] [-freq|-order|-rmin|-hmin|-hmax|-hgrad val] [-n nrm] [-o output[].sol]] [-p param[.mhes]]\n",prog);

  fprintf(stdout,"\nOptions and flags:\n\
  --help       show the syntax and exit.\n\
  --version    show the version and date of release and exit.\n\n\
  -e err       approximation error\n\
  -i           isotropic metric (default is aniso)\n\
  -l           metric for level sets\n\
  -n nrm       metric normalization method (default=1)\n\
  -g           create a map linked to the gradient (default=0)\n\
  -v           suppress any message (for use with function call).\n\
  +v           increase the verbosity level for output.\n\n\
  source.mesh    name of the source mesh\n\
  metric.sol     name of metric file\n\
  param.mhes     name of file containing metric specifications\n\
  output.sol     name of the output file\n");
  exit(1);
}


static int parsar(int argc,char *argv[],MSst *msst) {
  int      i;
  char    *ptr,*data;

  i = 1;
  while ( i < argc ) {
    if ( (*argv[i] == '-') || (*argv[i ]== '+') ) {
      switch(argv[i][1]) {
      case '-':  /* on-line help */
        if ( !strcmp(argv[i],"--help") )
          usage(argv[0]);
        else if ( !strcmp(argv[i],"--version") ) {
          fprintf(stdout,"%s: version: %s release: %s\n",argv[0],MS_VER,MS_REL);
          exit(1);
        }
        break;
      case 'h':
      case '?':
        if ( !strcmp(argv[i],"-hmin") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.hmin = strtod(argv[i],NULL);
        }
	else if ( !strcmp(argv[i],"-rmin") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.rmin = strtod(argv[i],NULL);
        }
	else if ( !strcmp(argv[i],"-freq") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.freq = strtod(argv[i],NULL);
        }
	else if ( !strcmp(argv[i],"-order") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.order = strtod(argv[i],NULL);
        }
        else if ( !strcmp(argv[i],"-hmax") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.hmax = strtod(argv[i],NULL);
        }
        else if ( !strcmp(argv[i],"-hgrad") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.hgrad = strtod(argv[i],NULL);
        }
        else
          usage(argv[0]);
        break;
      case 'e':
        if ( ++i < argc && isdigit(argv[i][0]) )
          msst->info.err = strtod(argv[i],NULL);
        else {
          fprintf(stderr,"%s: missing argument option\n",argv[0]);
          usage(argv[0]);
        }
        break;
      case 'i':
        msst->info.iso = 1;
        break;
      case 'l':
          msst->info.ls = 1;
        break;
      case 'f':
	if ( !strcmp(argv[i],"-freq") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.freq = strtod(argv[i],NULL);
	}
        break;
      case 'g':
	msst->info.grad = 1;
        break;
      case 'n':
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          msst->info.nrm = atoi(argv[i]);
          if ( msst->info.nrm < 0 || msst->info.nrm > 2 ) {
            fprintf(stdout,"%s: wrong parameter value: set default\n",argv[0]);
            msst->info.nrm = 1;
          }
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'p':
        if ( ++i < argc ) {
          msst->sol.namepar = argv[i];
          ptr = strstr(msst->sol.namepar,".mhes");
          if ( !ptr )  strcat(msst->sol.namepar,".mhes");
        }
        else {
          fprintf(stdout,"%s: missing parameter file\n", argv[0]);
          usage(argv[0]);
        }
        break;
      case 'r':
	if ( !strcmp(argv[i],"-rmin") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.rmin = strtod(argv[i],NULL);
	}
	break;
      case 'o':
	if ( !strcmp(argv[i],"-order") ) {
          if ( ++i < argc && isdigit(argv[i][0]) )
            msst->info.order = strtod(argv[i],NULL);
	}
        break;
      case 'v':
        if ( !strcmp(argv[i],"-v") )
          msst->info.verb = '0';
        else if ( !strcmp(argv[i],"+v") )
          msst->info.verb = '+';
        else {
          fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
          usage(argv[0]);
        }
        break;
      default:
        fprintf(stderr,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    else {
      if ( msst->mesh.name == NULL ) {
        data = (char*)calloc(strlen(argv[i])+10,sizeof(char));
        strcpy(data,argv[i]);
        ptr = strstr(data,".mesh");
        if ( !ptr )  strcat(data,".mesh");
        msst->mesh.name = data;
      }
      else {
        fprintf(stdout,"%s: illegal option %s\n",argv[0],argv[i]);
        usage(argv[0]);
      }
    }
    i++;
  }

  /* check params */
  if ( msst->mesh.name == NULL ) {
    if ( msst->info.verb != '0' )  fprintf(stderr,"%s: missing argument\n",argv[0]);
    usage(argv[0]);
  }
  if ( !msst->sol.name ) {
    msst->sol.name = (char *)calloc(128,sizeof(char));
    assert(msst->sol.name);
    strcpy(msst->sol.name,msst->mesh.name);
  }

  return(1);
}


int parsop(MSst *msst) {
  int      i,ret;
  char    *ptr,data[256];
  FILE    *in;

  /* check for parameter file */
  if ( !msst->sol.namepar ) {
    strcpy(data,msst->mesh.name);
    ptr = strstr(data,".mesh");
    if ( ptr )  *ptr = '\0';
    strcat(data,".mhes");
    in = fopen(data,"r");
    if ( !in ) {
      sprintf(data,"%s","DEFAULT.mhes");
      in = fopen(data,"r");
    }
  }
  else {
    strcpy(data,msst->sol.namepar);
    ptr = strstr(data,".mhes");
    if ( !ptr )  strcat(data,".mhes");
    in = fopen(data,"r");
  }
  if ( !in )  return(-1);
  if ( msst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* read metric specifications */
  while ( !feof(in) ) {
    /* scan line */
    ret = fscanf(in,"%s",data);
    if ( !ret || feof(in) )  break;
    for (i=0; i<strlen(data); i++) data[i] = tolower(data[i]);

    /* check for keywords */
    if ( !strcmp(data,"hmin") )
      fscanf(in,"%lf",&msst->info.hmin);
    else if ( !strcmp(data,"rmin") )
      fscanf(in,"%lf",&msst->info.rmin);
    else if ( !strcmp(data,"order") )
      fscanf(in,"%d",&msst->info.order);
    else if ( !strcmp(data,"freq") )
      fscanf(in,"%lf",&msst->info.freq);
    else if ( !strcmp(data,"hmax") )
      fscanf(in,"%lf",&msst->info.hmax);
    else if ( !strcmp(data,"err") )
      fscanf(in,"%lf",&msst->info.err);
    else if ( !strcmp(data,"iso") )
      msst->info.iso = 1;
    else if ( data[0] == '#' ) {
      fgets(data,255,in);
    }
  }
  fclose(in);

  return(1);
}


int main(int argc,char **argv) {
  MSst       msst;
  int        ier;
  char       stim[32];

  memset(&msst,0,sizeof(MSst));
  tminit(msst.info.ctim,TIMEMAX);
  chrono(ON,&msst.info.ctim[0]);

  /* trap exceptions */
  signal(SIGABRT,excfun);
  signal(SIGFPE,excfun);
  signal(SIGILL,excfun);
  signal(SIGSEGV,excfun);
  signal(SIGTERM,excfun);
  signal(SIGINT,excfun);
  signal(SIGBUS,excfun);

  /* init structure */
  memset(&msst.mesh,0,sizeof(Mesh));
  memset(&msst.sol,0,sizeof(Sol));

  /* default values */
	msst.info.dim   = 3;
  msst.info.err   = 0.01;
  msst.info.hmin  = 0.01;
  msst.info.hmax  = 1.0;
  msst.info.iso   = 0;
  msst.info.ls    = 0;
  msst.info.grad  = 0;
  msst.info.nrm   = 1;
  msst.info.ver   = 1;
  msst.info.rmin  = CTE_GEOPHY_HMIN;
  msst.info.freq  = CTE_GEOPHY_FREQ;
  msst.info.order = CTE_GEOPHY_ORDER;
  msst.info.verb  = '1';

  /* parse command line */
  if ( !parsar(argc,argv,&msst) )  return(1);

  /* loading data */
  chrono(ON,&msst.info.ctim[1]);

  if ( msst.info.verb != '0' ) {
    fprintf(stdout," - MSHMET, Release %s, %s\n   %s\n\n",MS_VER,MS_REL,MS_CPY);
    fprintf(stdout," - LOADING DATA\n");
  }

  /* loading mesh and solution */
  ier = loadMesh(&msst);
  if ( ier <= 0 )  return(1);
  ier = loadSol(&msst);
  if ( ier <= 0 )  return(1);
  if ( msst.sol.umax < MS_EPSD ) {
    if ( msst.info.verb != '0' )  fprintf(stdout," # Warning: solution modulus too small\n");
    return(1);
  }
  if ( !parsop(&msst) )  return(1);

  /* setting adjacencies */
  msst.info.dim == 2 ? hashel_2d(&msst) : hashel_3d(&msst);

  chrono(OFF,&msst.info.ctim[1]);
	printim(msst.info.ctim[1].gdif,stim);
  if ( msst.info.verb != '0' )  fprintf(stdout," - COMPLETED: %s\n",stim);

  /* build metric */
  chrono(ON,&msst.info.ctim[2]);
  if ( msst.info.verb != '0' )
    fprintf(stdout,"\n ** MODULE MSHMET: %s\n",MS_VER);
  ier = MS_mshmet(&msst);
  chrono(OFF,&msst.info.ctim[2]);
  if ( msst.info.verb != '0' ) {
		printim(msst.info.ctim[2].gdif,stim);
    if ( ier )
      fprintf(stdout," ** COMPLETED: %s\n\n",stim);
    else
      fprintf(stdout," ** NOT COMPLETED!: %s\n\n",stim);
  }

  /* save file */
  if ( msst.info.verb != '0' )  fprintf(stdout," - WRITING DATA\n");
  chrono(ON,&msst.info.ctim[3]);
  ier = saveMet(&msst);
  if ( !ier )  return(1);
  chrono(OFF,&msst.info.ctim[3]);
  if ( msst.info.verb != '0' ) {
    printim(msst.info.ctim[3].gdif,stim);
    fprintf(stdout," - COMPLETED: %s\n",stim);
  }

  /* free mem */
  free(msst.mesh.point);
  if ( msst.info.nt > 0 )  free(msst.mesh.tria);
  if ( msst.info.ne > 0 )  free(msst.mesh.tetra);
  free(msst.sol.u);
  free(msst.sol.m);

  chrono(OFF,&msst.info.ctim[0]);
  if ( msst.info.verb != '0' ) {
	  printim(msst.info.ctim[0].gdif,stim);
    fprintf(stdout,"\n ** Cumulative time: %s.\n",stim);
  }

  return(0);
}
