#include "mshmet.h"


/* read mesh data */
int loadMesh(MSst *msst) {
  pPoint       ppt;
  pTria        pt1;
  pTetra       pt;
  float        fp1,fp2,fp3;
  int          i,k,inm,ref;
  char        *ptr,data[256];

  strcpy(data,msst->mesh.name);
  ptr = strstr(data,".mesh");
  if ( !ptr ) {
    strcat(data,".meshb");
    if( !(inm = GmfOpenMesh(data,GmfRead,&msst->info.ver,&msst->info.dim)) ) {
      ptr = strstr(data,".mesh");
      *ptr = '\0';
      strcat(data,".mesh");
      if ( !(inm = GmfOpenMesh(data,GmfRead,&msst->info.ver,&msst->info.dim)) ) {
        fprintf(stderr," # %s: file not found.\n",data);
        return(0);
      }
    }
  }
  else if ( !(inm = GmfOpenMesh(data,GmfRead,&msst->info.ver,&msst->info.dim)) ) {
    fprintf(stderr," # %s: file not found.\n",data);
    return(0);
  }

  if ( msst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  msst->info.np = GmfStatKwd(inm,GmfVertices);
  msst->info.nt = GmfStatKwd(inm,GmfTriangles);
  msst->info.ne = GmfStatKwd(inm,GmfTetrahedra);
  if ( !msst->info.np ) {
    if ( msst->info.verb != '0' )  fprintf(stdout,"\n # missing data\n");
    return(0);
  }

  /* memory allocation */
  msst->mesh.point = (pPoint)calloc(msst->info.np+1,sizeof(Point));
  assert(msst->mesh.point);
  GmfGotoKwd(inm,GmfVertices);
  if ( msst->info.dim == 2 ) {
    /* 2d mesh */
    for (k=1; k<=msst->info.np; k++) {
      ppt = &msst->mesh.point[k];
      if ( msst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ref);
    }
  }
  else {
    /* 3d mesh */
    for (k=1; k<=msst->info.np; k++) {
      ppt = &msst->mesh.point[k];
      if ( msst->info.ver == GmfFloat ) {
        GmfGetLin(inm,GmfVertices,&fp1,&fp2,&fp3,&ref);
        ppt->c[0] = fp1;
        ppt->c[1] = fp2;
        ppt->c[2] = fp3;
      }
      else
        GmfGetLin(inm,GmfVertices,&ppt->c[0],&ppt->c[1],&ppt->c[2],&ref);
    }
  }
  if ( msst->info.nt > 0 ) {
    msst->mesh.tria  = (pTria)calloc(msst->info.nt+1,sizeof(Tria));
    assert(msst->mesh.tria);
    /* read mesh triangles */
    GmfGotoKwd(inm,GmfTriangles);
    for (k=1; k<=msst->info.nt; k++) {
      pt1 = &msst->mesh.tria[k];
      GmfGetLin(inm,GmfTriangles,&pt1->v[0],&pt1->v[1],&pt1->v[2],&ref);
    }
  }
  if ( msst->info.ne > 0 ) {
    msst->mesh.tetra  = (pTetra)calloc(msst->info.ne+1,sizeof(Tetra));
    assert(msst->mesh.tetra);
    /* read tetrahedra */
    GmfGotoKwd(inm,GmfTetrahedra);
    for (k=1; k<=msst->info.ne; k++) {
      pt = &msst->mesh.tetra[k];
      GmfGetLin(inm,GmfTetrahedra,&pt->v[0],&pt->v[1],&pt->v[2],&pt->v[3],&ref);
    }
  }
  GmfCloseMesh(inm);

  if ( msst->info.verb != '0' ) {
    fprintf(stdout," %d vertices",msst->info.np);
    if ( msst->info.nt )  fprintf(stdout,", %d triangles",msst->info.nt);
    if ( msst->info.ne )  fprintf(stdout,", %d tetrahedra",msst->info.ne);
    fprintf(stdout,"\n");
  }

  return(1);
}


/* load one solution */
int loadSol(MSst *msst) {
  double       dbuf[GmfMaxTyp];
  float        fbuf[GmfMaxTyp];
  int          k,i,inm,np,dim,ver,type,size,typtab[GmfMaxTyp];
  char        *ptr,data[128];

  if ( !msst->sol.name )  return(-1);
  strcpy(data,msst->sol.name);

  /* remove .mesh extension */
  ptr = strstr(data,".mesh");
  if ( ptr )  *ptr = '\0';

  /* look for data file */
  ptr = strstr(data,".sol");
  if ( ptr ) {
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
  }
  else {
    /* first try to read binary file */
    strcat(data,".solb");
    inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    if ( !inm ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      inm = GmfOpenMesh(data,GmfRead,&ver,&dim);
    }
  }
  if ( !inm )  return(-1);

  if ( dim != msst->info.dim )  return(-1);
  np = GmfStatKwd(inm,GmfSolAtVertices,&type,&size,&typtab);
  if ( !np || np != msst->info.np )  return(-1);

  if ( msst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* memory allocation */
  msst->sol.u = (double*)calloc(msst->info.np,sizeof(double));
  assert(msst->sol.u);

  /* look at first solution only (for now) */
  GmfGotoKwd(inm,GmfSolAtVertices);
  msst->sol.umax = -FLT_MAX;
  switch(typtab[0]) {
    case GmfSca:
      /* read scalar field */
      for (k=0; k<msst->info.np; k++) {
        if ( ver == GmfFloat ) {
          GmfGetLin(inm,GmfSolAtVertices,fbuf);
          msst->sol.u[k] = fbuf[0];
        }
        else {
          GmfGetLin(inm,GmfSolAtVertices,dbuf);
          msst->sol.u[k] = dbuf[0];
        }
        msst->sol.umax = MS_MAX(msst->sol.umax,msst->sol.u[k]);
      }
      break;
    case GmfVec:
      /* read vector field */
      for (k=0; k<msst->info.np; k++) {
        if ( ver == GmfFloat ) {
          GmfGetLin(inm,GmfSolAtVertices,fbuf);
          msst->sol.u[k] = sqrt(fbuf[0]*fbuf[0] + fbuf[1]*fbuf[1]);
        }
        else {
          GmfGetLin(inm,GmfSolAtVertices,dbuf);
          msst->sol.u[k] = sqrt(dbuf[0]*dbuf[0] + dbuf[1]*dbuf[1]);
        }
        msst->sol.umax = MS_MAX(msst->sol.umax,msst->sol.u[k]);
      }
      break;
    case GmfSymMat:
      /* read tensor field */
      for (k=0; k<msst->info.np; k++) {
        msst->sol.umax = MS_MAX(msst->sol.umax,msst->sol.u[k]);
      }
      break;
    default:
      break;
  }
  GmfCloseMesh(inm);

  if ( msst->info.verb != '0' ) {
    if ( typtab[0] == GmfSca)
      fprintf(stdout," %d scalar data\n",msst->info.np);
    else if ( typtab[0] == GmfVec )
      fprintf(stdout," %d vector data\n",msst->info.np);
    else
      fprintf(stdout," %d tensor data\n",msst->info.np);
  }

  return(1);  
}


/* save metric file */
int saveMet(MSst *msst) {
  double     db[6];
  int        i,k,outm,type,typtab[GmfMaxTyp];
  char      *ptr,data[128];

  strcpy(data,msst->sol.name);
  ptr = strstr(data,".mesh");
  if ( ptr )  {
    *ptr = '\0';
    strcat(data,".met.sol");
  }
  else {
    ptr = strstr(data,".met.sol");
    if ( !ptr )  strcat(data,".met.solb");
  }

  msst->info.ver = GmfDouble;
  if ( !(outm = GmfOpenMesh(data,GmfWrite,msst->info.ver,msst->info.dim)) ) {
    if ( msst->info.verb != '0' )  fprintf(stderr," # unable to open %s\n",data);
    return(0);
  }
  if ( msst->info.verb != '0' )  fprintf(stdout,"    %s:",data);

  /* write metric file */
  type = 1;
  if ( msst->info.iso ) {
    typtab[0] = GmfSca;
    GmfSetKwd(outm,GmfSolAtVertices,msst->info.np,type,typtab);
    for (k=0; k<msst->info.np; k++) {
      GmfSetLin(outm,GmfSolAtVertices,&msst->sol.m[k]);
    }
  }
  else {
    typtab[0] = GmfSymMat;
    GmfSetKwd(outm,GmfSolAtVertices,msst->info.np,type,typtab);
    if ( msst->info.dim == 2 ) {
      for (k=0; k<msst->info.np; k++) {
        for (i=0; i<3; i++)  db[i] = msst->sol.m[3*k+i];
        GmfSetLin(outm,GmfSolAtVertices,db);
      }
    }
  }
  GmfCloseMesh(outm);

  if ( msst->info.verb != '0' )
    fprintf(stdout," %d %s data\n",msst->info.np,msst->info.iso ? "scalar" : "tensor");

  return(1);
}


