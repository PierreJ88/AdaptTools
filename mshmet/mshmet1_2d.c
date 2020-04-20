#include "mshmet.h"
#include "eigen.h"

#define CTE2D    2.0 / 9.0

#define CTE_GEOPHY_HMIN  50.0
#define CTE_GEOPHY_HMAX 200.0


/* define 2d metric tensor field */
int defmet_2d(MSst *msst) {
  double   *m,*h,l[2],vp[2][2],dd,hmin,hmax;
  int       k;

  hmin = 1.0 / (msst->info.hmax*msst->info.hmax);
  hmax = 1.0 / (msst->info.hmin*msst->info.hmin);

  for (k=0; k<msst->info.np; k++) {
    h = &msst->sol.h[3*k];

    if ( !eigen_2d(h,l,vp) ) {
      if ( msst->info.verb != '0' )  fprintf(stdout," # Error: eigenvalue problem.\n");
      return(0);
    }
    /* truncation */
    l[0] = MS_MIN(MS_MAX(fabs(l[0]),hmin),hmax);
    l[1] = MS_MIN(MS_MAX(fabs(l[1]),hmin),hmax);

    if ( msst->info.iso ) {
      msst->sol.m[k] = (l[0] < l[1]) ? 1.0 / sqrt(l[1]) : 1.0 / sqrt(l[0]);
    }
    else {
      /* compute Mat = B M Bt */
      m = &msst->sol.m[3*k];
      m[0] = l[0]*vp[0][0]*vp[0][0] + l[1]*vp[1][0]*vp[1][0];
      m[1] = l[0]*vp[0][0]*vp[0][1] + l[1]*vp[1][0]*vp[1][1];
      m[2] = l[0]*vp[0][1]*vp[0][1] + l[1]*vp[1][1]*vp[1][1];
    }
  }

  return(1);
}

/* define 2d metric tensor field from a gradient vector */
int defgrad2met(MSst *msst) {
  pPoint    ppt;
  double   *m,*h,dd,rmin,hmin,hmax;
  int       k,i;
  double    l[2],vp[2][2],nrm_grad;
  double    *sol, sol_norm[msst->info.np], max_sol, min_sol;
  double    grad_norm[msst->info.np], max_grad, min_grad;

  if ( !msst->info.hmin ) {
    printf("No HMIN given chosen arbitrary to be :%lf\n",CTE_GEOPHY_HMIN);
    hmin = CTE_GEOPHY_HMIN;
  }
  else {
    hmin = msst->info.hmin;
    rmin = msst->info.hmin;
  }

  if ( !msst->info.rmin ) {
    printf("No RMIN given chosen arbitrary to be hmin :%lf\n", hmin);
    rmin = hmin;
  }
  else {
    rmin = msst->info.rmin;
  }

  if ( !msst->info.hmax ) {
    printf("No HMAX given chosen arbitrary to be :%lf\n",CTE_GEOPHY_HMAX);
    hmin = CTE_GEOPHY_HMIN;
    hmax = CTE_GEOPHY_HMAX;
  }
  else {
    hmax = msst->info.hmax;
  }

  /* Determine sol(min max) and grad(min max) */
    max_sol = 0.0;
    min_sol = 1.e200;
    max_grad = 0.0;
    min_grad = 1.e200;

    for (k=0; k<msst-> info.np; k++) {
      h = &msst->sol.g[msst->info.dim*k];
      sol = &msst->sol.u[k];

      max_sol = MS_MAX(*sol, max_sol);
      min_sol = MS_MIN(*sol, min_sol);

      nrm_grad = sqrt(h[0]*h[0] + h[1]*h[1]);

      max_grad = MS_MAX(nrm_grad, max_grad);
      min_grad = MS_MIN(nrm_grad, min_grad);
    }

    printf("max_sol =%lf\n",max_sol);
    printf("min_sol =%lf\n",min_sol);
    printf("max_grad =%lf\n",max_grad);
    printf("min_grad =%lf\n",min_grad);

    /* Normalize sol and grad */
    for (k=0; k<msst-> info.np; k++) {
      h = &msst->sol.g[msst->info.dim*k];
      sol = &msst->sol.u[k];

      sol_norm[k] = (*sol - min_sol) / (max_sol - min_sol);

      nrm_grad = sqrt(h[0]*h[0] + h[1]*h[1]);

      grad_norm[k] = (nrm_grad - min_grad) / (max_grad - min_grad);

      grad_norm[k] = pow(grad_norm[k]*0.7,msst->info.err);

      printf("grad k =%lf\n",grad_norm[k]);

    }

    if (msst->info.iso) {
      for (k=0; k<msst-> info.np; k++) {

	l[1] = hmax*sol_norm[k] + hmin*(1-sol_norm[k]);
	l[0] = rmin*grad_norm[k] + l[1]*(1-grad_norm[k]);

	/* l[0] = 1.0/(l[0]*l[0]); */
	/* l[1] = 1.0/(l[1]*l[1]); */

	msst->sol.m[k] = MS_MIN(l[0],l[1]);

      }
    }

    else {

      for (k=0; k<msst-> info.np; k++) {

	l[1] = hmax*sol_norm[k] + hmin*(1-sol_norm[k]);
	l[0] = rmin*grad_norm[k] + l[1]*(1-grad_norm[k]);

	l[0] = 1.0/(l[0]*l[0]);
	l[1] = 1.0/(l[1]*l[1]);


	h = &msst->sol.g[msst->info.dim*k];
	dd = sqrt(h[0]*h[0]+h[1]*h[1]);

	m = &msst->sol.m[3*k];
	if ( dd < MS_EPS ) {
	  m[0] = l[0];
	  m[1] = 0.;
	  m[2] = l[1];
	} else {
	  vp[0][0] = h[0]/dd;
	  vp[0][1] = h[1]/dd;

	  vp[1][0] = - vp[0][1];
	  vp[1][1] = vp[0][0];

	  /* compute Mat = B M Bt */
	  m[0] = l[0]*vp[0][0]*vp[0][0] + l[1]*vp[1][0]*vp[1][0];
	  m[1] = l[0]*vp[0][0]*vp[0][1] + l[1]*vp[1][0]*vp[1][1];
	  m[2] = l[0]*vp[0][1]*vp[0][1] + l[1]*vp[1][1]*vp[1][1];
	}
      }
    }

    return 1;
}


int nrmhes_2d(MSst *msst) {
  double   u,err;
  int      i,k;

  switch(msst->info.nrm) {
  case 0:  /* no normalization */
    for (k=0; k<msst->info.np; k++) {
      for (i=0; i<3; i++)  msst->sol.h[3*k+i] *= CTE2D / msst->info.err;
    }
    break;
  case 1:  /* relative value: M(u)= |H(u)| / err*||u||_inf */
  default:
    for (k=0; k<msst->info.np; k++) {
      for (i=0; i<3; i++)  msst->sol.h[3*k+i] *= CTE2D / (msst->info.err*msst->sol.umax);
    }
    break;
  case 2:  /* local norm: M(u)= |H(u)| / err*|u| */
    err = msst->sol.umax > 0.0 ? msst->sol.umax*0.01 : 0.01;
    for (k=0; k<msst->info.np; k++) {
      u = fabs(msst->sol.u[k]);
      for (i=0; i<3; i++)  msst->sol.h[3*k+i] *= CTE2D / MS_MAX(err,u);
    }
    break;
  }

  return(1);
}

int nrmgrad(MSst *msst) {
  double   u,err;
  int      i,k;

  switch(msst->info.nrm) {
  case 0:  /* no normalization */
  default:
    for (k=0; k<msst->info.np; k++) {
      for (i=0; i<msst->info.dim; i++)  msst->sol.g[msst->info.dim*k+i] /= msst->info.err;
    }
    break;
  }

  return(1);
}


/* least square approximation of hessian */
int hessls_2d(MSst *msst) {
  pPoint    p0,p1;
  pTria     pt;
  double   *ha,*hb,*ma,*mb,*g0,*h,a[6],u0,u1,dd;
  int       k,i;

  /* memory allocation */
  ha = (double*)calloc(msst->info.np*6,sizeof(double));
  assert(ha);
  hb = (double*)calloc(msst->info.np*3,sizeof(double));
  assert(hb);

  for (k=1; k<=msst->info.nt; k++) {
    pt = &msst->mesh.tria[k];

    for (i=0; i<3; i++) {
      p0 = &msst->mesh.point[pt->v[i]];
      u0 = msst->sol.u[pt->v[i]-1];
      g0 = &msst->sol.g[msst->info.dim*(pt->v[i]-1)];
      ma = &ha[6*(pt->v[i]-1)];
      mb = &hb[3*(pt->v[i]-1)];

      if ( pt->adj[(i+1)%3]/3 < k ) {
        memset(a,0,6*sizeof(double));
        p1 = &msst->mesh.point[pt->v[(i+2)%3]];
        u1 = msst->sol.u[pt->v[(i+2)%3]-1];

        a[0] += 0.5 * (p1->c[0]-p0->c[0]) * (p1->c[0]-p0->c[0]);
        a[1] += (p1->c[0]-p0->c[0]) * (p1->c[1]-p0->c[1]);
        a[2] += 0.5 * (p1->c[1]-p0->c[1]) * (p1->c[1]-p0->c[1]);
        dd = (u1-u0) - (p1->c[0]-p0->c[0])*g0[0] - (p1->c[1]-p0->c[1])*g0[1];

        /* M = At*A symmetric definite positive */
        ma[0] += a[0]*a[0];
        ma[1] += a[0]*a[1];
        ma[2] += a[1]*a[1];
        ma[3] += a[0]*a[2];
        ma[4] += a[1]*a[2];
        ma[5] += a[2]*a[2];

        /* c = At*b */
        mb[0] += a[0]*dd;
        mb[1] += a[1]*dd;
        mb[2] += a[2]*dd;
      }
      if ( pt->adj[(i+2)%3]/3 < k ) {
        memset(a,0,6*sizeof(double));
        p1 = &msst->mesh.point[pt->v[(i+1)%3]];
        u1 = msst->sol.u[pt->v[(i+1)%3]-1];

        a[0] += 0.5 * (p1->c[0]-p0->c[0]) * (p1->c[0]-p0->c[0]);
        a[1] += (p1->c[0]-p0->c[0]) * (p1->c[1]-p0->c[1]);
        a[2] += 0.5 * (p1->c[1]-p0->c[1]) * (p1->c[1]-p0->c[1]);
        dd = (u1-u0) - (p1->c[0]-p0->c[0])*g0[0] - (p1->c[1]-p0->c[1])*g0[1];

        /* M = At*A symmetric definite positive */
        ma[0] += a[0]*a[0];
        ma[1] += a[0]*a[1];
        ma[2] += a[1]*a[1];
        ma[3] += a[0]*a[2];
        ma[4] += a[1]*a[2];
        ma[5] += a[2]*a[2];

        /* c = At*b */
        mb[0] += a[0]*dd;
        mb[1] += a[1]*dd;
        mb[2] += a[2]*dd;
      }
    }
  }

  /* hessian evaluation */
  for (k=0; k<msst->info.np; k++) {
    ma = &ha[6*k];
    mb = &hb[3*k];

    /* direct solving */
    a[0] = ma[2]*ma[5] - ma[4]*ma[4];
    a[1] = ma[4]*ma[3] - ma[1]*ma[5];
    a[2] = ma[1]*ma[4] - ma[2]*ma[3];
    dd   = ma[0]*a[0] + ma[1]*a[1] + ma[3]*a[2];

    /* singular matrix */
    if ( fabs(dd) > MS_EPSD ) {
      h = &msst->sol.h[3*k];

      a[3] = ma[0]*ma[5] - ma[3]*ma[3];
      a[4] = ma[1]*ma[3] - ma[0]*ma[4];
      a[5] = ma[0]*ma[2] - ma[1]*ma[1];
      h[0] = (mb[0]*a[0] + mb[1]*a[1] + mb[2]*a[2]) / dd;
      h[1] = (mb[0]*a[1] + mb[1]*a[3] + mb[2]*a[4]) / dd;
      h[2] = (mb[0]*a[2] + mb[1]*a[4] + mb[2]*a[5]) / dd;
    }
    else {
	printf("Pt %d : SINGULAR hessian\n",k);
    }
  }
  free(ha);
  free(hb);

  return(1);
}


/* compute gradients at mesh vertices (least-square approx.) */
int gradls_2d(MSst *msst) {
  pTria     pt;
  pPoint    p0,p1;
  double   *ga,*gb,*a,*b,*g,u0,u1,dd;
  int       k,i;

  /* memory allocation */
  ga = (double*)calloc(msst->info.np*3,sizeof(double));
  assert(ga);
  gb = (double*)calloc(msst->info.np*2,sizeof(double));
  assert(gb);

  /* triangle contribution */
  for (k=1; k<=msst->info.nt; k++) {
    pt = &msst->mesh.tria[k];

    for (i=0; i<3; i++) {
      p0 = &msst->mesh.point[pt->v[i]];
      u0 = msst->sol.u[pt->v[i]-1];
      a  = &ga[3*(pt->v[i]-1)];
      b  = &gb[2*(pt->v[i]-1)];
      if ( pt->adj[(i+1)%3]/3 < k ) {
        p1 = &msst->mesh.point[pt->v[(i+2)%3]];
        u1 = msst->sol.u[pt->v[(i+2)%3]-1];
        /* M = At*A symmetric definite positive */
        a[0] += (p1->c[0]-p0->c[0]) * (p1->c[0]-p0->c[0]);
        a[1] += (p1->c[0]-p0->c[0]) * (p1->c[1]-p0->c[1]);
        a[2] += (p1->c[1]-p0->c[1]) * (p1->c[1]-p0->c[1]);
        //printf(" tria %d: point %d:  %e %e %e\n",k,pt->v[i],a[0],a[1],a[2]);
        /* b = A^t*du */
        b[0] += (p1->c[0]-p0->c[0]) * (u1-u0);
        b[1] += (p1->c[1]-p0->c[1]) * (u1-u0);
      }
      if ( pt->adj[(i+2)%3]/3 < k ) {
        p1 = &msst->mesh.point[pt->v[(i+1)%3]];
        u1 = msst->sol.u[pt->v[(i+1)%3]-1];

        /* M = At*A symmetric definite positive */
        a[0] += (p1->c[0]-p0->c[0]) * (p1->c[0]-p0->c[0]);
        a[1] += (p1->c[0]-p0->c[0]) * (p1->c[1]-p0->c[1]);
        a[2] += (p1->c[1]-p0->c[1]) * (p1->c[1]-p0->c[1]);

        /* b = A^t*du */
        b[0] += (p1->c[0]-p0->c[0]) * (u1-u0);
        b[1] += (p1->c[1]-p0->c[1]) * (u1-u0);
      }
    }
  }

  /* gradient evaluation */
  for (k=0; k<msst->info.np; k++) {
    /* solution of A(2,2)*grad(1,2) = b(1,2) */
    a  = &ga[3*k];
    b  = &gb[2*k];
    dd = a[0]*a[2] - a[1]*a[1];
    if ( fabs(dd) > MS_EPSD ) {
      g = &msst->sol.g[msst->info.dim*k];
      g[0] = (a[2]*b[0] - a[1]*b[1]) / dd;
      g[1] = (a[0]*b[1] - a[1]*b[0]) / dd;
    }
    else {
	printf("Singular grad pt %d\n",k);
    }
  }
  free(ga);
  free(gb);

  return(1);
}


int mshmet1_2d(MSst *msst) {
  int   k,ier;

  /* compute nodal gradients */
  msst->sol.g = (double*)calloc(msst->info.np,msst->info.dim*sizeof(double));
  assert(msst->sol.g);
  ier = gradls_2d(msst);
  if ( !ier ) {
    if ( msst->info.verb != '0' )  fprintf(stdout," # Error: unable to evaluate gradients\n");
    free(msst->sol.g);
    return(0);
  }

  if ( !msst->info.grad ) {
    /* compute nodal hessian matrix */
    msst->sol.h = (double*)calloc(msst->info.np,3*sizeof(double));
    assert(msst->sol.h);
    ier = hessls_2d(msst);
    free(msst->sol.g);
    if ( !ier ) {
      if ( msst->info.verb != '0' )  fprintf(stdout," # Error: unable to evaluate hessian\n");
      free(msst->sol.h);
      return(0);
    }


    /* normalize hessian */
    ier = nrmhes_2d(msst);
    if ( !ier ) {
      if ( msst->info.verb != '0' )  fprintf(stdout," # Error: unable to normalize hessian\n");
      free(msst->sol.h);
      return(0);
    }
    /* compute metric */
  if ( msst->info.iso )
    msst->sol.m = (double*)calloc(msst->info.np,sizeof(double));
  else
    msst->sol.m = (double*)calloc(msst->info.np,3*sizeof(double));
  assert(msst->sol.m);
  ier = defmet_2d(msst);
  free(msst->sol.h);
  if ( !ier ) {
    if ( msst->info.verb != '0' )  fprintf(stdout," # Error: unable to define metric\n");
    return(0);
  }

  }
  else {
    /* sensitivity to varations (gradient) */
    /*ier = nrmgrad(msst); */
    if ( !ier ) {
      if ( msst->info.verb != '0' )  fprintf(stdout," # Error: unable to compute sensitive gradiant\n");
      free(msst->sol.g);
      return(0);
    }
    if ( msst->info.iso )
      msst->sol.m = (double*)calloc(msst->info.np,sizeof(double));
    else
      msst->sol.m = (double*)calloc(msst->info.np,3*sizeof(double));
    assert(msst->sol.m);
    ier = defgrad2met(msst);
    free(msst->sol.g);
    if ( !ier ) {
      if ( msst->info.verb != '0' )  fprintf(stdout," # Error: unable to define metric\n");
      return(0);
    }
  }

  return(1);
}
