#include <vector>
#include <map>
#include <Rcpp.h>
using namespace std;



// http://docs.rexamine.com/R-devel/BLAS_8h.html#ab76361ff6b07fd3d58f2044f43a1cb93
// https://github.com/Microsoft/microsoft-r-open/blob/master/source/src/include/R_ext/Applic.h
//
#include <R_ext/Applic.h>



///////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018  Bowen Liu                                         //
//                                                                       //
// This program is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published by  //
// the Free Software Foundation, either version 3 of the License, or     //
// (at your option) any later version.                                   //
//                                                                       //
// This program is distributed in the hope that it will be useful,       //
// but WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
// GNU General Public License for more details.                          //
//                                                                       //
// You should have received a copy of the GNU General Public License     //
// along with this program.  If not, see <http://www.gnu.org/licenses/>. //
///////////////////////////////////////////////////////////////////////////






// https://github.com/Microsoft/microsoft-r-open/blob/master/source/src/appl/lbfgsb.c
/*
*  R : A Computer Language for Statistical Data Analysis
*  Copyright (C) 2000-2015 The R Core Team
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, a copy is available at
*  https://www.R-project.org/Licenses/
*/
/* l-bfgs-b.f -- translated by f2c (version 19991025).

From ?optim:
The code for method ‘"L-BFGS-B"’ is based on Fortran code by Zhu,
Byrd, Lu-Chen and Nocedal obtained from Netlib (file 'opt/lbfgs_bcm.shar')

The Fortran files contained no copyright information.

Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C.  (1995) A limited
memory algorithm for bound constrained optimization.
\emph{SIAM J. Scientific Computing}, \bold{16}, 1190--1208.
*/

/* <UTF8> all char uses here are ASCII */

#include <math.h>
#include <float.h> /* for DBL_EPSILON */
#include <string.h>
#include <R_ext/RS.h> /* for F77_CALL */
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Print.h> /* Rprintf */

#define FALSE_ 0
#define TRUE_ 1
#ifndef max
# define max(a, b) (a < b)?(b):(a)
# define min(a, b) (a > b)?(b):(a)
#endif


static int c__1 = 1;
static int c__11 = 11;

static void active(int, double *, double *, int *, double *, int *,
                   int, int *, int *, int *);
static void bmv(int, double *, double *, int *, double *, double *, int *);
static void cauchy(int, double *, double *,
                   double *, int *, double *, int *, int *,
                   double *, double *, double *, int, double *,
                   double *, double *, double *, double *, int * ,
                   int *, double *, double *, double *, double *,
                   int *, int, double *, int *, double *);
static void cmprlb(int, int, double *,
                   double *, double *, double *, double *,
                   double *, double *, double *, double *, int *,
                   double *, int *, int *, int *, int *,
                   int *);
static void dcsrch(double *, double *, double *,
                   double, double, double,
                   double, double, char *);
static void dcstep(double *, double *,
                   double *, double *, double *, double *,
                   double *, double *, double *, int *, double *,
                   double *);
static void errclb(int, int, double,
                   double *, double *, int *, char *, int *, int *);
static void formk(int, int *, int *, int *, int *, int *, int *,
                  int *, double *, double *, int, double *,
                  double *, double *, double *, int *, int *, int *);
static void formt(int, double *, double *,
                  double *, int *, double *, int *);
static void freev(int, int *, int *,
                  int *, int *, int *, int *, int *, int *,
                  int *, int, int *);
static void hpsolb(int, double *, int *, int);
static void lnsrlb(int, double *, double *,
                   int *, double *, double *, double *, double *,
                   double *, double *, double *, double *,
                   double *, double *, double *, double *,
                   double *, double *, double *, int *, int *,
                   int *, int *, int *, char *, int *, int *,
                   char *);
static void mainlb(int, int, double *,
                   double *, double *, int *, double *, double *,
                   double, double *, double *, double *,
                   double *, double *, double *, double *,
                   double *, double *, double *, double *,
                   double *, double *, int *, int *, int *, char *,
                   int, char *, int *);
static void matupd(int, int, double *, double *, double *,
                   double *, double *, double *, int *, int *,
                   int *, int *, double *, double *, double *,
                   double *, double *);
static void projgr(int, double *, double *,
                   int *, double *, double *, double *);
static void subsm(int, int, int *, int *, double *, double *,
                  int *, double *, double *, double *, double *,
                  double *, int *, int *, int *, double *,
                  double *, int, int *);

static void prn1lb(int n, int m, double *l, double *u, double *x,
                   int iprint, double epsmch);
static void prn2lb(int n, double *x, double *f, double *g, int iprint,
                   int iter, int nfgv, int nact, double sbgnrm,
                   int nint, char *word, int iword, int iback,
                   double stp, double xstep);
static void prn3lb(int n, double *x, double *f, char *task, int iprint,
                   int info, int iter, int nfgv, int nintol, int nskip,
                   int nact, double sbgnrm, int nint,
                   char *word, int iback, double stp, double xstep,
                   int k);


/* ================    L-BFGS-B (version 2.3)	========================== */
static
  void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
              double *f, double *g, double factr, double *pgtol,
              double *wa, int * iwa, char *task, int iprint, int *isave)
  {
 
    
    char csave[60];
    
    /* Local variables */
    static int lsnd, ld, lr, lt, lz, lwa, lwn, lss, lws, lwt, lsy, lwy;
    
    /* make sure csave is initialized */
    csave[0] = '\0';
    
    /* Parameter adjustments */
    --wa;
    
    if (strncmp(task, "START", 5) == 0) {
      lws = 1;
      lwy = lws + m * n;
      lsy = lwy + m * n;
      lss = lsy + m * m;
      lwt = lss + m * m;
      lwn = lwt + m * m;
      lsnd = lwn + (m * m << 2);
      lz = lsnd + (m * m << 2);
      lr = lz + n;
      ld = lr + n;
      lt = ld + n;
      lwa = lt + n;
    }
    
    mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol,
           &wa[lws], &wa[lwy], &wa[lsy],&wa[lss], &wa[lwt],&wa[lwn],
                                                              &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa],
                                                              iwa, &iwa[n], &iwa[n << 1], task, iprint,
                                                                                      csave, isave);
    return;
  } /* setulb */
    /* ======================= The end of setulb ============================= */
    
    static void mainlb(int n, int m, double *x,
                       double *l, double *u, int *nbd, double *f, double *g,
                       double factr, double *pgtol, double *ws, double * wy,
                       double *sy, double *ss, double *wt, double *wn,
                       double *snd, double *z, double *r, double *d,
                       double *t, double *wa, int *indx, int *iwhere,
                       int *indx2, char *task, int iprint,
                       char *csave, int *isave)
    {

      
      /* System generated locals */
      int ws_offset=0, wy_offset=0, sy_offset=0, ss_offset=0, wt_offset=0,
        wn_offset=0, snd_offset=0;
      double d__1, d__2;
      
      /* Local variables */
      double ddum;
      char word[4]; /* allow for terminator */
      int k = 0; /* -Wall */
      double xstep = 0.0; /* printed before being used */
      double dr, rr;
      int wrk;
      
      
      /* Parameter adjustments */
      --indx2;
      --iwhere;
      --indx;
      --t;
      --d;
      --r;
      --z;
      --g;
      --nbd;
      --u;
      --l;
      --x;
      --wa;
      // formerly lsave
      static int prjctd, cnstnd, boxed, updatd;
      // in isave
      static int nintol, iback, nskip, head, col, itail, iter, iupdat,
      nint, nfgv, info, ifun, iword, nfree, nact, ileave, nenter;
      // formerly dsave
      static double theta, fold, tol, dnorm, epsmch, gd, stpmx, sbgnrm,
      stp, gdold, dtd;
      
      /* Function Body */
      if (strncmp(task, "START", 5) == 0) {
        /*	  Generate the current machine precision. */
        epsmch = DBL_EPSILON;
        fold = 0.;
        dnorm = 0.;
        gd = 0.;
        sbgnrm = 0.;
        stp = 0.;
        xstep = 0.;
        stpmx = 0.;
        gdold = 0.;
        dtd = 0.;
        /*	  Initialize counters and scalars when task='START'. */
        /*	     for the limited memory BFGS matrices: */
        col = 0;
        head = 1;
        theta = 1.;
        iupdat = 0;
        updatd = FALSE_;
        iback = 0;
        itail = 0;
        ifun = 0;
        iword = 0;
        nact = 0;
        ileave = 0;
        nenter = 0;
        /*	     for operation counts: */
        iter = 0;
        nfgv = 0;
        nint = 0;
        nintol = 0;
        nskip = 0;
        nfree = n;
        /*	     for stopping tolerance: */
        tol = factr * epsmch;
        /*	     'word' records the status of subspace solutions. */
        strcpy(word, "---");
        /*	     'info' records the termination information. */
        info = 0;
        /*	  Check the input arguments for errors. */
        errclb(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k);
        if (strncmp(task, "ERROR", 5) == 0) {
          prn3lb(n, x+1, f, task, iprint, info, iter, nfgv, nintol, nskip,
                 nact, sbgnrm, nint, word, iback, stp, xstep, k);
          return;
        }
        
        prn1lb(n, m, l+1, u+1, x+1, iprint, epsmch);
        
        /*	  Initialize iwhere & project x onto the feasible set. */
        active(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], iprint, &prjctd,
               &cnstnd, &boxed);
        /*	  The end of the initialization. */
      } else {
        /*	After returning from the driver go to the point where execution */
        /*	is to resume. */
        if (strncmp(task, "FG_LN", 5) == 0)	goto L666;
        if (strncmp(task, "NEW_X", 5) == 0)     goto L777;
        if (strncmp(task, "FG_ST", 5) == 0)     goto L111;
        
        if (strncmp(task, "STOP", 4) == 0) {
          if (strncmp(task + 6, "CPU", 3) == 0) {
            /* restore the previous iterate. */
            F77_CALL(dcopy)(&n, &t[1], &c__1, &x[1], &c__1);
            F77_CALL(dcopy)(&n, &r[1], &c__1, &g[1], &c__1);
            *f = fold;
          }
          goto L999;
        }
      }
      /*     Compute f0 and g0. */
      strcpy(task, "FG_START");
      /*	    return to the driver to calculate f and g; reenter at 111. */
      goto L1000;
      L111:
        nfgv = 1;
      /*     Compute the infinity norm of the (-) projected gradient. */
      projgr(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
      
      if (iprint >= 1)
        Rprintf("At iterate %5d  f= %12.5g  |proj g|= %12.5g\n",
                iter, *f, sbgnrm);
      
      if (sbgnrm <= *pgtol) {
        /*				  terminate the algorithm. */
        strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
        goto L999;
      }
      /* ----------------- the beginning of the loop -------------------------- */
      L222:
        if (iprint >= 99) Rprintf("Iteration %5d\n", iter);
        iword = -1;
        
        if (! cnstnd && col > 0) {
          /*					      skip the search for GCP. */
          F77_CALL(dcopy)(&n, &x[1], &c__1, &z[1], &c__1);
          wrk = updatd;
          nint = 0;
          goto L333;
        }
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        
        /*     Compute the Generalized Cauchy Point (GCP). */
        
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        cauchy(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], &t[
                 1], &d[1], &z[1], m, &wy[wy_offset], &ws[ws_offset], &sy[
                 sy_offset], &wt[wt_offset], &theta, &col, &head, &wa[1], &wa[(m
                                                                                 << 1) + 1], &wa[(m << 2) + 1], &wa[m * 6 + 1], &nint, iprint, &
                 sbgnrm, &info, &epsmch);
        if (info != 0) {
          /*	   singular triangular system detected; refresh the lbfgs memory. */
          if (iprint >= 1)
            Rprintf("%s\n%s\n", "Singular triangular system detected;",
                    "   refresh the lbfgs memory and restart the iteration.");
          info = 0;
          col = 0;
          head = 1;
          theta = 1.;
          iupdat = 0;
          updatd = FALSE_;
          goto L222;
        }
        nintol += nint;
        /*     Count the entering and leaving variables for iter > 0; */
        /*     find the index set of free and active variables at the GCP. */
        freev(n, &nfree, &indx[1], &nenter, &ileave, &indx2[1], &iwhere[1], &
        wrk, &updatd, &cnstnd, iprint, &iter);
        nact = n - nfree;
        L333:
          /*     If there are no free variables or B=theta*I, then */
          /*					  skip the subspace minimization. */
          if (nfree == 0 || col == 0) {
            goto L555;
          }
          /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
          
          /*     Subspace minimization. */
          
          /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
          /*     Form  the LEL^T factorization of the indefinite */
          /*	 matrix	   K = [-D -Y'ZZ'Y/theta     L_a'-R_z'	] */
          /*		       [L_a -R_z	   theta*S'AA'S ] */
          /*	 where	   E = [-I  0] */
          /*		       [ 0  I] */
          if (wrk)
            formk(n, &nfree, &indx[1], &nenter, &ileave, &indx2[1], &iupdat, &
              updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset], &
              wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
          if (info != 0) {
            /*	    nonpositive definiteness in Cholesky factorization; */
            /*	    refresh the lbfgs memory and restart the iteration. */
            if (iprint >= 0)
              Rprintf("%s\n%s\n",
                      "Nonpositive definiteness in Cholesky factorization in formk;",
                      "   refresh the lbfgs memory and restart the iteration.");
            info = 0;
            col = 0;
            head = 1;
            theta = 1.;
            iupdat = 0;
            updatd = FALSE_;
            goto L222;
          }
          /*	  compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) */
          /*						     from 'cauchy'). */
          cmprlb(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset],
                 &sy[sy_offset], &wt[wt_offset], &z[1], &r[1], &wa[1], &indx[1],
                                                                            &theta, &col, &head, &nfree, &cnstnd, &info);
                                                                            if (info != 0)
                                                                              goto L444;
                                                                            /*	 call the direct method. */
                                                                            subsm(n, m, &nfree, &indx[1], &l[1], &u[1], &nbd[1], &z[1], &r[1], &
                                                                            ws[ws_offset], &wy[wy_offset], &theta, &col, &head, &iword, &wa[1]
                                                                                    , &wn[wn_offset], iprint, &info);
                                                                            L444:
                                                                              if (info != 0) {
                                                                                /*	    singular triangular system detected; */
                                                                                /*	    refresh the lbfgs memory and restart the iteration. */
                                                                                if (iprint >= 1)
                                                                                  Rprintf("%s\n%s\n", "Singular triangular system detected;",
                                                                                          "   refresh the lbfgs memory and restart the iteration.");
                                                                                info = 0;
                                                                                col = 0;
                                                                                head = 1;
                                                                                theta = 1.;
                                                                                iupdat = 0;
                                                                                updatd = FALSE_;
                                                                                goto L222;
                                                                              }
                                                                              L555:
                                                                                /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
                                                                                
                                                                                /*     Line search and optimality tests. */
                                                                                
                                                                                /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
                                                                                /*     Generate the search direction d:=z-x. */
                                                                                for (int i = 1; i <= n; ++i)
                                                                                  d[i] = z[i] - x[i];
                                                                              L666:
                                                                                lnsrlb(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1], &
                                                                                  d[1], &r[1], &t[1], &z[1], &stp, &dnorm, &dtd, &xstep, &
                                                                                  stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed, &cnstnd,
                                                                                  csave);
                                                                              if (info != 0 || iback >= 20) {
                                                                                /*	    restore the previous iterate. */
                                                                                F77_CALL(dcopy)(&n, &t[1], &c__1, &x[1], &c__1);
                                                                                F77_CALL(dcopy)(&n, &r[1], &c__1, &g[1], &c__1);
                                                                                *f = fold;
                                                                                if (col == 0) {
                                                                                  /*	       abnormal termination. */
                                                                                  if (info == 0) {
                                                                                    info = -9;
                                                                                    /*		  restore the actual number of f and g evaluations etc. */
                                                                                    --nfgv;
                                                                                    --ifun;
                                                                                    --iback;
                                                                                  }
                                                                                  strcpy(task, "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH");
                                                                                  ++iter;
                                                                                  goto L999;
                                                                                } else {
                                                                                  /*	       refresh the lbfgs memory and restart the iteration. */
                                                                                  if (iprint >= 1)
                                                                                    Rprintf("%s\n%s\n", "Bad direction in the line search;",
                                                                                            "   refresh the lbfgs memory and restart the iteration.");
                                                                                  if (info == 0)
                                                                                    --nfgv;
                                                                                  info = 0;
                                                                                  col = 0;
                                                                                  head = 1;
                                                                                  theta = 1.;
                                                                                  iupdat = 0;
                                                                                  updatd = FALSE_;
                                                                                  strcpy(task, "RESTART_FROM_LNSRCH");
                                                                                  goto L222;
                                                                                }
                                                                              } else if (strncmp(task, "FG_LN", 5) == 0) {
                                                                                /*	    return to the driver for calculating f and g; reenter at 666. */
                                                                                goto L1000;
                                                                              } else {
                                                                                /*	    calculate and print out the quantities related to the new X. */
                                                                                ++iter;
                                                                                /*	  Compute the infinity norm of the projected (-)gradient. */
                                                                                projgr(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
                                                                                /*	  Print iteration information. */
                                                                                prn2lb(n, x+1, f, g+1, iprint, iter, nfgv, nact,
                                                                                       sbgnrm, nint, word, iword, iback, stp, xstep);
                                                                                goto L1000;
                                                                              }
                                                                              L777:
                                                                                /*     Test for termination. */
                                                                                if (sbgnrm <= *pgtol) {
                                                                                  /*				  terminate the algorithm. */
                                                                                  strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
                                                                                  goto L999;
                                                                                }
                                                                                /* Computing MAX */
                                                                                d__1 = fabs(fold), d__2 = fabs(*f), d__1 = max(d__1, d__2);
                                                                                ddum = max(d__1, 1.);
                                                                                if (fold - *f <= tol * ddum) {
                                                                                  /*					  terminate the algorithm. */
                                                                                  strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
                                                                                  if (iback >= 10) info = -5;
                                                                                  /*	     i.e., to issue a warning if iback>10 in the line search. */
                                                                                  goto L999;
                                                                                }
                                                                                /*     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's. */
                                                                                for (int i = 1; i <= n; ++i)
                                                                                  r[i] = g[i] - r[i];
                                                                                rr = F77_CALL(ddot)(&n, &r[1], &c__1, &r[1], &c__1);
                                                                                if (stp == 1.) {
                                                                                  dr = gd - gdold;
                                                                                  ddum = -gdold;
                                                                                } else {
                                                                                  dr = (gd - gdold) * stp;
                                                                                  F77_CALL(dscal)(&n, &stp, &d[1], &c__1);
                                                                                  ddum = -gdold * stp;
                                                                                }
                                                                                if (dr <= epsmch * ddum) {
                                                                                  /*			      skip the L-BFGS update. */
                                                                                  ++nskip;
                                                                                  updatd = FALSE_;
                                                                                  if (iprint >= 1)
                                                                                    Rprintf("ys=%10.3e  -gs=%10.3e, BFGS update SKIPPED\n", dr, ddum);
                                                                                  goto L888;
                                                                                }
                                                                                /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
                                                                                
                                                                                /*     Update the L-BFGS matrix. */
                                                                                
                                                                                /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
                                                                                updatd = TRUE_;
                                                                                ++iupdat;
                                                                                /*     Update matrices WS and WY and form the middle matrix in B. */
                                                                                matupd(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], &ss[
                                                                                         ss_offset], &d[1], &r[1], &itail, &iupdat, &col, &head, &
                                                                                           theta, &rr, &dr, &stp, &dtd);
                                                                                /*     Form the upper half of the pds T = theta*SS + L*D^(-1)*L'; */
                                                                                /*	  Store T in the upper triangular of the array wt; */
                                                                                /*	  Cholesky factorize T to J*J' with */
                                                                                /*	     J' stored in the upper triangular of wt. */
                                                                                formt(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col, &theta, &
                                                                                info);
                                                                                if (info != 0) {
                                                                                  /*	    nonpositive definiteness in Cholesky factorization; */
                                                                                  /*	    refresh the lbfgs memory and restart the iteration. */
                                                                                  if (iprint >= 0)
                                                                                    Rprintf("%s\n%s\n",
                                                                                            "Nonpositive definiteness in Cholesky factorization in formk;",
                                                                                            "   refresh the lbfgs memory and restart the iteration.");
                                                                                  info = 0;
                                                                                  col = 0;
                                                                                  head = 1;
                                                                                  theta = 1.;
                                                                                  iupdat = 0;
                                                                                  updatd = FALSE_;
                                                                                  goto L222;
                                                                                }
                                                                                /*     Now the inverse of the middle matrix in B is */
                                                                                /*	 [  D^(1/2)	 O ] [ -D^(1/2)	 D^(-1/2)*L' ] */
                                                                                /*	 [ -L*D^(-1/2)	 J ] [	0	 J'	     ] */
                                                                                L888:
                                                                                  /* -------------------- the end of the loop ----------------------------- */
                                                                                  goto L222;
                                                                                L999:
                                                                                  L1000:
                                                                                  /*     Save local variables. */
                                                                                  //    isave[1-1] = nintol;
                                                                                  //    isave[3-1] = itfile;
                                                                                  //    isave[4-1] = iback;
                                                                                  //    isave[5-1] = nskip;
                                                                                  //    isave[6-1] = head;
                                                                                  //    isave[7-1] = col;
                                                                                  //    isave[8-1] = itail;
                                                                                  //    isave[9-1] = iter;
                                                                                  //    isave[10-1] = iupdat;
                                                                                  //    isave[12-1] = nint;
                                                                                  isave[13-1] = nfgv;
                                                                                //    isave[14-1] = info;
                                                                                //    isave[15-1] = ifun;
                                                                                //    isave[16-1] = iword;
                                                                                //    isave[17-1] = nfree;
                                                                                //    isave[18-1] = nact;
                                                                                //    isave[19-1] = ileave;
                                                                                //    isave[20-1] = nenter;
                                                                                
                                                                                prn3lb(n, x+1, f, task, iprint, info, iter, nfgv, nintol, nskip, nact,
                                                                                       sbgnrm, nint, word, iback, stp, xstep, k);
                                                                                return;
    } /* mainlb */
                                                                                  /* ======================= The end of mainlb ============================= */
                                                                                  
                                                                                  static void active(int n, double *l, double *u,
                                                                                                     int *nbd, double *x, int *iwhere, int iprint,
                                                                                                     int *prjctd, int *cnstnd, int *boxed)
                                                                                  {

                                                                                    /* Local variables */
                                                                                    int nbdd;
                                                                                    
                                                                                    /* Parameter adjustments */
                                                                                    --iwhere;
                                                                                    --x;
                                                                                    --nbd;
                                                                                    --u;
                                                                                    --l;
                                                                                    
                                                                                    /* Function Body */
                                                                                    
                                                                                    /*Initialize nbdd, prjctd, cnstnd and boxed. */
                                                                                    nbdd = 0;
                                                                                    *prjctd = FALSE_;
                                                                                    *cnstnd = FALSE_;
                                                                                    *boxed = TRUE_;
                                                                                    /*     Project the initial x to the easible set if necessary. */
                                                                                    for (int i = 1; i <= n; ++i) {
                                                                                      if (nbd[i] > 0) {
                                                                                        if (nbd[i] <= 2 && x[i] <= l[i]) {
                                                                                          if (x[i] < l[i]) {
                                                                                            *prjctd = TRUE_;
                                                                                            x[i] = l[i];
                                                                                          }
                                                                                          ++nbdd;
                                                                                        } else if (nbd[i] >= 2 && x[i] >= u[i]) {
                                                                                          if (x[i] > u[i]) {
                                                                                            *prjctd = TRUE_;
                                                                                            x[i] = u[i];
                                                                                          }
                                                                                          ++nbdd;
                                                                                        }
                                                                                      }
                                                                                    }
                                                                                    
                                                                                    /*     Initialize iwhere and assign values to cnstnd and boxed. */
                                                                                    for (int i = 1; i <= n; ++i) {
                                                                                      if (nbd[i] != 2)
                                                                                        *boxed = FALSE_;
                                                                                      if (nbd[i] == 0) {
                                                                                        /*				  this variable is always free */
                                                                                        iwhere[i] = -1;
                                                                                        /*	     otherwise set x(i)=mid(x(i), u(i), l(i)). */
                                                                                      } else {
                                                                                        *cnstnd = TRUE_;
                                                                                        if (nbd[i] == 2 && u[i] - l[i] <= 0.) {
                                                                                          /*		     this variable is always fixed */
                                                                                          iwhere[i] = 3;
                                                                                        } else
                                                                                          iwhere[i] = 0;
                                                                                      }
                                                                                    }
                                                                                    if (iprint >= 0) {
                                                                                      if (*prjctd)
                                                                                        Rprintf("The initial X is infeasible.  Restart with its projection.\n");
                                                                                      if (!*cnstnd) Rprintf("This problem is unconstrained.\n");
                                                                                    }
                                                                                    if (iprint > 0)
                                                                                      Rprintf("At X0, %d variables are exactly at the bounds\n", nbdd);
                                                                                    
                                                                                    return;
                                                                                  } /* active */
                                                                                          /* ======================= The end of active ============================= */
                                                                                          
                                                                                          static void bmv(int m, double *sy, double *wt,
                                                                                                          int *col, double *v, double *p, int *info)
                                                                                          {

                                                                                            /* System generated locals */
                                                                                            int sy_dim1, sy_offset, wt_dim1, wt_offset, Col;
                                                                                            
                                                                                            
                                                                                            /* Local variables */
                                                                                            int i2, k;
                                                                                            double sum;
                                                                                            
                                                                                            /* Parameter adjustments */
                                                                                            wt_dim1 = m;
                                                                                            wt_offset = 1 + wt_dim1 * 1;
                                                                                            wt -= wt_offset;
                                                                                            sy_dim1 = m;
                                                                                            sy_offset = 1 + sy_dim1 * 1;
                                                                                            sy -= sy_offset;
                                                                                            --p;
                                                                                            --v;
                                                                                            
                                                                                            /* Function Body */
                                                                                            if (*col == 0) {
                                                                                              return;
                                                                                            }
                                                                                            /*	PART I: solve [	 D^(1/2)      O ] [ p1 ] = [ v1 ]
                                                                                            *		      [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
                                                                                            *	solve Jp2=v2+LD^(-1)v1.
                                                                                            */
                                                                                            Col = *col;
                                                                                            p[*col + 1] = v[*col + 1];
                                                                                            for (int i = 2; i <= Col; ++i) {
                                                                                              i2 = *col + i;
                                                                                              sum = 0.;
                                                                                              for (int k = 1; k <= i - 1; ++k)
                                                                                                sum += sy[i + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
                                                                                              p[i2] = v[i2] + sum;
                                                                                            }
                                                                                            /*     Solve the triangular system */
                                                                                            F77_CALL(dtrsl)(&wt[wt_offset], &m, col, &p[*col + 1], &c__11, info);
                                                                                            if (*info != 0) {
                                                                                              return;
                                                                                            }
                                                                                            /*	 solve D^(1/2)p1=v1. */
                                                                                            for (int i = 1; i <= Col; ++i)
                                                                                              p[i] = v[i] / sqrt(sy[i + i * sy_dim1]);
                                                                                            
                                                                                            /*	PART II: solve [ -D^(1/2)   D^(-1/2)*L'	 ] [ p1 ] = [ p1 ]
                                                                                            *		       [  0	    J'		 ] [ p2 ]   [ p2 ].
                                                                                            *	solve J^Tp2=p2.
                                                                                            */
                                                                                            F77_CALL(dtrsl)(&wt[wt_offset], &m, col, &p[*col + 1], &c__1, info);
                                                                                            if (*info != 0)
                                                                                              return;
                                                                                            
                                                                                            /*	 compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2) */
                                                                                            /*		   =-D^(-1/2)p1 + D^(-1)L'p2. */
                                                                                            for (int i = 1; i <= Col; ++i) {
                                                                                              p[i] = -p[i] / sqrt(sy[i + i * sy_dim1]);
                                                                                            }
                                                                                            for (int i = 1; i <= Col; ++i) {
                                                                                              sum = 0.;
                                                                                              for (k = i + 1; k <= Col; ++k) {
                                                                                                sum += sy[k + i * sy_dim1] * p[*col + k] / sy[i + i * sy_dim1];
                                                                                              }
                                                                                              p[i] += sum;
                                                                                            }
                                                                                            return;
                                                                                          } /* bmv */
                                                                                            /* ======================== The end of bmv =============================== */
                                                                                            
                                                                                            static void cauchy(int n, double *x, double *l, double *u, int *nbd,
                                                                                                               double *g, int *iorder, int * iwhere, double *t,
                                                                                                               double *d, double *xcp, int m,
                                                                                                               double *wy, double *ws, double *sy, double *wt,
                                                                                                               double *theta, int *col, int *head, double *p,
                                                                                                               double *c, double *wbp, double *v, int *nint,
                                                                                                               int iprint, double *sbgnrm, int *info, double * epsmch)
                                                                                            {

                                                                                              /* System generated locals */
                                                                                              int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset,
                                                                                              wt_dim1, wt_offset, i__2;
                                                                                              double d__1;
                                                                                              
                                                                                              /* Local variables */
                                                                                              double bkmin, dibp, dibp2, zibp, neggi, tsum;
                                                                                              double f1, f2, f2_org__, dt, tj, tj0, tl= 0.0, tu=0.0, dtm, wmc, wmp, wmw;
                                                                                              
                                                                                              int ibp, iter, bnded, nfree, nleft, nbreak, ibkmin, pointr;
                                                                                              int xlower, xupper, col2;
                                                                                              
                                                                                              /* Parameter adjustments */
                                                                                              --xcp;
                                                                                              --d;
                                                                                              --t;
                                                                                              --iwhere;
                                                                                              --iorder;
                                                                                              --g;
                                                                                              --nbd;
                                                                                              --u;
                                                                                              --l;
                                                                                              --x;
                                                                                              --v;
                                                                                              --wbp;
                                                                                              --c;
                                                                                              --p;
                                                                                              wt_dim1 = m;    wt_offset = 1 + wt_dim1 * 1;    wt -= wt_offset;
                                                                                              sy_dim1 = m;    sy_offset = 1 + sy_dim1 * 1;    sy -= sy_offset;
                                                                                              ws_dim1 = n;    ws_offset = 1 + ws_dim1 * 1;    ws -= ws_offset;
                                                                                              wy_dim1 = n;    wy_offset = 1 + wy_dim1 * 1;    wy -= wy_offset;
                                                                                              
                                                                                              /* Function Body */
                                                                                              
                                                                                              /*     Check the status of the variables, reset iwhere(i) if necessary;
                                                                                              *	 compute the Cauchy direction d and the breakpoints t; initialize
                                                                                              *	 the derivative f1 and the vector p = W'd (for theta = 1).
                                                                                              */
                                                                                              
                                                                                              if (*sbgnrm <= 0.) {
                                                                                                if (iprint >= 0) Rprintf("Subgnorm = 0.  GCP = X.\n");
                                                                                                F77_CALL(dcopy)(&n, &x[1], &c__1, &xcp[1], &c__1);
                                                                                                return;
                                                                                              }
                                                                                              bnded = TRUE_;
                                                                                              nfree = n + 1;
                                                                                              nbreak = 0;
                                                                                              ibkmin = 0;
                                                                                              bkmin = 0.;
                                                                                              col2 = *col << 1;
                                                                                              f1 = 0.;
                                                                                              if (iprint >= 99)
                                                                                                Rprintf("\n---------------- CAUCHY entered-------------------\n\n");
                                                                                              
                                                                                              /*     We set p to zero and build it up as we determine d. */
                                                                                              for (int i = 1; i <= col2; ++i)
                                                                                                p[i] = 0.;
                                                                                              
                                                                                              /*     In the following loop we determine for each variable its bound */
                                                                                              /*	  status and its breakpoint, and update p accordingly. */
                                                                                              /*	  Smallest breakpoint is identified. */
                                                                                              
                                                                                              for (int i = 1; i <= n; ++i) {
                                                                                                neggi = -g[i];
                                                                                                if (iwhere[i] != 3 && iwhere[i] != -1) {
                                                                                                  /*	       if x(i) is not a constant and has bounds, */
                                                                                                  /*	       compute the difference between x(i) and its bounds. */
                                                                                                  if (nbd[i] <= 2) {
                                                                                                    tl = x[i] - l[i];
                                                                                                  }
                                                                                                  if (nbd[i] >= 2) {
                                                                                                    tu = u[i] - x[i];
                                                                                                  }
                                                                                                  /*	     If a variable is close enough to a bound */
                                                                                                  /*	       we treat it as at bound. */
                                                                                                  xlower = nbd[i] <= 2 && tl <= 0.;
                                                                                                  xupper = nbd[i] >= 2 && tu <= 0.;
                                                                                                  /*		reset iwhere(i). */
                                                                                                  iwhere[i] = 0;
                                                                                                  if (xlower) {
                                                                                                    if (neggi <= 0.) {
                                                                                                      iwhere[i] = 1;
                                                                                                    }
                                                                                                  } else if (xupper) {
                                                                                                    if (neggi >= 0.) {
                                                                                                      iwhere[i] = 2;
                                                                                                    }
                                                                                                  } else {
                                                                                                    if (fabs(neggi) <= 0.) {
                                                                                                      iwhere[i] = -3;
                                                                                                    }
                                                                                                  }
                                                                                                }
                                                                                                pointr = *head;
                                                                                                if (iwhere[i] != 0 && iwhere[i] != -1) {
                                                                                                  d[i] = 0.;
                                                                                                } else {
                                                                                                  d[i] = neggi;
                                                                                                  f1 -= neggi * neggi;
                                                                                                  /*	       calculate p := p - W'e_i* (g_i). */
                                                                                                  i__2 = *col;
                                                                                                  for (int j = 1; j <= i__2; ++j) {
                                                                                                    p[j] += wy[i + pointr * wy_dim1] * neggi;
                                                                                                    p[*col + j] += ws[i + pointr * ws_dim1] * neggi;
                                                                                                    pointr = pointr % m + 1;
                                                                                                  }
                                                                                                  if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.) {
                                                                                                    /*				   x(i) + d(i) is bounded; compute t(i). */
                                                                                                    ++nbreak;
                                                                                                    iorder[nbreak] = i;
                                                                                                    t[nbreak] = tl / (-neggi);
                                                                                                    if (nbreak == 1 || t[nbreak] < bkmin) {
                                                                                                      bkmin = t[nbreak];
                                                                                                      ibkmin = nbreak;
                                                                                                    }
                                                                                                  } else if (nbd[i] >= 2 && neggi > 0.) {
                                                                                                    /*				   x(i) + d(i) is bounded; compute t(i). */
                                                                                                    ++nbreak;
                                                                                                    iorder[nbreak] = i;
                                                                                                    t[nbreak] = tu / neggi;
                                                                                                    if (nbreak == 1 || t[nbreak] < bkmin) {
                                                                                                      bkmin = t[nbreak];
                                                                                                      ibkmin = nbreak;
                                                                                                    }
                                                                                                  } else {/*		  x(i) + d(i) is not bounded. */
                                                                                                    --nfree;
                                                                                                    iorder[nfree] = i;
                                                                                                    if (fabs(neggi) > 0.)
                                                                                                      bnded = FALSE_;
                                                                                                  }
                                                                                                }
                                                                                              } /* for(i = 1:n) */
                                                                                                    
                                                                                                    /*     The indices of the nonzero components of d are now stored */
                                                                                                    /*	 in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
                                                                                                    /*	 The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
                                                                                                    if (*theta != 1.) {
                                                                                                      /*		     complete the initialization of p for theta not= one. */
                                                                                                      F77_CALL(dscal)(col, theta, &p[*col + 1], &c__1);
                                                                                                    }
                                                                                                    /*     Initialize GCP xcp = x. */
                                                                                                    F77_CALL(dcopy)(&n, &x[1], &c__1, &xcp[1], &c__1);
                                                                                                    if (nbreak == 0 && nfree == n + 1) {
                                                                                                      /*		    is a zero vector, return with the initial xcp as GCP. */
                                                                                                      if (iprint > 100) {
                                                                                                        Rprintf("Cauchy X =  ");
                                                                                                        for (int i = 1; i <= n; i++) Rprintf("%g ", xcp[i]);
                                                                                                        Rprintf("\n");
                                                                                                      }
                                                                                                      return;
                                                                                                    }
                                                                                                    /*     Initialize c = W'(xcp - x) = 0. */
                                                                                                    for (int j = 1; j <= col2; ++j)
                                                                                                      c[j] = 0.;
                                                                                                    
                                                                                                    /*     Initialize derivative f2. */
                                                                                                    f2 = -(*theta) * f1;
                                                                                                    f2_org__ = f2;
                                                                                                    if (*col > 0) {
                                                                                                      bmv(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
                                                                                                      if (*info != 0) {
                                                                                                        return;
                                                                                                      }
                                                                                                      f2 -= F77_CALL(ddot)(&col2, &v[1], &c__1, &p[1], &c__1);
                                                                                                    }
                                                                                                    dtm = -f1 / f2;
                                                                                                    tsum = 0.;
                                                                                                    *nint = 1;
                                                                                                    if (iprint >= 99) Rprintf("There are %d  breakpoints\n", nbreak);
                                                                                                    
                                                                                                    /*     If there are no breakpoints, locate the GCP and return. */
                                                                                                    if (nbreak == 0) {
                                                                                                      goto L888;
                                                                                                    }
                                                                                                    nleft = nbreak;
                                                                                                    iter = 1;
                                                                                                    tj = 0.;
                                                                                                    /* ------------------- the beginning of the loop ------------------------- */
                                                                                                    L777:
                                                                                                      /*     Find the next smallest breakpoint; */
                                                                                                      /*	 compute dt = t(nleft) - t(nleft + 1). */
                                                                                                      tj0 = tj;
                                                                                                    if (iter == 1) {
                                                                                                      /*	   Since we already have the smallest breakpoint we need not do */
                                                                                                      /*	   heapsort yet. Often only one breakpoint is used and the */
                                                                                                      /*	   cost of heapsort is avoided. */
                                                                                                      tj = bkmin;
                                                                                                      ibp = iorder[ibkmin];
                                                                                                    } else {
                                                                                                      if (iter == 2) {
                                                                                                        /* Replace the already used smallest breakpoint with the */
                                                                                                        /* breakpoint numbered nbreak > nlast, before heapsort call. */
                                                                                                        if (ibkmin != nbreak) {
                                                                                                          t[ibkmin] = t[nbreak];
                                                                                                          iorder[ibkmin] = iorder[nbreak];
                                                                                                        }
                                                                                                      }
                                                                                                      /* Update heap structure of breakpoints */
                                                                                                      /* (if iter=2, initialize heap). */
                                                                                                      hpsolb(nleft, &t[1], &iorder[1], iter - 2);
                                                                                                      tj = t[nleft];
                                                                                                      ibp = iorder[nleft];
                                                                                                    }
                                                                                                    dt = tj - tj0;
                                                                                                    
                                                                                                    if (dt != 0 && iprint >=  100) {
                                                                                                      Rprintf("\nPiece    %3i f1, f2 at start point %11.4e %11.4e\n",
                                                                                                              *nint, f1, f2);
                                                                                                      Rprintf("Distance to the next break point =  %11.4e\n", dt);
                                                                                                      Rprintf("Distance to the stationary point =  %11.4e\n", dtm);
                                                                                                    }
                                                                                                    
                                                                                                    /*     If a minimizer is within this interval, */
                                                                                                    /*	 locate the GCP and return. */
                                                                                                    if (dtm < dt) {
                                                                                                      goto L888;
                                                                                                    }
                                                                                                    /*     Otherwise fix one variable and */
                                                                                                    /*	 reset the corresponding component of d to zero. */
                                                                                                    tsum += dt;
                                                                                                    --nleft;
                                                                                                    ++iter;
                                                                                                    dibp = d[ibp];
                                                                                                    d[ibp] = 0.;
                                                                                                    if (dibp > 0.) {
                                                                                                      zibp = u[ibp] - x[ibp];
                                                                                                      xcp[ibp] = u[ibp];
                                                                                                      iwhere[ibp] = 2;
                                                                                                    } else {
                                                                                                      zibp = l[ibp] - x[ibp];
                                                                                                      xcp[ibp] = l[ibp];
                                                                                                      iwhere[ibp] = 1;
                                                                                                    }
                                                                                                    if (iprint >= 100) Rprintf("Variable  %d  is fixed.\n", ibp);
                                                                                                    if (nleft == 0 && nbreak == n) {
                                                                                                      /*					       all n variables are fixed, */
                                                                                                      /*						  return with xcp as GCP. */
                                                                                                      dtm = dt;
                                                                                                      goto L999;
                                                                                                    }
                                                                                                    /*     Update the derivative information. */
                                                                                                    ++(*nint);
                                                                                                    dibp2 = dibp * dibp;
                                                                                                    /*     Update f1 and f2. */
                                                                                                    /*	  temporarily set f1 and f2 for col=0. */
                                                                                                    f1 += dt * f2 + dibp2 - *theta * dibp * zibp;
                                                                                                    f2 -= *theta * dibp2;
                                                                                                    if (*col > 0) {
                                                                                                      /*			    update c = c + dt*p. */
                                                                                                      F77_CALL(daxpy)(&col2, &dt, &p[1], &c__1, &c[1], &c__1);
                                                                                                      /*	     choose wbp, */
                                                                                                      /*	     the row of W corresponding to the breakpoint encountered. */
                                                                                                      pointr = *head;
                                                                                                      for (int j = 1; j <= *col; ++j) {
                                                                                                        wbp[j] = wy[ibp + pointr * wy_dim1];
                                                                                                        wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
                                                                                                        pointr = pointr % m + 1;
                                                                                                      }
                                                                                                      /*	     compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
                                                                                                      bmv(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
                                                                                                      if (*info != 0) {
                                                                                                        return;
                                                                                                      }
                                                                                                      wmc = F77_CALL(ddot)(&col2,  &c[1], &c__1, &v[1], &c__1);
                                                                                                      wmp = F77_CALL(ddot)(&col2,  &p[1], &c__1, &v[1], &c__1);
                                                                                                      wmw = F77_CALL(ddot)(&col2,&wbp[1], &c__1, &v[1], &c__1);
                                                                                                      /*	     update p = p - dibp*wbp. */
                                                                                                      d__1 = -dibp;
                                                                                                      F77_CALL(daxpy)(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
                                                                                                      /*	     complete updating f1 and f2 while col > 0. */
                                                                                                      f1 += dibp * wmc;
                                                                                                      f2 += (2. * dibp * wmp - dibp2 * wmw);
                                                                                                    }
                                                                                                    if(f2 < (d__1 = *epsmch * f2_org__)) f2 = d__1;
                                                                                                    if (nleft > 0) {
                                                                                                      dtm = -f1 / f2;
                                                                                                      goto L777;
                                                                                                      /*		   to repeat the loop for unsearched intervals. */
                                                                                                    } else if (bnded) {
                                                                                                      f1 = 0.;
                                                                                                      f2 = 0.;
                                                                                                      dtm = 0.;
                                                                                                    } else {
                                                                                                      dtm = -f1 / f2;
                                                                                                    }
                                                                                                    /* ------------------- the end of the loop ------------------------------- */
                                                                                                    L888:
                                                                                                      if (iprint >= 99) {
                                                                                                        Rprintf("\nGCP found in this segment\n");
                                                                                                        Rprintf("Piece    %3i f1, f2 at start point %11.4e %11.4e\n",
                                                                                                                *nint,f1,f2);
                                                                                                        Rprintf("Distance to the stationary point =  %11.4e\n", dtm);
                                                                                                      }
                                                                                                      
                                                                                                      if (dtm <= 0.) {
                                                                                                        dtm = 0.;
                                                                                                      }
                                                                                                      tsum += dtm;
                                                                                                      /*     Move free variables (i.e., the ones w/o breakpoints) and */
                                                                                                      /*	 the variables whose breakpoints haven't been reached. */
                                                                                                      F77_CALL(daxpy)(&n, &tsum, &d[1], &c__1, &xcp[1], &c__1);
                                                                                                      L999:
                                                                                                        /*     Update c = c + dtm*p = W'(x^c - x) */
                                                                                                        /*	 which will be used in computing r = Z'(B(x^c - x) + g). */
                                                                                                        if (*col > 0) {
                                                                                                          F77_CALL(daxpy)(&col2, &dtm, &p[1], &c__1, &c[1], &c__1);
                                                                                                        }
                                                                                                        if (iprint >= 100) {
                                                                                                          Rprintf("Cauchy X =  ");
                                                                                                          for (int i = 1; i <= n; i++) Rprintf("%g ", xcp[i]);
                                                                                                          Rprintf("\n");
                                                                                                        }
                                                                                                        
                                                                                                        if (iprint >= 99)
                                                                                                          Rprintf("\n---------------- exit CAUCHY----------------------\n\n");
                                                                                                        return;
                                                                                            } /* cauchy */
                                                                                                        /* ====================== The end of cauchy ============================== */
                                                                                                        
                                                                                                        static void cmprlb(int n, int m, double *x,
                                                                                                                           double *g, double *ws, double *wy, double *sy,
                                                                                                                           double *wt, double *z, double *r, double *wa,
                                                                                                                           int *indx, double *theta, int *col, int *head,
                                                                                                                           int *nfree, int *cnstnd, int *info)
                                                                                                        {

                                                                                                          /* System generated locals */
                                                                                                          int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset,
                                                                                                          wt_dim1, wt_offset, Col, n_f;
                                                                                                          
                                                                                                          /* Local variables */
                                                                                                          int k;
                                                                                                          double a1, a2;
                                                                                                          int pointr;
                                                                                                          
                                                                                                          /* Parameter adjustments */
                                                                                                          --indx;
                                                                                                          --r;
                                                                                                          --z;
                                                                                                          --g;
                                                                                                          --x;
                                                                                                          --wa;
                                                                                                          wt_dim1 = m;
                                                                                                          wt_offset = 1 + wt_dim1 * 1;
                                                                                                          wt -= wt_offset;
                                                                                                          sy_dim1 = m;
                                                                                                          sy_offset = 1 + sy_dim1 * 1;
                                                                                                          sy -= sy_offset;
                                                                                                          wy_dim1 = n;
                                                                                                          wy_offset = 1 + wy_dim1 * 1;
                                                                                                          wy -= wy_offset;
                                                                                                          ws_dim1 = n;
                                                                                                          ws_offset = 1 + ws_dim1 * 1;
                                                                                                          ws -= ws_offset;
                                                                                                          
                                                                                                          /* Function Body */
                                                                                                          Col = *col;
                                                                                                          if (! (*cnstnd) && Col > 0) {
                                                                                                            for (int i = 1; i <= n; ++i)
                                                                                                              r[i] = -g[i];
                                                                                                          }
                                                                                                          else {
                                                                                                            n_f = *nfree;
                                                                                                            for (int i = 1; i <= n_f; ++i) {
                                                                                                              k = indx[i];
                                                                                                              r[i] = -(*theta) * (z[k] - x[k]) - g[k];
                                                                                                            }
                                                                                                            bmv(m, &sy[sy_offset], &wt[wt_offset], col,
                                                                                                                &wa[(m << 1) + 1], &wa[1], info);
                                                                                                            if (*info != 0) {
                                                                                                              *info = -8;
                                                                                                              return;
                                                                                                            }
                                                                                                            pointr = *head;
                                                                                                            for (int j = 1; j <= Col; ++j) {
                                                                                                              a1 = wa[j];
                                                                                                              a2 = *theta * wa[Col + j];
                                                                                                              for (int i = 1; i <= n_f; ++i) {
                                                                                                                k = indx[i];
                                                                                                                r[i] += wy[k + pointr * wy_dim1] * a1 +
                                                                                                                  ws[k + pointr * ws_dim1] * a2;
                                                                                                              }
                                                                                                              pointr = pointr % m + 1;
                                                                                                            }
                                                                                                          }
                                                                                                          return;
                                                                                                        } /* cmprlb */
                                                                                                          /* ======================= The end of cmprlb ============================= */
                                                                                                          
                                                                                                          static void errclb(int n, int m, double factr, double *l, double *u,
                                                                                                                             int *nbd, char *task, int *info, int *k)
                                                                                                          {

                                                                                                            /* Parameter adjustments */
                                                                                                            --nbd;
                                                                                                            --u;
                                                                                                            --l;
                                                                                                            
                                                                                                            /* Function Body */
                                                                                                            /*     Check the input arguments for errors. */
                                                                                                            if (n <= 0)
                                                                                                              strcpy(task, "ERROR: N .LE. 0");
                                                                                                            if (m <= 0)
                                                                                                              strcpy(task, "ERROR: M .LE. 0");
                                                                                                            if (factr < 0.)
                                                                                                              strcpy(task, "ERROR: FACTR .LT. 0");
                                                                                                            
                                                                                                            /*     Check the validity of the arrays nbd(i), u(i), and l(i). */
                                                                                                            for (int i = 1; i <= n; ++i) {
                                                                                                              if (nbd[i] < 0 || nbd[i] > 3) {
                                                                                                                /*						     return */
                                                                                                                strcpy(task, "ERROR: INVALID NBD");
                                                                                                                *info = -6;
                                                                                                                *k = i;
                                                                                                              }
                                                                                                              if (nbd[i] == 2) {
                                                                                                                if (l[i] > u[i]) {
                                                                                                                  /*				      return */
                                                                                                                  strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
                                                                                                                  *info = -7;
                                                                                                                  *k = i;
                                                                                                                }
                                                                                                              }
                                                                                                            }
                                                                                                            return;
                                                                                                          } /* errclb */
                                                                                                                  /* ======================= The end of errclb ============================= */
                                                                                                                  
                                                                                                                  static void formk(int n, int *nsub, int *ind, int * nenter, int *ileave,
                                                                                                                                    int *indx2, int *iupdat, int * updatd, double *wn,
                                                                                                                                    double *wn1, int m, double *ws, double *wy, double *sy,
                                                                                                                                    double *theta, int *col, int *head, int *info)
                                                                                                                  {

                                                                                                                    /* System generated locals */
                                                                                                                    int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset,
                                                                                                                    wy_dim1, wy_offset, sy_dim1, sy_offset;
                                                                                                                    
                                                                                                                    /* Local variables */
                                                                                                                    int dend, pend;
                                                                                                                    int upcl;
                                                                                                                    double temp1, temp2, temp3, temp4;
                                                                                                                    //    int i, k;
                                                                                                                    int ipntr, jpntr, k1, m2, dbegin, is, js, iy, pbegin, is1, js1,
                                                                                                                    col2;
                                                                                                                    
                                                                                                                    /* Parameter adjustments */
                                                                                                                    --indx2;
                                                                                                                    --ind;
                                                                                                                    sy_dim1 = m;
                                                                                                                    sy_offset = 1 + sy_dim1 * 1;
                                                                                                                    sy -= sy_offset;
                                                                                                                    wy_dim1 = n;
                                                                                                                    wy_offset = 1 + wy_dim1 * 1;
                                                                                                                    wy -= wy_offset;
                                                                                                                    ws_dim1 = n;
                                                                                                                    ws_offset = 1 + ws_dim1 * 1;
                                                                                                                    ws -= ws_offset;
                                                                                                                    wn1_dim1 = 2 * m;
                                                                                                                    wn1_offset = 1 + wn1_dim1 * 1;
                                                                                                                    wn1 -= wn1_offset;
                                                                                                                    wn_dim1 = 2 * m;
                                                                                                                    wn_offset = 1 + wn_dim1 * 1;
                                                                                                                    wn -= wn_offset;
                                                                                                                    
                                                                                                                    /* Function Body */
                                                                                                                    
                                                                                                                    /*     Form the lower triangular part of */
                                                                                                                    /*		 WN1 = [Y' ZZ'Y	  L_a'+R_z'] */
                                                                                                                    /*		       [L_a+R_z	  S'AA'S   ] */
                                                                                                                    /*	  where L_a is the strictly lower triangular part of S'AA'Y */
                                                                                                                    /*		R_z is the upper triangular part of S'ZZ'Y. */
                                                                                                                    
                                                                                                                    if (*updatd) {
                                                                                                                      if (*iupdat > m) {/*		shift old part of WN1. */
                                                                                                                    for (int jy = 1; jy <= m - 1; ++jy) {
                                                                                                                      js = m + jy;
                                                                                                                      int i__2 = m - jy;
                                                                                                                      F77_CALL(dcopy)(&i__2, &wn1[jy + 1 + (jy + 1)* wn1_dim1], &c__1,
                                                                                                                               &wn1[jy + jy * wn1_dim1], &c__1);
                                                                                                                      F77_CALL(dcopy)(&i__2, &wn1[js + 1 + (js + 1)* wn1_dim1], &c__1,
                                                                                                                               &wn1[js + js * wn1_dim1], &c__1);
                                                                                                                      i__2 = m - 1;
                                                                                                                      F77_CALL(dcopy)(&i__2, &wn1[m + 2 + (jy + 1) * wn1_dim1], &c__1,
                                                                                                                               &wn1[m + 1 + jy * wn1_dim1], &c__1);
                                                                                                                    }
                                                                                                                      }
                                                                                                                      /*	    put new rows in blocks (1,1), (2,1) and (2,2). */
                                                                                                                      pbegin = 1;
                                                                                                                      pend = *nsub;
                                                                                                                      dbegin = *nsub + 1;
                                                                                                                      dend = n;
                                                                                                                      iy = *col;
                                                                                                                      is = m + *col;
                                                                                                                      ipntr = *head + *col - 1;
                                                                                                                      if (ipntr > m) {
                                                                                                                        ipntr -= m;
                                                                                                                      }
                                                                                                                      jpntr = *head;
                                                                                                                      for (int jy = 1; jy <= *col; ++jy) {
                                                                                                                        js = m + jy;
                                                                                                                        temp1 = 0.;
                                                                                                                        temp2 = 0.;
                                                                                                                        temp3 = 0.;
                                                                                                                        /*	       compute element jy of row 'col' of Y'ZZ'Y */
                                                                                                                        for (int k = pbegin; k <= pend; ++k) {
                                                                                                                          k1 = ind[k];
                                                                                                                          temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                        }
                                                                                                                        /*	       compute elements jy of row 'col' of L_a and S'AA'S */
                                                                                                                        for (int k = dbegin; k <= dend; ++k) {
                                                                                                                          k1 = ind[k];
                                                                                                                          temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
                                                                                                                          temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                        }
                                                                                                                        wn1[iy + jy * wn1_dim1] = temp1;
                                                                                                                        wn1[is + js * wn1_dim1] = temp2;
                                                                                                                        wn1[is + jy * wn1_dim1] = temp3;
                                                                                                                        jpntr = jpntr % m + 1;
                                                                                                                      }
                                                                                                                      /*	    put new column in block (2,1). */
                                                                                                                      int jy = *col;
                                                                                                                      jpntr = *head + *col - 1;
                                                                                                                      if (jpntr > m) {
                                                                                                                        jpntr -= m;
                                                                                                                      }
                                                                                                                      ipntr = *head;
                                                                                                                      for (int i = 1; i <= *col; ++i) {
                                                                                                                        is = m + i;
                                                                                                                        temp3 = 0.;
                                                                                                                        /*	       compute element i of column 'col' of R_z */
                                                                                                                        for (int k = pbegin; k <= pend; ++k) {
                                                                                                                          k1 = ind[k];
                                                                                                                          temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                        }
                                                                                                                        ipntr = ipntr % m + 1;
                                                                                                                        wn1[is + jy * wn1_dim1] = temp3;
                                                                                                                      }
                                                                                                                      upcl = *col - 1;
                                                                                                                    } else {
                                                                                                                      upcl = *col;
                                                                                                                    }
                                                                                                                    /*	 modify the old parts in blocks (1,1) and (2,2) due to changes */
                                                                                                                    /*	 in the set of free variables. */
                                                                                                                    ipntr = *head;
                                                                                                                    for (int iy = 1; iy <= upcl; ++iy) {
                                                                                                                      is = m + iy;
                                                                                                                      jpntr = *head;
                                                                                                                      for (int jy = 1; jy <= iy; ++jy) {
                                                                                                                        js = m + jy;
                                                                                                                        temp1 = 0.;
                                                                                                                        temp2 = 0.;
                                                                                                                        temp3 = 0.;
                                                                                                                        temp4 = 0.;
                                                                                                                        for (int k = 1; k <= *nenter; ++k) {
                                                                                                                          k1 = indx2[k];
                                                                                                                          temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                          temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
                                                                                                                        }
                                                                                                                        for (int k = *ileave; k <= n; ++k) {
                                                                                                                          k1 = indx2[k];
                                                                                                                          temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                          temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
                                                                                                                        }
                                                                                                                        wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
                                                                                                                        wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
                                                                                                                        jpntr = jpntr % m + 1;
                                                                                                                      }
                                                                                                                      ipntr = ipntr % m + 1;
                                                                                                                    }
                                                                                                                    /*	 modify the old parts in block (2,1). */
                                                                                                                    ipntr = *head;
                                                                                                                    for (int is = m + 1; is <= m + upcl; ++is) {
                                                                                                                      jpntr = *head;
                                                                                                                      for (int jy = 1; jy <= upcl; ++jy) {
                                                                                                                        temp1 = 0.;
                                                                                                                        temp3 = 0.;
                                                                                                                        for (int k = 1; k <= *nenter; ++k) {
                                                                                                                          k1 = indx2[k];
                                                                                                                          temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                        }
                                                                                                                        for (int k = *ileave; k <= n; ++k) {
                                                                                                                          k1 = indx2[k];
                                                                                                                          temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
                                                                                                                        }
                                                                                                                        if (is <= jy + m) {
                                                                                                                          wn1[is + jy * wn1_dim1] +=  temp1 - temp3;
                                                                                                                        } else {
                                                                                                                          wn1[is + jy * wn1_dim1] += -temp1 + temp3;
                                                                                                                        }
                                                                                                                        jpntr = jpntr % m + 1;
                                                                                                                      }
                                                                                                                      ipntr = ipntr % m + 1;
                                                                                                                    }
                                                                                                                    /*     Form the upper triangle of WN = [D+Y' ZZ'Y/theta	  -L_a'+R_z' ] */
                                                                                                                    /*				       [-L_a +R_z	 S'AA'S*theta] */
                                                                                                                    m2 = m << 1;
                                                                                                                    for (int iy = 1; iy <= *col; ++iy) {
                                                                                                                      is = *col + iy;
                                                                                                                      is1 = m + iy;
                                                                                                                      for (int jy = 1; jy <= iy;  ++jy) {
                                                                                                                        js = *col + jy;
                                                                                                                        js1 = m + jy;
                                                                                                                        wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
                                                                                                                        wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
                                                                                                                      }
                                                                                                                      for (int jy = 1; jy <= iy - 1; ++jy) {
                                                                                                                        wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
                                                                                                                      }
                                                                                                                      for (int  jy = iy; jy <= *col; ++jy) {
                                                                                                                        wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
                                                                                                                      }
                                                                                                                      wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
                                                                                                                    }
                                                                                                                    /*     Form the upper triangle of */
                                                                                                                    /*	    WN= [  LL'		  L^-1(-L_a'+R_z')] */
                                                                                                                    /*		[(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
                                                                                                                    /*	  first Cholesky factor (1,1) block of wn to get LL' */
                                                                                                                    /*			    with L' stored in the upper triangle of wn. */
                                                                                                                    F77_CALL(dpofa)(&wn[wn_offset], &m2, col, info);
                                                                                                                    if (*info != 0) {
                                                                                                                      *info = -1;
                                                                                                                      return;
                                                                                                                    }
                                                                                                                    /*	  then form L^-1(-L_a'+R_z') in the (1,2) block. */
                                                                                                                    col2 = *col << 1;
                                                                                                                    for (int js = *col + 1; js <= col2; ++js) {
                                                                                                                      F77_CALL(dtrsl)(&wn[wn_offset], &m2, col,
                                                                                                                               &wn[js * wn_dim1 + 1], &c__11, info);
                                                                                                                    }
                                                                                                                    /*     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
                                                                                                                    /*	  upper triangle of (2,2) block of wn. */
                                                                                                                    for (int is = *col + 1; is <= col2; ++is) {
                                                                                                                      for (int js = is; js <= col2; ++js) {
                                                                                                                        wn[is + js * wn_dim1] +=
                                                                                                                          F77_CALL(ddot)(col, &wn[is * wn_dim1 + 1], &c__1,
                                                                                                                                   &wn[js * wn_dim1 + 1], &c__1);
                                                                                                                      }
                                                                                                                    }
                                                                                                                    /*     Cholesky factorization of (2,2) block of wn. */
                                                                                                                    F77_CALL(dpofa)(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
                                                                                                                    if (*info != 0) {
                                                                                                                      *info = -2;
                                                                                                                      return;
                                                                                                                    }
                                                                                                                    return;
                                                                                                                  } /* formk */
                                                                                                                    /* ======================= The end of formk ============================== */
                                                                                                                    
                                                                                                                    static void formt(int m, double *wt, double *sy, double *ss,
                                                                                                                                      int *col, double *theta, int *info)
                                                                                                                    {

                                                                                                                      /* System generated locals */
                                                                                                                      int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset;
                                                                                                                      
                                                                                                                      /* Local variables */
                                                                                                                      double ddum;
                                                                                                                      int k, k1;
                                                                                                                      
                                                                                                                      /* Parameter adjustments */
                                                                                                                      ss_dim1 = m;
                                                                                                                      ss_offset = 1 + ss_dim1 * 1;
                                                                                                                      ss -= ss_offset;
                                                                                                                      sy_dim1 = m;
                                                                                                                      sy_offset = 1 + sy_dim1 * 1;
                                                                                                                      sy -= sy_offset;
                                                                                                                      wt_dim1 = m;
                                                                                                                      wt_offset = 1 + wt_dim1 * 1;
                                                                                                                      wt -= wt_offset;
                                                                                                                      
                                                                                                                      /* Function Body */
                                                                                                                      
                                                                                                                      /*     Form the upper half of  T = theta*SS + L*D^(-1)*L', */
                                                                                                                      /*	  store T in the upper triangle of the array wt. */
                                                                                                                      for (int j = 1; j <= *col; ++j)
                                                                                                                        wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
                                                                                                                      for (int i = 2; i <= *col; ++i) {
                                                                                                                        for (int j = i; j <= *col; ++j) {
                                                                                                                          k1 = min(i,j) - 1;
                                                                                                                          ddum = 0.;
                                                                                                                          for (k = 1; k <= k1; ++k) {
                                                                                                                            ddum += sy[i + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k +
                                                                                                                              k * sy_dim1];
                                                                                                                          }
                                                                                                                          wt[i + j * wt_dim1] = ddum + *theta * ss[i + j * ss_dim1];
                                                                                                                        }
                                                                                                                      }
                                                                                                                      /*     Cholesky factorize T to J*J' with */
                                                                                                                      /*	  J' stored in the upper triangle of wt. */
                                                                                                                      F77_CALL(dpofa)(&wt[wt_offset], &m, col, info);
                                                                                                                      if (*info != 0) {
                                                                                                                        *info = -3;
                                                                                                                      }
                                                                                                                      return;
                                                                                                                    } /* formt */
                                                                                                                      /* ======================= The end of formt ============================== */
                                                                                                                      
                                                                                                                      static void freev(int n, int *nfree, int *indx,
                                                                                                                                        int *nenter, int *ileave, int *indx2, int *iwhere,
                                                                                                                                        int *wrk, int *updatd, int *cnstnd, int iprint,
                                                                                                                                        int *iter)
                                                                                                                      {

                                                                                                                        /* Local variables */
                                                                                                                        int iact, k;
                                                                                                                        
                                                                                                                        /* Parameter adjustments */
                                                                                                                        --iwhere;
                                                                                                                        --indx2;
                                                                                                                        --indx;
                                                                                                                        
                                                                                                                        /* Function Body */
                                                                                                                        *nenter = 0;
                                                                                                                        *ileave = n + 1;
                                                                                                                        if (*iter > 0 && *cnstnd) {/* count the entering and leaving variables. */
                                                                                                                        for (int i = 1; i <= *nfree; ++i) {
                                                                                                                          k = indx[i];
                                                                                                                          if (iwhere[k] > 0) {
                                                                                                                            --(*ileave);
                                                                                                                            indx2[*ileave] = k;
                                                                                                                            if (iprint >= 100)
                                                                                                                              Rprintf("Variable %d leaves the set of free variables\n",
                                                                                                                                      k);
                                                                                                                          }
                                                                                                                        }
                                                                                                                        for (int i = *nfree + 1; i <= n; ++i) {
                                                                                                                          k = indx[i];
                                                                                                                          if (iwhere[k] <= 0) {
                                                                                                                            ++(*nenter);
                                                                                                                            indx2[*nenter] = k;
                                                                                                                            if (iprint >= 100)
                                                                                                                              Rprintf("Variable %d enters the set of free variables\n",
                                                                                                                                      k);
                                                                                                                          }
                                                                                                                          if (iprint >= 100)
                                                                                                                            Rprintf("%d variables leave; %d variables enter\n",
                                                                                                                                    n + 1 - *ileave, *nenter);
                                                                                                                        }
                                                                                                                        }
                                                                                                                        *wrk = *ileave < n + 1 || *nenter > 0 || *updatd;
                                                                                                                        /*     Find the index set of free and active variables at the GCP. */
                                                                                                                        *nfree = 0;
                                                                                                                        iact = n + 1;
                                                                                                                        for (int i = 1; i <= n; ++i) {
                                                                                                                          if (iwhere[i] <= 0) {
                                                                                                                            ++(*nfree);
                                                                                                                            indx[*nfree] = i;
                                                                                                                          } else {
                                                                                                                            --iact;
                                                                                                                            indx[iact] = i;
                                                                                                                          }
                                                                                                                        }
                                                                                                                        if (iprint >= 99)
                                                                                                                          Rprintf("%d  variables are free at GCP on iteration %d\n",
                                                                                                                                  *nfree, *iter + 1);
                                                                                                                        return;
                                                                                                                      } /* freev */
                                                                                                                        /* ======================= The end of freev ============================== */
                                                                                                                        
                                                                                                                        static void hpsolb(int n, double *t, int *iorder, int iheap)
                                                                                                                        {

                                                                                                                          /* Local variables */
                                                                                                                          double ddum;
                                                                                                                          int i, j, indxin, indxou;
                                                                                                                          double out;
                                                                                                                          
                                                                                                                          /* Parameter adjustments */
                                                                                                                          --iorder;
                                                                                                                          --t;
                                                                                                                          
                                                                                                                          /* Function Body */
                                                                                                                          if (iheap == 0) {
                                                                                                                            /*	  Rearrange the elements t(1) to t(n) to form a heap. */
                                                                                                                            for (int k = 2; k <= n; ++k) {
                                                                                                                              ddum = t[k];
                                                                                                                              indxin = iorder[k];
                                                                                                                              /*	     Add ddum to the heap. */
                                                                                                                              i = k;
                                                                                                                              h_loop:
                                                                                                                                if (i > 1) {
                                                                                                                                  j = i / 2;
                                                                                                                                  if (ddum < t[j]) {
                                                                                                                                    t[i] = t[j];
                                                                                                                                    iorder[i] = iorder[j];
                                                                                                                                    i = j;
                                                                                                                                    goto h_loop;
                                                                                                                                  }
                                                                                                                                }
                                                                                                                                t[i] = ddum;
                                                                                                                                iorder[i] = indxin;
                                                                                                                            }
                                                                                                                          }
                                                                                                                          /*     Assign to 'out' the value of t(1), the least member of the heap, */
                                                                                                                          /*	  and rearrange the remaining members to form a heap as */
                                                                                                                          /*	  elements 1 to n-1 of t. */
                                                                                                                          if (n > 1) {
                                                                                                                            i = 1;
                                                                                                                            out = t[1];
                                                                                                                            indxou = iorder[1];
                                                                                                                            ddum = t[n];
                                                                                                                            indxin = iorder[n];
                                                                                                                            /*	  Restore the heap */
                                                                                                                            Loop:
                                                                                                                              j = i + i;
                                                                                                                            if (j <= n - 1) {
                                                                                                                              if (t[j + 1] < t[j]) {
                                                                                                                                ++j;
                                                                                                                              }
                                                                                                                              if (t[j] < ddum) {
                                                                                                                                t[i] = t[j];
                                                                                                                                iorder[i] = iorder[j];
                                                                                                                                i = j;
                                                                                                                                goto Loop;
                                                                                                                              }
                                                                                                                            }
                                                                                                                            t[i] = ddum;
                                                                                                                            iorder[i] = indxin;
                                                                                                                            /*     Put the least member in t(n). */
                                                                                                                            t[n] = out;
                                                                                                                            iorder[n] = indxou;
                                                                                                                          }
                                                                                                                          return;
                                                                                                                        } /* hpsolb */
                                                                                                                            /* ====================== The end of hpsolb ============================== */
                                                                                                                            
                                                                                                                            static void lnsrlb(int n, double *l, double *u,
                                                                                                                                               int *nbd, double *x, double *f, double *fold,
                                                                                                                                               double *gd, double *gdold, double *g, double *d,
                                                                                                                                               double *r, double *t, double *z, double *stp,
                                                                                                                                               double *dnorm, double *dtd, double *xstep,
                                                                                                                                               double *stpmx, int *iter, int *ifun, int *iback, int *nfgv,
                                                                                                                                               int *info, char *task, int *boxed, int *cnstnd,
                                                                                                                                               char *csave)
                                                                                                                            {

                                                                                                                              /* For dcsrch(): */
                                                                                                                              const double stpmin = 0.;
                                                                                                                              const double ftol = .001;
                                                                                                                              const double gtol = .9;
                                                                                                                              const double xtol = .1;
                                                                                                                              
                                                                                                                              /* System generated locals */
                                                                                                                              double d1;
                                                                                                                              
                                                                                                                              /* Local variables */
                                                                                                                              double a1, a2;
                                                                                                                              
                                                                                                                              /* Parameter adjustments */
                                                                                                                              --z;
                                                                                                                              --t;
                                                                                                                              --r;
                                                                                                                              --d;
                                                                                                                              --g;
                                                                                                                              --x;
                                                                                                                              --nbd;
                                                                                                                              --u;
                                                                                                                              --l;
                                                                                                                              
                                                                                                                              /* Function Body */
                                                                                                                              if (strncmp(task, "FG_LN", 5) == 0) {
                                                                                                                                goto L556;
                                                                                                                              }
                                                                                                                              *dtd = F77_CALL(ddot)(&n, &d[1], &c__1, &d[1], &c__1);
                                                                                                                              *dnorm = sqrt(*dtd);
                                                                                                                              /*     Determine the maximum step length. */
                                                                                                                              *stpmx = 1e10;
                                                                                                                              if (*cnstnd) {
                                                                                                                                if (*iter == 0) {
                                                                                                                                  *stpmx = 1.;
                                                                                                                                } else {
                                                                                                                                  for (int i = 1; i <= n; ++i) {
                                                                                                                                    a1 = d[i];
                                                                                                                                    if (nbd[i] != 0) {
                                                                                                                                      if (a1 < 0. && nbd[i] <= 2) {
                                                                                                                                        a2 = l[i] - x[i];
                                                                                                                                        if (a2 >= 0.) {
                                                                                                                                          *stpmx = 0.;
                                                                                                                                        } else if (a1 * *stpmx < a2) {
                                                                                                                                          *stpmx = a2 / a1;
                                                                                                                                        }
                                                                                                                                      } else if (a1 > 0. && nbd[i] >= 2) {
                                                                                                                                        a2 = u[i] - x[i];
                                                                                                                                        if (a2 <= 0.) {
                                                                                                                                          *stpmx = 0.;
                                                                                                                                        } else if (a1 * *stpmx > a2) {
                                                                                                                                          *stpmx = a2 / a1;
                                                                                                                                        }
                                                                                                                                      }
                                                                                                                                    }
                                                                                                                                  }
                                                                                                                                }
                                                                                                                              }
                                                                                                                              if (*iter == 0 && ! (*boxed)) {
                                                                                                                                d1 = 1. / *dnorm;
                                                                                                                                *stp = min(d1,*stpmx);
                                                                                                                              } else {
                                                                                                                                *stp = 1.;
                                                                                                                              }
                                                                                                                              F77_CALL(dcopy)(&n, &x[1], &c__1, &t[1], &c__1);
                                                                                                                              F77_CALL(dcopy)(&n, &g[1], &c__1, &r[1], &c__1);
                                                                                                                              *fold = *f;
                                                                                                                              *ifun = 0;
                                                                                                                              *iback = 0;
                                                                                                                              strcpy(csave, "START");
                                                                                                                              L556:
                                                                                                                                *gd = F77_CALL(ddot)(&n, &g[1], &c__1, &d[1], &c__1);
                                                                                                                              if (*ifun == 0) {
                                                                                                                                *gdold = *gd;
                                                                                                                                if (*gd >= 0.) {
                                                                                                                                  /*				 the directional derivative >=0. */
                                                                                                                                  /*				 Line search is impossible. */
                                                                                                                                  *info = -4;
                                                                                                                                  return;
                                                                                                                                }
                                                                                                                              }
                                                                                                                              dcsrch(f, gd, stp, ftol, gtol, xtol, stpmin, *stpmx, csave);
                                                                                                                              *xstep = *stp * *dnorm;
                                                                                                                              if (strncmp(csave, "CONV", 4) != 0 && strncmp(csave, "WARN", 4) != 0) {
                                                                                                                                strcpy(task, "FG_LNSRCH");
                                                                                                                                ++(*ifun);
                                                                                                                                ++(*nfgv);
                                                                                                                                *iback = *ifun - 1;
                                                                                                                                if (*stp == 1.) {
                                                                                                                                  F77_CALL(dcopy)(&n, &z[1], &c__1, &x[1], &c__1);
                                                                                                                                } else {
                                                                                                                                  for (int i = 1; i <= n; ++i) {
                                                                                                                                    x[i] = *stp * d[i] + t[i];
                                                                                                                                  }
                                                                                                                                }
                                                                                                                              } else {
                                                                                                                                strcpy(task, "NEW_X");
                                                                                                                              }
                                                                                                                              return;
                                                                                                                            } /* lnsrlb */
                                                                                                                                  /* ======================= The end of lnsrlb ============================= */
                                                                                                                                  
                                                                                                                                  static void matupd(int n, int m, double *ws,
                                                                                                                                                     double *wy, double *sy, double *ss, double *d,
                                                                                                                                                     double *r, int *itail, int *iupdat, int *col,
                                                                                                                                                     int *head, double *theta, double *rr, double *dr,
                                                                                                                                                     double *stp, double *dtd)
                                                                                                                                  {

                                                                                                                                    /* System generated locals */
                                                                                                                                    int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset,
                                                                                                                                    ss_dim1, ss_offset;
                                                                                                                                    
                                                                                                                                    /* Local variables */
                                                                                                                                    int pointr;
                                                                                                                                    
                                                                                                                                    /* Parameter adjustments */
                                                                                                                                    --r;
                                                                                                                                    --d;
                                                                                                                                    ss_dim1 = m;
                                                                                                                                    ss_offset = 1 + ss_dim1 * 1;
                                                                                                                                    ss -= ss_offset;
                                                                                                                                    sy_dim1 = m;
                                                                                                                                    sy_offset = 1 + sy_dim1 * 1;
                                                                                                                                    sy -= sy_offset;
                                                                                                                                    wy_dim1 = n;
                                                                                                                                    wy_offset = 1 + wy_dim1 * 1;
                                                                                                                                    wy -= wy_offset;
                                                                                                                                    ws_dim1 = n;
                                                                                                                                    ws_offset = 1 + ws_dim1 * 1;
                                                                                                                                    ws -= ws_offset;
                                                                                                                                    
                                                                                                                                    /* Function Body */
                                                                                                                                    
                                                                                                                                    /*     Set pointers for matrices WS and WY. */
                                                                                                                                    if (*iupdat <= m) {
                                                                                                                                      *col = *iupdat;
                                                                                                                                      *itail = (*head + *iupdat - 2) % m + 1;
                                                                                                                                    } else {
                                                                                                                                      *itail = *itail % m + 1;
                                                                                                                                      *head = *head % m + 1;
                                                                                                                                    }
                                                                                                                                    /*     Update matrices WS and WY. */
                                                                                                                                    F77_CALL(dcopy)(&n, &d[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
                                                                                                                                    F77_CALL(dcopy)(&n, &r[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
                                                                                                                                    /*     Set theta=yy/ys. */
                                                                                                                                    *theta = *rr / *dr;
                                                                                                                                    /*     Form the middle matrix in B. */
                                                                                                                                    /*	  update the upper triangle of SS, */
                                                                                                                                    /*					   and the lower triangle of SY: */
                                                                                                                                    if (*iupdat > m) {
                                                                                                                                      /*				move old information */
                                                                                                                                      for (int j = 1; j <= *col - 1; ++j) {
                                                                                                                                        F77_CALL(dcopy)(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1,
                                                                                                                                                 &ss[j * ss_dim1 + 1], &c__1);
                                                                                                                                        int i__2 = *col - j;
                                                                                                                                        F77_CALL(dcopy)(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1,
                                                                                                                                                 &sy[j + j * sy_dim1], &c__1);
                                                                                                                                      }
                                                                                                                                    }
                                                                                                                                    /*	  add new information: the last row of SY */
                                                                                                                                    /*					       and the last column of SS: */
                                                                                                                                    pointr = *head;
                                                                                                                                    for (int j = 1; j <= *col - 1; ++j) {
                                                                                                                                      sy[*col + j * sy_dim1] =
                                                                                                                                        F77_CALL(ddot)(&n, &d[1], &c__1, &wy[pointr * wy_dim1 + 1], &c__1);
                                                                                                                                      ss[j + *col * ss_dim1] =
                                                                                                                                        F77_CALL(ddot)(&n, &ws[pointr * ws_dim1 + 1], &c__1, &d[1], &c__1);
                                                                                                                                      pointr = pointr % m + 1;
                                                                                                                                    }
                                                                                                                                    if (*stp == 1.) {
                                                                                                                                      ss[*col + *col * ss_dim1] = *dtd;
                                                                                                                                    } else {
                                                                                                                                      ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
                                                                                                                                    }
                                                                                                                                    sy[*col + *col * sy_dim1] = *dr;
                                                                                                                                    return;
                                                                                                                                  } /* matupd */
                                                                                                                                    /* ======================= The end of matupd ============================= */
                                                                                                                                    
                                                                                                                                    static void projgr(int n, double *l, double *u,
                                                                                                                                                       int *nbd, double *x, double *g, double *sbgnrm)
                                                                                                                                    {
                                                                                                                                      double gi, d__1;
                                                                                                                                      
                                                                                                                                      *sbgnrm = 0.;
                                                                                                                                      for (int i = 0; i < n; ++i) {
                                                                                                                                        gi = g[i];
                                                                                                                                        if (nbd[i] != 0) {
                                                                                                                                          if (gi < 0.) {
                                                                                                                                            if (nbd[i] >= 2) {
                                                                                                                                              if(gi < (d__1 = x[i] - u[i]))
                                                                                                                                                gi = d__1;
                                                                                                                                            }
                                                                                                                                          } else {
                                                                                                                                            if (nbd[i] <= 2) {
                                                                                                                                              if(gi > (d__1 = x[i] - l[i]))
                                                                                                                                                gi = d__1;
                                                                                                                                            }
                                                                                                                                          }
                                                                                                                                        }
                                                                                                                                        if(*sbgnrm < (d__1 = fabs(gi))) *sbgnrm = d__1;
                                                                                                                                      }
                                                                                                                                      return;
                                                                                                                                    } /* projgr */
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      static void subsm(int n, int m, int *nsub, int *ind,
                                                                                                                                                        double *l, double *u, int *nbd, double *x,
                                                                                                                                                        double *d, double *ws, double *wy, double *theta,
                                                                                                                                                        int *col, int *head, int *iword, double *wv,
                                                                                                                                                        double *wn, int iprint, int *info)
                                                                                                                                      {

                                                                                                                                        /* System generated locals */
                                                                                                                                        int ws_offset, wn_dim1, wn_offset;
                                                                                                                                        
                                                                                                                                        /* Local variables */
                                                                                                                                        double alpha, dk, temp1, temp2;
                                                                                                                                        int k, m2, js, pointr, ibd = 0, col2, ns;
                                                                                                                                        
                                                                                                                                        /* Parameter adjustments */
                                                                                                                                        --d;
                                                                                                                                        --u;
                                                                                                                                        --l;
                                                                                                                                        --x;
                                                                                                                                        --ind;
                                                                                                                                        --nbd;
                                                                                                                                        --wv;
                                                                                                                                        wn_dim1 = 2 * m;
                                                                                                                                        wn_offset = 1 + wn_dim1 * 1;
                                                                                                                                        wn -= wn_offset;
                                                                                                                                        /* ws[] and wy[] are both  [n x m ] :*/
                                                                                                                                        ws_offset = 1 + n * 1;
                                                                                                                                        ws -= ws_offset;
                                                                                                                                        wy -= ws_offset;
                                                                                                                                        
                                                                                                                                        ns = *nsub;
                                                                                                                                        if (ns <= 0)
                                                                                                                                          return;
                                                                                                                                        
                                                                                                                                        /*     Compute wv = W'Zd. */
                                                                                                                                        pointr = *head;
                                                                                                                                        for (int i = 1; i <= *col; ++i) {
                                                                                                                                          temp1 = 0.;
                                                                                                                                          temp2 = 0.;
                                                                                                                                          for (int j = 1; j <= ns; ++j) {
                                                                                                                                            k = ind[j];
                                                                                                                                            temp1 += wy[k + pointr * n] * d[j];
                                                                                                                                            temp2 += ws[k + pointr * n] * d[j];
                                                                                                                                          }
                                                                                                                                          wv[i] = temp1;
                                                                                                                                          wv[*col + i] = *theta * temp2;
                                                                                                                                          pointr = pointr % m + 1;
                                                                                                                                        }
                                                                                                                                        /*     Compute wv:=K^(-1)wv. */
                                                                                                                                        m2 = m << 1;
                                                                                                                                        col2 = *col << 1;
                                                                                                                                        F77_CALL(dtrsl)(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
                                                                                                                                        if (*info != 0) {
                                                                                                                                          return;
                                                                                                                                        }
                                                                                                                                        for (int i = 1; i <= *col; ++i)
                                                                                                                                          wv[i] = -wv[i];
                                                                                                                                        
                                                                                                                                        F77_CALL(dtrsl)(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
                                                                                                                                        if (*info != 0) {
                                                                                                                                          return;
                                                                                                                                        }
                                                                                                                                        /*     Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
                                                                                                                                        pointr = *head;
                                                                                                                                        for (int jy = 1; jy <= *col; ++jy) {
                                                                                                                                          js = *col + jy;
                                                                                                                                          for (int i = 1; i <= ns; ++i) {
                                                                                                                                            k = ind[i];
                                                                                                                                            d[i] += (wy[k + pointr * n] * wv[jy] / *theta +
                                                                                                                                              ws[k + pointr * n] * wv[js]);
                                                                                                                                          }
                                                                                                                                          pointr = pointr % m + 1;
                                                                                                                                        }
                                                                                                                                        
                                                                                                                                        for (int i = 1; i <= ns; ++i)
                                                                                                                                          d[i] /= *theta;
                                                                                                                                        
                                                                                                                                        /*     Backtrack to the feasible region. */
                                                                                                                                        alpha = 1.;
                                                                                                                                        temp1 = alpha;
                                                                                                                                        for (int i = 1; i <= ns; ++i) {
                                                                                                                                          k = ind[i];
                                                                                                                                          dk = d[i];
                                                                                                                                          if (nbd[k] != 0) {
                                                                                                                                            if (dk < 0. && nbd[k] <= 2) {
                                                                                                                                              temp2 = l[k] - x[k];
                                                                                                                                              if (temp2 >= 0.) {
                                                                                                                                                temp1 = 0.;
                                                                                                                                              } else if (dk * alpha < temp2) {
                                                                                                                                                temp1 = temp2 / dk;
                                                                                                                                              }
                                                                                                                                            } else if (dk > 0. && nbd[k] >= 2) {
                                                                                                                                              temp2 = u[k] - x[k];
                                                                                                                                              if (temp2 <= 0.) {
                                                                                                                                                temp1 = 0.;
                                                                                                                                              } else if (dk * alpha > temp2) {
                                                                                                                                                temp1 = temp2 / dk;
                                                                                                                                              }
                                                                                                                                            }
                                                                                                                                            if (temp1 < alpha) {
                                                                                                                                              alpha = temp1;
                                                                                                                                              ibd = i;
                                                                                                                                            }
                                                                                                                                          }
                                                                                                                                        }
                                                                                                                                        if (alpha < 1.) {
                                                                                                                                          dk = d[ibd];
                                                                                                                                          k = ind[ibd];
                                                                                                                                          if (dk > 0.) {
                                                                                                                                            x[k] = u[k];
                                                                                                                                            d[ibd] = 0.;
                                                                                                                                          } else if (dk < 0.) {
                                                                                                                                            x[k] = l[k];
                                                                                                                                            d[ibd] = 0.;
                                                                                                                                          }
                                                                                                                                        }
                                                                                                                                        for (int i = 1; i <= ns; ++i)
                                                                                                                                          x[ind[i]] += alpha * d[i];
                                                                                                                                        
                                                                                                                                        *iword = (alpha < 1.) ? 1 : 0;
                                                                                                                                        
                                                                                                                                        return;
                                                                                                                                      } /* subsm */
                                                                                                                                        /* ====================== The end of subsm =============================== */
                                                                                                                                        
                                                                                                                                        static void dcsrch(double *f, double *g, double *stp,
                                                                                                                                                           /*Chgd: the next five are no longer pointers:*/
                                                                                                                                                           double ftol, double gtol, double xtol,
                                                                                                                                                           double stpmin, double stpmax, char *task)
                                                                                                                                        {

                                                                                                                                          /* Local variables */
                                                                                                                                          static int stage, brackt;
                                                                                                                                          static double ginit, gtest, gx, gy, finit, fx, fy, stx, sty,
                                                                                                                                          stmin, stmax, width, width1;
                                                                                                                                          double ftest, fm, gm, fxm, fym, gxm, gym;
                                                                                                                                          
                                                                                                                                          /* Function Body */
                                                                                                                                          
                                                                                                                                          /*     Initialization block. */
                                                                                                                                          if (strncmp(task, "START", 5) == 0) {
                                                                                                                                            /*	  Check the input arguments for errors. */
                                                                                                                                            if (*stp < stpmin)	strcpy(task, "ERROR: STP .LT. STPMIN");
                                                                                                                                            if (*stp > stpmax)	strcpy(task, "ERROR: STP .GT. STPMAX");
                                                                                                                                            if (*g >= 0.)		strcpy(task, "ERROR: INITIAL G .GE. ZERO");
                                                                                                                                            if (ftol < 0.)		strcpy(task, "ERROR: FTOL .LT. ZERO");
                                                                                                                                            if (gtol < 0.)		strcpy(task, "ERROR: GTOL .LT. ZERO");
                                                                                                                                            if (xtol < 0.)		strcpy(task, "ERROR: XTOL .LT. ZERO");
                                                                                                                                            if (stpmin < 0.)	strcpy(task, "ERROR: STPMIN .LT. ZERO");
                                                                                                                                            if (stpmax < stpmin)	strcpy(task, "ERROR: STPMAX .LT. STPMIN");
                                                                                                                                            
                                                                                                                                            /*	  Exit if there are errors on input. */
                                                                                                                                            if (strncmp(task, "ERROR", 5) == 0) {
                                                                                                                                              return;
                                                                                                                                            }
                                                                                                                                            /*	  Initialize local variables. */
                                                                                                                                            brackt = FALSE_;
                                                                                                                                            stage = 1;
                                                                                                                                            finit = *f;
                                                                                                                                            ginit = *g;
                                                                                                                                            gtest = ftol * ginit;
                                                                                                                                            width = stpmax - stpmin;
                                                                                                                                            width1 = width / .5;
                                                                                                                                            /*	  The variables stx, fx, gx contain the values of the step, */
                                                                                                                                            /*	  function, and derivative at the best step. */
                                                                                                                                            /*	  The variables sty, fy, gy contain the value of the step, */
                                                                                                                                            /*	  function, and derivative at sty. */
                                                                                                                                            /*	  The variables stp, f, g contain the values of the step, */
                                                                                                                                            /*	  function, and derivative at stp. */
                                                                                                                                            stx = 0.;	fx = finit;	gx = ginit;
                                                                                                                                            sty = 0.;	fy = finit;	gy = ginit;
                                                                                                                                            stmin = 0.;
                                                                                                                                            stmax = *stp + *stp * 4.;
                                                                                                                                            strcpy(task, "FG");
                                                                                                                                            return;
                                                                                                                                          }
                                                                                                                                          /*      If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
                                                                                                                                          /*      algorithm enters the second stage. */
                                                                                                                                          ftest = finit + *stp * gtest;
                                                                                                                                          if (stage == 1 && *f <= ftest && *g >= 0.)
                                                                                                                                            stage = 2;
                                                                                                                                          /*	Test for warnings. */
                                                                                                                                          if (brackt && (*stp <= stmin || *stp >= stmax))
                                                                                                                                            strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
                                                                                                                                          if (brackt && stmax - stmin <= xtol * stmax)
                                                                                                                                            strcpy(task, "WARNING: XTOL TEST SATISFIED");
                                                                                                                                          if (*stp == stpmax && *f <= ftest && *g <= gtest)
                                                                                                                                            strcpy(task, "WARNING: STP = STPMAX");
                                                                                                                                          if (*stp == stpmin && (*f > ftest || *g >= gtest))
                                                                                                                                            strcpy(task, "WARNING: STP = STPMIN");
                                                                                                                                          /*	Test for convergence. */
                                                                                                                                          if (*f <= ftest && fabs(*g) <= gtol * (-ginit))
                                                                                                                                            strcpy(task, "CONVERGENCE");
                                                                                                                                          /*	Test for termination. */
                                                                                                                                          if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4) == 0)
                                                                                                                                            return;
                                                                                                                                          
                                                                                                                                          /*     A modified function is used to predict the step during the */
                                                                                                                                          /*     first stage if a lower function value has been obtained but */
                                                                                                                                          /*     the decrease is not sufficient. */
                                                                                                                                          if (stage == 1 && *f <= fx && *f > ftest) {
                                                                                                                                            /*	  Define the modified function and derivative values. */
                                                                                                                                            fm = *f - *stp * gtest;
                                                                                                                                            fxm = fx - stx * gtest;
                                                                                                                                            fym = fy - sty * gtest;
                                                                                                                                            gm = *g - gtest;
                                                                                                                                            gxm = gx - gtest;
                                                                                                                                            gym = gy - gtest;
                                                                                                                                            /*	  Call dcstep to update stx, sty, and to compute the new step. */
                                                                                                                                            dcstep(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
                                                                                                                                            stmin, &stmax);
                                                                                                                                            /*	  Reset the function and derivative values for f. */
                                                                                                                                            fx = fxm + stx * gtest;
                                                                                                                                            fy = fym + sty * gtest;
                                                                                                                                            gx = gxm + gtest;
                                                                                                                                            gy = gym + gtest;
                                                                                                                                          } else {
                                                                                                                                            /*	 Call dcstep to update stx, sty, and to compute the new step. */
                                                                                                                                            dcstep(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
                                                                                                                                            stmax);
                                                                                                                                          }
                                                                                                                                          /*     Decide if a bisection step is needed. */
                                                                                                                                          if (brackt) {
                                                                                                                                            if (fabs(sty - stx) >= width1 * .66) {
                                                                                                                                              *stp = stx + (sty - stx) * .5;
                                                                                                                                            }
                                                                                                                                            width1 = width;
                                                                                                                                            width = fabs(sty - stx);
                                                                                                                                          }
                                                                                                                                          /*     Set the minimum and maximum steps allowed for stp. */
                                                                                                                                          if (brackt) {
                                                                                                                                            stmin = min(stx,sty);
                                                                                                                                            stmax = max(stx,sty);
                                                                                                                                          } else {
                                                                                                                                            stmin = *stp + (*stp - stx) * 1.1;
                                                                                                                                            stmax = *stp + (*stp - stx) * 4.;
                                                                                                                                          }
                                                                                                                                          /*     Force the step to be within the bounds stpmax and stpmin. */
                                                                                                                                          if(*stp < stpmin) *stp = stpmin;
                                                                                                                                          if(*stp > stpmax) *stp = stpmax;
                                                                                                                                          
                                                                                                                                          /*     If further progress is not possible, let stp be the best */
                                                                                                                                          /*     point obtained during the search. */
                                                                                                                                          if ((brackt && (*stp <= stmin || *stp >= stmax)) ||
                                                                                                                                          (brackt && (stmax - stmin <= xtol * stmax))) {
                                                                                                                                            *stp = stx;
                                                                                                                                          }
                                                                                                                                          /*     Obtain another function and derivative. */
                                                                                                                                          strcpy(task, "FG");
                                                                                                                                          return;
                                                                                                                                        } /* dcsrch */
                                                                                                                                          /* ====================== The end of dcsrch ============================== */
                                                                                                                                          
                                                                                                                                          static void dcstep(double *stx, double *fx, double *dx,
                                                                                                                                                             double *sty, double *fy, double *dy, double *stp,
                                                                                                                                                             double *fp, double *dp, int *brackt, double *stpmin,
                                                                                                                                                             double *stpmax)
                                                                                                                                          {

                                                                                                                                            /* System generated locals */
                                                                                                                                            double d__1, d__2;
                                                                                                                                            
                                                                                                                                            /* Local variables */
                                                                                                                                            double sgnd, stpc, stpf, stpq, p, q, gamm, r__, s, theta;
                                                                                                                                            
                                                                                                                                            sgnd = *dp * (*dx / fabs(*dx));
                                                                                                                                            /*     First case: A higher function value. The minimum is bracketed. */
                                                                                                                                            /*     If the cubic step is closer to stx than the quadratic step, the */
                                                                                                                                            /*     cubic step is taken, otherwise the average of the cubic and */
                                                                                                                                            /*     quadratic steps is taken. */
                                                                                                                                            if (*fp > *fx) {
                                                                                                                                              theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
                                                                                                                                              /* Computing MAX */
                                                                                                                                              d__1 = fabs(theta), d__2 = fabs(*dx),
                                                                                                                                                d__1 = max(d__1,d__2), d__2 = fabs(*dp);
                                                                                                                                              s = max(d__1,d__2);
                                                                                                                                              /* Computing 2nd power */
                                                                                                                                              d__1 = theta / s;
                                                                                                                                              gamm = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
                                                                                                                                              if (*stp < *stx) {
                                                                                                                                                gamm = -gamm;
                                                                                                                                              }
                                                                                                                                              p = gamm - *dx + theta;
                                                                                                                                              q = gamm - *dx + gamm + *dp;
                                                                                                                                              r__ = p / q;
                                                                                                                                              stpc = *stx + r__ * (*stp - *stx);
                                                                                                                                              stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp
                                                                                                                                                                                                                - *stx);
                                                                                                                                              if (fabs(stpc - *stx) < fabs(stpq - *stx)) {
                                                                                                                                                stpf = stpc;
                                                                                                                                              } else {
                                                                                                                                                stpf = stpc + (stpq - stpc) / 2.;
                                                                                                                                              }
                                                                                                                                              *brackt = TRUE_;
                                                                                                                                              /*     Second case: A lower function value and derivatives of opposite */
                                                                                                                                              /*     sign. The minimum is bracketed. If the cubic step is farther from */
                                                                                                                                              /*     stp than the secant step, the cubic step is taken, otherwise the */
                                                                                                                                              /*     secant step is taken. */
                                                                                                                                            } else if (sgnd < 0.) {
                                                                                                                                              theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
                                                                                                                                              /* Computing MAX */
                                                                                                                                              d__1 = fabs(theta), d__2 = fabs(*dx),
                                                                                                                                                d__1 = max(d__1,d__2), d__2 = fabs(*dp);
                                                                                                                                              s = max(d__1,d__2);
                                                                                                                                              /* Computing 2nd power */
                                                                                                                                              d__1 = theta / s;
                                                                                                                                              gamm = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
                                                                                                                                              if (*stp > *stx) {
                                                                                                                                                gamm = -gamm;
                                                                                                                                              }
                                                                                                                                              p = gamm - *dp + theta;
                                                                                                                                              q = gamm - *dp + gamm + *dx;
                                                                                                                                              r__ = p / q;
                                                                                                                                              stpc = *stp + r__ * (*stx - *stp);
                                                                                                                                              stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
                                                                                                                                              if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
                                                                                                                                                stpf = stpc;
                                                                                                                                              } else {
                                                                                                                                                stpf = stpq;
                                                                                                                                              }
                                                                                                                                              *brackt = TRUE_;
                                                                                                                                              /*     Third case: A lower function value, derivatives of the same sign, */
                                                                                                                                              /*     and the magnitude of the derivative decreases. */
                                                                                                                                            } else if (fabs(*dp) < fabs(*dx)) {
                                                                                                                                              /*	  The cubic step is computed only if the cubic tends to infinity */
                                                                                                                                              /*	  in the direction of the step or if the minimum of the cubic */
                                                                                                                                              /*	  is beyond stp. Otherwise the cubic step is defined to be the */
                                                                                                                                              /*	  secant step. */
                                                                                                                                              theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
                                                                                                                                              /* Computing MAX */
                                                                                                                                              d__1 = fabs(theta), d__2 = fabs(*dx),
                                                                                                                                                d__1 = max(d__1,d__2), d__2 = fabs(*dp);
                                                                                                                                              s = max(d__1,d__2);
                                                                                                                                              /*	  The case gamm = 0 only arises if the cubic does not tend */
                                                                                                                                              /*	  to infinity in the direction of the step. */
                                                                                                                                              /* Computing MAX */
                                                                                                                                              /* Computing 2nd power */
                                                                                                                                              d__1 = theta / s;
                                                                                                                                              d__1 = d__1 * d__1 - *dx / s * (*dp / s);
                                                                                                                                              gamm = d__1 < 0 ? 0. : s * sqrt(d__1);
                                                                                                                                              if (*stp > *stx) {
                                                                                                                                                gamm = -gamm;
                                                                                                                                              }
                                                                                                                                              p = gamm - *dp + theta;
                                                                                                                                              q = gamm + (*dx - *dp) + gamm;
                                                                                                                                              r__ = p / q;
                                                                                                                                              if (r__ < 0. && gamm != 0.) {
                                                                                                                                                stpc = *stp + r__ * (*stx - *stp);
                                                                                                                                              } else if (*stp > *stx) {
                                                                                                                                                stpc = *stpmax;
                                                                                                                                              } else {
                                                                                                                                                stpc = *stpmin;
                                                                                                                                              }
                                                                                                                                              stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
                                                                                                                                              if (*brackt) {
                                                                                                                                                /*	     A minimizer has been bracketed. If the cubic step is */
                                                                                                                                                /*	     closer to stp than the secant step, the cubic step is */
                                                                                                                                                /*	     taken, otherwise the secant step is taken. */
                                                                                                                                                if (fabs(stpc - *stp) < fabs(stpq - *stp)) {
                                                                                                                                                  stpf = stpc;
                                                                                                                                                } else {
                                                                                                                                                  stpf = stpq;
                                                                                                                                                }
                                                                                                                                                d__1 = *stp + (*sty - *stp) * .66;
                                                                                                                                                if (*stp > *stx) {
                                                                                                                                                  stpf = min(d__1,stpf);
                                                                                                                                                } else {
                                                                                                                                                  stpf = max(d__1,stpf);
                                                                                                                                                }
                                                                                                                                              } else {
                                                                                                                                                /*	     A minimizer has not been bracketed. If the cubic step is */
                                                                                                                                                /*	     farther from stp than the secant step, the cubic step is */
                                                                                                                                                /*	     taken, otherwise the secant step is taken. */
                                                                                                                                                if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
                                                                                                                                                  stpf = stpc;
                                                                                                                                                } else {
                                                                                                                                                  stpf = stpq;
                                                                                                                                                }
                                                                                                                                                stpf = min(*stpmax,stpf);
                                                                                                                                                stpf = max(*stpmin,stpf);
                                                                                                                                              }
                                                                                                                                              /*     Fourth case: A lower function value, derivatives of the */
                                                                                                                                              /*     same sign, and the magnitude of the derivative does not */
                                                                                                                                              /*     decrease. If the minimum is not bracketed, the step is either */
                                                                                                                                              /*     stpmin or stpmax, otherwise the cubic step is taken. */
                                                                                                                                            } else {
                                                                                                                                              if (*brackt) {
                                                                                                                                                theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
                                                                                                                                                /* Computing MAX */
                                                                                                                                                d__1 = fabs(theta), d__2 = fabs(*dy), d__1 = max(d__1,d__2), d__2 =
                                                                                                                                                fabs(*dp);
                                                                                                                                                s = max(d__1,d__2);
                                                                                                                                                /* Computing 2nd power */
                                                                                                                                                d__1 = theta / s;
                                                                                                                                                gamm = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
                                                                                                                                                if (*stp > *sty) {
                                                                                                                                                  gamm = -gamm;
                                                                                                                                                }
                                                                                                                                                p = gamm - *dp + theta;
                                                                                                                                                q = gamm - *dp + gamm + *dy;
                                                                                                                                                r__ = p / q;
                                                                                                                                                stpc = *stp + r__ * (*sty - *stp);
                                                                                                                                                stpf = stpc;
                                                                                                                                              } else if (*stp > *stx) {
                                                                                                                                                stpf = *stpmax;
                                                                                                                                              } else {
                                                                                                                                                stpf = *stpmin;
                                                                                                                                              }
                                                                                                                                            }
                                                                                                                                            /*     Update the interval which contains a minimizer. */
                                                                                                                                            if (*fp > *fx) {
                                                                                                                                              *sty = *stp;
                                                                                                                                              *fy = *fp;
                                                                                                                                              *dy = *dp;
                                                                                                                                            } else {
                                                                                                                                              if (sgnd < 0.) {
                                                                                                                                                *sty = *stx;
                                                                                                                                                *fy = *fx;
                                                                                                                                                *dy = *dx;
                                                                                                                                              }
                                                                                                                                              *stx = *stp;
                                                                                                                                              *fx = *fp;
                                                                                                                                              *dx = *dp;
                                                                                                                                            }
                                                                                                                                            /*     Compute the new step. */
                                                                                                                                            *stp = stpf;
                                                                                                                                            return;
                                                                                                                                          } /* dcstep */
                                                                                                                                            /* ====================== The end of dcstep ============================== */
                                                                                                                                            
                                                                                                                                            static void pvector(const char *title, double *x, int n)
                                                                                                                                            {
                                                                                                                                              Rprintf("%s ", title);
                                                                                                                                              for (int i = 0; i < n; i++) Rprintf("%g ", x[i]);
                                                                                                                                              Rprintf("\n");
                                                                                                                                            }


static void prn1lb(int n, int m, double *l, double *u, double *x,
                   int iprint, double epsmch)
{
  if (iprint >=  0) {
    Rprintf("N = %d, M = %d machine precision = %g\n", n, m, epsmch);
    if (iprint >= 100){
      pvector("L =", l, n);
      pvector("X0 =",x, n);
      pvector("U =", u, n);
    }
  }
}


static void prn2lb(int n, double *x, double *f, double *g, int iprint,
                   int iter, int nfgv, int nact, double sbgnrm,
                   int nint, char *word, int iword, int iback,
                   double stp, double xstep)
{
  if (iprint >=  99) {
    Rprintf("LINE SEARCH %d times; norm of step = %g\n", iback, xstep);
    if (iprint > 100) {
      pvector("X =", x, n);
      pvector("G =", g, n);
    }
  } else if (iprint > 0 && iter%iprint == 0) {
    Rprintf("At iterate %5d  f = %12.5g  |proj g|=  %12.5g\n",
            iter, *f, sbgnrm);
  }
}

static void prn3lb(int n, double *x, double *f, char *task, int iprint,
                   int info, int iter, int nfgv, int nintol, int nskip,
                   int nact, double sbgnrm, int nint,
                   char *word, int iback, double stp, double xstep,
                   int k)
{
  if(strncmp(task, "CONV", 4) == 0) {
    if (iprint >= 0) {
      Rprintf("\niterations %d\nfunction evaluations %d\nsegments explored during Cauchy searches %d\nBFGS updates skipped %d\nactive bounds at final generalized Cauchy point %d\nnorm of the final projected gradient %g\nfinal function value %g\n\n", iter, nfgv, nintol, nskip, nact, sbgnrm, *f);
    }
    if (iprint >= 100) pvector("X =", x, n);
    if (iprint >= 1) Rprintf("F = %g\n", *f);
  }
  if (iprint >= 0) {
    switch(info) {
    case -1: Rprintf("Matrix in 1st Cholesky factorization in formk is not Pos. Def."); break;
    case -2: Rprintf("Matrix in 2st Cholesky factorization in formk is not Pos. Def."); break;
    case -3: Rprintf("Matrix in the Cholesky factorization in formt is not Pos. Def."); break;
    case -4: Rprintf("Derivative >= 0, backtracking line search impossible."); break;
    case -5: Rprintf("l(%d) > u(%d).  No feasible solution", k, k); break;
    case -6: Rprintf("Input nbd(%d) is invalid", k); break;
    case -7: Rprintf("Warning:  more than 10 function and gradient evaluations\n   in the last line search\n"); break;
    case -8: Rprintf("The triangular system is singular."); break;
    case -9: Rprintf("%s\n%s\n", "Line search cannot locate an adequate point after 20 function", "and gradient evaluations"); break;
    default: break;
    }
  }
}









/* 
 AnDarl.c
$Revision: 1.1 $  $Date: 2014/06/09 07:53:21 $
Original C code by G. and J. Marsaglia
R interface by Adrian Baddeley
 https://github.com/cran/goftest
*/
static double adinf(double z) { 
  if(z<2.) return (
      exp(-1.2337141/z)/sqrt(z)
  )*(
      2.00012+(.247105-	
      (.0649821-
      (.0347962-
      (.011672-.00168691*z)
         *z)*z)*z)*z);
  /* max |error| < .000002 for z<2, (p=.90816...) */
  return exp(
    -exp(1.0776-(2.30695-(.43424-(.082433-(.008056 -.0003146*z)
                                    *z)*z)*z)*z));
                                    /* max |error|<.0000008 for 4<z<infinity */
}



static double AD(int n,double z){
  double c,v,x;
  x=adinf(z);
  /* now x=adinf(z). Next, get v=errfix(n,x) and return x+v; */
  if(x>.8) {
    v=(-130.2137+(745.2337-(1705.091-(1950.646-(1116.360-255.7844*x)
                                        *x)*x)*x)*x)/n;
    return x+v;
  }
  c=.01265+.1757/n;
  if(x<c){ 
    v=x/c;
    v=sqrt(v)*(1.-v)*(49*v-102);
    return(x+v*(.0037/(n*n)+.00078/n+.00006)/n);
  }
  v=(x-c)/(.8-c);
  v=-.00022633+(6.54034-(14.6538-(14.458-(8.259-1.91864*v)*v)*v)*v)*v;
  return (x+v*(.04213+.01365/n)/n);
}



static double ADtest(int n, double *x)
{ int i;
  double t,z=0;
  for(i=0;i<n;i++)   {
    t=x[i]*(1.-x[n-1-i]);
    z=z-(i+i+1)*log(t);
  }
  return AD(n,-n+z/n);
}






///////////////////////
class emp_mode_class {
public:
  emp_mode_class(std::vector<double> & data_){
    data = data_;
    n = data_.size();
    // std::sort(data.begin(), data.end());
    subsample_size = 1000;
    do_choose_m = true;
    do_density = false;
    n_points = 512;
  }

  ///////////////////

  
  void set_m(int m_){
    m = m_;
  }
  
  void set_do_choose_m(bool TorF){
    do_choose_m = TorF;
  }
  
  void set_do_density(bool TorF){
    do_density = TorF;
  }  
  
  void set_n_points(int np){
    n_points = np;
  }
  
  void set_sample_size(int t_){
    subsample_size = t_;
  }
  
  double mix_cdf_c(double u) {
    if(u < 0){
      return 0;
    }else if(u >= 1){
      return 1;
    }else if(u == 0){
      int cum_index = 0;
      while(cum_index < n && data[cum_index] <= 0){
        cum_index ++ ;
      }
      double Fn_value = cum_index / (double) n;
      return Fn_value;
    }else if ( coef_map_cdf_bern.find(m) == coef_map_cdf_bern.end() ) {
      // not found
      // coef.vec.cdf.bern <- Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1) / (m + 1)
      // coef.vec.cdf.bern <- log(Fn(0 : m / m) / (m + 1)) + dbeta_r(u, 0 : m + 1, m - (0 : m) + 1, log = TRUE)
      std::vector<double> coef_vec_cdf_bern(m + 1);
      std::vector<double> Fn_vec_cdf_bern(m + 1);
      double summation = 0;
      int cum_index = 0;
      double dbeta_log = m * log(1 - u);
      for(int i = 0; i < (m + 1); i++) {
        // double dbeta_log = R::dbeta(u, i + 1, m - i + 1, 1);
        if(i > 0){
          dbeta_log += log((m - i + 1) / (double) i * u / (1 - u));
        }
        double criterion = i / (double) m;
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value = cum_index / (double) n; // / (double) (m + 1);
        Fn_vec_cdf_bern[i] = Fn_value;
        coef_vec_cdf_bern[i] = log(Fn_value) + dbeta_log;
        summation += Fn_value * exp(dbeta_log);
      }
      // m.vec.cdf.bern <<- c(m.vec.cdf.bern, m)
      coef_map_cdf_bern[m] = coef_vec_cdf_bern;
      Fn_map_cdf_bern[m] = Fn_vec_cdf_bern;
      u_map_cdf_bern[m] = u;
      // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
      return summation;
    } else {
      // found
      double & u_post = u_map_cdf_bern[m];
      std::vector<double> & coef_vec_cdf_bern = coef_map_cdf_bern[m];
      // sum.cdf.bern <- crossprod(coef.ls.cdf.bern[[which(m.vec.cdf.bern == m)]], 
      // (u / u.post) ^ (0 : m) * ((1 - u) / (1 - u.post)) ^ (m : 0))
      double summation = 0;
      for(int i = 0; i < (m + 1); i++){
        summation += exp(coef_vec_cdf_bern[i] + 
          log(u / u_post) * i + log((1 - u) / (1 - u_post)) * (m - i));
      }

      if(isnan(summation) || isinf(summation) || summation < 0 || summation > 1){
        // coef.vec.cdf.bern <- Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1) / (m + 1)
        // coef.vec.cdf.bern <- log(Fn(0 : m / m) / (m + 1)) + dbeta(u, 0 : m + 1, m - (0 : m) + 1, log = TRUE)
        // u.vec.cdf.bern[which(m.vec.cdf.bern == m)] <<- u
        // coef.ls.cdf.bern[[which(m.vec.cdf.bern == m)]] <<- coef.vec.cdf.bern
        // return(sum(exp(coef.vec.cdf.bern)))
        summation = 0;
        std::vector<double> & Fn_vec_cdf_bern = Fn_map_cdf_bern[m];
        double dbeta_log = m * log(1 - u);
        for(int i = 0; i < (m + 1); i++){
          // double dbeta_log = R::dbeta(u, i + 1, m - i + 1, 1);
          if(i > 0){
            dbeta_log += log((m - i + 1) / (double) i * u / (1 - u));
          }
          double Fn_value = Fn_vec_cdf_bern[i];
          coef_vec_cdf_bern[i] = log(Fn_value) + dbeta_log;
          summation += Fn_value * exp(dbeta_log);
        }        
        // coef_map_cdf_bern[m] = coef_vec_cdf_bern;
        u_post = u;
        // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
        return summation;          
        
      }else{
        // return(sum.cdf.bern)
        return summation;
      }
    }
  }
  
  std::vector<double> cdf_bern_c(std::vector<double> & u){
    int n = u.size();
    std::vector<double> out(n);
    for(int i = 0; i < n; i++){
      out[i] = mix_cdf_c(u[i]);
    }
    return out;
  }
  
  
  double mix_pdf_c(double u) {
    if(u * (1 - u) < 0){
      return 0;
    }else if (u == 0){
      int cum_index = 0;
      while(cum_index < n && data[cum_index] <= 0){
        cum_index ++ ;
      }
      double Fn_value_0 = cum_index / (double) n;
      while(cum_index < n && data[cum_index] <= 1 / (double) m){
        cum_index ++ ;
      }
      double Fn_value_1 = cum_index / (double) n;
      return (Fn_value_1 - Fn_value_0) * m;
    }else if (u == 1){
      int cum_index = n;
      while(cum_index > 0 && data[cum_index - 1] > 1 - 1 / (double) m){
        cum_index -- ;
      }
      double Fn_value = cum_index / (double) n;
      return (1 - Fn_value) * m;
    }else if ( coef_map_pdf_bern.find(m) == coef_map_pdf_bern.end() ) {
      // not found
      // log(Fn((1 : m) / m) - Fn((0 : (m - 1)) / m)) + dbeta(u, 1 : m, m - (1 : m) + 1, log = TRUE)
      std::vector<double> coef_vec_pdf_bern(m);
      std::vector<double> Fn_vec_pdf_bern(m);
      double summation = 0;
      int cum_index = 0;
      double criterion = 0; // m
      while(cum_index < n && data[cum_index] <= criterion){
        cum_index ++ ;
      }
      double Fn_value_0 = cum_index / (double) n;
      double dbeta_log = log(m) + (m - 1) * log(1 - u);
      
      for(int i = 0; i < m; i++) {
        // double dbeta_log = R::dbeta(u, i + 1, m - i, 1);
        if(i > 0){
          dbeta_log += log((m - i) / (double) i * u / (1 - u));
        }
        criterion = (i + 1) / (double) m;
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value = cum_index / (double) n;
        Fn_vec_pdf_bern[i] = Fn_value - Fn_value_0;
        coef_vec_pdf_bern[i] = log(Fn_value - Fn_value_0) + dbeta_log;
        summation += (Fn_value - Fn_value_0) * exp(dbeta_log);
        Fn_value_0 = Fn_value;
      }
      // m.vec.pdf.bern <<- c(m.vec.pdf.bern, m)
      coef_map_pdf_bern[m] = coef_vec_pdf_bern;
      Fn_map_pdf_bern[m] = Fn_vec_pdf_bern;
      u_map_pdf_bern[m] = u;
      // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
      return summation;
    } else {
      // found
      double & u_post = u_map_pdf_bern[m];
      std::vector<double> & coef_vec_pdf_bern = coef_map_pdf_bern[m];
      // sum.pdf.bern <- crossprod(coef.ls.pdf.bern[[which(m.vec.pdf.bern == m)]], 
      // (u / u.post) ^ (0 : m) * ((1 - u) / (1 - u.post)) ^ (m : 0))
      double summation = 0;
      for(int i = 0; i < m; i++){
        summation += exp(coef_vec_pdf_bern[i] + 
          log(u / u_post) * i + log((1 - u) / (1 - u_post)) * (m - 1 - i));
      }
      
      if(isnan(summation) || summation < 0){
        // coef.vec.pdf.bern <- Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1) / (m + 1)
        // coef.vec.pdf.bern <- log(Fn(0 : m / m) / (m + 1)) + dbeta(u, 0 : m + 1, m - (0 : m) + 1, log = TRUE)
        // u.vec.pdf.bern[which(m.vec.pdf.bern == m)] <<- u
        // coef.ls.pdf.bern[[which(m.vec.pdf.bern == m)]] <<- coef.vec.pdf.bern
        // return(sum(exp(coef.vec.pdf.bern)))
        summation = 0;
        std::vector<double> & Fn_vec_pdf_bern = Fn_map_pdf_bern[m];
        double dbeta_log = log(m) + (m - 1) * log(1 - u);
        
        for(int i = 0; i < m; i++){
          // double dbeta_log = R::dbeta(u, i + 1, m - i, 1);
          if(i > 0){
            dbeta_log += log((m - i) / (double) i * u / (1 - u));
          }
          double Fn_value = Fn_vec_pdf_bern[i];
          coef_vec_pdf_bern[i] = log(Fn_value) + dbeta_log;
          summation += Fn_value * exp(dbeta_log);
        }        
        // coef_map_pdf_bern[m] = coef_vec_pdf_bern;
        u_post = u;
        // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
        return summation;          
        
      }else{
        // return(sum.pdf.bern)
        return summation;
      }
    }
  }
  
  std::vector<double> pdf_bern_c(std::vector<double> & u){
    int n = u.size();
    std::vector<double> out(n);
    for(int i = 0; i < n; i++){
      out[i] = mix_pdf_c(u[i]);
    }
    return out;
  }
  
  
  double mix_dpdf_c(double u) {
      // Obtain environment containing function
      // Rcpp::Environment base("package:stats"); 
      // Make function callable from C++
      // Rcpp::Function dbeta_r = base["dbeta"];   
      // Rcpp::Function sum_r = base["sum"];   
      
      if(u * (1 - u) < 0){
        return 0;
      }else if (u == 0){
        int cum_index = 0;
        while(cum_index < n && data[cum_index] <= 0){
          cum_index ++ ;
        }
        double Fn_value_0 = cum_index / (double) n;
        while(cum_index < n && data[cum_index] <= 1 / (double) m){
          cum_index ++ ;
        }
        double Fn_value_1 = cum_index / (double) n;
        while(cum_index < n && data[cum_index] <= 2 / (double) m){
          cum_index ++ ;
        }
        double Fn_value_2 = cum_index / (double) n;
        return (Fn_value_2 - 2 * Fn_value_1 + Fn_value_0) * m * (m - 1);
      }else if (u == 1){
        int cum_index = n;
        while(cum_index > 0 && data[cum_index - 1] > 1 - 1 / (double) m){
          cum_index -- ;
        }
        double Fn_value = cum_index / (double) n;
        while(cum_index > 0 && data[cum_index - 1] > 1 - 2 / (double) m){
          cum_index -- ;
        }
        double Fn_value_ = cum_index / (double) n;
        return (1 - 2 * Fn_value + Fn_value_) * m * (m - 1);
      }else if ( coef_map_dpdf_bern.find(m) == coef_map_dpdf_bern.end() ) {
        // not found
        // (Fn(2 : m / m) - 2 * Fn((2 : m - 1) / m) + Fn((2 : m - 2) / m)) * m * dbeta(u, 2 : m - 1, m - (2 : m) + 1)
        std::vector<double> coef_vec_dpdf_bern(m - 1);
        std::vector<double> Fn_vec_dpdf_bern(m - 1);
        std::vector<bool> sign_vec_dpdf_bern(m - 1);
        double summation = 0;
        int cum_index = 0;
        double criterion = 0; // m
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value_0 = cum_index / (double) n;
        
        criterion = 1 / (double) m;
        while(cum_index < n && data[cum_index] <= criterion){
          cum_index ++ ;
        }
        double Fn_value_1 = cum_index / (double) n;
        double dbeta_log = log(m * (m - 1)) + (m - 2) * log(1 - u);
        
        for(int i = 0; i < (m - 1); i++) {
          // double dbeta_log = R::dbeta(u, i + 1, m - i - 1, 1);
          if(i > 0){
            dbeta_log += log((m - i - 1) / (double) i * u / (1 - u));
          }
          criterion = (i + 2) / (double) m;
          while(cum_index < n && data[cum_index] <= criterion){
            cum_index ++ ;
          }
          double Fn_value = cum_index / (double) n;
          Fn_vec_dpdf_bern[i] = (Fn_value - 2 * Fn_value_1 + Fn_value_0); // * m;
          sign_vec_dpdf_bern[i] = Fn_value - 2 * Fn_value_1 + Fn_value_0 > 0;
          coef_vec_dpdf_bern[i] = log(std::abs(Fn_value - 2 * Fn_value_1 + Fn_value_0)) + dbeta_log; // * m;
          summation += (Fn_value - 2 * Fn_value_1 + Fn_value_0) * exp(dbeta_log); // * m;
          Fn_value_0 = Fn_value_1;
          Fn_value_1 = Fn_value;
        }
        // m.vec.dpdf.bern <<- c(m.vec.dpdf.bern, m)
        coef_map_dpdf_bern[m] = coef_vec_dpdf_bern;
        sign_map_dpdf_bern[m] = sign_vec_dpdf_bern;
        Fn_map_dpdf_bern[m] = Fn_vec_dpdf_bern;
        u_map_dpdf_bern[m] = u;
        // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
        return summation;
      } else {
        // found
        double & u_post = u_map_dpdf_bern[m];
        std::vector<double> & coef_vec_dpdf_bern = coef_map_dpdf_bern[m];
        std::vector<bool> & sign_vec_dpdf_bern = sign_map_dpdf_bern[m];
        // sum.dpdf.bern <- crossprod(coef.ls.dpdf.bern[[which(m.vec.dpdf.bern == m)]], 
        // (u / u.post) ^ (0 : m) * ((1 - u) / (1 - u.post)) ^ (m : 0))
        double summation = 0;
        for(int i = 0; i < (m - 1); i++){
          summation += (2 * sign_vec_dpdf_bern[i] - 1) * exp(coef_vec_dpdf_bern[i] + 
            log(u / u_post) * i + log((1 - u) / (1 - u_post)) * (m - 2 - i));
        }
        
        if(isnan(summation)){
          summation = 0;
          std::vector<double> & Fn_vec_dpdf_bern = Fn_map_dpdf_bern[m];
          double dbeta_log = log(m * (m - 1)) + (m - 2) * log(1 - u);
          
          for(int i = 0; i < (m - 1); i++){
            // double dbeta_log = R::dbeta(u, i + 1, m - i - 1, 1);
            if(i > 0){
              dbeta_log += log((m - i - 1) / (double) i * u / (1 - u));
            }
            double Fn_value = Fn_vec_dpdf_bern[i];
            coef_vec_dpdf_bern[i] = log(std::abs(Fn_value)) + dbeta_log;
            summation += Fn_value * exp(dbeta_log);
          }        
          // coef_map_dpdf_bern[m] = coef_vec_dpdf_bern;
          u_post = u;
          // return(sum(Fn(0 : m / m) * dbeta(u, 0 : m + 1, m - (0 : m) + 1)) / (m + 1)) # 1 / (m+2) 
          return summation;          
          
        }else{
          // return(sum.dpdf.bern)
          return summation;
        }
      }
    }
  
  
  std::vector<double> dpdf_bern_c(std::vector<double> & u){
    int n = u.size();
    std::vector<double> out(n);
    for(int i = 0; i < n; i++){
      out[i] = mix_dpdf_c(u[i]);
    }
    return out;
  }
  
  
  double fminfn(double * x){return -mix_pdf_c(*x); }
  
  void fmingr(double * x, double * res){res[0] = -mix_dpdf_c(*x); }
  
  
  //////////////////////////////
  // https://github.com/Microsoft/microsoft-r-open/blob/master/source/src/appl/optim.c
  double * vect(int n)
  {
    return (double *)R_alloc(n, sizeof(double));
  }
  
  void lbfgsb(double *x, double *l, double *u, int *nbd, char *msg, int *fncount, int *grcount, 
              double *Fmin, int *fail,
              int n = 1, int m = 5, double factr = 1e7, double pgtol = 0.,
              int maxit = 100, int trace = 0, int nREPORT = 10)
  {
    char task[60];
    double f, *g,  *wa;
    int tr = -1, iter = 0, *iwa, isave[21];
    isave[12] = 0; // -Wall
    tr = -1; 
    *fail = 0;
    g = vect(n);
    /* this needs to be zeroed for snd in mainlb to be zeroed */ 
    // std::vector<double> wa_vec (2*m*n+4*n+11*m*m+8*m);
    // wa = &wa_vec[0];
    // std::vector<int> iwa_vec (3*n);
    // iwa = &iwa_vec[0];
    wa = (double *) S_alloc(2*m*n+4*n+11*m*m+8*m, sizeof(double));
    iwa = (int *) R_alloc(3*n, sizeof(int));
    strcpy(task, "START");
    while(1) {
      setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task,
             tr, isave);
      /*	Rprintf("in lbfgsb - %s\n", task);*/
      if (strncmp(task, "FG", 2) == 0) {
        f = fminfn(x);
        if (!R_FINITE(f))
          break; // error(("L-BFGS-B needs finite values of 'fn'"));
        fmingr(x, g);
      } else if (strncmp(task, "NEW_X", 5) == 0) {
        iter++;
        if(trace == 1 && (iter % nREPORT == 0)) {
          Rprintf("iter %4d value %f\n", iter, f);
        }
        if (iter > maxit) {
          *fail = 1;
          break;
        }
      } else if (strncmp(task, "WARN", 4) == 0) {
        *fail = 51;
        break;
      } else if (strncmp(task, "CONV", 4) == 0) {
        break;
      } else if (strncmp(task, "ERROR", 5) == 0) {
        *fail = 52;
        break;
      } else { /* some other condition that is not supposed to happen */
      *fail = 52;
        break;
      }
    }
    *Fmin = f;
    *fncount = *grcount = isave[12];
    strcpy(msg, task);
  }
  
  
  // https://github.com/SurajGupta/r-source/blob/master/src/library/stats/src/optimize.c
  /* Formerly in src/appl/fmim.c 
  fmin.f -- translated by f2c (version 19990503).
  R's  optimize() :   function	fmin(ax,bx,f,tol)
  =    ==========		~~~~~~~~~~~~~~~~~
  an approximation  x  to the point where  f  attains a minimum  on
  the interval  (ax,bx)  is determined.
  INPUT..
  ax    left endpoint of initial interval
  bx    right endpoint of initial interval
  f     function which evaluates  f(x, info)  for any  x
  in the interval  (ax,bx)
  tol   desired length of the interval of uncertainty of the final
  result ( >= 0.)
  OUTPUT..
  fmin  abcissa approximating the point where  f  attains a minimum */
  
  double Brent_fmin(double ax, double bx, double tol = pow(DBL_EPSILON, 0.25))
  {
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;
    
    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    
    /*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);
    
    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    
    d = 0.;/* -Wall */
    e = 0.;
    fx = fminfn(&x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    
    /*  main loop starts here ----------------------------------- */
    
    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
      
      /* check stopping criterion */
      
      if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) { /* fit parabola */
      
      r = (x - w) * (fx - fv);
        q = (x - v) * (fx - fw);
        p = (x - v) * q - (x - w) * r;
        q = (q - r) * 2.;
        if (q > 0.) p = -p; else q = -q;
        r = e;
        e = d;
      }
      
      if (fabs(p) >= fabs(q * .5 * r) ||
          p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */
      
      if (x < xm) e = b - x; else e = a - x;
      d = c * e;
      }
      else { /* a parabolic-interpolation step */
      
      d = p / q;
        u = x + d;
        
        /* f must not be evaluated too close to ax or bx */
        
        if (u - a < t2 || b - u < t2) {
          d = tol1;
          if (x >= xm) d = -d;
        }
      }
      
      /* f must not be evaluated too close to x */
      
      if (fabs(d) >= tol1)
        u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;
      
      fu = fminfn(&u);
      
      /*  update  a, b, v, w, and x */
      
      if (fu <= fx) {
        if (u < x) b = x; else a = x;
        v = w;    w = x;   x = u;
        fv = fw; fw = fx; fx = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w == x) {
          v = w; fv = fw;
          w = u; fw = fu;
        } else if (fu <= fv || v == x || v == w) {
          v = u; fv = fu;
        }
      }
    }
    /* end of main loop */
    
    return x;
  }
  
  
  
  
  Rcpp::List emp_mode_c(double & fix_lower, double & fix_upper){
    Rcpp::List res;
    double min_data = *data.begin();
    // double max_data = data.back();
    
    if(data.size() == 1){
      res["optimal alpha"] = R_NilValue;
      res["optimal p-value"] = R_NilValue;
      res["mode"] = min_data;
      Rcpp::Rcout << "The sample size is 1 and no smoothing is applied." << std::endl;
      return res;
    }

    // double delta = 1 / (double) n; 
    lower = fix_lower;
    upper = fix_upper;
      
    for(int i = 0; i < n; i++){
      data[i] = (data[i] - lower) / (upper - lower);
    }
    
    
    
    std::vector<double> newdata;
    if(n > subsample_size){
      std::vector<double> newdata_(subsample_size);
      for(int i = 0 ; i < subsample_size ; i++){
        newdata_[i] = data[n * (i + 1) / (subsample_size + 1)];
      }
      newdata = newdata_;
    }else{
      newdata = data;
    }
    
    
    double max_pvalue = -1;
    double max_alpha = -1; 
    if(do_choose_m){
    for (double alpha = 0.3; alpha <= 0.95; alpha += 0.05) {
      // Rcpp::Rcout << alpha;
      set_m((int) ceil(pow(n, alpha)));
      std::vector<double> v = cdf_bern_c(newdata);
      double* x = &v[0];
      double pvalue_curr = 1. - ADtest(newdata.size(), x);
      
      // int len = newdata.size();
      // // ADtest
      // double t, z = 0;
      // for(int i = 0; i < len; i++)   {
      //   t = mix_cdf_c(newdata[i]) * (1. - mix_cdf_c(newdata[len - 1 - i]));
      //   z = z - (i + i + 1) * log(t);
      // }
      // double pvalue_curr = AD(len, -len + z / len);
      
      
      // double adstat = ADstat(newdata.size(), x);
      if(pvalue_curr > max_pvalue){
        max_pvalue = pvalue_curr;
        max_alpha = alpha;
      }
      if(max_pvalue > 0.95){
        break;
      }
    }
    set_m((int) ceil(pow(n, max_alpha)));
    }else{
      std::vector<double> v = cdf_bern_c(newdata);
      double* x = &v[0];
      double pvalue_curr = 1. - ADtest(newdata.size(), x);
      max_pvalue = pvalue_curr;
    }
    
    
    // gradient descent
    double max_x = -1;
    double max_density = -1;
    
    double pdfpar_vec[] = {0, 0.5, 1};
    for(int i = 0; (unsigned) i < sizeof(pdfpar_vec) / sizeof(pdfpar_vec[0]); ++i){
    try{
      // https://github.com/Microsoft/microsoft-r-open/blob/master/source/src/library/stats/src/optim.c
      double dpar = pdfpar_vec[i]; // 0.5, 1
      double lower = 0;
      double upper = 1;
      int nbd = 2;
      double val = 0.0;
      int ifail = 0;
      char msg[60];
      int fncount = 0, grcount = 0;
      lbfgsb(&dpar, &lower, &upper, &nbd, msg, &fncount, &grcount, &val, &ifail); 
      if(-val > max_density){
        max_x = dpar;
        max_density = -val;
      }
    }catch(std::exception e){
      Rcpp::Rcout << e.what() << std::endl;
    }
    }
    
    double pdfilon_vec[] = {0, 0.1, 0.5};
    for(int i = 0; (unsigned) i < sizeof(pdfilon_vec) / sizeof(pdfilon_vec[0]); ++i){
      try{
        double dpar = Brent_fmin(0 - pdfilon_vec[i], 1 + pdfilon_vec[i]);
        double val = mix_pdf_c(dpar);
        if(val > max_density){
          max_x = dpar;
          max_density = val;
        }
      }catch(std::exception e){
        Rcpp::Rcout << e.what() << std::endl;
      }
    }
          
    std::vector<double> & Fn_vec_cdf_bern = Fn_map_cdf_bern[m];
    double prob_left = Fn_vec_cdf_bern[2];
    double prob_right = 1 - Fn_vec_cdf_bern[Fn_vec_cdf_bern.size() - 3];
    int cum_index = 0;
    while(cum_index < n && data[cum_index] <= max_x - 1 / (double) m){
      cum_index ++ ;
    }
    double Fn_value_0 = cum_index / (double) n;
    while(cum_index < n && data[cum_index] <= max_x + 1 / (double) m){
      cum_index ++ ;
    }
    double Fn_value_1 = cum_index / (double) n;
    double prob_mid = Fn_value_1 - Fn_value_0;
    
    if(max_pvalue > 0.05 || (prob_mid > prob_left && prob_mid > prob_right)){
      Rcpp::Rcout << "The mode estimate is choosen as the maximum point of Bernstein polynomials density estimate." << std::endl;
      res["mode"] = max_x * (upper - lower) + lower;
    }else if(prob_left > prob_right){
      Rcpp::Rcout << "The mode estimate is choosen as a boundary point and Anderson-Darling criterion is not satisfied with p-value less than 0.05." << std::endl;
      res["mode"] = 0 * (upper - lower) + lower;
    }else if(prob_left < prob_right){
      Rcpp::Rcout << "The mode estimate is choosen as a boundary point and Anderson-Darling criterion is not satisfied with p-value less than 0.05." << std::endl;
      res["mode"] = 1 * (upper - lower) + lower;
    }else{
      Rcpp::Rcout << "The mode estimate is choosen as a boundary point and Anderson-Darling criterion is not satisfied with p-value less than 0.05." << std::endl;
      res["mode"] = floor(rand() / double(RAND_MAX) * 2) * (upper - lower) + lower;
    }
    
    if(do_choose_m){
      res["optimal alpha"] = max_alpha;
      res["optimal p-value"] = max_pvalue;
    }else{
      res["given degree"] = m;
      res["p-value"] = max_pvalue;
    }
    
    if(do_density){
      std::vector<double> x_cord(n_points);
      for(int i = 0; i < n_points; i++){
        x_cord[i] = i / (double) n_points;
      }
      res["y"] = pdf_bern_c(x_cord);
      for(int i = 0; i < n_points; i++){
        x_cord[i] = i / (double) n_points * (upper - lower) + lower;
      }
      res["x"] = x_cord;
    }
    
  return res;
  }
  
  
// private:
  //////////////////////
  std::vector<double> data;
  int n;
  int m;
  int subsample_size;
  double lower, upper;
  bool do_choose_m;
  bool do_density;
  int n_points;

  std::map<int, std::vector<double> > coef_map_cdf_bern;
  std::map<int, std::vector<double> > Fn_map_cdf_bern;
  std::map<int, double> u_map_cdf_bern;
  
  std::map<int, std::vector<double> > coef_map_pdf_bern;
  std::map<int, std::vector<double> > Fn_map_pdf_bern;
  std::map<int, double> u_map_pdf_bern;
  
  std::map<int, std::vector<double> > coef_map_dpdf_bern;
  std::map<int, std::vector<bool> > sign_map_dpdf_bern;
  std::map<int, std::vector<double> > Fn_map_dpdf_bern;
  std::map<int, double> u_map_dpdf_bern;
};








