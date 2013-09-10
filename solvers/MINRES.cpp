/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include<cmath>
#include "blas.h"
#include "matrix.h"

// Preconditioner for A has to be positive definite
unsigned MinRes(const Matrix<double>& A, double* const b, double* const x,
                double& tol, unsigned& nsteps)
{
  unsigned n = A.n, k = 0;
  double *u = new double[n], *uold = new double[n];
  double *w = new double[n], *wold = new double[n], *woold = new double[n];
  double *v = new double[n], *vold = new double[n];
  double *r = new double[n], *z = new double[n];
  double alpha, beta, eta, np, np0, c=1.0, s=0.0, cold=1.0, sold=0.0;

  blas::setzero(n, uold);                 // u_old := 0
  blas::setzero(n, vold);                 // v_old := 0
  blas::setzero(n, w);                    // w     := 0
  blas::setzero(n, wold);                 // w_old := 0

  blas::copy(n, b, r);                    // r := b - A*x
  A.amux(D_MONE, x, r);
  blas::copy(n, r, z);                    // z := C*r
  A.precond_apply(z);

  beta = sqrt(blas::scpr(n, r, z));       // beta := sqrt(r'*z)
  if (beta<1e-10) goto LEAVE;             // happy breakdown
  eta = beta;

  blas::scal(n, 1.0/beta, r, v);          // v  := r / beta
  blas::scal(n, 1.0/beta, z, u);          // u  := z / beta
  np0 = np = blas::nrm2(n, z);            // np := ||z||

  do {
    // Lanczos
    blas::setzero(n, r);
    A.amux(D_ONE, u, r);                 // r     := A*u
    alpha = blas::scpr(n, r, u);         // alpha := r'*u
    blas::copy(n, r, z);                 // z     := B*r
    A.precond_apply(z);

    blas::axpy(n, -alpha, v, r);         // r := r - alpha v
    blas::axpy(n, -beta, vold, r);       // r := r - beta v_old
    blas::axpy(n, -alpha, u, z);         // z := z - alpha u
    blas::axpy(n, -beta, uold, z);       // z := z - beta u_old

    double betaold = beta;

    beta = sqrt(blas::scpr(n, r, z));    // beta := sqrt(r'*z)
    if (beta<1e-14) goto LEAVE;          // happy breakdown

    // QR factorisation

    double coold = cold, soold = sold;
    cold = c;
    sold = s;

    double rho0 = cold * alpha - coold * sold * betaold;
    double rho1 = sqrt(rho0*rho0 + beta*beta);
    double rho2 = sold * alpha + coold * cold * betaold;
    double rho3 = soold * betaold;

    // Givens rotation

    c = rho0 / rho1;
    s = beta / rho1;

    // Update

    blas::copy(n, wold, woold);          // w_oold := w_old
    blas::copy(n, w, wold);              // w_old  := w

    blas::copy(n, u, w);                 // w := u
    blas::axpy(n, -rho2, wold, w);       // w := w - rho2 w_old
    blas::axpy(n, -rho3, woold, w);      // w := w - rho3 w_oold
    blas::scal(n, 1.0/rho1, w);          // w := w / rho1

    blas::axpy(n, c*eta, w, x);          // x := x + c eta w
    eta *= - s;

    blas::copy(n, v, vold);
    blas::scal(n, 1.0/beta, r, v);       // v := r / beta
    blas::copy(n, u, uold);
    blas::scal(n, 1.0/beta, z, u);       // u := z / beta

    np *= fabs(s);

  } while (np>tol*np0 && ++k<nsteps);

  tol = np/np0;

LEAVE:
  tol = np/np0;
  delete [] z;
  delete [] r;
  delete [] vold;
  delete [] v;
  delete [] woold;
  delete [] wold;
  delete [] w;
  delete [] uold;
  delete [] u;
  if (k<nsteps) {
    nsteps = k;
    return 0;
  } else return 1;
}

