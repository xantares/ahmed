/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "matrix.h"
#include "blas.h"

// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// A and C are symmetric H-matrices
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//      x  --  approximate solution to Ax = b
// nsteps  --  the number of iterations performed before the
//             tolerance was reached
//    eps  --  the residual after the final iteration


unsigned CG(const Matrix<double>& A, double* const b, double* const x,
            double& eps, unsigned& nsteps)
{
  unsigned N = A.n;
  double *p, *q, *r, *rhat, rho, rho1 = 0.0;

  p = new double[4*N];
  q = p + N;
  r = q + N;
  rhat = r + N;

  double nrmb = blas::nrm2(N, b);
  if (nrmb<D_PREC) {
    blas::setzero(N, x);
    eps = 0.0;
    nsteps = 0;
    delete [] p;
    return 0;
  }

  // r0 = b - Ax0
  blas::copy(N, b, r);
  A.amux(D_MONE, x, r);

  double resid = blas::nrm2(N, r);
  if (resid<=eps*nrmb) {
    eps = resid/nrmb;
    nsteps = 0;
    delete [] p;
    return 0;
  }

  for (unsigned l=1; l<=nsteps; ++l) {

#ifndef NDEBUG
    std::cout << "Step " << l << ", resid=" << resid/nrmb << std::endl;
#endif

    // r^ = C r
    blas::copy(N, r, rhat);
    A.precond_apply(rhat);

    // rho = r * r^;
    rho = blas::scpr(N, r, rhat);

    if (l>1) {
      double beta = rho / rho1;
      // p = r^ + beta * p
      blas::scal(N, beta, p);
      blas::axpy(N, D_ONE, rhat, p);
    }
    else
      blas::copy(N, rhat, p);

    // q = Ap
    blas::setzero(N, q);
    A.amux(D_ONE, p, q);

    // alpha = rho / p*q
    double alpha = rho / blas::scpr(N, p, q);

    // x += alpha * p
    blas::axpy(N, alpha, p, x);

    // r -= alpha * q
    blas::axpy(N, -alpha, q, r);

    resid = blas::nrm2(N, r);

    if (resid<=eps*nrmb) {
      eps = resid/nrmb;
      nsteps = l;
      delete [] p;
      return 0;
    }

    rho1 = rho;
  }

  eps = resid / nrmb;
  delete [] p;
  return 1;
}


