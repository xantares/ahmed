/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "matrix.h"
#include "blas.h"

// solves Ax=b for x using the stabilized bicg algorithm
// input:  x starting vector, eps D_PRECision, N number of unknowns
//         bl a block cluster, A sequence of Matrix<double> blocks,
//         b right hand side,
//         nsteps max. number of steps, C the D_PREConditioner
// output: x the (approximate) solution vector, eps attained D_PRECision,
//         nsteps the number of steps
//         if the algorithm converges 0 is returned, otherwise a nonzero number

unsigned BiCGStab(const Matrix<double>& A, double* const b, double* const x,
                  double& eps, unsigned& nsteps)
{
  double *v, *p, *phat, *s, *shat, *t, *r, *r0;
  unsigned N = A.n;

  v = new double[8*N];
  p = v + N;
  phat = p + N;
  s = phat + N;
  shat = s + N;

  t = shat + N;
  r = t + N;
  r0 = r + N;

  double resid;

  // normb = |b|
  double nrmb = blas::nrm2(N, b);
  if (nrmb<D_PREC) nrmb = D_ONE;

  // r = r0 = b - A x0
  blas::copy(N, b, r0);
  A.amux(D_MONE, x, r0);
  blas::copy(N, r0, r);

  resid = blas::nrm2(N, r) / nrmb;

  if (resid<eps) {
    eps = resid;
    nsteps = 0;
    delete [] v;
    return 0;
  }

  double alpha = D_ZERO, omega = D_ZERO, rho2 = D_ZERO;

  for (unsigned l=1; l<=nsteps; ++l) {

    // rho1 = r0 * r
    const double rho1 = blas::scpr(N, r0, r);
    if (fabs(rho1)<D_PREC) {
      eps = blas::nrm2(N, r) / nrmb;
      delete [] v;
      return 2;
    }

    if (l==1)
      blas::copy(N, r, p);                // p = r
    else {
      blas::axpy(N, -omega, v, p);        // p = (p-omega v)*beta+r
      const double beta = rho1*alpha/(rho2*omega);
      blas::scal(N, beta, p);
      blas::axpy(N, D_ONE, r, p);
    }

    // p^ = C p
    blas::copy(N, p, phat);
    A.precond_apply(phat);

    // v = A p^
    blas::setzero(N, v);
    A.amux(D_ONE, phat, v);

    alpha = rho1 / blas::scpr(N, r0, v);

    // s = r - alpha v
    blas::copy(N, r, s);
    blas::axpy(N, -alpha, v, s);

    resid = blas::nrm2(N, s) / nrmb;
#ifndef NDEBUG
    std::cout << "Step " << l << ", resid=" << resid << std::endl;
#endif
    if (resid<eps) {
      // x += alpha p^
      blas::axpy(N, alpha, phat, x);
      eps = resid;
      nsteps = l;
      delete [] v;
      return 0;
    }

    // s^ = C s
    blas::copy(N, s, shat);
    A.precond_apply(shat);

    // t = A s^
    blas::setzero(N, t);
    A.amux(D_ONE, shat, t);

    // omega = t*s / t*t
    omega = blas::scpr(N, t, s) / blas::scpr(N, t, t);

    // x += alpha p^ + omega s^
    blas::axpy(N, alpha, phat, x);
    blas::axpy(N, omega, shat, x);

    // r = s - omega t
    blas::copy(N, s, r);
    blas::axpy(N, -omega, t, r);

    rho2 = rho1;

    resid = blas::nrm2(N, r) / nrmb;

    if (resid<eps) {
      eps = resid;
      nsteps = l;
      delete [] v;
      return 0;
    }

    if (fabs(omega)<D_PREC) {
      eps = resid;
      delete [] v;
      return 3;
    }
  }

  eps = resid;
  delete [] v;
  return 1;
}


unsigned BiCGStab(const Matrix<dcomp>& A, dcomp* const b, dcomp* const x,
                  double& eps, unsigned& nsteps)
{
  dcomp *v, *p, *phat, *s, *shat, *t, *r, *r0;
  unsigned N = A.n;

  v = new dcomp[8*N];
  p = v + N;
  phat = p + N;
  s = phat + N;
  shat = s + N;

  t = shat + N;
  r = t + N;
  r0 = r + N;

  double resid;

  // normb = |b|
  double nrmb = blas::nrm2(N, b);
  if (nrmb<D_PREC) nrmb = D_ONE;

  // r = r0 = b - A x0
  blas::copy(N, b, r0);
  A.amux(Z_MONE, x, r0);
  blas::copy(N, r0, r);

  resid = blas::nrm2(N, r) / nrmb;

  if (resid<eps) {
    eps = resid;
    nsteps = 0;
    delete [] v;
    return 0;
  }

  dcomp alpha = Z_ZERO, omega = Z_ONE, rho2 = Z_ONE;

  for (unsigned l=1; l<=nsteps; ++l) {

    // rho1 = r0 * r
    const dcomp rho1 = blas::scpr(N, r0, r);
    if (abs(rho1)<D_PREC) {
      eps = blas::nrm2(N, r) / nrmb;
      delete [] v;
      return 2;
    }

    if (l==1)
      blas::copy(N, r, p);                // p = r
    else {
      blas::axpy(N, D_MONE*omega, v, p);        // p = (p-omega v)*beta+r
      const dcomp beta = rho1*alpha/(rho2*omega);
      blas::scal(N, beta, p);
      blas::axpy(N, Z_ONE, r, p);
    }

    // p^ = C p
    blas::copy(N, p, phat);
    A.precond_apply(phat);

    // v = A p^
    blas::setzero(N, v);
    A.amux(Z_ONE, phat, v);

    alpha = rho1 / blas::scpr(N, r0, v);

    // s = r - alpha v
    blas::copy(N, r, s);
    blas::axpy(N, D_MONE*alpha, v, s);

    resid = blas::nrm2(N, s) / nrmb;
#ifndef NDEBUG
    std::cout << "Step " << l << ", resid=" << resid << std::endl;
#endif
    if (resid<eps) {
      // x += alpha p^
      blas::axpy(N, alpha, phat, x);
      eps = resid;
      nsteps = l;
      delete [] v;
      return 0;
    }

    // s^ = C s
    blas::copy(N, s, shat);
    A.precond_apply(shat);

    // t = A s^
    blas::setzero(N, t);
    A.amux(Z_ONE, shat, t);

    // omega = t*s / t*t
    omega = blas::scpr(N, t, s) / blas::scpr(N, t, t);

    // x += alpha p^ + omega s^
    blas::axpy(N, alpha, phat, x);
    blas::axpy(N, omega, shat, x);

    // r = s - omega t
    blas::copy(N, s, r);
    blas::axpy(N, D_MONE*omega, t, r);

    rho2 = rho1;

    resid = blas::nrm2(N, r) / nrmb;

    if (resid<eps) {
      eps = resid;
      nsteps = l;
      delete [] v;
      return 0;
    }

    if (abs(omega)<D_PREC) {
      eps = resid;
      delete [] v;
      return 3;
    }
  }

  eps = resid;
  delete [] v;
  return 1;
}

