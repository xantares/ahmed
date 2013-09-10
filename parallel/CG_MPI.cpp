/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <mpi.h>
#include "matrix.h"
#include "blas.h"

extern MPI::Intracomm COMM_AHMED;

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


unsigned CG_MPI(const Matrix<double>& A, double* const b, double* const x,
                double& eps, unsigned& nsteps)
{
  unsigned rank  = COMM_AHMED.Get_rank(),
           nproc = COMM_AHMED.Get_size();

  double nrmb, resid;
  unsigned N = A.n;
  double *p, *q, *r, *rhat, rho, rho1 = 0.0;

  if (rank==0) {
    p = new double[4*N];
    q = p + N;
    r = q + N;
    rhat = r + N;

    nrmb = blas::nrm2(N, b);
    // direct other processors whether to continue of to leave
    unsigned inf = (nrmb<D_PREC);
    for (unsigned j=1; j<nproc; ++j)
      COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, j, 9);

    if (inf) {
      blas::setzero(N, x);
      eps = 0.0;
      nsteps = 0;
      delete [] p;
      return 0;
    }

    // r0 = b - Ax0
    blas::copy(N, b, r);
  }
  else {
    unsigned inf;
    COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
    if (inf) return 0;
  }

  A.amux(D_MONE, x, r);

  if (rank==0) {
    resid = blas::nrm2(N, r);

    // direct other processors whether to continue of to leave
    unsigned inf = (resid<=eps*nrmb);
    for (unsigned j=1; j<nproc; ++j)
      COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, j, 9);

    if (inf) {
      eps = resid/nrmb;
      nsteps = 0;
      delete [] p;
      return 0;
    }
  } else {
    unsigned inf;
    COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
    if (inf) return 0;
  }

  for (unsigned l=1; l<=nsteps; ++l) {

    if (rank==0) {
#ifndef NDEBUG
      std::cout << "Step " << l << ", resid=" << resid/nrmb << std::endl;
#endif

      // r^ = Cr r
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
      else blas::copy(N, rhat, p);

      // q = Ap
      blas::setzero(N, q);
    }

    A.amux(D_ONE, p, q);

    if (rank==0) {

      // alpha = rho / p*q
      double alpha = rho / blas::scpr(N, p, q);

      // x += alpha * p
      blas::axpy(N, alpha, p, x);

      // r -= alpha * q
      blas::axpy(N, -alpha, q, r);

      resid = blas::nrm2(N, r);

      // direct other processors whether to continue of to leave
      unsigned inf = (resid<=eps*nrmb);
      for (unsigned j=1; j<nproc; ++j)
        COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, j, 9);

      if (inf) {
        eps = resid/nrmb;
        nsteps = l;
        delete [] p;
        return 0;
      }

      rho1 = rho;
    } else {
      unsigned inf;
      COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
      if (inf) return 0;
    }
  }

  if (rank==0)   {
    eps = resid;
    delete [] p;
    return 1;
  } else return 1;
}


