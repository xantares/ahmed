/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method
//
// The return value indicates convergence within nsteps (input)
// iterations (0), or no convergence within nsteps iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//      x  --  approximate solution to Ax = b
// nsteps  --  the number of iterations performed before the
//             tolerance was reached
//    eps  --  the residual after the final iteration
//
//*****************************************************************
#include <mpi.h>
#include "blas.h"
#include "matrix.h"

extern MPI::Intracomm COMM_AHMED;

static void genPlRot(double dx, double dy, double& cs, double& sn)
{
  if (dy==0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (fabs(dy)>fabs(dx)) {
    double tmp = dx / dy;
    sn = 1.0 / sqrt(1.0+tmp*tmp);
    cs = tmp * sn;
  } else {
    double tmp = dy / dx;
    cs = 1.0 / sqrt(1.0+tmp*tmp);
    sn = tmp * cs;
  }
}


inline void applPlRot(double& dx, double& dy, double cs, double sn)
{
  double tmp = cs*dx + sn*dy;
  dy = cs*dy - sn*dx;
  dx = tmp;
}


// solve H y = s and update x += MVy
static void update(const Matrix<double>& A, unsigned k, double* H,
                   unsigned ldH, double* s, double* V, double* x)
{
  double *y = new double[k];
  double *xh = new double[A.n];
  blas::copy(k, s, y);
  int inf;

  dtrtrs_(JOB_STR+5, JOB_STR, JOB_STR, &k, &N_ONE, H, &ldH, y, &k, &inf);
  assert(inf==0);

  // x += M V y
  blas::setzero(A.n, xh);
  blas::gemva(A.n, k, D_ONE, V, y, xh);
  A.precond_apply(xh);
  blas::add(A.n, xh, x);

  delete [] xh;
  delete [] y;
}

unsigned GMRes_MPI(const Matrix<double>& A, double* const b, double* const x,
                   double& eps, const unsigned m, unsigned& nsteps)
{
  unsigned rank  = COMM_AHMED.Get_rank(),
           nproc = COMM_AHMED.Get_size();

  double resid, normb, beta;
  unsigned j = 1;
  double *r, *V, *H, *cs, *sn, *s, *xh;

  if (rank==0) {
    r = new double[2*A.n+(A.n+m+4)*(m+1)];   // n
    V = r + A.n;                         // n x (m+1)
    H = V + A.n*(m+1);                   // m+1 x m
    cs = H + (m+1)*m;                    // m+1
    sn = cs + m+1;                       // m+1
    s = sn + m+1;                        // m+1
    xh = s + m+1;                        // m+1

    // normb = norm(b)
    normb = blas::nrm2(A.n, b);

    bool inf = (normb<D_PREC);
    for (unsigned p=1; p<nproc; ++p)
      COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

    if (inf) {
      blas::setzero(A.n, x);
      eps = 0.0;
      nsteps = 0;
      delete [] r;
      return 0;
    }

    // r = b - Ax
    blas::copy(A.n, b, r);

  } else {
    unsigned inf;
    COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
    if (inf) return 0;
  }

  A.amux(D_MONE, x, r);

  if (rank==0) {
    beta = blas::nrm2(A.n, r);

    bool inf = (resid=beta/normb<=eps);
    for (unsigned p=1; p<nproc; ++p)
      COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

    if (inf) {
      eps = resid;
      nsteps = 0;
      delete [] r;
      return 0;
    }
  } else {
    unsigned inf;
    COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
    if (inf) return 0;
  }

  while (j<=nsteps) {

    if (rank==0) {
      blas::copy(A.n, r, V);                // v0 first orthonormal vector
      blas::scal(A.n, 1.0/beta, V);

      s[0] = beta;
      blas::setzero(m, s+1);
    }

    for (unsigned i=0; i<m && j<=nsteps; i++, j++) {

      // w = A M * v[i];
      if (rank==0) 
	blas::copy(A.n, V+i*A.n, xh);

      A.precond_apply(xh);

      if (rank==0)
	blas::setzero(A.n, V+(i+1)*A.n);

      A.amux(D_ONE, xh, V+(i+1)*A.n);

      if (rank==0) {
	for (unsigned k=0; k<=i; k++) {
	  H[k+i*(m+1)] = blas::scpr(A.n, V+(i+1)*A.n, V+k*A.n);
	  blas::axpy(A.n, -H[k+i*(m+1)], V+k*A.n, V+(i+1)*A.n);
	}

	H[i*(m+2)+1] = blas::nrm2(A.n, V+(i+1)*A.n);
	blas::scal(A.n, 1.0/H[i*(m+2)+1], V+(i+1)*A.n);

	// apply old Givens rotations to the last column in H
	for (unsigned k=0; k<i; k++)
	  applPlRot(H[k+i*(m+1)], H[k+1+i*(m+1)], cs[k], sn[k]);

	// generate new Givens rotation which eleminates H[i*(m+2)+1]
	genPlRot(H[i*(m+2)], H[i*(m+2)+1], cs[i], sn[i]);

	// apply it to H and s
	applPlRot(H[i*(m+2)], H[i*(m+2)+1], cs[i], sn[i]);
	applPlRot(s[i], s[i+1], cs[i], sn[i]);
      
	unsigned inf = ((resid=fabs(s[i+1]/normb)) < eps);
	for (unsigned p=1; p<nproc; ++p)
	  COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

	if (inf) {
	  update(A, i+1, H, m+1, s, V, x);
	  eps = resid;
	  nsteps = j;
	  delete [] r;
	  return 0;
	}
#ifndef NDEBUG
	std::cout << "Step " << j << ", resid=" << resid << std::endl;
#endif
      } else {
	unsigned inf;
	COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
	if (inf) return 0;
      }
    }

    if (rank==0) {
      update(A, m, H, m+1, s, V, x);

      // r = b - A x;
      blas::copy(A.n, b, r);
    }

    A.amux(D_MONE, x, r);

    if (rank==0) {
      beta = blas::nrm2(A.n, r);

      unsigned inf = ((resid=beta/normb) < eps);
      for (unsigned p=1; p<nproc; ++p)
        COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

      if (inf) {
	eps = resid;
	nsteps = j;
	delete [] r;
	return 0;
      }
    } else {
      unsigned inf;
      COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
      if (inf) return 0;
    }
  }

  if (rank==0) {
    eps = resid;
    delete [] r;
    return 1;
  }
  else return 1;
}


// ----------------------------------------------------------------------------
// complex version

static void genPlRot(dcomp a, dcomp b, double &cs, dcomp &sn)
{
  double k1, k2;
  if (Re(b)==0 && Im(b)==0) {
    cs = 1.0;
    sn = dcomp(0., 0.);
  } else if (Re(a)==0 && Im(a) == 0) {
    cs = 0.0;
    sn = dcomp(1., 0.);
  } else {
    if (Re(a)>=Im(a)) {
      k1 = Re(a) * Re(b) - Im(a) * Im(b);
      k1 /= abs2(a);
      k2 = (k1 * Im(a) - Im(b)) / Re(a);
    } else {
      k1 = Re(a) * Re(b) + Im(a) * Im(b);
      k1 /= abs2(a);
      k2 = (Re(b) - k1 * Re(a)) / Im(a);
    }
    cs = 1.0 / sqrt(1 + k1*k1 + k2*k2);
    sn = dcomp(cs * k1, cs * k2);
  }
}

static void applPlRot(dcomp &a, dcomp &b, double cs, dcomp sn)
{
  double ra, ia, rb, ib;
  ra = cs * Re(a) + Re(sn) * Re(b) - Im(sn) * Im(b);
  ia = cs * Im(a) + Re(sn) * Im(b) + Im(sn) * Re(b);
  rb = cs * Re(b) - Re(sn) * Re(a) - Im(sn) * Im(a);
  ib = cs * Im(b) - Re(sn) * Im(a) + Im(sn) * Re(a);

  a = dcomp(ra, ia);
  b = dcomp(rb, ib);
}

// H y = s, x += Vy
static void update(const Matrix<dcomp>& A, unsigned k, dcomp* H,
                   unsigned ldH, dcomp* s, dcomp* V, dcomp* x)
{
  dcomp *xh = new dcomp[A.n];
  dcomp *y = new dcomp[k];
  blas::copy(k, s, y);

  int inf;
  ztrtrs_(JOB_STR+5, JOB_STR, JOB_STR, &k, &N_ONE, H, &ldH, y, &k, &inf);
  assert(inf==0);


  // x += M V y
  blas::setzero(A.n, xh);
  blas::gemva(A.n, k, Z_ONE, V, y, xh);
  A.precond_apply(xh);
  blas::add(A.n, xh, x);

  delete [] y;
  delete [] xh;
}

unsigned GMRes_MPI(const Matrix<dcomp>& A, dcomp* const b, dcomp* const x,
                   double& eps, const unsigned m, unsigned& nsteps)
{
  unsigned rank  = COMM_AHMED.Get_rank(),
           nproc = COMM_AHMED.Get_size();

  double resid, normb, beta;
  unsigned j = 1;

  dcomp *r, *V, *H, *s, *sn, *xh;
  double *cs;

  if (rank==0) {
    r   = new dcomp[A.n];		// nx1, Residuum r
    V   = new dcomp[A.n * (m+1)];	// nx(m+1), orthonormale Matrix
    H   = new dcomp[(m+1)*m];	// (m+1)xm, Hessenbergmatrix
    s   = new dcomp[m+1];		// (m+1)x1, rechte Seite des GS Hy=s
    sn  = new dcomp[m+1];		// m+1, Sinus-Teil der Rotation
    cs = new double[m+1];		// m+1, Kosinus-Teil der Rotation
    xh = new dcomp[A.n];

    // normb = norm(b)
    normb = blas::nrm2(A.n, b);
    unsigned inf = (normb<D_PREC);
    for (unsigned p=1; p<nproc; ++p)
      COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

    if (inf) {
      blas::setzero(A.n, x);
      eps = 0.0;
      nsteps = 0;
      delete [] xh;
      delete [] cs;
      delete [] sn;
      delete [] s;
      delete [] H;
      delete [] V;
      delete [] r;
      return 0;
    }
    // r = b - Ax
    blas::copy(A.n, b, r);
  }
  else {
    unsigned inf;
    COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
    if (inf) return 0;
  }

  A.amux(Z_MONE, x, r);

  if (rank==0) {
    beta = blas::nrm2(A.n, r);
    unsigned inf = ((resid=beta/normb)<=eps);
    for (unsigned p=1; p<nproc; ++p)
      COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

    if (inf) {
      eps = resid;
      nsteps = 0;
      delete [] cs;
      delete [] sn;
      delete [] xh;
      delete [] s;
      delete [] H;
      delete [] V;
      delete [] r;
      return 0;
    }
  } else {
    unsigned inf;
    COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
    if (inf) return 0;
  }

  while (j<=nsteps) {

    if (rank==0) {
      blas::copy(A.n, r, V);                 // v0 first orthonormal vector
      blas::scal(A.n, 1.0/beta, V);

      s[0] = beta;
      blas::setzero(m, s+1);
    }

    for (unsigned i=0; i<m && j<=nsteps; i++, j++) {

	// w = A M * v[i];
      if (rank==0) 
	blas::copy(A.n, V+i*A.n, xh);

      A.precond_apply(xh);

      if (rank==0)
	blas::setzero(A.n, V+(i+1)*A.n);

      A.amux(Z_ONE, xh, V+(i+1)*A.n);

      if (rank==0) {
	unsigned k;

	for (k=0; k<=i; k++) {
	  H[k+i*(m+1)] = blas::scpr(A.n, V+k*A.n, V+(i+1)*A.n);
	  blas::axpy(A.n, Z_MONE * H[k+i*(m+1)], V+k*A.n, V+(i+1)*A.n);
	}

	H[i*(m+2)+1] = blas::nrm2(A.n, V+(i+1)*A.n);
	blas::scal(A.n, 1.0/H[i*(m+2)+1], V+(i+1)*A.n);

	// apply old Givens rotations to the last column in H
	for (k=0; k<i; k++)
	  applPlRot(H[k+i*(m+1)], H[k+1+i*(m+1)], cs[k], sn[k]);
	
	// generate new Givens rotation which eleminates H[i*(m+2)+1]
	genPlRot(H[i*(m+2)], H[i*(m+2)+1], cs[i], sn[i]);

	// apply it to H and s
	applPlRot(H[i*(m+2)], H[i*(m+2)+1], cs[i], sn[i]);
	applPlRot(s[i], s[i+1], cs[i], sn[i]);

	unsigned inf = ((resid=fabs(abs(s[i+1])/normb)) < eps);
	for (unsigned p=1; p<nproc; ++p)
	  COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

	if (inf) {
	  update(A, i+1, H, m+1, s, V, x);
	  eps = resid;
	  nsteps = j;
	  delete [] cs;
	  delete [] sn;
	  delete [] xh;
	  delete [] s;
	  delete [] H;
	  delete [] V;
	  delete [] r;
	  return 0;
	}
#ifndef NDEBUG
      std::cout << "Step " << j << ", resid=" << resid << std::endl;
#endif
      } else {
	unsigned inf;
	COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
	if (inf) return 0;
      }

    }


    if (rank==0) {
      update(A, m, H, m+1, s, V, x);

      // r = b - A * x);
      blas::copy(A.n, b, r);
    }

    A.amux((dcomp)-1.0, xh, r);

    if (rank==0) {
      beta = blas::nrm2(A.n, r);

      unsigned inf = ((resid=beta/normb) < eps);
      for (unsigned p=1; p<nproc; ++p)
        COMM_AHMED.Send(&inf, 1, MPI::UNSIGNED, p, 9);

      if (inf) {
	eps = resid;
	nsteps = j;
	delete [] r;
	return 0;
      }
    } else {
      unsigned inf;
      COMM_AHMED.Recv(&inf, 1, MPI::UNSIGNED, 0, 9);
      if (inf) return 0;
    }
  }

  if (rank==0) {
    eps = resid;
    delete [] cs;
    delete [] sn;
    delete [] xh;
    delete [] s;
    delete [] H;
    delete [] V;
    delete [] r;
    return 1;
  }
  else return 1;
}
