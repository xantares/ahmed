/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "cmplx.h"

/*
  -----------------------------
  A times a vector
  -----------------------------
  multiplies a matrix by a vector using the dot product form
  Matrix A is stored in compressed sparse row storage.

  on entry:
  ----------
  n     = row dimension of A
  x     = real array of length equal to the column dimension of
          the A matrix.
  A, jA,
  iA    = input matrix in compressed sparse row storage format.

  on return:
  -----------
  y     = real array of length n, containing the product y += d * Ax

*/

template<class T> static
void amuxCRS_(const unsigned n, const T d, const T* const x,
              T* const y, const unsigned* const iA, const unsigned* const jA,
              const T* const A)
{
  for (unsigned i=0; i<n; ++i) {

    // compute the inner product of row i with vector x
    T t = 0.0;
    for (unsigned k=iA[i]; k<iA[i+1]; ++k)
      t += A[k] * x[jA[k]];

    // store result in y(i)
    y[i] += d*t;
  }
}

void amuxCRS(unsigned n, double d, double* x, double* y,
             unsigned* iA, unsigned* jA, double* A)
{
  amuxCRS_(n, d, x, y, iA, jA, A);
}



void amuxCRS(unsigned n, float d, float* x, float* y,
             unsigned* iA, unsigned* jA, float* A)
{
  amuxCRS_(n, d, x, y, iA, jA, A);
}



void amuxCRS(unsigned n, dcomp d, dcomp* x, dcomp* y,
             unsigned* iA, unsigned* jA, dcomp* A)
{
  amuxCRS_(n, d, x, y, iA, jA, A);
}

void amuxCRS(unsigned n, scomp d, scomp* x, scomp* y,
             unsigned* iA, unsigned* jA, scomp* A)
{
  amuxCRS_(n, d, x, y, iA, jA, A);
}


