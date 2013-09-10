/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "cmplx.h"

/*
  --------------------------------
  transp( A ) times a vector
  --------------------------------
  multiplies the transpose of a matrix by a vector when the original
  matrix is stored in compressed sparse row storage. Can also be
  viewed as the product of a matrix by a vector when the original
  matrix is stored in the compressed sparse column format.
  --------------------------------

  on entry:
  ----------
  n     = row dimension of A
  x     = real array of length equal to the column dimension of
             the A matrix.
  A, jA,
  iA    = input matrix in compressed sparse row format.

  on return:
  -----------
  y     = real array of length n, containing the product y += d*transp(A)*x

*/

template<class T> static
void atmuxCRS_(const unsigned n, const T d, T* x, T* y,
               unsigned *iA, unsigned *jA, T* A)
{
  for (unsigned i=0; i<n; ++i)
    for (unsigned k=iA[i]; k<iA[i+1]; ++k)
      y[jA[k]] += d*x[i]*A[k];
}

void atmuxCRS(const unsigned n, const double d, double* x, double* y,
              unsigned* iA, unsigned *jA, double* A)
{
  atmuxCRS_(n, d, x, y, iA, jA, A);
}

void atmuxCRS(const unsigned n, const float d, float* x, float* y,
              unsigned* iA, unsigned *jA, float* A)
{
  atmuxCRS_(n, d, x, y, iA, jA, A);
}

void atmuxCRS(const unsigned n, const dcomp d, dcomp* x, dcomp* y,
              unsigned* iA, unsigned *jA, dcomp* A)
{
  atmuxCRS_(n, d, x, y, iA, jA, A);
}

void atmuxCRS(const unsigned n, const scomp d, scomp* x, scomp* y,
              unsigned* iA, unsigned *jA, scomp* A)
{
  atmuxCRS_(n, d, x, y, iA, jA, A);
}
