/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*
  -----------------------------
  A times a vector
  -----------------------------
  multiplies a symmetric matrix by a vector
  Matrix A is stored in compressed sparse row storage.

  on entry:
  ----------
  n     = row dimension of A
  x     = real array of length equal to the column dimension of
          the A matrix.
  A, jA,
  iA    = input matrix in compressed sparse row format.

  on return:
  -----------
  y     = real array of length n, containing the product y += d * Ax

*/

void amuxSymmCRS(unsigned n, double d, double* x, double* y,
                 unsigned *iA, unsigned* jA, double* A)
{
  for (unsigned i=0; i<n; ++i) {
    unsigned idx = iA[i+1];
    for (unsigned k=iA[i]; k<idx; ++k) {
      unsigned j = jA[k];
      if (i!=j) {
        y[i] += d * A[k] * x[j];
        y[j] += d * A[k] * x[i];
      } else
        y[i] += d * x[i] * A[k];
    }
  }
}

// assumes that a_{ij}=1 and a_{ii}=0
void amuxSymmCRS(unsigned n, double d, double* x, double* y,
                 unsigned *iA, unsigned* jA)
{
  for (unsigned i=0; i<n; ++i) {
    unsigned idx = iA[i+1];
    for (unsigned k=iA[i]; k<idx; ++k) {
      unsigned j = jA[k];
      y[i] += d*x[j];
      y[j] += d*x[i];
    }
  }
}
