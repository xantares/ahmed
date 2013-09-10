/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


// remove lower triangular part from a CRS matrix
void CRS2CRSSym(unsigned n, double* A, unsigned* jA, unsigned* iA,
                double* &A_new, unsigned* &jA_new, unsigned* &iA_new)
{
  unsigned i, j;
  unsigned nnz = 0;

  // count number of non-zeros in the upper triangular part
  for (i=0; i<n; i++) {
    unsigned idx = iA[i+1];
    for (j=iA[i]; j<idx; j++) if (jA[j] >= i) ++nnz;
  }

  A_new = new double[nnz];
  jA_new = new unsigned[nnz];
  iA_new = new unsigned[n+1];

  iA[0] = nnz = 0;

  for (i=0; i<n; i++) {
    unsigned idx = iA[i+1];
    for (j=iA[i]; j<idx; j++) {
      if (jA[j] >= i) {
        A_new[nnz] = A[j];
        jA_new[nnz++] = jA[j];
      }
    }
    iA_new[i+1] = nnz;
  }
}
