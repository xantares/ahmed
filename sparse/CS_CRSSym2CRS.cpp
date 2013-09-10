/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "basmod.h"
// adds (symmetric) lower triangular part to upper triangular matrix

void CRSSym2CRS(unsigned n, unsigned* &iA, unsigned* &jA, double* &A)
{
  unsigned i;
  // count entries of each column
  unsigned* cnt = new unsigned[n];
  for (i=0; i<n; ++i) cnt[i] = 0;

  for (i=0; i<n; ++i) {
    unsigned j(iA[i]), idx(iA[i+1]);
    while (j < idx) {
      if (jA[j] != i) cnt[jA[j]]++;
      j++;
    }
  }

  // summing up entries
  for (i=2; i<n; ++i) cnt[i] += cnt[i-1];

  unsigned* iAn (new unsigned[n+1]);
  iAn[0] = 0;
  for (i=1; i<=n; ++i) iAn[i] = iA[i] + cnt[i-1];

  double *An (new double [iAn[n]]);
  unsigned *jAn (new unsigned [iAn[n]]);
  for (unsigned k=0; k<n; k++) cnt[k] = iAn[k];

  for (i=0; i<n; ++i) {
    unsigned j(iA[i]), idx(iA[i+1]);
    while (j < idx) {
      if (i == jA[j]) {
        An[cnt[i]] = A[j];
        jAn[cnt[i]++] = jA[j];
      } else {
        An[cnt[i]] = A[j];
        An[cnt[jA[j]]] = A[j];
        jAn[cnt[i]++] = jA[j];
        jAn[cnt[jA[j]]++] = i;
      }
      j++;
    }
  }

  swap(iA, iAn);
  swap(jA, jAn);
  swap(A, An);

  delete [] jAn;
  delete [] iAn;
  delete [] An;
  delete [] cnt;
}
