/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <iostream>
#include "blas.h"

bool generateDiagPrecond (const unsigned n, const double* const A, 
			  const unsigned* const jA, const unsigned* const iA,
                          double* diag)
{
  unsigned idx; // first idx of next row
  unsigned c; // column
  unsigned j;
  bool has_no_diag;

  for (unsigned r(0); r<n; ++r) {
    idx=iA[r+1];
    has_no_diag=true;
    for (j=iA[r]; j<idx && has_no_diag; ++j) {
      c=jA[j];
      if (c==r) {
        has_no_diag=false;
        diag[r] = 1.0/A[j];
      }
    }
    if (j==idx && has_no_diag) {
      std::cout << "row " << r << " has no diagonal element " << std::endl;
      return false;
    }
  }
  return true;
}

void scaleCRSSym(const unsigned n, const double* const A, const unsigned* const jA,
		 const unsigned* const iA, double*& diag, double* A_New)
{
  blas::copy(iA[n], A, A_New);
  diag = new double[n];
  generateDiagPrecond(n, A, jA, iA, diag);

  //blas::load(n, 1.0, diag);

  for(unsigned i=0;i<n;i++)
    for(unsigned k=iA[i];k<iA[i+1];k++)
      A_New[k] *= diag[jA[k]]*diag[i];
}
