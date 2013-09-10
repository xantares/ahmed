/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef VECTORITER_H
#define VECTORITER_H

#include "blas.h"

template<class T>
double vectorIter(const unsigned counter, const unsigned m, const T* A)
{
  T* x = new T[m];
  T* y = new T[m];
  double nrm2;

  blas::load(m, (T) 1.0, y);
  for (unsigned i=0; i<counter; i++) {
    nrm2 = 1.0/blas::nrm2(m,y);
    blas::scal(m, nrm2, y);
    swap(x, y);
    blas::gemv(m, m, (T)1.0, A, x, y);
  }

  nrm2 = abs(blas::scpr(m, x, y));

  delete [] x;
  delete [] y;
  return nrm2;
}

#endif
