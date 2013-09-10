/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "mblock.h"
#include "blcluster.h"
#include <cmath>

template<class T> static
void Hrownrms2_(blcluster* bl, mblock<T>** A, T* rownrm)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->rownrms2(rownrm+bl->getb1());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) Hrownrms2_(son, A, rownrm);
      }
  }
}

template<class T> static
void Hcolnrms2_(blcluster* bl, mblock<T>** A, T* colnrm)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->colnrms2(colnrm+bl->getb2());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) Hcolnrms2_(son, A, colnrm);
      }
  }
}

template<class T> static
void Hnrms2_sym(blcluster* bl, mblock<T>** A, T* nrm)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->symnrms2(nrm+bl->getb1(), nrm+bl->getb2());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) Hnrms2_sym(son, A, nrm);
      }
  }
}

void Hcolnrms(unsigned n, blcluster* bl, mblock<double>** A, double* D)
{
  blas::setzero(n, D);
  Hcolnrms2_(bl, A, D);
  for (unsigned i=0; i<n; ++i) D[i] = sqrt(D[i]);
}

void Hrownrms(unsigned n, blcluster* bl, mblock<double>** A, double* D)
{
  blas::setzero(n, D);
  Hrownrms2_(bl, A, D);
  for (unsigned i=0; i<n; ++i) D[i] = sqrt(D[i]);
}

// compute the column/row norm of the symmetric matrix
void Hnrms_sym(unsigned n, blcluster* bl, mblock<double>** A, double* D)
{
  blas::setzero(n, D);
  Hnrms2_sym(bl, A, D);
  for (unsigned i=0; i<n; ++i) D[i] = sqrt(D[i]);
}

template<class T> static
void Hdiag_sym(blcluster* bl, mblock<T>** A, T* diag)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->diag_sym(diag+bl->getb1());
  } else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) Hdiag_sym(bl->getson(i, i), A, diag);
  }
}

template<class T> static
void Hscale_rows_(blcluster* bl, mblock<T>** A, T* D)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->scale_rows(D+bl->getb1());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) Hscale_rows_(son, A, D);
      }
  }
}

template<class T> static
void Hscale_cols_(blcluster* bl, mblock<T>** A, T* D)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->scale_cols(D+bl->getb2());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) Hscale_cols_(son, A, D);
      }
  }
}

template<class T> static
void Hscale_sym(blcluster* bl, mblock<T>** A, T* D)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    A[idx]->scale_sym(D+bl->getb1(), D+bl->getb2());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) Hscale_sym(son, A, D);
      }
  }
}

// A -> DA, where D_i is the inverse Euclidean norm of the i-th row of A
void Hequilib_rows(unsigned n, blcluster* bl, mblock<double>** A, double* D)
{
  blas::setzero(n, D);
  Hrownrms2_(bl, A, D);
  for (unsigned i=0; i<n; ++i) {
    double e = D[i];
    if (e!=0.0) D[i] = 1.0/sqrt(e);
  }

  Hscale_rows_(bl, A, D);
}


// A -> AD, where D_i is the inverse Euclidean norm of the i-th column of A
void Hequilib_cols(unsigned n, blcluster* bl, mblock<double>** A, double* D)
{
  blas::setzero(n, D);
  Hcolnrms2_(bl, A, D);
  for (unsigned i=0; i<n; ++i) {
    double e = D[i];
    if (e!=0.0) D[i] = 1.0/sqrt(e);
  }

  Hscale_cols_(bl, A, D);
}

// A -> DAD, where D_i is the square root of the inverse i-th diag entry
// of the symmetric matrix A
void Hequilib_sym(unsigned n, blcluster* bl, mblock<double>** A, double* D)
{
  blas::setzero(n, D);
  Hdiag_sym(bl, A, D);

  for (unsigned i=0; i<n; ++i) {
    double e = abs(D[i]);
    if (e!=0.0) D[i] = 1.0/sqrt(e);
  }

  Hscale_sym(bl, A, D);
}


