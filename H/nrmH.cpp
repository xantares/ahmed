/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <cmath>
#include "blcluster.h"
#include "H.h"

//frobeniusnorm of H-matrix
template<class T> static
double nrmF2GeH_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    return A[idx]->nrmF2();
  } else {
    double nrm2 = 0.0;
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) nrm2 += nrmF2GeH_(son, A);
      }
    }
    return nrm2;
  }
}
//symmetric case:
template<class T> static
double nrmF2HeH_(blcluster* bl, mblock<T>** A)
{
   if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    return A[idx]->nrmF2();
  } else {
    double nrm2 = 0.0;
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      nrm2 += nrmF2HeH_(bl->getson(i, i), A);
      for (unsigned j=i+1; j<ns; ++j) 
        nrm2 += 2. * nrmF2GeH_(bl->getson(i,j), A);      
    }
    return nrm2;
  }
}
// Symmetric A
template<class T> static
double spctrlnrm_sym_(blcluster* bl, mblock<T>** A)
{
  const unsigned n = bl->getn1();
  T* const x = new T[n];
  T* const y = new T[n];
  blas::load(n, (T) 1.0, x);             // starting vector
  double d;

  for (unsigned l=0; l<10; ++l) {
    blas::setzero(n, y);
    mltaHeHVec((T) 1.0, bl, A, x, y);

    d = blas::nrm2(n, y);                    // normalize
    const double e = 1.0/d;
    blas::scal(n, (T) e, y);
    blas::copy(n, y, x);
  }

  delete [] y;
  delete [] x;

  return d;
}

template<class T> static
double spctrlnrm_(unsigned n, blcluster* blA, mblock<T>** A,
                  blcluster* blB, mblock<T>** B)
{
  T* const x = new T[n];
  T* const y = new T[n];
  blas::load(n, (T) 1.0, x);   // starting vector
  double d;

  for (unsigned l=0; l<10; ++l) {
    // compute (I-AB)^T (I-AB) x
    blas::setzero(n, y);
    mltaGeHVec((T) 1.0, blB, B, x, y);
    mltaGeHVec((T) -1.0, blA, A, y, x);
    blas::setzero(n, y);
    mltaGeHhVec((T) 1.0, blA, A, x, y);
    mltaGeHhVec((T) -1.0, blB, B, y, x);

    d = blas::nrm2(n, x);                    // normalize
    double e = 1.0/d;
    blas::scal(n, (T) e, x);
  }

  delete[] y;
  delete[] x;

  return sqrt(d);
}


// Symmetric A, B
template<class T> static
double spctrlnrm_sym_(unsigned n, blcluster* blA, mblock<T>** A,
                      blcluster* blB, mblock<T>** B)
{
  T* const x = new T[n];
  T* const y = new T[n];
  blas::load(n, (T) 1.0, x);             // starting vector
  double d;

  for (unsigned l=0; l<10; ++l) {
    // compute (I-AB)^T (I-AB) x
    blas::setzero(n, y);
    mltaHeHVec((T) 1.0, blB, B, x, y);
    mltaHeHVec((T) -1.0, blA, A, y, x);
    blas::setzero(n, y);
    mltaHeHVec((T) 1.0, blA, A, x, y);
    mltaHeHVec((T) -1.0, blB, B, y, x);

    d = blas::nrm2(n, x);                    // normalize
    const double e = 1.0/d;
    blas::scal(n, (T) e, x);
  }

  delete [] y;
  delete [] x;

  return sqrt(d);
}

///////////////////////////////////////////////////////////////////////////////
// Instanzen

double nrmF2GeH(blcluster* bl, mblock<double>** A)
{
  return nrmF2GeH_(bl, A);
}

double nrmF2GeH(blcluster* bl, mblock<float>** A)
{
  return nrmF2GeH_(bl, A);
}
double nrmF2GeH(blcluster* bl, mblock<dcomp>** A)
{
  return nrmF2GeH_(bl, A);
}

double nrmF2GeH(blcluster* bl, mblock<scomp>** A)
{
  return nrmF2GeH_(bl, A);
}

double nrmF2HeH(blcluster* bl, mblock<double>** A)
{
  return nrmF2HeH_(bl, A);
}

double nrmF2HeH(blcluster* bl, mblock<float>** A)
{
  return nrmF2HeH_(bl, A);
}
double nrmF2HeH(blcluster* bl, mblock<dcomp>** A)
{
  return nrmF2HeH_(bl, A);
}
double nrmF2HeH(blcluster* bl, mblock<scomp>** A)
{
  return nrmF2HeH_(bl, A);
}

double spctrlnrm_sym(blcluster* blA, mblock<double>** A)
{
  return spctrlnrm_sym_(blA, A);
}
double spctrlnrm_sym(blcluster* blA, mblock<float>** A)
{
  return spctrlnrm_sym_(blA, A);
}
double spctrlnrm_sym(blcluster* blA, mblock<dcomp>** A)
{
  return spctrlnrm_sym_(blA, A);
}
double spctrlnrm_sym(blcluster* blA, mblock<scomp>** A)
{
  return spctrlnrm_sym_(blA, A);
}

double spctrlnrm(unsigned n, blcluster* blA, mblock<double>** A,
                 blcluster* blB, mblock<double>** B)
{
  return spctrlnrm_(n, blA, A, blB, B);
}
double spctrlnrm(unsigned n, blcluster* blA, mblock<float>** A,
                 blcluster* blB, mblock<float>** B)
{
  return spctrlnrm_(n, blA, A, blB, B);
}
double spctrlnrm(unsigned n, blcluster* blA, mblock<dcomp>** A,
                 blcluster* blB, mblock<dcomp>** B)
{
  return spctrlnrm_(n, blA, A, blB, B);
}
double spctrlnrm(unsigned n, blcluster* blA, mblock<scomp>** A,
                 blcluster* blB, mblock<scomp>** B)
{
  return spctrlnrm_(n, blA, A, blB, B);
}

double spctrlnrm_sym(unsigned n, blcluster* blA, mblock<double>** A,
                     blcluster* blB, mblock<double>** B)
{
  return spctrlnrm_sym_(n, blA, A, blB, B);
}
double spctrlnrm_sym(unsigned n, blcluster* blA, mblock<float>** A,
                     blcluster* blB, mblock<float>** B)
{
  return spctrlnrm_sym_(n, blA, A, blB, B);
}
double spctrlnrm_sym(unsigned n, blcluster* blA, mblock<dcomp>** A,
                     blcluster* blB, mblock<dcomp>** B)
{
  return spctrlnrm_sym_(n, blA, A, blB, B);
}
double spctrlnrm_sym(unsigned n, blcluster* blA, mblock<scomp>** A,
                     blcluster* blB, mblock<scomp>** B)
{
  return spctrlnrm_sym_(n, blA, A, blB, B);
}



