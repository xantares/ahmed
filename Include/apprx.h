/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef APPRX_H
#define APPRX_H

//#define CHECK_ACA_ERROR

#include <fstream>
#include "cmplx.h"
#include "blas.h"
#include "bemblcluster.h"
#include "ACA.h"

#ifdef CHECK_ACA_ERROR

double ACA_error_max = 0.0;

template<class T,class T1,class T2, class MATGEN_T> static
void check_error(MATGEN_T& MatGen, bemblcluster<T1,T2>* bl, double eps,
		 mblock<T>* mbl, unsigned i0)
{
  std::cout << "checking ACA error ... " << std::flush;
  unsigned n1 = bl->getn1(), n2 = bl->getn2();
  T *A = new T[n1*n2];
  mbl->convLrM_toGeM(A, n1);
  T *B = new T[n1*n2];

  MatGen.cmpbl(bl->getb1(), n1, bl->getb2(), n2, B);

  double dnrm = 0.0, dmax = 0.0;
  for (unsigned j=0; j<n2; ++j) {
    for (unsigned i=0; i<n1; ++i) {
      dnrm += abs2(B[i+j*n1]);
      dmax += abs2(A[i+j*n1]-B[i+j*n1]);
    }
  }

  if ((dnrm>0.0 && sqrt(dmax/dnrm)>20*eps) || (dnrm==0.0 && dmax>eps)) {
    double maxerr = sqrt(dmax/dnrm);
    std::cout << "eps=" << eps << ", maxerr=" << maxerr << std::endl;

    if (ACA_error_max<maxerr) ACA_error_max = maxerr;
    /*
    std::ofstream os("out.dat");
    os << n1 << ' ' << n2 << ' ' << i0 << std::endl;
    os.precision(16);
    for (unsigned i=0; i<n1; ++i) {
      for (unsigned j=0; j<n2; ++j) os << B[i+j*n1] << ' ';
      os << std::endl;
    }
    os.close();
    // check whether it can be approximated using SVD
    double* S = new double[MIN(n1, n2)];
    unsigned nwk = 5*(n1+n2);
    T* wk = new T[nwk];
    blas::svals(n1, n2, B, S, nwk, wk);
    unsigned kt = MIN(n1, n2);
    while (S[kt-1]<S[0]*eps) --kt;
    std::cout << "matrix is " << n1 << 'x' << n2
        << ", required rank " << kt << std::endl;
    for (unsigned i=0; i<MIN(n1, n2); ++i)
      std::cout << "sval[" << i << "]="<< S[i] << std::endl;
    delete [] wk;
    delete [] S;
    
    exit(1);
    */
  } else std::cout << "ok." << std::endl;

  delete [] B;
  delete [] A;
}

#endif


template<class T,class T1, class MATGEN_T>
void apprx_sym(MATGEN_T& MatGen, mblock<T>* &mbl, bemblcluster<T1,T1>* bl,
	       double eps, unsigned rankmax, const bool& cmplx_sym)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2();
  unsigned b1 = bl->getb1(), b2 = bl->getb2();

  // admit a slightly higher rank
  unsigned maxk = n1*n2/(n1+n2)+4;
  if (maxk>rankmax) maxk = rankmax;

  mbl = new mblock<T>(n1, n2);
  assert(mbl!=NULL);

  bool succ = false;
  
  if (bl->isadm() && maxk) {
    unsigned k;

    T *U, *V;

    unsigned i0 = bl->getcl1()->geticom() - b1;

    succ = ::ACAr(MatGen, b1, n1, b2, n2, eps, maxk, i0, k, U, V);

    if (succ) {
      mbl->cpyLrM_cmpr(k, U, n1, V, n2, eps, k);      
#ifdef CHECK_ACA_ERROR
      check_error(MatGen, bl, eps, mbl, i0);
#endif
    }
    delete [] V;
    delete [] U;
  }

  if (!succ) {
    if (bl->isdbl()) {
      if(cmplx_sym) mbl->setSyM();
      else mbl->setHeM();
      
      MatGen.cmpblsym(b1, n1, mbl->getdata());
    } else {
      mbl->setGeM();
      MatGen.cmpbl(b1, n1, b2, n2, mbl->getdata());
    }
  }

}

// matrix generation procedure for unsymmetric matrices
template<class T,class T1,class T2, class MATGEN_T>
void apprx_unsym(MATGEN_T& MatGen, mblock<T>* &mbl, bemblcluster<T1,T2>* bl,
                 double eps, unsigned rankmax)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2();
  unsigned b1 = bl->getb1(), b2 = bl->getb2();

  unsigned maxk = n1*n2/(n1+n2)+4;
  if (maxk>rankmax) maxk = rankmax;

  mbl = new mblock<T>(n1, n2);
  assert(mbl!=NULL);

  bool succ = false;

  if (bl->isadm() && maxk) {
    unsigned k;
    T *U, *V;

    unsigned i0 = bl->getcl1()->geticom() - b1;
    succ = ::ACA(MatGen, b1, n1, b2, n2, eps, maxk, i0, k, U, V);

    if (succ) {
      mbl->cpyLrM_cmpr(k, U, n1, V, n2, eps, k);
#ifdef CHECK_ACA_ERROR
      check_error(MatGen, bl, eps, mbl, i0);
#endif
    }
    delete [] V;
    delete [] U;
  }

  if (!succ) {
    mbl->setGeM();
    MatGen.cmpbl(b1, n1, b2, n2, mbl->getdata());
  }
}

// matrix generation procedure for unsymmetric matrices
template<class T,class T1,class T2, class MATGEN_T>
void apprxSVD_unsym(MATGEN_T& MatGen, mblock<T>* &mbl, bemblcluster<T1,T2>* bl,
                    double eps, unsigned rankmax)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2();
  unsigned b1 = bl->getb1(), b2 = bl->getb2();
  unsigned maxk = n1*n2/(n1+n2);
  if (maxk>rankmax) maxk = rankmax;

  mbl = new mblock<T>(n1, n2);
  assert(mbl!=NULL);

  bool succ = false;

  if (bl->isadm() && maxk) {
    unsigned min12 = MIN(n1, n2);
    double* S = new double[min12];
    T *U = new T[n1*n2], *VT = new T[n1*n2];
    MatGen(b1, n1, b2, n2, U);

    unsigned nwk = 5*(n1+n2);
    T* wk = new T[nwk];
    blas::gesvd(n1, n2, U, S, U, n1, VT, MIN(n1,n2), nwk, wk);
    unsigned kt = min12;
    while (kt>0 && S[kt-1]<S[0]*eps) --kt;

    if (kt>0 && kt*(n1+n2)<n1*n2) {
      succ = true;
      mbl->setrank(kt);
      blas::copy(kt*n1, U, mbl->getdata());
      blas::transpose(kt, n2, VT, mbl->getdata()+kt*n1);
    } else
      succ = false;

    delete [] wk;
    delete [] VT;
    delete [] U;
    delete [] S;
  }

  if (!succ) {
    mbl->setGeM();
    MatGen(b1, n1, b2, n2, mbl->getdata());
  }
}

#endif
