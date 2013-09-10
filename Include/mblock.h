/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file  cluster.h
  \brief Include file for the class mblock
  \author M. Bebendorf
*/

#ifndef MBLOCK_H
#define MBLOCK_H

#include <assert.h>
#include <iostream>
#include <fstream>
#include "blas.h"
#include "cmplx.h"
#include "preserveVec.h"
#include "cluster.h"
#include "basmod.h"

template<class T> class mblock
{
  friend class blcluster;
  friend void copymbl_(mblock<double>*, mblock<float>*);
  friend void copymbl_(mblock<dcomp>*, mblock<scomp>*);

protected:
  T *data;                // if low-rank-repr. UV^H (twice MAX_RANK columns)
  unsigned n1, n2;        // n1 number of rows, n2 number of columns
  unsigned bl_rank;       // the rank of this block

  struct st_info {
    unsigned is_LrM : 1;        // low rank representation ? else dense
    unsigned is_HeM : 1;        // for dense matrices: is hermitian ?
    unsigned is_SyM : 1;        // for dense matrices: is symmetric ?
    unsigned is_LtM : 1;        // for dense matrices: is lower triangular ?
    unsigned is_UtM : 1;        // for dense matrices: is upper triangular ?
  } info;

  // a low-rank matrix U V^H is stored columnwise : (U,V)
  // a dense matrix is stored column by column
  // a dense symmetric/hermitian matrix is stored as an upper triangular matrix
  // an upper triangular matrix is stored columnwise (a11,a12,a22,...)
  // a lower triangular matrix is stored columnwise (a11,a21,a31,...)
  //         on the diagonal of it is a permutation P of the rows of L defined
  //         by (Pl)_{L_{ii}} = l_i

  // add low-rank matrix to this low-rank matrix with truncation
  void addtrll(unsigned k, T* U, unsigned ldU, T* V, unsigned ldV,
               double delta, unsigned kgoal, contLowLevel<T>* haarInfo=NULL,
               T* X=NULL, unsigned ldX=0, T* Y_=NULL, unsigned ldY_=0);

  // append low-rank matrix to this low-rank matrix
  void append(unsigned k, T* U, unsigned ldU, T* V, unsigned ldV);

  // add low-rank matrix to this low-rank matrix with truncation
  // returns remainder
  void addtrll_rmnd(unsigned, T*, unsigned, T*, unsigned,
                    double, unsigned, unsigned&, T*&, T*&);

  // compute singular values of this low-rank matrix
  void get_svals_LrM(double*) const;
  void get_svals_LrM(float*) const;

  // unify two low-rank matrices
  void unify_cols_LrMLrM(double, unsigned, const mblock&, const mblock&,
			 contLowLevel<T>* haarInfo=NULL, T* X=NULL, unsigned ldX=0,
			 T* Y_=NULL, unsigned ldY_=0);
  void unify_rows_LrMLrM(double, unsigned, const mblock&, const mblock&,
			 contLowLevel<T>* haarInfo=NULL, T* X=NULL, unsigned ldX=0,
			 T* Y_=NULL, unsigned ldY_=0);

  // *************************************************************************
  // PostScript output for low-rank and dense matrices
  //
  void psout_LrM(std::ostream&, unsigned, unsigned, unsigned, bool,
                 unsigned level=NULL, cluster* cl=NULL,
                 void (*f)(cluster* cl, const unsigned lvl,
                           double* X, const double* const basis)=NULL) const;

  void psout_GeM(std::ostream&, unsigned, unsigned, unsigned, bool) const;

  // *************************************************************************
  // multiplication routines
  //




  // multiply dense by vector: y += d A x  (A is a dense matrix)
  void mltaGeMVec(T d, T* x, T* y) const {
    assert(isGeM());
    blas::gemva(n1, n2, d, data, x, y);
  }

  // multiply dense by dense: Y += d A X  (A is a dense matrix)
  void mltaGeMGeM(T d, unsigned p, T* X, unsigned ldX,
                     T* Y, unsigned ldY) const {
    assert(isGeM());
    blas::gemma(n1, n2, p, d, data, n1, X, ldX, Y, ldY);
  }

  // multiply vector by herm transposed dense: y += A^H x  (A is a dense matrix)
  void mltaGeMhVec(T d, T* x, T* y) const {
    assert(isGeM());
    blas::gemhva(n1, n2, d, data, x, y);
  }

  // multiply dense by herm. transposed dense: Y += d A^H X  (A is a dense matrix)
  void mltaGeMhGeM(T d, unsigned p, T* X, unsigned ldX,
                      T* Y, unsigned ldY) const {
    assert(isGeM());
    blas::gemhma(n1, n2, p, d, data, n1, X, ldX, Y, ldY);
  }

  // multiply vector by transposed dense: y += A^T x  (A is a dense matrix)
  void mltaGeMtVec(T d, T* x, T* y) const {
    assert(isGeM());
    blas::gemtva(n1, n2, d, data, x, y);
  }

  // multiply dense by transposed dense: Y += d A^T X  (A is a dense matrix)
  void mltaGeMtGeM(T d, unsigned p, T* X, unsigned ldX,
                      T* Y, unsigned ldY) const {
    assert(isGeM());
    blas::gemtma(n1, n2, p, d, data, n1, X, ldX, Y, ldY);
  }

  // multiply hermitian packed by vector: y += d A x (A is herm. dense)
  void mltaHeMVec(T d, T* x, T* y) const {
    assert(isGeM() && isHeM() && n1==n2);
    blas::hemva(n1, d, data, x, y);
  }  

  // multiply hermitian packed by dense: Y += d A X (A is herm. dense)
  void mltaHeMGeM(T d, unsigned p, T* X, unsigned ldX,
                        T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaHeMVec(d, X+ldX*l, Y+ldY*l);
  }

  // multiply sym packed by vector: y += d A x (A is sym. dense)
  void mltaSyMVec(T d, T* x, T* y) const {
    assert(isGeM() && isSyM() && n1==n2);
    blas::symva(n1, d, data, x, y);
  }

  // multiply sym packed by dense: Y += d A X (A is sym. dense)
  void mltaSyMGeM(T d, unsigned p, T* X, unsigned ldX,
                        T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaSyMVec(d, X+ldX*l, Y+ldY*l);
  }

  // multiply transposed of herm. packed by vector: y += d A^T x (A is herm. dense)
  void mltaHeMtVec(T d, T* x, T* y) const {
    assert(isGeM() && isHeM() && n1==n2);
    blas::conj(n1, x);
    blas::conj(n1, y);
    mltaHeMVec(conj(d), x, y);
    blas::conj(n1, x);
    blas::conj(n1, y);
  }  

  // multiply transposed of herm packed by dense: Y += d A^T X (A is herm. dense)
  void mltaHeMtGeM(T d, unsigned p, T* X, unsigned ldX,
                        T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaHeMtVec(d, X+ldX*l, Y+ldY*l);
  }
  
  // multiply herm transp of sym. packed by vector: y += d A^H x (A is sym. dense)
  void mltaSyMhVec(T d, T* x, T* y) const {
    assert(isGeM() && isSyM() && n1==n2);
    blas::conj(n1, x);
    blas::conj(n1, y);
    mltaSyMVec(conj(d), x, y);
    blas::conj(n1, x);
    blas::conj(n1, y);
  }

  // multiply herm trans of sym packed by dense: y += d A^H x (A is sym. dense)
  void mltaSyMhGeM(T d, unsigned p, T* X, unsigned ldX,
		   T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaSyMhVec(d, X+ldX*l, Y+ldY*l);
  }
  
  

  // multiply lower triangular by vector: y += d PL x
  void mltaLtMVec(T d, T* x, T* y) const;

  // multiply lower triangular by dense: Y += d PL X
  void mltaLtMGeM(T d, unsigned p, T* X, unsigned ldX,
                        T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaLtMVec(d, X+ldX*l, Y+ldY*l);
  }

  // multiply upper triangular by vector: y += d A x (A is utr)
  void mltaUtMVec(T d, T* x, T* y) const;

  // multiply upper triangular by dense: Y += d A X (A is utr)
  void mltaUtMGeM(T d, unsigned p, T* X, unsigned ldX,
                        T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaUtMVec(d, X+ldX*l, Y+ldY*l);
  }

  // multiply herm transp of lower triangular by vector: y += d (PL)^H x = d L^H P^{-1} x
  void mltaLtMhVec(T d, T* x, T* y) const;

  // multiply herm transp of lower triangular by dense: Y += d (PL)^H X = d L^H P^{-1} X
  void mltaLtMhGeM(T d, unsigned p, T* X, unsigned ldX,
		   T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaLtMhVec(d, X+ldX*l, Y+ldY*l);
  }

  // multiply herm transp of upper triangular by vector: y += d A^H x (A is utr)
  void mltaUtMhVec(T d, T* x, T* y) const;

  // multiply herm transposed of upper triang by dense: Y += d A^H X (A is utr)
  void mltaUtMhGeM(T d, unsigned p, T* X, unsigned ldX,
		   T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaUtMhVec(d, X+ldX*l, Y+ldY*l);
  }

  // multiply transposed of lower triang by vector: y += d (PL)^T x = d L^T P^{-1} x
  void mltaLtMtVec(T d, T* x, T* y) const {
    assert(isGeM() && isLtM() && n1==n2);
    blas::conj(n1, x);
    blas::conj(n1, y);
    mltaLtMhVec(conj(d), x, y);
    blas::conj(n1, x);
    blas::conj(n1, y);
  }

  // multiply transposed of  lower triang by dense: Y += d (PL)^T X = d L^T P^{-1} X
  void mltaLtMtGeM(T d, unsigned p, T* X, unsigned ldX,
		   T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaLtMtVec(d, X+ldX*l, Y+ldY*l);
  }

  
  // multiply transposed of upper triangular by vector: y += d A^T x (A is utr)
  void mltaUtMtVec(T d, T* x, T* y) const {
    assert(isGeM() && isUtM() && n1==n2);
    blas::conj(n1, x);
    blas::conj(n1, y);
    mltaUtMhVec(conj(d), x, y);
    blas::conj(n1, x);
    blas::conj(n1, y);
  }

  // multiply transposed of upper triangular by dense: Y += d A^T X (A is utr)
  void mltaUtMtGeM(T d, unsigned p, T* X, unsigned ldX,
		   T* Y, unsigned ldY) const {
    for (unsigned l=0; l<p; ++l) mltaUtMtVec(d, X+ldX*l, Y+ldY*l);
  }

  


  // *************************************************************************
  // PUBLIC routines
  //

public:

  // Konstruktor
  mblock(unsigned m, unsigned n) : n1(m), n2(n) {
    bl_rank = info.is_HeM = info.is_SyM = info.is_LtM = info.is_UtM = 0;
    info.is_LrM = 1;
    data = NULL;
  }

  // Destruktor
  ~mblock() {
    delete[] data;
  }

  // number of values in data
  unsigned long nvals() const {
    if (isGeM()) {
      if (isHeM() || isSyM() || isLtM() || isUtM()) return n1*(n1+1)/2;
      else return n1*n2;
    } else return (n1+n2)*bl_rank;
  };

  // return size of this mblock
  unsigned long size() const {
    return sizeof(T)*nvals() + sizeof(mblock<T>);
  };

  // *************************************************************************
  // set/read attributes

  unsigned getn1() {
    return n1;
  }
  unsigned getn2() {
    return n2;
  }

  bool isLrM() const {
    return info.is_LrM;
  }
  bool isGeM() const {
    return !isLrM();
  }
  bool isHeM() const {
    assert(isGeM());
    return info.is_HeM;
  }
  bool isSyM() const {
    assert(isGeM());
    return info.is_SyM;
  }
  bool isLtM() const {
    assert(isGeM());
    return info.is_LtM;
  }
  bool isUtM() const {
    assert(isGeM());
    return info.is_UtM;
  }

  unsigned rank() const {
    assert(isLrM());
    return bl_rank;
  }

  void freedata() {
    delete[] data;
    data = NULL;
    info.is_LrM = 1;
    bl_rank = 0;
  }

  T* getdata() const {
    return data;
  }

  void setrank(const unsigned k) {
    delete[] data;
    bl_rank = k;
    info.is_LrM = 1;
    if (k) {
      data = new T[k*(n1+n2)];
      assert(data!=NULL);
    } else data = NULL;
  }

  void setGeM() {
    delete[] data;
    info.is_LrM = info.is_HeM = info.is_SyM = info.is_LtM = info.is_UtM = 0;
    data = new T[n1*n2];
    assert(data!=NULL);
  }

  void setHeM() {
    assert(n1==n2);
    delete[] data;
    info.is_LrM = info.is_SyM = info.is_UtM = info.is_LtM = 0;
    info.is_HeM = 1;
    data = new T[n1*(n1+1)/2];
    assert(data!=NULL);
  }

  void setSyM() {
    assert(n1==n2);
    delete[] data;
    info.is_LrM = info.is_HeM = info.is_UtM = info.is_LtM = 0;
    info.is_SyM = 1;
    data = new T[n1*(n1+1)/2];
    assert(data!=NULL);
  }

  void setLtM() {
    assert(n1==n2);
    delete[] data;
    info.is_LrM = info.is_HeM = info.is_SyM = info.is_UtM = 0;
    info.is_LtM = 1;
    data = new T[n1*(n1+1)/2];
    assert(data!=NULL);
  }

  void setUtM() {
    assert(n1==n2);
    delete[] data;
    info.is_LrM = info.is_HeM = info.is_SyM = info.is_LtM = 0;
    info.is_UtM = 1;
    data = new T[n1*(n1+1)/2];
    assert(data!=NULL);
  }

  void get_prop(unsigned* inf) {
    if ((inf[1] = isLrM())) inf[0] = bl_rank;
    else {
      inf[2] = isHeM();
      inf[3] = isLtM();
      inf[4] = isUtM();
    }
  }


  // *************************************************************************
  // Initialization routines
  //

  // init with 0

  void init0_GeM(unsigned m, unsigned n) {
    n1 = m;
    n2 = n;
    setGeM();
    blas::setzero(n1*n2, data);
  }

  void init0_HeM(unsigned n) {
    n1 = n2 = n;
    setHeM();
    blas::setzero(n1*(n1+1)/2, data);
  }

  void init0_SyM(unsigned n) {
    n1 = n2 = n;
    setSyM();
    blas::setzero(n1*(n1+1)/2, data);
  }

  void init0_UtM(unsigned n) {
    n1 = n2 = n;
    setUtM();
    blas::setzero(n1*(n1+1)/2, data);
  }

  void init0_LtM(unsigned n) {
    n1 = n2 = n;
    setLtM();
    blas::fill0_ltr(n, data);
  }

  // init with identity
  void initId_GeM(unsigned n) {
    n1 = n2 = n;
    setGeM();
    blas::fillId(n, data);
  }

  void initId_HeM(unsigned n) {
    n1 = n2 = n;
    setHeM();
    blas::fillId_utr(n, data);
  }

  void initId_SyM(unsigned n) {
    n1 = n2 = n;
    setSyM();
    blas::fillId_utr(n, data);
  }

  void initId_LtM(unsigned n) {
    n1 = n2 = n;
    setLtM();
    blas::fillId_ltr(n, data);
  }

  void initId_UtM(unsigned n) {
    n1 = n2 = n;
    setUtM();
    blas::fillId_utr(n, data);
  }

  // *************************************************************************
  // Addition routines
  //

  // add dns to dns
  void addGeM_toGeM(T* const A, const unsigned ldA) {
    assert(ldA>=n1 && isGeM());
    for (unsigned j=0; j<n2; ++j) blas::add(n1, A+j*ldA, data+j*n1);
  }

  void addHeM_toHeM(T* const A) {
    assert(isGeM() && (isHeM() || isSyM()) );
    blas::axpy(n1*(n1+1)/2, (T) 1.0, A, data);
  }

  void addGeM_toHeM(T* const A, const unsigned ldA);

  // add lwr to dns
  void addLrM_toGeM(unsigned k, T* U, unsigned ldU,
                      T* V, unsigned ldV) {
    assert(isGeM());
    if (k>0) blas::gemmha(n1, k, n2, U, ldU, V, ldV, data, n1);
  }

  // add lwr to herm. dns
  void addLrM_toHeM(unsigned mult, unsigned k, T* U, unsigned ldU,
			   T* V, unsigned ldV);

  // add lwr to utr. dns
  void addLrM_toUtM(unsigned mult, unsigned k, T* U, unsigned ldU,
                         T* V, unsigned ldV);

  // add lwr to mblock and truncate to min(eps, kgoal)
  void addLrM(unsigned k, T* U, unsigned ldU, T* V, unsigned ldV,
	      double eps, unsigned kgoal, contLowLevel<T>* haarInfo=NULL, T* X=NULL,
	      unsigned ldX=0, T* Y_=NULL, unsigned ldY_=0);

  // add lwr to mblock and truncate to min(eps, kgoal)
  void addLrM_Exact(unsigned k, T* U, unsigned ldU, T* V, unsigned ldV);

  // add lwr to mblock and truncate to min(eps, kgoal), returns remainder
  void addLrM_rmnd(unsigned k, T* U, unsigned ldU, T* V,
		   unsigned ldV, double eps, unsigned kgoal,
		   unsigned& kR, T* &UR, T* &VR);

  // add dns to mblock and truncate to min(eps, MAX_RANK)
  void addGeM(T* A, unsigned ldA, double eps, unsigned rankmax,
	      contLowLevel<T>* haarInfo=NULL, T* X=NULL, unsigned ldX=0,
	      T* Y_=NULL, unsigned ldY_=0);

  void unify_rows(double, unsigned, const mblock&, const mblock&,
                  contLowLevel<T>* haarInfo=NULL, T* X=NULL, unsigned ldX=0,
                  T* Y_=NULL, unsigned ldY_=0);
  void unify_cols(double, unsigned, const mblock&, const mblock&,
                  contLowLevel<T>* haarInfo=NULL, T* X=NULL, unsigned ldX=0,
                  T* Y_=NULL, unsigned ldY_=0);

  void addMbl(double eps, unsigned rankmax, mblock* mbl, contLowLevel<T>* haarInfo=NULL,
	      T* X=NULL, unsigned ldX=0, T* Y_=NULL, unsigned ldY_=0) {
    if (mbl->isLrM()) addLrM(mbl->rank(), mbl->data, mbl->n1,
			     mbl->data+mbl->rank()*mbl->n1, mbl->n2,
			     eps, rankmax, haarInfo, X, ldX, Y_, ldY_);
    else {
      // if (mbl->isUtM()) add_utr(mbl->data, mbl->n1);
      //      else if (mbl->isLtM()) add_ltr(mbl->data, mbl->n1);
      if (mbl->isHeM()|| mbl->isSyM()) addHeM_toHeM(mbl->data);
      else addGeM(mbl->data, mbl->n1, eps, rankmax, haarInfo, X, ldX, Y_, ldY_);
    }
  }

  void add_dbl_Id() {
    assert(!isLrM() && !isUtM() && !isLtM());
    if (isHeM() || isSyM())
      for (unsigned i=0; i<n1; ++i) data[i*(i+3)/2] += (T) 1.0;
    else
      for (unsigned i=0; i<n1; ++i) data[i*(n1+1)] += (T) 1.0;
  }

  // *************************************************************************
  // mutiply mblock or its transposed by vector


  bool mltaVec(T d, T* x, T* y) const;

  bool mltaGeM(T d, unsigned p, T* X, unsigned ldX, T* Y, unsigned ldY) const;

  bool mltahVec(T d, T* x, T* y) const;  

  bool mltahGeM(T d, unsigned p, T* X, unsigned ldX, T* Y, unsigned ldY) const;

  bool mltatVec(T d, T* x, T* y) const;

  bool mltatGeM(T d, unsigned p, T* X, unsigned ldX, T* Y, unsigned ldY) const;

  // multiply upper triangular by transposed upper triangular Y += d UU^H
  void mltaUtMUtMh_toHeM(T d, T* U) {
    assert(isHeM());
    blas::putrmmh(d, n1, U, data);
  }

  // multiply transposed upper triangular by upper triangular Y += d U^HU
  void mltaUtMhUtM_toHeM(T d, T* U) {
    assert(isHeM());
    blas::putrmhm(d, n1, U, data);
  }  

  // multiply upper triangular by lower triangular Y += d UL
  void mltaUtMLtM(T d, T* U, T* L) {
    assert(n1==n2);
    blas::putrltrmm(d, n1, U, L, data);
  }

  // multiply mbl by transposed upper (packed) triangular Y += d AU^H
  void mltaMblUtM(T d, mblock* mbl, T* U, double eps, unsigned rankmax) {
    assert(n1==mbl->n1 && n2==mbl->n2);
    T* tmp;
    if (mbl->isGeM()) {
      tmp = new T[n1*n2];
      blas::setzero(n1*n2, tmp);
      blas::geputrmmh(d, n1, n2, mbl->getdata(), U, tmp);
      addGeM(tmp, n1, eps, rankmax);
    } else {
      tmp = new T[n2*mbl->rank()];
      blas::setzero(n2*mbl->rank(), tmp);
      blas::putrgemm(d, n2, U, mbl->rank(),
                     mbl->getdata()+n1*mbl->rank(), tmp);
      addLrM(mbl->rank(), mbl->getdata(), n1, tmp, n2, eps, rankmax);
    }
    delete [] tmp;
  }


  // *************************************************************************
  // Inversion routines
  //

  void invertGeM(T* dataC) {
    assert(isGeM() && n1==n2);
    blas::copy(n1*n1, data, dataC);
    lapack::geinv(n1, dataC);
  }

  void invertHeM(T* dataC) {
    assert(isGeM() && isHeM() && n1==n2);
    blas::copy(n1*(n1+1)/2, data, dataC);
    lapack::geinv_herm(n1, dataC);
  }

  void invertSyM(T* dataC) {
    assert(isGeM() && isSyM() && n1==n2);
    blas::copy(n1*(n1+1)/2, data, dataC);
    lapack::geinv_sym(n1, dataC);
  }

  int decomp_Cholesky() {
    assert(isGeM() && isHeM() && n1==n2);
    int inf = lapack::pptrf(n1, data);
    info.is_LtM = info.is_HeM = 0;
    info.is_UtM = 1;
    return inf;
  }

  int decomp_LU(mblock*, mblock*); // generates LU decomp., destroys data

  int decomp_UhDU(int* piv) { // generates U^H D U decomp. of herm, destroys data
    assert(isGeM() && isHeM() && n1==n2);

    // copy columnwise packed matrix (data) to lower triang. part of tmp
    T* tmp = new T[n1*n1];
    unsigned idx1, idx2=0;
    for (unsigned i=0; i<n1; ++i) {
      idx1 = i;
      for (unsigned j=0; j<=i; ++j) {
        tmp[idx1] = data[idx2++];
        idx1 += n1;
      }
    }

    int lwork = -1;
    T* work, optlwork;
    lapack::hetrf(n1, tmp, piv, &optlwork, lwork); // calculate optimal lwork
    lwork = (int) Re(optlwork);
    work = new T[lwork];
    int inf = lapack::hetrf(n1, tmp, piv, work, lwork);
    delete[] work;

    // write back lower triangular part of tmp to data
    idx2 = 0;
    for (unsigned i=0; i<n1; ++i) {
      idx1 = i;
      for (unsigned j=0; j<=i; ++j) {
        data[idx2++] = tmp[idx1];
        idx1 += n1;
      }
    }

    delete[] tmp;
    return inf;
  }

  // *************************************************************************
  // forward/backward substitution
  //

  // solves PL X = B for X, L is unit lower triangular, X is stored in B
  void ltr_solve(unsigned, T*, unsigned) const;

  // solves (PL)^H X = B for X, L is unit lower triangular, X is stored in B
  void ltrh_solve(unsigned, T*, unsigned) const;

  /* solves X (PL)^H = B for X, L is unit lower triangular, X is stored in B
  void ltrh_solve_left(unsigned, T*, unsigned) const;
  */

  // solves U X = B for X
  // U is upper triangular, X is stored in B
  void utr_solve(unsigned m, T* B, unsigned ldB) const {
    assert(isGeM() && isUtM());
    lapack::utrs(n1, data, m, B, ldB);
  }

  // solves U^H X = B for X, U is upper triangular, X is stored in B
  void utrh_solve(unsigned m, T* B, unsigned ldB) const {
    assert(isGeM() && isUtM());
    lapack::utrhs(n1, data, m, B, ldB);
  }

  // solves X U = B for X, U is upper triangular, X is stored in B
  void utr_solve_left(unsigned, T*, unsigned, T*, unsigned) const;


  // *************************************************************************
  // Misc routines
  //

  // Frobenius norm
  double nrmF2() const;

  // truncate this block to absolute accuracy eps
  void trunc_abs(double);

  // truncate this block to given rank
  void trunc_rank(unsigned);

  void outp(std::ostream& os) const {
    for (unsigned i=0; i<n1; ++i) {
      for (unsigned j=0; j<n2; ++j) os << data[i+j*n1] << ' ';
      os << std::endl;
    }
  }

  void psout(std::ostream& os, unsigned nmax, unsigned b1, unsigned b2,
             bool refl, unsigned level=0, cluster* cl1=NULL,
             void (*f)(cluster*, const unsigned, double*,
                       const double* const)=NULL) const {
    if (isLrM()) psout_LrM(os, nmax, b1, b2, refl, level, cl1, f);
    else psout_GeM(os, nmax, b1, b2, refl);
  }

  void save(std::ofstream& os) {
    // 0 isLrM, 1 isGeM, 2 isLtM, 3 isUtM, 4 dense
    os.write((char*) &n1, sizeof(unsigned));
    os.write((char*) &n2, sizeof(unsigned));
    unsigned status, rankt;
    if (isLrM()) {
      status = 0;
      os.write((char*) &status, sizeof(unsigned));
      rankt = rank();
      os.write((char*) &rankt, sizeof(unsigned));
      if (rankt>0)
        os.write((char*) data, rankt*(n1+n2)*sizeof(T));
    } else { // mbl is dense
      assert(isGeM());
      if (isHeM()) {
        status = 1;
        os.write((char*) &status, sizeof(unsigned));
        os.write((char*) data, n1*(n1+1)/2*sizeof(T));
      } else if (isLtM()) {
        status = 2;
        os.write((char*) &status, sizeof(unsigned));
        os.write((char*) data, n1*(n1+1)/2*sizeof(T));
      } else if (isUtM()) {
        status = 3;
        os.write((char*) &status, sizeof(unsigned));
        os.write((char*) data, n1*(n1+1)/2*sizeof(T));
      } if (isSyM()) {
        status = 5;
        os.write((char*) &status, sizeof(unsigned));
        os.write((char*) data, n1*(n1+1)/2*sizeof(T));
      }else {
        status = 4;
        os.write((char*) &status, sizeof(unsigned));
        os.write((char*) data, n1*n2*sizeof(T));
      }
    }
  }

  void load(std::ifstream& is) {
    // 0 isLrM, 1 isGeM, 2 isLtM, 3 isUtM, 4 dense
    unsigned status, rank;
    is.read((char*) &n1, sizeof(unsigned));
    is.read((char*) &n2, sizeof(unsigned));
    is.read((char*) &status, sizeof(unsigned));
    if (status == 0) {
      is.read((char*) &rank, sizeof(unsigned));
      setrank(rank);
      is.read((char*) data, rank*(n1+n2)*sizeof(T));
    } else if (status == 1) {
      setHeM();
      is.read((char*) data, n1*(n1+1)/2*sizeof(T));
    } else if (status == 2) {
      setLtM();
      is.read((char*) data, n1*(n1+1)/2*sizeof(T));
    } else if (status == 3) {
      setUtM();
      is.read((char*) data, n1*(n1+1)/2*sizeof(T));
    } else if (status == 4) {
      setGeM();
      is.read((char*) data, n1*n2*sizeof(T));
    } else if (status == 5) {
      setSyM();
      is.read((char*) data, n1*(n1+1)/2*sizeof(T));
    }else {
      std::cerr<<"Error while loading mbls"<<std::endl;
      exit(1);
    }
  }

  // copy mblock
  void copy(mblock&);

  // copy transposed of off-diagonal mblocks
  void copyH(mblock&);

  // copy low-rank matrix to this block
  void cpyLrM(unsigned k, T* U, T* V) {
    setrank(k);
    if (k>0) {
      blas::copy(k*n1, U, data);
      blas::copy(k*n2, V, data+k*n1);
    }
  }

  // copy low-rank matrix to this block with recompression
  void cpyLrM_cmpr(unsigned k, T* U, unsigned ldU, T* V, unsigned ldV,
                   double eps, unsigned rankmax) {
    freedata();
    addtrll(k, U, ldU, V, ldV, eps, rankmax);
  }

  void cpyGeM(T* A) {
    setGeM();
    blas::copy(n1*n2, A, data);
  }

  void cpyHeM(T* A) {
    setHeM();
    blas::copy(n1*(n1+1)/2, A, data);
  }

  void cpySyM(T* A) {
    setSyM();
    blas::copy(n1*(n1+1)/2, A, data);
  }

  void cpyUtM(T* A) {
    setUtM();
    blas::copy(n1*(n1+1)/2, A, data);
  }

  void cpyLtM(T* A) {
    setLtM();
    blas::copy(n1*(n1+1)/2, A, data);
  }

  void cpy_mbl(unsigned inf[5], T* A) {
    if (inf[1]) cpyLrM(inf[0], A, A+inf[0]*n1);
    else {
      if (inf[2]) cpyHeM(A);
      else if (inf[3]) cpyLtM(A);
      else if (inf[4]) cpyUtM(A);
      else cpyGeM(A);
    }
  }

  void convLrM_toGeM();  // internal conversion
  void convLrM_toGeM(unsigned, T*, unsigned, T*, unsigned);
  void convGeM_toLrM(double);   // internal conversion

  void convHeM_toGeM(T*, unsigned) const;
  void convSyM_toGeM(T*, unsigned) const;
  void convLrM_toGeM(T*, unsigned) const;
  void convGeM_toGeM(T*, unsigned) const;  // different leading dim

  // compute minimum diagonal entry of upper triangular matrix
  T get_min_dentr() const;

  // compute row/column norms of block
  void rownrms2(T* sum)
  {
    if (isGeM()) {
      assert(!isHeM() && !isSyM() && !isLtM() && !isUtM());
      for (unsigned j=0; j<n2; ++j)
	for (unsigned i=0; i<n1; ++i) sum[i] += abs2(data[i+j*n1]);
    } else {

      T* skpV = new T[bl_rank*(bl_rank+1)/2];
      for (unsigned l2=0; l2<bl_rank; ++l2) {
	for (unsigned l1=0; l1<=l2; ++l1)
	  skpV[UTR(l1, l2)] = blas::scpr(n2, data+bl_rank*n1 + l2*n2,
					 data+bl_rank*n1 + l1*n2);
      }

      for (unsigned l2=0; l2<bl_rank; ++l2) {
	for (unsigned l1=0; l1<=l2; ++l1) {
	  const T e = skpV[UTR(l1, l2)];
	  for (unsigned i=0; i<n1; ++i)
	    sum[i] += Re(data[i+l1*n1]*conj(data[i+l2*n1])*e);
	}

	for (unsigned l1=l2+1; l1<bl_rank; ++l1) {
	  const T e = skpV[UTR(l2, l1)];
	  for (unsigned i=0; i<n1; ++i)
	    sum[i] += Re(data[i+l1*n1]*conj(data[i+l2*n1])*e);
	}
      }

      delete [] skpV;
    }
  }

  void colnrms2(T* sum)
  {
    if (isGeM()) {
      assert(!isHeM() && !isSyM() &&!isLtM() && !isUtM());
      for (unsigned j=0; j<n2; ++j)
	for (unsigned i=0; i<n1; ++i) sum[j] += abs2(data[i+j*n1]);
    } else {
      T* skpU = new T[bl_rank*(bl_rank+1)/2];
      for (unsigned l2=0; l2<bl_rank; ++l2) {
	for (unsigned l1=0; l1<=l2; ++l1)
	  skpU[UTR(l1, l2)] = blas::scpr(n1, data + l2*n1, data + l1*n1);
      }

      T* const V = data + bl_rank*n1;
      for (unsigned l2=0; l2<bl_rank; ++l2) {
	for (unsigned l1=0; l1<=l2; ++l1) {
	  const T e = skpU[UTR(l1, l2)];
	  for (unsigned j=0; j<n2; ++j)
	    sum[j] += Re(V[j+l1*n2]*conj(V[j+l2*n2])*e);
	}

	for (unsigned l1=l2+1; l1<bl_rank; ++l1) {
	  const T e = skpU[UTR(l2, l1)];
	  for (unsigned j=0; j<n2; ++j)
	    sum[j] += Re(V[j+l1*n2]*conj(V[j+l2*n2])*e);
	}
      }
      delete [] skpU;
    }
  }

 void symnrms2(T* sum1, T* sum2)
  {
    if (isGeM()) {
      if (isHeM()||isSyM()) {
        for (unsigned j=0; j<n2; ++j) {
          for (unsigned i=0; i<j; ++i) {
            const T e = abs2(data[UTR(i, j)]);
            sum2[i] += e;
            sum2[j] += e;
          }
          sum2[j] += abs2(data[UTR(j, j)]);
        }
      } else {
        assert(!isLtM() && !isUtM());
        for (unsigned j=0; j<n2; ++j)
          for (unsigned i=0; i<n1; ++i) {
            const T d = abs2(data[i+j*n1]);
            sum2[j] += d;
            sum1[i] += d;
          }
      }
    } else {
      T* skp = new T[bl_rank*(bl_rank+1)/2];
      for (unsigned l2=0; l2<bl_rank; ++l2) {
        for (unsigned l1=0; l1<=l2; ++l1)
          skp[UTR(l1, l2)] = blas::scpr(n1, data + l2*n1, data + l1*n1);
      }
      
      T* const V = data+bl_rank*n1;

      for (unsigned l2=0; l2<bl_rank; ++l2) {
        for (unsigned l1=0; l1<=l2; ++l1) {
          const T e = skp[UTR(l1, l2)];
          for (unsigned j=0; j<n2; ++j)
            sum2[j] += Re(V[j+l1*n2]*conj(V[j+l2*n2])*e);
       }
        
        for (unsigned l1=l2+1; l1<bl_rank; ++l1) {
          const T e = skp[UTR(l2, l1)];
          for (unsigned j=0; j<n2; ++j)
            sum2[j] += Re(V[j+l1*n2]*conj(V[j+l2*n2])*e);
        }
      }

      for (unsigned l2=0; l2<bl_rank; ++l2) {
        for (unsigned l1=0; l1<=l2; ++l1)
          skp[UTR(l1, l2)] = blas::scpr(n2, V + l2*n2, V + l1*n2);
      }
      
      for (unsigned l2=0; l2<bl_rank; ++l2) {
        for (unsigned l1=0; l1<=l2; ++l1) {
          const T e = skp[UTR(l1, l2)];
          for (unsigned i=0; i<n1; ++i)
            sum1[i] += Re(data[i+l1*n1]*conj(data[i+l2*n1])*e);
        }

        for (unsigned l1=l2+1; l1<bl_rank; ++l1) {
          const T e = skp[UTR(l2, l1)];
          for (unsigned i=0; i<n1; ++i)
            sum1[i] += Re(data[i+l1*n1]*conj(data[i+l2*n1])*e);
        }
      }

      delete [] skp;
    }
  }


  void diag_sym(T* diag)
  {
    assert(isGeM() && (isHeM()||isSyM()) );
    for (unsigned j=0; j<n2; ++j) diag[j] = data[UTR(j, j)];
  }

  void scale_cols(T* D)
  {
    if (isGeM()) {
      assert(!isHeM() && !isSyM() && !isLtM() && !isUtM());
      for (unsigned j=0; j<n2; ++j) blas::scal(n1, D[j], data+j*n1);
    } else {
      for (unsigned l=0; l<bl_rank; ++l)
	blas::scal(n2, data+bl_rank*n1+l*n2, D);
    }
  }

  void scale_rows(T* D)
  {
    if (isGeM()) {
      assert(!isHeM() && !isSyM() && !isLtM() && !isUtM());
      for (unsigned j=0; j<n2; ++j)
	blas::scal(n1, data+j*n1, D);
    } else {
      for (unsigned l=0; l<bl_rank; ++l)
	blas::scal(n1, data+l*n1, D);
    }
  }

  // A -> D_1 A D_2
  void scale_sym(T* D1, T* D2)
  {
    if (isGeM()) {
      if (isHeM() || isSyM() ) {
	for (unsigned j=0; j<n2; ++j)
	  for (unsigned i=0; i<=j; ++i) data[UTR(i, j)] *= D1[i]*D2[j];
      } else {
	assert(!isLtM() && !isUtM());
	for (unsigned j=0; j<n2; ++j)
	  for (unsigned i=0; i<n1; ++i) data[i+j*n1] *= D1[i]*D2[j];
      }
    } else {
      for (unsigned l=0; l<bl_rank; ++l)
	blas::scal(n1, data+l*n1, D1);
      for (unsigned l=0; l<bl_rank; ++l)
	blas::scal(n2, data+bl_rank*n1+l*n2, D2);
    }
  }



private:
  void copy_(mblock<T>&);
  void copyH_(mblock<T>&);

  void initrnd_(mblock<T>&);

  // multiply lwr by vector y += d U V^H x, returns 1 if y has changed
  bool mltaLrMVec_(T d, T* x, T* y) const;

  // multiply herm transposed lwr by vector y += d V U^H x, returns 1 if y changed
  bool mltaLrMhVec_(T d, T* x, T* y) const;

 // multiply transposed lwr by vector y += d (U V^H)^T x, returns 1 if y changed
  bool mltaLrMtVec_(T d, T* x, T* y) const;

  // multiply lwr by dense Y += d U V^H X, returns 1 if X has changed
  bool mltaLrMGeM_(T d, unsigned p, T* X, unsigned ldX,
                      T* Y, unsigned ldY) const;

  // multiply transposed lwr by vector Y += d V U^H X, returns 1 if X changed
  bool mltaLrMhGeM_(T d, unsigned p, T* X, unsigned ldX,
                       T* Y, unsigned ldY) const;

  // multiply herm. transposed lwr by vector Y += d (U V^H)^T X, returns 1 if X changed
  bool mltaLrMtGeM_(T d, unsigned p, T* X, unsigned ldX,
		    T* Y, unsigned ldY) const;

  bool mltaVec_(T d, T* x, T* y) const;
  bool mltaGeM_(T d, unsigned p, T* X, unsigned ldX,
                  T* Y, unsigned ldY) const;
  bool mltahVec_(T d, T* x, T* y) const;
  bool mltahGeM_(T d, unsigned p, T* X, unsigned ldX,
                   T* Y, unsigned ldY) const;
  bool mltatVec_(T d, T* x, T* y) const;
  bool mltatGeM_(T d, unsigned p, T* X, unsigned ldX,
                   T* Y, unsigned ldY) const;
  double nrmF2_() const;
  void ltrh_solve_left_(unsigned, T*, unsigned) const;
};

template<class T>
void mblock<T>::psout_GeM(std::ostream& os, unsigned nmax,
                          unsigned b1, unsigned b2, bool refl) const
{
  assert(isGeM());

  double no2, mo1, mo2;
  mo1 = nmax-b1;

  if (refl) {
    no2 = b2+n1;
    mo2 = mo1-n2;
  } else {
    no2 = b2+n2;
    mo2 = mo1-n1;
  }

  os << "0.1 setlinewidth" << std::endl;
  os << no2 << ' ' << mo2 << " m" << std::endl;
  if ((!(isHeM()||isSyM()) || refl) && !isUtM())
    os << b2 << ' ' << mo2 << " l" << std::endl;
  os << b2 << ' ' << mo1 << " l" << std::endl;
  if ((!(isHeM()||isSyM()) || !refl) && !isLtM())
    os << no2 << ' ' << mo1 << " l" << std::endl;
  os << no2 << ' ' << mo2 << " l" << std::endl;
  if (refl) os << "flr";
  else os << "fr";
  os << std::endl << "cs" << std::endl;

  os << "bl" << std::endl;
  os << no2 << ' ' << mo2 << " m" << std::endl;
  if (!(isHeM()||isSyM()) && !isUtM()) {
    os << b2 << ' ' << mo2 << " l" << std::endl;
    os << b2 << ' ' << mo1 << " l" << std::endl;
  }
  if (!isLtM()) os << no2 << ' ' << mo1 << " l" << std::endl;
  os << no2 << ' ' << mo2 << " l" << std::endl;
  os << "cs" << std::endl;
}

#include <iomanip>

template<class T>
void mblock<T>::psout_LrM(std::ostream& os, unsigned nmax, unsigned b1,
                          unsigned b2, bool refl, unsigned level, cluster* cl1,
                          void (*f)(cluster* cl, const unsigned lvl,
                                    double* X, const double* const basis)) const
{
  assert(isLrM());
  double mo1 = nmax-b1, no2, mo2;

  if (refl) {
    no2 = b2+n1;
    mo2 = mo1-n2;
  } else {
    no2 = b2+n2;
    mo2 = mo1-n1;
  }

  const unsigned k = rank();
  unsigned cols = 0;
  double* val = NULL;
  double* basis = NULL;
  if (f!=NULL) {
    cols = pow2(level);
    val = new double[2*cols*cl1->size()];
    f(cl1, level, val, basis);
  }
  // parameter to adjust
  const double scal = 0.4/cols;
  const double scalf = 0.3;
  const double max_jump = 0.5;

  if (k) {
    os << "0.1 setlinewidth" << std::endl;
    os << b2 << ' ' << mo2 << " m" << std::endl;
    os << b2 << ' ' << mo1 << " l" << std::endl;
    os << no2 << ' ' << mo1 << " l" << std::endl;
    os << no2 << ' ' << mo2 << " l" << std::endl;
    os << b2 << ' ' << mo2 << " l" << std::endl;
    unsigned t;
    if (refl) {
      t = n2;
      os << "flg";
    } else {
      t = n1;
      os << "fg";
    }
    os << std::endl << "cs" << std::endl;
    if (f!=NULL) {
      for (unsigned i=0; i<cols; i++) {
        os << "0.7 setgray" << std::endl;
        os << b2 << ' ' << mo1-(i+1)*scal*t << " m" << std::endl;
        os << no2 << ' ' << mo1-(i+1)*scal*t << " l" << std::endl;
        os << "stroke" << std::endl;
        os << "fblue" <<std::endl;
        os << "0.8 setlinewidth" << std::endl;
        for (unsigned j=0; j<cl1->size()-1; j++) {
          if (val[j+i*cl1->size()]!=0 && val[j+1+i*cl1->size()]!=0 &&
              std::abs(val[j+i*cl1->size()]-val[j+1+i*cl1->size()])<max_jump) {
            os << b2+j+1 << ' ' << mo1-((i+1)-val[j+i*cl1->size()]*scalf)*t*scal
            << " m" <<std::endl;
            os << b2+j+2 << ' ' << mo1-((i+1)-val[j+1+i*cl1->size()]*scalf)*t*scal
            << " l" <<std::endl;
            os << "stroke" << std::endl;
          }
        }
        os << "0.1 setlinewidth" << std::endl;
      }
    }
    os << "0 setgray" << std::endl;
    os << b2 << ' ' << mo2 << " m" << std::endl;
    os << b2 << ' ' << mo1 << " l" << std::endl;
    os << no2 << ' ' << mo1 << " l" << std::endl;
    os << no2 << ' ' << mo2 << " l" << std::endl;
    os << b2 << ' ' << mo2 << " l" << std::endl;
    os << "cs" << std::endl << "hf " << 0.5*MIN(n1,n2) << " sf" << std::endl;
    //os << "cs" << std::endl << "hf " << 0.3*MIN(n1,n2) << " sf" << std::endl;
    os << b2+0.1*n2 << ' ' << mo2+0.1*n1 << " m ("
    << k << " " << ") show" << std::endl;
  }
  delete [] val;
}

#endif  // MBLOCK_H
