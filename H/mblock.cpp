/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "mblock.h"

/* solves X (PL)^H = X L^H P^{-1} = B for X
// L is unit lower triangular, X is stored in B
template<class T>
void mblock<T>::ltrh_solve_left_(unsigned m, T* B, unsigned ldB) const
{
  assert(isGeM() && isLtM());
  T *p = data, *tmp = new T[n1*m];

  for (unsigned j=0; j<n1; ++j) {
    unsigned ip = (unsigned) Re(*p);
    blas::copy(n1, B+j*ldB, tmp+ip*n1);
    p+= n1-j;
  }
  
  for (unsigned j=0; j<m; ++j) {
    for (unsigned k=0; k<j; ++k)
      blas::axpy(n1, data[LTR(j,k,n1)], tmp+k*n1, tmp+j*n1);
  }

  delete [] tmp;
}
*/

template<class T>
double mblock<T>::nrmF2_() const
{
  double nrm2 = 0.0;
  if (isLrM()) {
    T *U = data, *V = data+bl_rank*n1;
    for (unsigned k=0; k<bl_rank; ++k) {
      T s = 0.0;
      for (unsigned l=0; l<k; ++l)
	s += blas::scpr(n1, U+l*n1, U+k*n1) * blas::scpr(n2, V+k*n2, V+l*n2);
      nrm2 += 2*Re(s)+blas::sqrsum(n1, U+k*n1)*blas::sqrsum(n2, V+k*n2);
    }
  } else {
    if (isHeM() || isSyM()) {
      for (unsigned j=0; j<n2; ++j) {
	for (unsigned i=0; i<j; ++i) nrm2 += 2*abs2(data[UTR(i, j)]);
	nrm2 += abs2(data[UTR(j, j)]);
      }
    } else if (isLtM()) {
      for (unsigned j=0; j<n2; ++j) {
	for (unsigned i=j+1; i<n1; ++i) nrm2 += abs2(data[LTR(i,j,n1)]);
      }
      nrm2 += n2;
    } else if (isUtM()) {
      for (unsigned j=0; j<n2; ++j) {
	for (unsigned i=0; i<=j; ++i) nrm2 += abs2(data[UTR(i,j)]);
      }
    } else nrm2 = blas::sqrsum(n1*n2, data);
  }
  return nrm2;
}


// copies mbl to this
template<class T> void mblock<T>::copy_(mblock<T>& mbl)
{
  n1 = mbl.n1;
  n2 = mbl.n2;

  if (mbl.isLrM()) {
    const unsigned k = mbl.rank();
    setrank(k);
    if (k>0) blas::copy(k*(n1+n2), mbl.data, data);
  } else { // mbl is dense
    assert(mbl.isGeM());
    if (mbl.isHeM()) {
      setHeM();
      blas::copy(n1*(n1+1)/2, mbl.data, data);
    }
    else if (mbl.isSyM()) {
      setSyM();
      blas::copy(n1*(n1+1)/2, mbl.data, data);
    }
    else if (mbl.isLtM()) {
      setLtM();
      blas::copy(n1*(n1+1)/2, mbl.data, data);
    }
    else if (mbl.isUtM()) {
      setUtM();
      blas::copy(n1*(n1+1)/2, mbl.data, data);
    }
    else {
      setGeM();
      blas::copy(n1*n2, mbl.data, data);
    }
  }
}


// copies transposed of mbl to this, mbl is off-diagonal block!
template<class T> void mblock<T>::copyH_(mblock<T>& mbl)
{
  n1 = mbl.n2;
  n2 = mbl.n1;

  if (mbl.isLrM()) {
    const unsigned k = mbl.rank();
    setrank(k);
    if (k>0) {
      blas::copy(k*n1, mbl.data+k*n2, data);
      blas::copy(k*n2, mbl.data, data+k*n1);
    }
  } else { // mbl is dense
    setGeM();
    blas::transpose(n2, n1, mbl.data, data);
  }
}

// multiply lwr by vector y += d U V^H x
template<class T>
bool mblock<T>::mltaLrMVec_(T d, T* x, T* y) const
{
  assert(isLrM());
  const unsigned k = rank();
  if (k>0) {
    for (unsigned l=0; l<k; ++l) {
      const T e = d * blas::scpr(n2, data+k*n1+l*n2, x);
      blas::axpy(n1, e, data+l*n1, y);
    }
    return true;
  } else return false;
}

// multiply lwr by vector Y += d U V^H X
template<class T>
bool mblock<T>::mltaLrMGeM_(T d, unsigned p, T* X, unsigned ldX, T* Y,
				    unsigned ldY) const
{
  assert(isLrM());
  const unsigned k = rank();
  if (k>0) {
    T* tmp = new T[k*p];
    blas::gemhm(n2, k, p, (T) 1.0, data+k*n1, n2, X, ldX, tmp, k);
    blas::gemma(n1, k, p, d, data, n1, tmp, k, Y, ldY);
    delete [] tmp;
    return true;
  } else return false;
}


// multiply herm transposed lwr by vector y += d V U^H x
template<class T>
bool mblock<T>::mltaLrMhVec_(T d, T* x, T* y) const
{
  assert(isLrM());
  const unsigned k = rank();
  if (k>0) {
    for (unsigned l=0; l<k; ++l) {
      const T e = d * blas::scpr(n1, data+l*n1, x);
      blas::axpy(n2, e, data+k*n1+l*n2, y);
    }
    return true;
  } else return false;
}

// multiply herm transposed lwr by vector Y += d V U^H X
template<class T>
bool mblock<T>::mltaLrMhGeM_(T d, unsigned p, T* X, unsigned ldX, T* Y,
				     unsigned ldY) const
{
  assert(isLrM());
  const unsigned k = rank();
  if (k>0) {
    T* tmp = new T[k*p];
    blas::gemhm(n1, k, p, (T) 1.0, data, n1, X, ldX, tmp, k);
    blas::gemma(n2, k, p, d, data+k*n1, n2, tmp, k, Y, ldY);
    delete [] tmp;
    return true;
  } else return false;
}


// multiply transposed lwr by vector y += d (U V^H)^T x
template<class T>
bool mblock<T>::mltaLrMtVec_(T d, T* x, T* y) const
{
  assert(isLrM());
  
  //blas::conj(n2, x);
  //blas::conj(n1, y);
  //bool succ = mltaLrMhVec_(conj(d), x, y);
  //blas::conj(n2, x);
  //blas::conj(n1, y);

  blas::conj(bl_rank*(n1+n2), data);
  bool succ = mltaLrMhVec_(d, x, y);
  blas::conj(bl_rank*(n1+n2), data);

  return succ;
}
// multiply transposed lwr by vector Y += d (U V^H)^T X
template<class T>
bool mblock<T>::mltaLrMtGeM_(T d, unsigned p, T* X, unsigned ldX, T* Y,
				     unsigned ldY) const
{
  bool succ = false;
  for(unsigned l=0; l<p; ++l) succ |= mltaLrMtVec_(d, X+ldX*l, Y+ldY*l);
  return succ;
}


template<class T>
bool mblock<T>::mltaVec_(T d, T* x, T* y) const
{
  if (isLrM()) return mltaLrMVec_(d, x, y);
  else {
    if (isHeM()) mltaHeMVec(d, x, y);    
    else if (isSyM()) mltaSyMVec(d, x, y);
    else if (isLtM()) mltaLtMVec(d, x, y);
    else if (isUtM()) mltaUtMVec(d, x, y);
    else mltaGeMVec(d, x, y);
    return true;
  }
}


template<class T>
bool mblock<T>::mltaGeM_(T d, unsigned p, T* X, unsigned ldX,
			   T* Y, unsigned ldY) const
{
  if (isLrM()) return mltaLrMGeM_(d, p, X, ldX, Y, ldY);
  else {
    if (isHeM()) mltaHeMGeM(d, p, X, ldX, Y, ldY);
    else if (isSyM()) mltaSyMGeM(d, p, X, ldX, Y, ldY);
    else if (isLtM()) mltaLtMGeM(d, p, X, ldX, Y, ldY);
    else if (isUtM()) mltaUtMGeM(d, p, X, ldX, Y, ldY);
    else mltaGeMGeM(d, p, X, ldX, Y, ldY);
    return true;
  }
}


template<class T>
bool mblock<T>::mltahVec_(T d, T* x, T* y) const
{
  if (isLrM()) return mltaLrMhVec_(d, x, y);
  else {
    if (isHeM()) mltaHeMVec(d, x, y);
    else if (isSyM()) mltaSyMhVec(d, x, y);
    else if (isLtM()) mltaLtMhVec(d, x, y);
    else if (isUtM()) mltaUtMhVec(d, x, y);
    else mltaGeMhVec(d, x, y);
    return true;
  }
}


template<class T>
bool mblock<T>::mltahGeM_(T d, unsigned p, T* X, unsigned ldX,
			  T* Y, unsigned ldY) const
{
  if (isLrM()) return mltaLrMhGeM_(d, p, X, ldX, Y, ldY);
  else {
    if (isHeM()) mltaHeMGeM(d, p, X, ldX, Y, ldY);
    else if (isSyM()) mltaSyMhGeM(d, p, X, ldX, Y, ldY);
    else if (isLtM()) mltaLtMhGeM(d, p, X, ldX, Y, ldY);
    else if (isUtM()) mltaUtMhGeM(d, p, X, ldX, Y, ldY);
    else mltaGeMhGeM(d, p, X, ldX, Y, ldY);
    return true;
  }
}

template<class T>
bool mblock<T>::mltatVec_(T d, T* x, T* y) const
{
  if (isLrM()) return mltaLrMtVec_(d, x, y);
  else {
    if (isHeM()) mltaHeMtVec(d, x, y);
    else if (isSyM()) mltaSyMVec(d, x, y);
    else if (isLtM()) mltaLtMtVec(d, x, y);
    else if (isUtM()) mltaUtMtVec(d, x, y);
    else mltaGeMtVec(d, x, y);
    return true;
  }
}

template<class T>
bool mblock<T>::mltatGeM_(T d, unsigned p, T* X, unsigned ldX,
			  T* Y, unsigned ldY) const
{
  if (isLrM()) return mltaLrMtGeM_(d, p, X, ldX, Y, ldY);
  else {
    if (isHeM()) mltaHeMtGeM(d, p, X, ldX, Y, ldY);
    else if (isSyM()) mltaSyMGeM(d, p, X, ldX, Y, ldY);
    else if (isLtM()) mltaLtMtGeM(d, p, X, ldX, Y, ldY);
    else if (isUtM()) mltaUtMtGeM(d, p, X, ldX, Y, ldY);
    else mltaGeMtGeM(d, p, X, ldX, Y, ldY);
    return true;
  }
}



// Instanzen
template<> void mblock<double>::copyH(mblock<double>& mbl) { copyH_(mbl); }
template<> void mblock<float>::copyH(mblock<float>& mbl) { copyH_(mbl); }
template<> void mblock<dcomp>::copyH(mblock<dcomp>& mbl) { copyH_(mbl); }
template<> void mblock<scomp>::copyH(mblock<scomp>& mbl) { copyH_(mbl); }

template<> void mblock<double>::copy(mblock<double>& mbl) { copy_(mbl); }
template<> void mblock<float>::copy(mblock<float>& mbl) { copy_(mbl); }
template<> void mblock<dcomp>::copy(mblock<dcomp>& mbl) { copy_(mbl); }
template<> void mblock<scomp>::copy(mblock<scomp>& mbl) { copy_(mbl); }

template<>
bool mblock<double>::mltaVec(double d, double* x, double* y) const
{ return mltaVec_(d, x, y); }

template<>
bool mblock<float>::mltaVec(float d, float* x, float* y) const
{ return mltaVec_(d, x, y); }

template<>
bool mblock<dcomp>::mltaVec(dcomp d, dcomp* x, dcomp* y) const
{ return mltaVec_(d, x, y); }

template<>
bool mblock<scomp>::mltaVec(scomp d, scomp* x, scomp* y) const
{ return mltaVec_(d, x, y); }



template<>
bool mblock<double>::mltaGeM(double d, unsigned p, double* X, unsigned ldX,
			     double* Y, unsigned ldY) const
{ return mltaGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<float>::mltaGeM(float d, unsigned p, float* X, unsigned ldX,
			    float* Y, unsigned ldY) const
{ return mltaGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<dcomp>::mltaGeM(dcomp d, unsigned p, dcomp* X, unsigned ldX,
			    dcomp* Y, unsigned ldY) const
{ return mltaGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<scomp>::mltaGeM(scomp d, unsigned p, scomp* X, unsigned ldX,
			    scomp* Y, unsigned ldY) const
{ return mltaGeM_(d, p, X, ldX, Y, ldY); }



template<>
bool mblock<double>::mltahVec(double d, double* x, double* y) const
{ return mltahVec_(d, x, y); }

template<> bool mblock<float>::mltahVec(float d, float* x, float* y) const
{ return mltahVec_(d, x, y); }

template<> bool mblock<dcomp>::mltahVec(dcomp d, dcomp* x, dcomp* y) const
{ return mltahVec_(d, x, y); }

template<> bool mblock<scomp>::mltahVec(scomp d, scomp* x, scomp* y) const
{ return mltahVec_(d, x, y); }



template<>
bool mblock<double>::mltahGeM(double d, unsigned p, double* X,
			      unsigned ldX, double* Y, unsigned ldY) const
{ return mltahGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<float>::mltahGeM(float d, unsigned p, float* X, unsigned ldX,
			     float* Y, unsigned ldY) const
{ return mltahGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<dcomp>::mltahGeM(dcomp d, unsigned p, dcomp* X, unsigned ldX,
			     dcomp* Y, unsigned ldY) const
{ return mltahGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<scomp>::mltahGeM(scomp d, unsigned p, scomp* X, unsigned ldX,
			     scomp* Y, unsigned ldY) const
{ return mltahGeM_(d, p, X, ldX, Y, ldY); }


template<>
bool mblock<double>::mltatVec(double d, double* x, double* y) const
{ return mltatVec_(d, x, y); }

template<> bool mblock<float>::mltatVec(float d, float* x, float* y) const
{ return mltatVec_(d, x, y); }

template<> bool mblock<dcomp>::mltatVec(dcomp d, dcomp* x, dcomp* y) const
{ return mltatVec_(d, x, y); }

template<> bool mblock<scomp>::mltatVec(scomp d, scomp* x, scomp* y) const
{ return mltatVec_(d, x, y); }



template<>
bool mblock<double>::mltatGeM(double d, unsigned p, double* X,
			      unsigned ldX, double* Y, unsigned ldY) const
{ return mltatGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<float>::mltatGeM(float d, unsigned p, float* X, unsigned ldX,
			     float* Y, unsigned ldY) const
{ return mltatGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<dcomp>::mltatGeM(dcomp d, unsigned p, dcomp* X, unsigned ldX,
			       dcomp* Y, unsigned ldY) const
{ return mltatGeM_(d, p, X, ldX, Y, ldY); }

template<>
bool mblock<scomp>::mltatGeM(scomp d, unsigned p, scomp* X, unsigned ldX,
			     scomp* Y, unsigned ldY) const
{ return mltatGeM_(d, p, X, ldX, Y, ldY); }


template<> double mblock<double>::nrmF2() const { return nrmF2_(); }
template<> double mblock<float>::nrmF2() const { return nrmF2_(); }
template<> double mblock<dcomp>::nrmF2() const { return nrmF2_(); }
template<> double mblock<scomp>::nrmF2() const { return nrmF2_(); }

/*
template<> void mblock<double>::ltrh_solve_left(unsigned m, double* B,
						unsigned ldB) const
{ ltrh_solve_left_(m, B, ldB); }

template<> void mblock<float>::ltrh_solve_left(unsigned m, float* B,
					       unsigned ldB) const
{ ltrh_solve_left_(m, B, ldB); }

template<> void mblock<dcomp>::ltrh_solve_left(unsigned m, dcomp* B,
						unsigned ldB) const
{ ltrh_solve_left_(m, B, ldB); }

template<> void mblock<scomp>::ltrh_solve_left(unsigned m, scomp* B,
						unsigned ldB) const
{ ltrh_solve_left_(m, B, ldB); }
*/
