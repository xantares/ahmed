/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <cmath>
#include "mblock.h"
#include "basmod.h"

#define EPS0 1e-64

template<>
void mblock<double>::trunc_abs(double eps)
{
  assert(isLrM());
  unsigned k = bl_rank;

  if (k) {
    unsigned LWORK = 10*k;
    unsigned size = k*(n1+n2)+LWORK+(2*k+3)*k;

    double *tmp1=new double[size], *tmp2=tmp1+k*n1; // tmp2 = new V
    assert(tmp1!=NULL);

    // copy data to tmp1
    blas::copy(k*(n1+n2), data, tmp1);


    freedata();

    // QR-Zerlegung von A und B
    int INFO;
    double* const WORK = tmp2 + n2*k; // LWORK
    double* const tau1 = WORK + LWORK;   // mmin
    double* const tau2 = tau1 + k;    // nmin

    INFO = blas::geqrf(n1, k, tmp1, tau1, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n2, k, tmp2, tau2, LWORK, WORK);
    assert(INFO==0);

    double* const R = tau2 + k;       // mmin*nmin

    // Berechne R1*R2^T
    blas::utrmmh(k, k, k, tmp1, n1, tmp2, n2, R);

    // eliminate small entries
    const double thresh = eps*blas::nrm2(k*k, R)/k;
    for (unsigned i=0; i<k; ++i)
      for (unsigned j=0; j<k; ++j)
        if (fabs(R[i+j*k])<thresh) R[i+j*k] = 0.0;

    double* const S = R + k*k;     // amin
    double* const VT = S + k;         // amin*nmin

    // SVD von R1*R2^T
    INFO = blas::gesvd(k, k, R, S, VT, k, LWORK, WORK);
    assert(INFO==0);

    unsigned kt = k;
    while (kt>0 && S[kt-1]<eps) --kt;

    if (kt>0) {
      setrank(kt);

      // copy R to data und scale the columns with S
      unsigned j;
      for (j=0; j<kt; ++j) {
        unsigned i = 0;
        for (; i<k; ++i) data[i+j*n1] = S[j] * R[i+j*k];
        for (; i<n1; ++i) data[i+j*n1] = D_ZERO;
      }

      double* dataV = data + kt*n1;

      // copy V to dataV
      for (j=0; j<kt; ++j) {
        unsigned i;
        for (i=0; i<k; ++i) dataV[i+j*n2] = VT[i*k+j];
        for (; i<n2; ++i) dataV[i+j*n2] = D_ZERO;
      }

      INFO = blas::ormqr(n1, kt, k, tmp1, tau1, data, LWORK, WORK);
      assert(INFO==0);

      INFO = blas::ormqr(n2, kt, k, tmp2, tau2, dataV, LWORK, WORK);
      assert(INFO==0);
    }

    delete [] tmp1;
  }
}

template<>
void mblock<double>::trunc_rank(unsigned rank_new)
{
  assert(isLrM());
  unsigned k = bl_rank;

  if (k) {
    unsigned LWORK = 10*k;
    unsigned size = k*(n1+n2)+LWORK+(2*k+3)*k;

    double *tmp1=new double[size], *tmp2=tmp1+k*n1; // tmp2 = new V
    assert(tmp1!=NULL);

    // copy data to tmp1
    blas::copy(k*(n1+n2), data, tmp1);


    freedata();

    // QR-Zerlegung von A und B
    int INFO;
    double* const WORK = tmp2 + n2*k; // LWORK
    double* const tau1 = WORK + LWORK;   // mmin
    double* const tau2 = tau1 + k;    // nmin

    INFO = blas::geqrf(n1, k, tmp1, tau1, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n2, k, tmp2, tau2, LWORK, WORK);
    assert(INFO==0);

    double* const R = tau2 + k;       // mmin*nmin

    // Berechne R1*R2^T
    blas::utrmmh(k, k, k, tmp1, n1, tmp2, n2, R);

    // eliminate small entries
    const double thresh = EPS0*blas::nrm2(k*k, R)/k;
    for (unsigned i=0; i<k; ++i)
      for (unsigned j=0; j<k; ++j)
        if (fabs(R[i+j*k])<thresh) R[i+j*k] = 0.0;

    double* const S = R + k*k;     // amin
    double* const VT = S + k;         // amin*nmin

    // SVD von R1*R2^T
    INFO = blas::gesvd(k, k, R, S, VT, k, LWORK, WORK);
    assert(INFO==0);

    unsigned kt = rank_new;

    if (kt>0) {
      setrank(kt);

      // copy R to data und scale the columns with S
      unsigned j;
      for (j=0; j<kt; ++j) {
        unsigned i = 0;
        for (; i<k; ++i) data[i+j*n1] = S[j] * R[i+j*k];
        for (; i<n1; ++i) data[i+j*n1] = D_ZERO;
      }

      double* dataV = data + kt*n1;

      // copy V to dataV
      for (j=0; j<kt; ++j) {
        unsigned i;
        for (i=0; i<k; ++i) dataV[i+j*n2] = VT[i*k+j];
        for (; i<n2; ++i) dataV[i+j*n2] = D_ZERO;
      }

      INFO = blas::ormqr(n1, kt, k, tmp1, tau1, data, LWORK, WORK);
      assert(INFO==0);

      INFO = blas::ormqr(n2, kt, k, tmp2, tau2, dataV, LWORK, WORK);
      assert(INFO==0);
    }

    delete [] tmp1;
  }
}

template<>
double mblock<double>::get_min_dentr() const
{
  double min = data[0];
  for (unsigned i=1; i<n1; ++i) {
    double d = data[i*(i+3)/2];
    if (d<min) min = d;
  }
  return min;
}

template<>
void mblock<double>::get_svals_LrM(double* sigma) const
{
  assert(isLrM() && bl_rank>0);

  const unsigned k = bl_rank;
  const unsigned LWORK = 10*k;
  const unsigned size = k*(n1+n2+k+2)+LWORK;

  double* const tmp1 = new double[size];
  assert(tmp1!=NULL);

  // copy data to tmp1
  blas::copy(k*(n1+n2), data, tmp1);

  double* const tmp2 = tmp1+k*n1;  // tmp2 points to V

  // QR-Zerlegung von U und V
  int INFO;
  double* const WORK = tmp2 + n2*k;    // LWORK
  double* const tau1 = WORK + LWORK;   // k
  double* const tau2 = tau1 + k;       // k

  INFO = blas::geqrf(n1, k, tmp1, tau1, LWORK, WORK);
  assert(INFO==0);
  INFO = blas::geqrf(n2, k, tmp2, tau2, LWORK, WORK);
  assert(INFO==0);

  double* const R = tau2 + k;       // k*k

  // Berechne R1*R2^T
  blas::utrmmh(k, k, k, tmp1, n1, tmp2, n2, R);

  // SVD von R1*R2^T
  INFO = blas::svals(k, k, R, sigma, LWORK, WORK);
  assert(INFO==0);

  delete [] tmp1;
}



template<>
void mblock<double>::convLrM_toGeM(unsigned k, double* U, unsigned ldU,
				   double* V, unsigned ldV)
{
  if (k>0) {
    setGeM();
    blas::gemmh(n1, k, n2, D_ONE, U, ldU, V, ldV, data, n1);
  } else init0_GeM(n1, n2);
}


template<>
void mblock<double>::convLrM_toGeM()
{
  assert(isLrM());

  if (bl_rank>0) {
    double* const new_data = new double[n1*n2];
    assert(new_data!=NULL);

    blas::gemmh(n1, bl_rank, n2, D_ONE, data, n1, data+bl_rank*n1, n2,
                new_data, n1);

    delete [] data;
    data = new_data;
    info.is_LrM = info.is_UtM = info.is_LtM = info.is_HeM = info.is_SyM = 0;
  } else init0_GeM(n1, n2);
}


template<>
void mblock<double>::convLrM_toGeM(double* A, unsigned ldA) const
{
  assert(isLrM());

  if (bl_rank>0)
    blas::gemmh(n1, bl_rank, n2, D_ONE, data, n1, data+bl_rank*n1, n2,
                A, ldA);
  else
    for (unsigned j=0; j<n2; ++j) blas::setzero(n1, A+j*ldA);
}


template<>
void mblock<double>::convGeM_toLrM(double eps)
{
  assert(isGeM());

  unsigned nmin = MIN(n1, n2), LWORK = n1+n2+5*nmin;
  int INFO;
  double *tmp = new double[n1*n2+nmin*(n2+1)+LWORK];
  double *S = tmp + n1*n2, *VT = S + nmin, *WORK = VT + nmin*n2;

  blas::copy(n1*n2, data, tmp);

  INFO = blas::gesvd(n1, n2, tmp, S, VT, nmin, LWORK, WORK);
  assert(INFO==0);

  unsigned kt = nmin;
  while (kt>0 && (S[kt-1]<=eps*S[0] || S[kt-1]<EPS0)) --kt;

  setrank(kt);
  blas::copy(kt*n1, tmp, data);

  for (unsigned l=0; l<kt; ++l)
    for (unsigned j=0; j<n2; ++j)
      data[j+n2*l+kt*n1] = S[l] * VT[nmin*j+l];

  delete [] tmp;
}


template<>
void mblock<double>::convGeM_toGeM(double* A, unsigned ldA) const
{
  assert(isGeM());
  for (unsigned j=0; j<n2; ++j) blas::copy(n1, data+j*n1, A+j*ldA);
}


template<>
void mblock<double>::convHeM_toGeM(double* A, unsigned ldA) const
{
  assert(isHeM()||isSyM());

  // dns herm. upper triangular part stored columnwise

  double *p = data;

  for (unsigned j=0; j<n2; ++j) {
    for (unsigned i=0; i<j; ++i) A[i+j*ldA] = A[j+i*ldA] = *p++;
    A[j*(ldA+1)] = *p++;
  }
}

template<>
void mblock<double>::convSyM_toGeM(double* A, unsigned ldA) const
{
  convHeM_toGeM(A, ldA);
}

// appends lwr to lwr
template<>
void mblock<double>::append(unsigned k, double* U, unsigned ldU, double* V,
                            unsigned ldV)
{
  assert(isLrM());

  unsigned rank_new = bl_rank+k;
  double *tmp=new double[rank_new*(n1+n2)];
  blas::copy(bl_rank*n1, data, tmp);
  for (unsigned i=0; i<k; i++)
    blas::copy(n1, U+i*ldU, tmp+bl_rank*n1+i*n1);
  blas::copy(bl_rank*n2, data+bl_rank*n1, tmp+rank_new*n1);
  for (unsigned i=0; i<k; i++)
    blas::copy(n2, V+i*ldV, tmp+rank_new*n1+bl_rank*n2+i*n2);
  delete [] data;
  data = tmp;
  bl_rank = rank_new;
}

// adds lwr to lwr and truncates to min(delta, kgoal)
template<>
void mblock<double>::addtrll(unsigned k, double* U, unsigned ldU, double* V,
                             unsigned ldV, double delta, unsigned kgoal,
                             contLowLevel<double>* haarInfo, double* X, unsigned ldX,
                             double* Y_, unsigned ldY_)
{
  assert(isLrM());

  if (k>0) {
    if (haarInfo==0) {
      unsigned ksum = bl_rank + k, LWORK = 10*ksum, l;
      unsigned mmin = MIN(n1, ksum), nmin = MIN(n2, ksum), amin = MIN(mmin, nmin);
      unsigned size = ksum*(n1+n2)+LWORK+mmin+(mmin+amin+1)*nmin+amin;

      double *tmp1=new double[size], *tmp2=tmp1+ksum*n1; // tmp2 = new V
      assert(tmp1!=NULL);

      // copy data(U) to tmp1
      const unsigned k1U = bl_rank*n1;
      if (k1U>0) blas::copy(k1U, data, tmp1);

      // copy U to tmp1+k1U
      for (l=0; l<k; ++l) blas::copy(n1, U+l*ldU, tmp1+k1U+l*n1);

      // copy data(V) to tmp+k1U+k2U
      const unsigned k1V = bl_rank*n2;
      if (k1V>0) blas::copy(k1V, data+k1U, tmp2);

      // copy V to tmp2+k1V
      for (l=0; l<k; ++l) blas::copy(n2, V+l*ldV, tmp2+k1V+l*n2);

      freedata();

      // QR-Zerlegung von A und B
      int INFO;
      double* const WORK = tmp2 + n2*ksum; // LWORK
      double* const tau1 = WORK + LWORK;   // mmin
      double* const tau2 = tau1 + mmin;    // nmin
      assert(tmp1!=NULL);assert(tmp2!=NULL);assert(tau1!=NULL);assert(tau2!=NULL);
      INFO = blas::geqrf(n1, ksum, tmp1, tau1, LWORK, WORK);
      assert(INFO==0);
      INFO = blas::geqrf(n2, ksum, tmp2, tau2, LWORK, WORK);
      assert(INFO==0);

      double* const R = tau2 + nmin;       // mmin*nmin

      // Berechne R1*R2^T
      blas::utrmmh(mmin, ksum, nmin, tmp1, n1, tmp2, n2, R);

      // eliminate small entries
      const double thresh = delta*blas::nrm2(mmin*nmin, R)/sqrt((float)mmin*nmin);
      for (unsigned i=0; i<mmin; ++i)
        for (unsigned j=0; j<nmin; ++j)
	  if (fabs(R[i+j*mmin])<thresh) R[i+j*mmin] = 0.0;

      double* const S = R + mmin*nmin;     // amin
      double* const VT = S + amin;         // amin*nmin

      // SVD von R1*R2^T
      INFO = blas::gesvd(mmin, nmin, R, S, VT, amin, LWORK, WORK);
      assert(INFO==0);

      unsigned kt = MIN(amin, kgoal);
      while (kt>0 && (S[kt-1]<=delta*S[0] || S[kt-1]<EPS0)) --kt;

      if (kt>0) {
        setrank(kt);

        // copy R to data und scale the columns with S
        unsigned j;
        for (j=0; j<kt; ++j) {
          unsigned i;
          for (i=0; i<mmin; ++i) data[i+j*n1] = S[j] * R[i+j*mmin];
          for (; i<n1; ++i) data[i+j*n1] = D_ZERO;
        }

        double* dataV = data + kt*n1;

        // copy V to dataV
        for (j=0; j<kt; ++j) {
          unsigned i;
          for (i=0; i<nmin; ++i) dataV[i+j*n2] = VT[i*amin+j];
          for (; i<n2; ++i) dataV[i+j*n2] = D_ZERO;
        }

        INFO = blas::ormqr(n1, kt, mmin, tmp1, tau1, data, LWORK, WORK);
        assert(INFO==0);

        INFO = blas::ormqr(n2, kt, nmin, tmp2, tau2, dataV, LWORK, WORK);
        assert(INFO==0);
      }

      delete [] tmp1;
    } else {
      //std::cout<<"bl_rank: "<<bl_rank<<std::endl;
      /*std::cout<<"U: "<<std::endl;
	for(unsigned i=0;i<n1*k;i++)
      	std::cout<<U[i]<<" ";
	std::cout<<"\nV: "<<std::endl;
	for(unsigned i=0;i<n2*k;i++)
      	std::cout<<V[i]<<" ";
	std::cout<<"\ndata: "<<std::endl;
	for(unsigned i=0;i<(n1+n2)*bl_rank;i++)
      	std::cout<<data[i]<<" ";*/
      addLowRankNB(delta, kgoal, n1, k, bl_rank, n2, U, ldU, V, ldV, haarInfo,
                   X, ldX, Y_, ldY_, bl_rank, data);
      //std::cout<<"Rang: "<<bl_rank;
    }
  }
}


// adds lwr to lwr and truncates to min(delta, kgoal), gives back remainder
// returns true if required rank satisfies rank*(n1+n2)<n1*n2
template<>
void mblock<double>::addtrll_rmnd(unsigned k, double* U, unsigned ldU,
                                  double* V, unsigned ldV, double delta,
                                  unsigned kgoal, unsigned& kR, double* &UR,
                                  double* &VR)
{
  assert(isLrM());
  kR = 0;
  UR = VR = NULL;

  if (k>0) {
    unsigned ksum = bl_rank + k, LWORK = 10*ksum, l;
    unsigned mmin=MIN(n1, ksum), nmin=MIN(n2, ksum), amin=MIN(mmin, nmin);
    unsigned size = ksum*(n1+n2)+LWORK+mmin+(mmin+amin+1)*nmin+amin;

    double *tmp1=new double[size], *tmp2=tmp1+ksum*n1; // tmp2 = new V
    assert(tmp1!=NULL);

    // copy data(U) to tmp1
    const unsigned k1U = bl_rank*n1;
    if (k1U>0) blas::copy(k1U, data, tmp1);

    // copy U to tmp1+k1U
    for (l=0; l<k; ++l) blas::copy(n1, U+l*ldU, tmp1+k1U+l*n1);

    // copy data(V) to tmp+k1U+k2U
    const unsigned k1V = bl_rank*n2;
    if (k1V>0) blas::copy(k1V, data+k1U, tmp2);

    // copy V to tmp2+k1V
    for (l=0; l<k; ++l) blas::copy(n2, V+l*ldV, tmp2+k1V+l*n2);

    freedata();

    // QR-Zerlegung von A und B
    int INFO;
    double* const WORK = tmp2 + n2*ksum; // LWORK
    double* const tau1 = WORK + LWORK;   // mmin
    double* const tau2 = tau1 + mmin;    // nmin

    INFO = blas::geqrf(n1, ksum, tmp1, tau1, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n2, ksum, tmp2, tau2, LWORK, WORK);
    assert(INFO==0);

    double* const R = tau2 + nmin;       // mmin*nmin

    // Berechne R1*R2^T
    blas::utrmmh(mmin, ksum, nmin, tmp1, n1, tmp2, n2, R);

    double* const S = R + mmin*nmin;     // amin
    double* const VT = S + amin;         // amin*nmin

    // SVD von R1*R2^T
    INFO = blas::gesvd(mmin, nmin, R, S, VT, amin, LWORK, WORK);
    assert(INFO==0);

    unsigned kt = MIN(amin, kgoal);
    while (kt>0 && (S[kt-1]<=delta*S[0] || S[kt-1]<EPS0)) --kt;

    if (kt>0) {
      setrank(kt);

      // copy R to data und scale the columns with S, copy V to dataV
      double* dataV = data + kt*n1;
      for (unsigned j=0; j<kt; ++j) {
        unsigned i;
        const double sS = sqrt(S[j]);
        for (i=0; i<mmin; ++i) data[i+j*n1] = sS * R[i+j*mmin];
        for (; i<n1; ++i) data[i+j*n1] = D_ZERO;

        for (i=0; i<nmin; ++i) dataV[i+j*n2] = sS * VT[i*amin+j];
        for (; i<n2; ++i) dataV[i+j*n2] = D_ZERO;
      }

      INFO = blas::ormqr(n1, kt, mmin, tmp1, tau1, data, LWORK, WORK);
      assert(INFO==0);

      INFO = blas::ormqr(n2, kt, nmin, tmp2, tau2, dataV, LWORK, WORK);
      assert(INFO==0);
    }

    // compute remainder

    unsigned kl = amin;
    while (kl>kt && (S[kl-1]<=1e-16*S[0] || S[kl-1]<EPS0)) --kl;
    kR = kl - kt;

    if (kR>0) {
      UR = new double[n1*kR];
      VR = new double[n2*kR];
      assert(UR!=NULL && VR!=NULL);

      // copy R to UR und scale the columns with S, copy V to VR
      for (unsigned j=0; j<kR; ++j) {
        unsigned i;
        double sS = sqrt(S[j+kt]);
        for (i=0; i<mmin; ++i) UR[i+j*n1] = sS * R[i+(j+kt)*mmin];
        for (; i<n1; ++i) UR[i+j*n1] = D_ZERO;

        for (i=0; i<nmin; ++i) VR[i+j*n2] = sS * VT[i*amin+j+kt];
        for (; i<n2; ++i) VR[i+j*n2] = D_ZERO;
      }

      INFO = blas::ormqr(n1, kR, mmin, tmp1, tau1, UR, LWORK, WORK);
      assert(INFO==0);

      INFO = blas::ormqr(n2, kR, nmin, tmp2, tau2, VR, LWORK, WORK);
      assert(INFO==0);

    }

    delete [] tmp1;
  }
}

template<>
void mblock<double>::addGeM_toHeM(double* const A, const unsigned ldA)
{
  assert(ldA>=n1 && isGeM() && (isHeM()||isSyM()) && n1==n2);
  double *p = data;
  for (unsigned j=0; j<n2; ++j)
    for (unsigned i=0; i<=j; ++i) *p++ += A[i+j*ldA];
}


// admissible block: add a dense block and truncate to min(eps_adm, MAX_RANK)
template<>
void mblock<double>::addGeM(double* A, unsigned ldA, double eps,
			    unsigned rankmax, contLowLevel<double>* haarInfo,
			    double* X, unsigned ldX,
			    double* Y_, unsigned ldY_)
{
  if (isLrM()) {
    const unsigned min12 = MIN(n1, n2);
    const unsigned mtn = n1*n2;
    const unsigned lwork = 5*(n1+n2);
    double* const tmp = new double[mtn+min12*(n2+1)+lwork];
    assert(tmp!=NULL);

    for (unsigned j=0; j<n2; ++j)
      blas::copy(n1, A+j*ldA, tmp+j*n1);

    if (bl_rank>0)     // convert rank-k matrix to dense matrix, add it to tmp
      blas::gemmha(n1, bl_rank, n2, data, n1, data+bl_rank*n1, n2, tmp, n1);

    freedata();  // tmp[1..mtn] contains the sum

    // compute SVD
    double* const S = tmp + mtn;            // min12
    double* const VT = S + min12;           // min12*n2
    double* const work = VT + min12*n2;     // lwork
    int INFO = blas::gesvd(n1, n2, tmp, S, VT, min12, lwork, work);
    assert(INFO==0);

    if (haarInfo==NULL) {
      unsigned kt = MIN(min12, rankmax);
      while (kt>0 && (S[kt-1]<=eps*S[0] || S[kt-1]<EPS0)) --kt;

      if (kt>0) {                      // generate rank-k approximant
        setrank(kt);

        double *p = data, *q = tmp;

        for (unsigned l=0; l<kt; ++l) {
          for (unsigned i=0; i<n1; ++i ) *p++ = S[l] * *q++;
          blas::copy(n2, VT+l, min12, data+kt*n1+l*n2, N_ONE);
        }
      }
    } else {
      double* const V = new double[min12*n2];
      for (unsigned i=0; i<n2; i++)
        for (unsigned j=0; j<min12; j++)
          V[i+j*n2] = VT[j+i*min12];
      for (unsigned i=0; i<min12; i++)
        blas::scal(n1,S[i],&tmp[i*n1]);
      unsigned khat;
      createLowRankMatHouseholderNB(eps, rankmax, n1, min12, n2, tmp, V, haarInfo,
                                    X, ldX, Y_, ldY_, khat, data);

      bl_rank = khat;
      delete [] V;
    }
    delete [] tmp;
  } else { // is dense
    if (isHeM() || isSyM()) addGeM_toHeM(A, ldA);
    else addGeM_toGeM(A, ldA);
  }
}

// unify (U1 V1^T, U2 V2^T) to U V^T with rank at most kgoal
template<>
void mblock<double>::unify_cols_LrMLrM(double delta, unsigned kgoal,
				       const mblock<double>& mbl1,
				       const mblock<double>& mbl2,
				       contLowLevel<double>* haarInfo, 
				       double* X, unsigned ldX,
				       double* Y_, unsigned ldY_)
{
  assert(mbl1.isLrM() && mbl2.isLrM());
  const unsigned k1 = mbl1.bl_rank, k2 = mbl2.bl_rank, ksum = k1 + k2;

  if (ksum>0) {
    const unsigned n2a = mbl1.n2, n2b = mbl2.n2, kU = MIN(n1, ksum);
    const unsigned size=kU*ksum, LWORK=10*ksum;
    const unsigned sizeU1=k1*n1, sizeU2=k2*n1, sizeV1=k1*n2a, sizeV2=k2*n2b;
    double* const tmp=new double[(n1+1)*ksum+2*(size+kU)+LWORK+sizeV1+sizeV2];
    assert(tmp!=NULL);

    double* const R1 = tmp;                                    // n1*k1
    double* const R2 = tmp + sizeU1;                           // n1*k2
    double* const tau = R2 + sizeU2;                           // kU
    double* const tau1 = tau + kU;                             // k1
    double* const tau2 = tau1 + k1;                            // k2
    double* const M = tau2 + k2;                               // size
    double* const S = M + size;                                // kU
    double* const VT = S + kU;                                 // size
    double* const WORK = VT + size;                            // LWORK
    double* const V1 = WORK + LWORK;                           // sizeV1
    double* const V2 = V1 + sizeV1;                            // sizeV2

    // copy U1 and U2 to tmp / V1 and V2
    blas::copy(sizeU1, mbl1.data, R1);
    blas::copy(sizeV1, mbl1.data+sizeU1, V1);
    blas::copy(sizeU2, mbl2.data, R2);
    blas::copy(sizeV2, mbl2.data+sizeU2, V2);

    int INFO;
    // QR decompositions : (U1,U2) = Q (R1,R2), V1 = Q1 T1, V2 = Q2 T2
    INFO = blas::geqrf(n1, ksum, tmp, tau, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n2a, k1, V1, tau1, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n2b, k2, V2, tau2, LWORK, WORK);
    assert(INFO==0);

    // compute M=(R1 T1^T, R2 T2^T)
    blas::utrmmh(kU, k1, k1, R1, n1, V1, n2a, M);

    for (unsigned j=0; j<k2; ++j) {
      for (unsigned i=0; i<kU; ++i) {
        double d = D_ZERO;
        for (unsigned l=MAX(i,j+k1)-k1; l<k2; ++l) d += R2[i+l*n1]*V2[j+l*n2b];
        M[i+kU*(j+k1)] = d;
      }
    }

    // SVD von M
    INFO = blas::gesvd(kU, ksum, M, S, VT, kU, LWORK, WORK);
    assert(INFO==0);

    if (haarInfo==NULL) {
      unsigned kt = MIN(kU, kgoal);
      while (kt>0 && (S[kt-1]<=delta*S[0] || S[kt-1]<EPS0)) --kt;

      bl_rank = kt;

      if (kt>0) {
        data = new double[kt*(n1+n2)];
        assert(data!=NULL);

        INFO = blas::orgqr(n1, kU, tmp, tau, LWORK, WORK);
        assert(INFO==0);
        for (unsigned l=0; l<kt; ++l) blas::scal(kU, S[l], M+l*kU);
        blas::gemm(n1, kU, kt, D_ONE, tmp, n1, M, kU, data, n1);

        INFO = blas::orgqr(n2a, k1, V1, tau1, LWORK, WORK);
        assert(INFO==0);
        blas::gemmh(n2a, k1, kt, D_ONE, V1, n2a, VT, kU, data+kt*n1, n2);

        INFO = blas::orgqr(n2b, k2, V2, tau2, LWORK, WORK);
        assert(INFO==0);
        blas::gemmh(n2b, k2, kt, D_ONE, V2, n2b, VT+k1*kU, kU,
                    data+kt*n1+n2a, n2);
      }
    } else {
      double* const V = new double[size];
      for (unsigned i=0; i<kU; i++) {
        for (unsigned j=0; j<ksum; j++) {
          V[j+i*ksum] = VT[i+j*kU];
        }
      }
      INFO = blas::orgqr(n1, kU, tmp, tau, LWORK, WORK);
      assert(INFO==0);
      INFO = blas::orgqr(n2a, k1, V1, tau1, LWORK, WORK);
      assert(INFO==0);
      INFO = blas::orgqr(n2b, k2, V2, tau2, LWORK, WORK);
      assert(INFO==0);
      
      for (unsigned i=0; i<kU; i++) blas::scal(kU,S[i],&M[i*kU]);
      double* datatemp;
      unsigned cols = haarInfo->col;
      double* Xnew = new double[cols*ksum];
      double* Y_new = new double[cols*kU];
      blas::gemhm(n2a,k1,cols,D_ONE,V1,n2a,X,ldX,Xnew,ksum);
      blas::gemhm(n2b,k2,cols,D_ONE,V2,n2b,X+n2a,ldX,Xnew+k1,ksum);
      blas::gemhm(n1,kU,cols,D_ONE,tmp,n1,Y_,ldY_,Y_new,kU);
      createLowRankMatHouseholderNB(delta, kgoal, kU, kU, ksum, M, V, haarInfo,
                                    Xnew, ksum, Y_new, kU, bl_rank, datatemp);

      data = new double[bl_rank*(n1+n2)];
      assert(data!=NULL);
      blas::gemm(n1, kU, bl_rank, D_ONE, tmp, n1, datatemp, kU, data, n1);

      blas::gemm(n2a, k1, bl_rank, D_ONE, V1, n2a, &datatemp[bl_rank*kU], ksum,
                 data+bl_rank*n1, n2);

      blas::gemm(n2b, k2, bl_rank, D_ONE, V2, n2b, &datatemp[bl_rank*kU+k1],
                 ksum, data+bl_rank*n1+n2a, n2);

      delete [] datatemp;
      delete [] V;
      delete [] Xnew;
      delete [] Y_new;
    }

    delete [] tmp;
  }
}

// unify mblocks (A, B) to U V^T with rank at most kgoal
template<>
void mblock<double>::unify_cols(double delta, unsigned kgoal,
                                const mblock<double>& mbl1,
                                const mblock<double>& mbl2,
                                contLowLevel<double>* haarInfo,
				double* X, unsigned ldX,
                                double* Y_, unsigned ldY_)
{
  assert(mbl1.n1==mbl2.n1 && n1==mbl1.n1 && n2==mbl1.n2+mbl2.n2);
  freedata();

  if (mbl1.isLrM() && mbl2.isLrM())
    unify_cols_LrMLrM(delta, kgoal, mbl1, mbl2, haarInfo, X, ldX, Y_, ldY_);
  else {
    unsigned n2a=mbl1.n2, n2b=mbl2.n2, nmin=MIN(n1, n2), LWORK=5*(n1+n2);
    double *tmp = new double[n1*n2+nmin*(1+n2)+LWORK];
    assert(tmp!=NULL);

    // convert (A,B) to dense matrix tmp
    if (mbl1.isLrM())
      blas::gemmh(n1, mbl1.bl_rank, n2a, D_ONE, mbl1.data, n1,
                  mbl1.data+mbl1.bl_rank*n1, n2a, tmp, n1);
    else
      blas::copy(n1*n2a, mbl1.data, tmp);

    if (mbl2.isLrM())
      blas::gemmh(n1, mbl2.bl_rank, n2b, D_ONE, mbl2.data, n1,
                  mbl2.data+mbl2.bl_rank*n1, n2b, tmp+n1*n2a, n1);
    else
      blas::copy(n1*n2b, mbl2.data, tmp+n1*n2a);

    // SVD von tmp

    double *S = tmp+n1*n2, *VT = S+nmin, *WORK = VT+nmin*n2;
    int INFO = blas::gesvd(n1, n2, tmp, S, VT, nmin, LWORK, WORK);
    assert(INFO==0);

    if (haarInfo==NULL) {
      unsigned kt = MIN(nmin, kgoal);
      while (kt>0 && (S[kt-1]<=delta*S[0] || S[kt-1]<EPS0)) --kt;

      bl_rank = kt;

      if (kt>0) {
        data = new double[kt*(n1+n2)];
        assert(data!=NULL);

        blas::copy(n1*kt, tmp, data);

        for (unsigned l=0; l<kt; ++l)
          for (unsigned j=0; j<n2; ++j)
            data[kt*n1+l*n2+j] = S[l] * VT[j*nmin+l];
      }
    } else {
      double* const V = new double[nmin*n2];
      for (unsigned i=0; i<n2; i++) {
        for (unsigned j=0; j<nmin; j++) {
          V[i+j*n2] = VT[j+i*nmin];
        }
      }
      for (unsigned i=0; i<nmin; i++)
        blas::scal(n1,S[i],&tmp[i*n1]);
      createLowRankMatHouseholderNB(delta, kgoal, n1, nmin, n2, tmp, V, haarInfo,
                                    X, ldX, Y_, ldY_, bl_rank, data);
      delete [] V;
    }
    delete [] tmp;
  }
}


// unify (U1 V1^T \\ U2 V2^T) to U V^T with rank at most kgoal
template<>
void mblock<double>::unify_rows_LrMLrM(double delta, unsigned kgoal,
				       const mblock<double>& mbl1,
				       const mblock<double>& mbl2,
				       contLowLevel<double>* haarInfo,
				       double* X, unsigned ldX,
				       double* Y_, unsigned ldY_)
{
  assert(mbl1.isLrM() && mbl2.isLrM());
  const unsigned k1 = mbl1.bl_rank, k2 = mbl2.bl_rank, ksum = k1 + k2;

  if (ksum>0) {
    const unsigned n1a = mbl1.n1, n1b = mbl2.n1, kV = MIN(n2, ksum);
    const unsigned size=kV*ksum, LWORK=10*ksum;
    const unsigned sizeU1=k1*n1a, sizeU2=k2*n1b, sizeV1=k1*n2, sizeV2=k2*n2;
    double* const tmp=new double[(n2+1)*ksum+2*(size+kV)+LWORK+sizeU1+sizeU2];
    assert(tmp!=NULL);

    double* const R1 = tmp;                                    // n2*k1
    double* const R2 = tmp + sizeV1;                           // n2*k2
    double* const tau = R2 + sizeV2;                           // kV
    double* const tau1 = tau + kV;                             // k1
    double* const tau2 = tau1 + k1;                            // k2
    double* const M = tau2 + k2;                               // size
    double* const S = M + size;                                // kV
    double* const VT = S + kV;                                 // size
    double* const WORK = VT + size;                            // LWORK
    double* const U1 = WORK + LWORK;                           // sizeU1
    double* const U2 = U1 + sizeU1;                            // sizeU2

    // copy V1 and V2 to tmp / U1 and U2
    blas::copy(sizeU1, mbl1.data, U1);
    blas::copy(sizeV1, mbl1.data+sizeU1, R1);
    blas::copy(sizeU2, mbl2.data, U2);
    blas::copy(sizeV2, mbl2.data+sizeU2, R2);

    int INFO;
    // QR decompositions : (V1,V2) = Q (R1,R2), U1 = Q1 T1, U2 = Q2 T2
    INFO = blas::geqrf(n2, ksum, R1, tau, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n1a, k1, U1, tau1, LWORK, WORK);
    assert(INFO==0);
    INFO = blas::geqrf(n1b, k2, U2, tau2, LWORK, WORK);
    assert(INFO==0);

    // compute M=(T1 R1^T \\ T2 R2^T)
    unsigned j;
    for (j=0; j<k1; ++j) {
      unsigned i;
      for (i=0; i<k1; ++i) {
        double d = D_ZERO;
        for (unsigned l=MAX(i,j); l<k1; ++l) d += U1[i+l*n1a]*R1[j+l*n2];
        M[i+ksum*j] = d;
      }
      for (i=0; i<k2; ++i) {
        double d = D_ZERO;
        for (unsigned l=i; l<k2; ++l) d += U2[i+l*n1b]*R2[j+l*n2];
        M[i+k1+ksum*j] = d;
      }
    }

    for (j=k1; j<kV; ++j) {
      blas::setzero(k1, M+ksum*j);
      for (unsigned i=0; i<k2; ++i) {
        double d = D_ZERO;
        for (unsigned l=MAX(i,j-k1); l<k2; ++l) d += U2[i+l*n1b]*R2[j+l*n2];
        M[i+k1+ksum*j] = d;
      }
    }

    // SVD von M
    INFO = blas::gesvd(ksum, kV, M, S, VT, kV, LWORK, WORK);
    assert(INFO==0);

    if (haarInfo==NULL) {
      unsigned kt = MIN(kV, kgoal);
      while (kt>0 && (S[kt-1]<=delta*S[0] || S[kt-1]<EPS0)) --kt;

      bl_rank = kt;
      if (kt>0) {
        for (unsigned l=0; l<kt; ++l) blas::scal(ksum, S[l], M+l*ksum);
        data = new double[kt*(n1+n2)];
        assert(data!=NULL);

        INFO = blas::orgqr(n1a, k1, U1, tau1, LWORK, WORK);
        assert(INFO==0);
        blas::gemm(n1a, k1, kt, D_ONE, U1, n1a, M, ksum, data, n1);

        INFO = blas::orgqr(n1b, k2, U2, tau2, LWORK, WORK);
        assert(INFO==0);
        blas::gemm(n1b, k2, kt, D_ONE, U2, n1b, M+k1, ksum, data+n1a, n1);

        INFO = blas::orgqr(n2, kV, tmp, tau, LWORK, WORK);
        assert(INFO==0);
        blas::gemmh(n2, kV, kt, D_ONE, tmp, n2, VT, kV, data+kt*n1, n2);
      }
    } else {
      double* const V = new double[size];
      for (unsigned i=0; i<kV; i++) {
        for (unsigned j=0; j<kV; j++) {
          V[j+i*kV] = VT[i+j*kV];
        }
      }
      INFO = blas::orgqr(n1a, k1, U1, tau1, LWORK, WORK);
      assert(INFO==0);
      INFO = blas::orgqr(n1b, k2, U2, tau2, LWORK, WORK);
      assert(INFO==0);
      INFO = blas::orgqr(n2, kV, tmp, tau, LWORK, WORK);
      assert(INFO==0);

      for (unsigned i=0; i<kV; i++) blas::scal(ksum,S[i],&M[i*ksum]);
      double* datatemp;
      const unsigned cols = haarInfo->col;
      double* Xnew = new double[cols*kV];
      double* Y_new = new double[cols*ksum];
      blas::gemhm(n2,kV,cols,D_ONE,tmp,n2,X,ldX,Xnew,kV);
      blas::gemhm(n1a,k1,cols,D_ONE,U1,n1a,Y_,ldY_,Y_new,ksum);
      blas::gemhm(n1b,k2,cols,D_ONE,U2,n1b,Y_+n1a,ldY_,Y_new+k1,ksum);
      createLowRankMatHouseholderNB(delta, kgoal, ksum, kV, kV, M, V, haarInfo,
                                    Xnew, kV, Y_new, ksum, bl_rank, datatemp);

      data = new double[bl_rank*(n1+n2)];
      assert(data!=NULL);
      blas::gemm(n1a, k1, bl_rank, D_ONE, U1, n1a, datatemp, ksum, data, n1);

      blas::gemm(n1b, k2, bl_rank, D_ONE, U2, n1b, datatemp+k1, ksum,
                 data+n1a, n1);

      blas::gemm(n2, kV, bl_rank, D_ONE, tmp, n2, datatemp+bl_rank*ksum,
                 kV, data+bl_rank*n1, n2);

      delete [] datatemp;
      delete [] V;
      delete [] Xnew;
      delete [] Y_new;
    }

    delete [] tmp;
  }
}


// unify mblocks (A \\ B) to U V^T with rank at most kgoal
template<>
void mblock<double>::unify_rows(double delta, unsigned kgoal,
                                const mblock<double>& mbl1,
                                const mblock<double>& mbl2,
                                contLowLevel<double>* haarInfo,
				double* X, unsigned ldX,
                                double* Y_, unsigned ldY_)
{
  assert(mbl1.n2==mbl2.n2 && n2==mbl1.n2 && n1==mbl1.n1+mbl2.n1);
  freedata();

  if (mbl1.isLrM() && mbl2.isLrM())
    unify_rows_LrMLrM(delta, kgoal, mbl1, mbl2, haarInfo, X, ldX, Y_, ldY_);
  else {
    unsigned n1a=mbl1.n1, n1b=mbl2.n1, nmin=MIN(n1, n2), LWORK=5*(n1+n2);
    double *tmp = new double[n1*n2+nmin*(1+n2)+LWORK];
    assert(tmp!=NULL);

    // convert (A \\ B) to dense matrix tmp
    if (mbl1.isLrM())
      blas::gemmh(n1a, mbl1.bl_rank, n2, D_ONE, mbl1.data, n1a,
                  mbl1.data+mbl1.bl_rank*n1a, n2, tmp, n1);
    else
      for (unsigned l=0; l<n2; ++l)
        blas::copy(n1a, mbl1.data+n1a*l, tmp+n1*l);

    if (mbl2.isLrM())
      blas::gemmh(n1b, mbl2.bl_rank, n2, D_ONE, mbl2.data, n1b,
                  mbl2.data+mbl2.bl_rank*n1b, n2, tmp+n1a, n1);
    else
      for (unsigned l=0; l<n2; ++l)
        blas::copy(n1b, mbl2.data+n1b*l, tmp+n1*l+n1a);


    // SVD von tmp
    double *S = tmp+n1*n2, *VT = S+nmin, *WORK = VT+nmin*n2;
    int INFO = blas::gesvd(n1, n2, tmp, S, VT, nmin, LWORK, WORK);
    assert(INFO==0);

    if (haarInfo==NULL) {
      unsigned kt = MIN(nmin, kgoal);
      while (kt>0 && (S[kt-1]<=delta*S[0] || S[kt-1]<EPS0)) --kt;

      bl_rank = kt;

      if (kt>0) {
        data = new double[kt*(n1+n2)];
        assert(data!=NULL);

        blas::copy(n1*kt, tmp, data);

        for (unsigned l=0; l<kt; ++l)
          for (unsigned j=0; j<n2; ++j)
            data[kt*n1+l*n2+j] = S[l] * VT[j*nmin+l];
      }
    } else {
      double* const V = new double[nmin*n2];
      for (unsigned i=0; i<n2; i++)
        for (unsigned j=0; j<nmin; j++)
          V[i+j*n2] = VT[j+i*nmin];
      for (unsigned i=0; i<nmin; i++)
        blas::scal(n1,S[i],&tmp[i*n1]);
      createLowRankMatHouseholderNB(delta, kgoal, n1, nmin, n2, tmp, V, haarInfo,
                                    X, ldX, Y_, ldY_, bl_rank, data);
      delete [] V;
    }
    delete [] tmp;
  }
}


template<>
int mblock<double>::decomp_LU(mblock<double>* L, mblock<double> *U)
// destroys data
{
  assert(n1==n2);

  unsigned j, *perm = new unsigned[2*n1], *ipiv = perm + n1;
  for (j=0; j<n1; ++j) perm[j] = j;

  int inf = blas::getrf(n1, data, ipiv);

  L->setLtM();
  U->setUtM();

  double *pA = data, *pL = L->data, *pU = U->data;

  for (j=0; j<n1; ++j) {

    unsigned itmp = perm[j];
    perm[j] = perm[ipiv[j]-1];
    perm[ipiv[j]-1] = itmp;

    unsigned i;
    for (i=0; i<=j; ++i) *pU++ = *pA++;
    *pL++ = perm[j];
    for (; i<n1; ++i) *pL++ = *pA++;
  }

  delete [] perm;
  return inf;
}


// multiply lower triangular by vector: y += d PL x
template<>
void mblock<double>::mltaLtMVec(double d, double* x, double* y) const
{
  assert(isGeM() && isLtM() && n1==n2);
  unsigned *ip = new unsigned[n1];
  double* z = new double[n1];
  blas::setzero(n1, z);

  double *p = data;
  for (unsigned j=0; j<n2; ++j) {
    ip[j] = (unsigned) *p++;    // diag element is permutation info
    double e = d * x[j];
    z[j] += e;
    for (unsigned i=j+1; i<n1; ++i) z[i] += e * *p++;
  }

  for (unsigned i=0; i<n1; ++i) y[ip[i]] += z[i];

  delete [] z;
  delete [] ip;
}


// multiply upper triangular by vector: y += d A x (A is utr)
template<>
void mblock<double>::mltaUtMVec(double d, double* x, double* y) const
{
  assert(isGeM() && isUtM() && n1==n2);

  double *p = data;
  for (unsigned j=0; j<n2; ++j) {
    double e = d * x[j];
    for (unsigned i=0; i<=j; ++i) y[i] += e * *p++;
  }
}


// multiply lower triangular by vector: y += d (PL)^T x = d L^T P^{-1} x
template<>
void mblock<double>::mltaLtMhVec(double d, double* x, double* y) const
{
  assert(isGeM() && isLtM() && n1==n2);
  double *p = data, *z = new double[n1];

  for (unsigned j=0; j<n2; ++j) {
    unsigned ip = (unsigned) *p;    // diag element is permutation info
    z[j] = x[ip];
    p += n1-j;
  }

  blas::ltrphv(n1, data, z);
  blas::axpy(n1, d, z, y);

  delete [] z;
}


// multiply upper triangular by vector: y += d A^T x (A is utr)
template<>
void mblock<double>::mltaUtMhVec(double d, double* x, double* y) const
{
  assert(isGeM() && isUtM());
  double *p = data;
  for (unsigned i=0; i<n1; ++i) {
    double e = D_ZERO;
    for (unsigned j=0; j<=i; ++j) e += x[j] * *p++;
    y[i] += d * e;
  }
}


template<>
void mblock<double>::ltr_solve(unsigned m, double* B, unsigned ldB) const
{
  assert(isGeM() && isLtM() && n1==n2);
  double *p = data, *Z = new double[n1*m];
  unsigned j, *ip = new unsigned[n1];

  for (j=0; j<n2; ++j) {
    ip[j] = (unsigned) *p;    // diag element is permutation info
    p += n1-j;
  }

  for (j=0; j<m; ++j)
    for (unsigned i=0; i<n1; ++i) Z[i+j*n1] = B[ip[i]+j*ldB];

  lapack::ltrs(n1, data, m, Z, n1);

  for (j=0; j<m; ++j) blas::copy(n1, Z+j*n1, B+j*ldB);

  delete [] ip;
  delete [] Z;
}


// solves (PL)^T X = L^T P^{-1} X = B for X
// L is unit lower triangular, X is stored in B
template<>
void mblock<double>::ltrh_solve(unsigned m, double* B, unsigned ldB) const
{
  assert(isGeM() && isLtM() && n1==n2);
  double *p = data, *Z = new double[n1*m];
  unsigned j, *ip = new unsigned[n1];

  for (j=0; j<n2; ++j) {
    ip[j] = (unsigned) *p;    // diag element is permutation info
    p += n1-j;
  }

  for (j=0; j<m; ++j) blas::copy(n1, B+j*ldB, Z+j*n1);

  lapack::ltrhs(n1, data, m, Z, n1);

  for (j=0; j<m; ++j)
    for (unsigned i=0; i<n1; ++i) B[ip[i]+j*ldB] = Z[i+j*n1];

  delete [] ip;
  delete [] Z;
}

// solves X U = B for X, U is upper triangular
template<>
void mblock<double>::utr_solve_left(unsigned m, double* B, unsigned ldB,
                                    double* X, unsigned ldX) const
{
  assert(isGeM() && isUtM());
  double *p = data;

  for (unsigned j=0; j<n1; ++j) {
    blas::copy(m, B+j*ldB, X+j*ldX);
    for (unsigned l=0; l<j; ++l) {
      double d = -*p++;
      blas::axpy(m, d, X+l*ldX, X+j*ldX);
    }
    double e = 1.0 / *p++;
    blas::scal(m, e, X+j*ldX);
  }
}

// add lwr to herm. dense
template<>
void mblock<double>::addLrM_toHeM(unsigned mult, unsigned k, double* U,
				  unsigned ldU, double* V, unsigned ldV)
{
  assert(isGeM() && (isHeM()||isSyM()) && n1==n2);
  for (unsigned l=0; l<k; ++l)
    for (unsigned j=0; j<n1; ++j)
      blas::axpy(j+1, mult*V[j+l*ldV], U+ldU*l, data+j*(j+1)/2);
}

// add lwr to herm. dense
template<>
void mblock<double>::addLrM_toUtM(unsigned mult, unsigned k, double* U,
				  unsigned ldU, double* V, unsigned ldV)
{
  assert(isGeM() && isUtM());
  for (unsigned l=0; l<k; ++l)
    for (unsigned j=0; j<n1; ++j)
      blas::axpy(j+1, mult*V[j+l*ldV], U+ldU*l, data+j*(j+1)/2);
}

// add lwr to mblock<double> and truncate to min(eps, kgoal)
template<>
void mblock<double>::addLrM(unsigned k, double* U, unsigned ldU, double* V,
			    unsigned ldV, double eps, unsigned kgoal,
			    contLowLevel<double>* haarInfo,
			    double* X, unsigned ldX,
			    double* Y_, unsigned ldY_)
{
  if (isLrM()) { // cannot be symmetric since it is lwr
    addtrll(k, U, ldU, V, ldV, eps, kgoal, haarInfo, X, ldX, Y_, ldY_);
    if (bl_rank*(n1+n2)>n1*n2) convLrM_toGeM();
  } else {
    if (isHeM() || isSyM()) addLrM_toHeM(1, k, U, ldU, V, ldV);
    else if (isUtM()) addLrM_toUtM(1, k, U, ldU, V, ldV);
    else addLrM_toGeM(k, U, ldU, V, ldV);
  }
}

// add lwr to mblock<double> and
template<>
void mblock<double>::addLrM_Exact(unsigned k, double* U, unsigned ldU,
				  double* V, unsigned ldV)
{
  if (isLrM()) { // cannot be symmetric since it is lwr
    append(k, U, ldU, V, ldV);
    if (bl_rank*(n1+n2)>n1*n2) convLrM_toGeM();
  } else {
    if (isHeM() || isSyM()) addLrM_toHeM(1, k, U, ldU, V, ldV);
    else addLrM_toGeM(k, U, ldU, V, ldV);
  }
}

// add lwr to mblock<double> and truncate to min(eps, kgoal)
template<>
void mblock<double>::addLrM_rmnd(unsigned k, double* U, unsigned ldU,
				 double* V, unsigned ldV, double eps,
				 unsigned kgoal, unsigned& kR, double* &UR,
				 double* &VR)
{
  if (isLrM()) {
    addtrll_rmnd(k, U, ldU, V, ldV, eps, kgoal, kR, UR, VR);
    if (bl_rank*(n1+n2)>n1*n2) {
      convLrM_toGeM();
      addLrM_toGeM(kR, UR, n1, VR, n2);
      delete [] UR;
      delete [] VR;
      kR = 0;
      UR = VR = NULL;
    }
  }
  else {
    if (isHeM() || isSyM()) addLrM_toHeM(1, k, U, ldU, V, ldV);
    else addLrM_toGeM(k, U, ldU, V, ldV);
    kR = 0;
    UR = VR = NULL;
  }
}
