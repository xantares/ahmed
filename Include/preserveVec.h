/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef LOWRANKMAT_H
#define LOWRANKMAT_H

#include "blas.h"
#include "vectorIter.h"

#define MVEC

class blcluster;
template<class T> class contBasis;

struct contConst
{
  // higher accuracy for non-admissible blocks
  const double c_adm;
  // threshold for vectors to be preserved
  const double eps2;
  
contConst(double c_admn, double eps2n)
:c_adm(c_admn), eps2(eps2n){}
};

// container for information passed to the low level routines
template<class T>
struct contLowLevel{
  const bool isadm;
  const unsigned col;
  contConst& constants;
  
contLowLevel(contBasis<T>* haar)
: col(haar->getcols()), isadm(haar->getbl()->isadm()), constants(haar->getconstants()){}
};

template<class T>
static void errorTest(T* A, T* M2, double* SR2, T* FT2,
			T* M, double* SR, T* FT, 
			double* e)
{
  T* t = new T[4];
  
  t[0] = M2[0]*SR2[0]*FT2[0] + M2[2]*SR2[1]*FT2[1];
  t[1] = M2[1]*SR2[0]*FT2[0] + M2[3]*SR2[1]*FT2[1];
  t[2] = M2[0]*SR2[0]*FT2[2] + M2[2]*SR2[1]*FT2[3];
  t[3] = M2[1]*SR2[0]*FT2[2] + M2[3]*SR2[1]*FT2[3];
  
  double nrm = blas::nrm2(4, t);
  
  t[0] -= M[0]*SR[0]*FT[0] + M[2]*SR[1]*FT[1];
  t[1] -= M[1]*SR[0]*FT[0] + M[3]*SR[1]*FT[1];
  t[2] -= M[0]*SR[0]*FT[2] + M[2]*SR[1]*FT[3];
  t[3] -= M[1]*SR[0]*FT[2] + M[3]*SR[1]*FT[3];
  
  double error = blas::nrm2(4,t);
  delete [] t;
  e[0] = error/nrm;

  A[0] -= M[0]*SR[0]*FT[0] + M[2]*SR[1]*FT[1];
  A[1] -= M[1]*SR[0]*FT[0] + M[3]*SR[1]*FT[1];
  A[2] -= M[0]*SR[0]*FT[2] + M[2]*SR[1]*FT[3];
  A[3] -= M[1]*SR[0]*FT[2] + M[3]*SR[1]*FT[3];

  error = blas::nrm2(4,A);
  e[1] = error/nrm;
}

// svd of a 2x2 matrix using Jacobi rotation
// U is stored in A
template<class T> static
void gesvd(T* A, double* S, T* FT, const unsigned nwk, T* wk)
{
  assert(nwk>8);
    
  // compute A^HA
  T* AHA = wk;
  AHA[0] = SQR(A[0]) + SQR(A[1]);
  AHA[1] = A[0]*A[2] + A[1]*A[3];
  AHA[2] = SQR(A[2]) + SQR(A[3]);

  // compute Jacobi transformation AHA
  wk[3] = (AHA[2]-AHA[0])/(2.0*AHA[1]);
  wk[4] = SIGN(wk[3])/(abs(wk[3])+ sqrt(1.0+SQR(wk[3])));
  // compute singular values
  S[0] = sqrt(AHA[0]-wk[4]*AHA[1]);
  S[1] = sqrt(AHA[2]+wk[4]*AHA[1]);

  if(S[0]>S[1]){
    //std::cout<<"upper:"<<SIGN(wk[3])<<" "<<SIGN(AHA[1]);
    FT[0] = sqrt(1.0/(1.0 + SQR(wk[4])));
    FT[3] = FT[0];
    FT[1] = wk[4]*FT[0];
    FT[2] = -FT[1];
  }else{
    //std::cout<<"lower:";
    FT[2] = SIGN(wk[4])*sqrt(1.0/(1.0 + SQR(wk[4])));
    FT[1] = -FT[2];
    FT[3] = wk[4]*FT[2];
    FT[0] = FT[3];
    
    double St = Re(S[0]);
    S[0] = S[1];
    S[1] = St;
  }
    
  // compute AA^H
  T* AAH = wk;
  AAH[0] = SQR(A[0]) + SQR(A[2]);
  AAH[1] = A[0]*A[1] + A[2]*A[3];
  AAH[2] = SQR(A[1]) + SQR(A[3]);
    
  // compute Jacobi transformation AAH
  wk[3] = (AAH[2]-AAH[0])/(2.0*AAH[1]);
  wk[4] = SIGN(wk[3])/(abs(wk[3]) + sqrt(1.0+SQR(wk[3])));

  T S1 = AAH[0]-wk[4]*AAH[1];
  T S2 = AAH[2]+wk[4]*AAH[1];
  
  wk[5] = A[0];
  wk[6] = A[1];
  wk[7] = A[2];
  wk[8] = A[3];

  if(S1>S2){
    //std::cout<<"upper";
    A[0] = sqrt(1.0/(1.0 + SQR(wk[4])));
    A[3] = A[0];
    A[2] = wk[4]*A[0];
    A[1] = -A[2];
  }else{
    //std::cout<<"lower";
    A[1] = SIGN(wk[4])*sqrt(1.0/(1.0 + SQR(wk[4])));
    A[2] = -A[1];
    A[0] = wk[4]*A[1];
    A[3] = A[0];
  }

  S1 = SIGN(FT[0]*(A[0]*wk[5] + A[1]*wk[6]) + FT[2]*(A[0]*wk[7] + A[1]*wk[8]));
  S2 = SIGN(FT[1]*(A[2]*wk[5] + A[3]*wk[6]) + FT[3]*(A[2]*wk[7] + A[3]*wk[8]));
  FT[0] *= S1;
  FT[1] *= S2;
  FT[2] *= S1;
  FT[3] *= S2;
}



// U is mxk and V is nxk
template<class T>
void createLowRankMatHouseholderNBSingle_(const double eps, const unsigned kgoal,
					  const unsigned m,
					  const unsigned k, const unsigned n,
					  T* const U, T* const V,
					  contLowLevel<T>* haarInfo,
					  const T* const X, const T* const Y_,
					  unsigned& khat, T*& data)
{
  const double eps2 = haarInfo->constants.eps2; //accuracy for preserved vectors
  assert(eps2<eps);

  unsigned nwk = MAX(10*k,2);
  int INFO;
  const unsigned mk = MIN(m,k);
  const unsigned nk = MIN(n,k);
  const unsigned mnk = MIN(mk,n);
  const unsigned mnkgoal = MIN(mnk,kgoal);

  T* const tmp = new T[nwk+2+2*k+m+n+mk*nk+mnk*nk];
  assert(tmp!=NULL);
  T* const wk = tmp;					//nwk
  T* const tauH1 = tmp+nwk;				//1
  T* const H1 = tauH1+1;				//m
  T* const tauH2 = H1+m;				//1
  T* const H2 = tauH2+1;				//n
  T* const tauH1U = H2+n;				//k
  T* const tauH2V = tauH1U+k;				//k
  T* const R = tauH2V+k;				//mk*nk
  T* const VT = R+mk*nk;				//mnk*nk

  // qr-decomposition
  // create H1
  blas::copy(m, Y_, H1);
#ifdef MVEC
  INFO = blas::geqrf(m,1,H1,tauH1,nwk,wk);
  assert(INFO==0);
#else
  blas::geqrfs(m, H1, tauH1);
#endif
  // create H2
  blas::copy(n, X, H2);
#ifdef MVEC
  INFO = blas::geqrf(n,1,H2,tauH2,nwk,wk);
#else
  blas::geqrfs(n, H2, tauH2);
#endif
  
  // qr-decomposition
  // of U
  INFO = blas::geqrf(m, k, U, tauH1U, nwk, wk);
  assert(INFO==0);
  // of V
  INFO = blas::geqrf(n, k, V, tauH2V, nwk, wk);
  assert(INFO==0);

  // Berechne R1*R2^T
  blas::utrmmh(mk, k, nk, U, m, V, n, R);

  double* const S = new double[mnk];

  // SVD von R1*R2^T
  INFO = blas::gesvd(mk, nk, R, S, VT, mnk, nwk, wk);
  assert(INFO==0);

  unsigned kt = mnkgoal;
  
  double c_adm = 1.0;
  if(!haarInfo->isadm)
    c_adm = haarInfo->constants.c_adm;
    
  while (kt>0 && S[kt-1]<=c_adm*eps*S[0]) --kt;
  khat = mnkgoal;
  while (khat>0 && S[khat-1]<=c_adm*eps2*S[0]) --khat;

  if (kt+2>=khat) {
    if (khat==0) {
      data = NULL;
    } else {
      data = new T[khat*(m+n)];

      // copy R to data und scale the columns with S
      for (unsigned j=0; j<khat; ++j) {
        unsigned i;
        for (i=0; i<mk; ++i) data[i+j*m] = S[j] * R[i+j*mk];
        for (; i<m; ++i) data[i+j*m] = (T)0.0;
      }

      T* dataV = data + khat*m;

      // copy V to dataV
      for (unsigned j=0; j<khat; ++j) {
        unsigned i;
        for (i=0; i<nk; ++i) dataV[i+j*n] = VT[i*mnk+j];
        for (; i<n; ++i) dataV[i+j*n] = (T)0.0;
      }

      INFO = blas::ormqr(m, khat, mk, U, tauH1U, data, nwk, wk);
      assert(INFO==0);

      INFO = blas::ormqr(n, khat, nk, V, tauH2V, dataV, nwk, wk);
      assert(INFO==0);
    }
  } else {
    const unsigned kr = khat-kt;
    T* const tmp2 = new T[10+kr*(m+n)+(m+n-2)];
    T* const UR = tmp2;			//kr*m
    T* const VR = UR + kr*m;		//kr*n
    T* const M = VR + kr*n;		//4
    T* const R1DH = M + 4;	        //(m-1)
    T* const tauR1 = R1DH + (m-1);	//1
    T* const R2CH = tauR1 + 1;	        //(n-1)
    T* const tauR2 = R2CH + n-1;	//1
    T* const FT = tauR2 + 1;		//4

    // copy rest of R to UR und scale the columns with S
    for (unsigned j=0; j<kr; ++j) {
      unsigned i;
      for (i=0; i<mk; ++i) UR[i+j*m] = S[j+kt] * R[i+(j+kt)*mk];
      for (; i<m; ++i) UR[i+j*m] = (T)0.0;
    }

    INFO = blas::ormqr(m,kr,mk,U,tauH1U,UR,nwk,wk);
    assert(INFO==0);

#ifdef MVEC
    INFO = blas::ormqrh(m,kr,1,H1,m,tauH1,UR,m,nwk,wk);
    assert(INFO==0);
#else
    blas::ormqrsh(m, kr, H1, tauH1, UR, m, nwk, wk);
#endif

    // copy rest of V to VR
    for (unsigned j=0; j<kr; ++j) {
      unsigned i;
      for (i=0; i<nk; ++i) VR[i+j*n] = VT[i*mnk+j+kt];
      for (; i<n; ++i) VR[i+j*n] = (T)0.0;
    }

    INFO = blas::ormqr(n,kr,nk,V,tauH2V,VR,nwk,wk);
    assert(INFO==0);

#ifdef MVEC
    INFO = blas::ormqrh(n,kr,1,H2,n,tauH2,VR,n,nwk,wk);
    assert(INFO==0);
#else
    blas::ormqrsh(n, kr, H2, tauH2, VR, n, nwk, wk);
#endif

    // create CDH
    //blas::gemmh(1,kr,1,1.0,UR,m,VR,n,M,2);
    M[0] = (T)0.0;
    for(unsigned i=0;i<kr;i++)
      M[0] += UR[i*m]*VR[i*n];

    // create R1DH and decompose
    blas::gemmh(m-1,kr,1,1.0,UR+1,m,VR,n,R1DH,m-1);

#ifdef MVEC
    INFO = blas::geqrf(m-1, 1, R1DH, tauR1, nwk, wk);
    assert(INFO==0);
#else
    blas::geqrfs(m-1, R1DH, tauR1);
#endif

    M[1] = R1DH[0];

    // create R2CH and decompose
    blas::gemmh(n-1, kr, 1, 1.0, VR+1, n, UR, m, R2CH, n-1);

#ifdef MVEC
    INFO = blas::geqrf(n-1, 1, R2CH, tauR2, nwk, wk);
    assert(INFO==0);
#else
    blas::geqrfs(n-1, R2CH, tauR2);
#endif
    
    M[2] = R2CH[0];
    M[3] = (T)0.0;

    /*T* M2 = new T[4];
    T* Atmp = new T[4];
    double* SR2 = new double[4];
    T* FT2 = new T[4];
    blas::copy(4, M, M2);
    blas::copy(4, M, Atmp);

    std::cout<<"A:"<<std::endl;
    std::cout<<"["<<M[0]<<","<<M[2]<<";"<<M[1]<<","<<M[3]<<"]"<<std::endl;*/

    double* SR = new double[2];
    //gesvd(M, SR, FT, nwk, wk);

    // SVD of M
    INFO = blas::gesvd(2, 2, M, SR, FT, 2, nwk, wk);
    assert(INFO==0);
    
    /*double* errorVec = new double[2];
    errorTest(Atmp,M,SR,FT,M2,SR2,FT2,errorVec);
    std::cout<<"Error: "<<errorVec[0]<<" "<<errorVec[1]<<std::endl;


    {      
      {
	std::cout<<"m: "<<m<<" n: "<<n<<std::endl;
	  
	std::cout<<"\n\nU:"<<std::endl;
	std::cout<<"["<<M2[0]<<","<<M2[2]<<";"<<M2[1]<<","<<M2[3]<<"]"<<std::endl;
	std::cout<<"S:"<<std::endl;
	std::cout<<"["<<SR2[0]<<",0.0;0.0,"<<SR2[1]<<"]"<<std::endl;      
	std::cout<<"VT:"<<std::endl;
	std::cout<<"["<<FT2[0]<<","<<FT2[2]<<";"<<FT2[1]<<","<<FT2[3]<<"]"<<std::endl;
	
	std::cout<<"U:"<<std::endl;
	std::cout<<"["<<M[0]<<","<<M[2]<<";"<<M[1]<<","<<M[3]<<"]"<<std::endl;
	std::cout<<"S:"<<std::endl;
	std::cout<<"["<<SR[0]<<",0.0;0.0,"<<SR[1]<<"]"<<std::endl;      
	std::cout<<"VT:"<<std::endl;
	std::cout<<"["<<FT[0]<<","<<FT[2]<<";"<<FT[1]<<","<<FT[3]<<"]"<<std::endl;
      }

      if(errorVec[0]>1e-4){
	exit(1);
      }

      delete [] M2;
      delete [] SR2;
      delete [] FT2;
    }
    delete [] errorVec;*/
    
    unsigned rank2=2;
    while (rank2>0 && SR[rank2-1]<=c_adm*eps2*S[0]) --rank2;

    khat = rank2+kt;
    data = new T[khat*(m+n)];
    T* const dataR = data+kt*m;

    // copy R to data und scale the columns with S
    for (unsigned j=0; j<kt; ++j) {
      unsigned i;
      for (i=0; i<mk; ++i) data[i+j*m] = S[j] * R[i+j*mk];
      for (; i<m; ++i) data[i+j*m] = (T)0.0;
    }
    for (unsigned j=0; j<rank2; ++j) {
      unsigned i;
      for (i=0; i<2; ++i) dataR[i+j*m] = SR[j] * M[i+j*2];
      for (; i<m; ++i) dataR[i+j*m] = (T)0.0;
    }

    INFO = blas::ormqr(m, kt, mk, U, tauH1U, data, nwk, wk);
    assert(INFO==0);

#ifdef MVEC
    INFO = blas::ormqr(m-1, rank2, 1, R1DH, m-1, tauR1, dataR+1, m, nwk, wk);
    assert(INFO==0);
#else
    blas::ormqrs(m-1, rank2, R1DH, tauR1, dataR+1, m, nwk, wk);
#endif

#ifdef MVEC
    INFO = blas::ormqr(m, rank2, 1, H1, tauH1, dataR, nwk, wk);
    assert(INFO==0);
#else
    blas::ormqrs(m, rank2, H1, tauH1, dataR, m, nwk, wk);
#endif

    T* const dataV = data + khat*m;
    T* const dataVR = dataV + kt*n;

    // copy V to dataV
    for (unsigned j=0; j<kt; ++j) {
      unsigned i;
      for (i=0; i<nk; ++i) dataV[i+j*n] = VT[i*mnk+j];
      for (; i<n; ++i) dataV[i+j*n] = (T)0.0;
    }
    for (unsigned j=0; j<rank2; ++j) {
      unsigned i;
      for (i=0; i<2; ++i) dataVR[i+j*n] = FT[i*2+j];
      for (; i<n; ++i) dataVR[i+j*n] = (T)0.0;
    }

    INFO = blas::ormqr(n, kt, nk, V, tauH2V, dataV, nwk, wk);
    assert(INFO==0);

#ifdef MVEC
    INFO = blas::ormqr(n-1, rank2, 1, R2CH, n-1, tauR2, dataVR+1, n, nwk, wk);
    assert(INFO==0);
#else
    blas::ormqrs(n-1, rank2, R2CH, tauR2, dataVR+1, n, nwk, wk);
#endif

#ifdef MVEC
    INFO = blas::ormqr(n, rank2, 1, H2, tauH2, dataVR, nwk, wk);
    assert(INFO==0);
#else
    blas::ormqrs(n, rank2, H2, tauH2, dataVR, n, nwk, wk);
#endif
    
    delete [] tmp2;
    delete [] SR;
  }
  delete [] S;
  delete [] tmp;  
}

// U is mxk and V is nxk
template<class T>
void createLowRankMatHouseholderNBMulti_(const double eps, const unsigned kgoal,
					 const unsigned m,
					 const unsigned k, const unsigned n,
					 T* const U, T* const V,
					 contLowLevel<T>* haarInfo,
					 const T* const X, const unsigned ldX,
					 const T* const Y_, const unsigned ldY_,
					 unsigned& khat, T*& data)
{
  const unsigned l = haarInfo->col;
  const double eps2 = 1e-14; //accuracy for preserved vectors
  assert(eps2<eps);

  unsigned nwk = MAX(10*k,2*l);
  int INFO;
  const unsigned mk = MIN(m,k);
  const unsigned nk = MIN(n,k);
  const unsigned mnk = MIN(mk,n);
  const unsigned mnkgoal = MIN(mnk,kgoal);
  const unsigned ll = 2*l;

  T* const tmp = new T[nwk+2*l+2*k+m*l+n*l+mk*nk+mnk*nk];
  assert(tmp!=NULL);
  T* const wk = tmp;					//nwk
  T* const tauH1 = tmp+nwk;				//l
  T* const H1 = tauH1+l;				//ml
  T* const tauH2 = H1+m*l;				//l
  T* const H2 = tauH2+l;				//nl
  T* const tauH1U = H2+n*l;				//k
  T* const tauH2V = tauH1U+k;				//k
  T* const R = tauH2V+k;				//mk*nk
  T* const VT = R+mk*nk;				//mnk*nk

  // qr-decomposition
  // create H1
  for(unsigned i=0;i<l;i++) blas::copy(m,&Y_[ldY_*i],&H1[m*i]);
  INFO = blas::geqrf(m,l,H1,tauH1,nwk,wk);
  assert(INFO==0);
  // create H2
  for (unsigned i=0; i<l; i++) blas::copy(n,&X[ldX*i],&H2[n*i]);
  INFO = blas::geqrf(n,l,H2,tauH2,nwk,wk);
  assert(INFO==0);
  
  // qr-decomposition
  // of U
  INFO = blas::geqrf(m, k, U, tauH1U, nwk, wk);
  assert(INFO==0);
  // of V
  INFO = blas::geqrf(n, k, V, tauH2V, nwk, wk);
  assert(INFO==0);

  // Berechne R1*R2^T
  blas::utrmmh(mk, k, nk, U, m, V, n, R);

  double* const S = new double[mnk];

  // SVD von R1*R2^T
  INFO = blas::gesvd(mk, nk, R, S, VT, mnk, nwk, wk);
  assert(INFO==0);

  unsigned kt = mnkgoal;
  while (kt>0 && S[kt-1]<=eps*S[0]) --kt;
  khat = mnkgoal;
  while (khat>0 && S[khat-1]<=eps2*S[0]) --khat;

  if (kt+ll>=khat) {
    if (khat==0) {
      data = NULL;
    } else {
      data = new T[khat*(m+n)];

      // copy R to data und scale the columns with S
      for (unsigned j=0; j<khat; ++j) {
        unsigned i;
        for (i=0; i<mk; ++i) data[i+j*m] = S[j] * R[i+j*mk];
        for (; i<m; ++i) data[i+j*m] = (T)0.0;
      }

      T* dataV = data + khat*m;

      // copy V to dataV
      for (unsigned j=0; j<khat; ++j) {
        unsigned i;
        for (i=0; i<nk; ++i) dataV[i+j*n] = VT[i*mnk+j];
        for (; i<n; ++i) dataV[i+j*n] = (T)0.0;
      }

      INFO = blas::ormqr(m, khat, mk, U, tauH1U, data, nwk, wk);
      assert(INFO==0);

      INFO = blas::ormqr(n, khat, nk, V, tauH2V, dataV, nwk, wk);
      assert(INFO==0);
    }
  } else {
    const unsigned kr = khat-kt;
    const unsigned mll = MIN(m-l,l);
    const unsigned nll = MIN(n-l,l);
    T* const tmp2 = new T[2*SQR(ll)+kr*(m+n)+(m+n-2*l)*l+mll+nll];
    T* const UR = tmp2;			//kr*m
    T* const VR = UR + kr*m;		//kr*n
    T* const M = VR + kr*n;		//SQR(ll)
    T* const R1DH = M + SQR(ll);	//(m-l)*l
    T* const tauR1 = R1DH + (m-l)*l;	//mll
    T* const R2CH = tauR1 + mll;	//(n-l)*l
    T* const tauR2 = R2CH + (n-l)*l;	//nll
    T* const FT = tauR2 + nll;		//SQR(ll)

    // copy rest of R to UR und scale the columns with S
    for (unsigned j=0; j<kr; ++j) {
      unsigned i;
      for (i=0; i<mk; ++i) UR[i+j*m] = S[j+kt] * R[i+(j+kt)*mk];
      for (; i<m; ++i) UR[i+j*m] = (T)0.0;
    }

    INFO = blas::ormqr(m,kr,mk,U,tauH1U,UR,nwk,wk);
    assert(INFO==0);

    INFO = blas::ormqrh(m,kr,l,H1,m,tauH1,UR,m,nwk,wk);
    assert(INFO==0);

    // copy rest of V to VR
    for (unsigned j=0; j<kr; ++j) {
      unsigned i;
      for (i=0; i<nk; ++i) VR[i+j*n] = VT[i*mnk+j+kt];
      for (; i<n; ++i) VR[i+j*n] = (T)0.0;
    }

    INFO = blas::ormqr(n,kr,nk,V,tauH2V,VR,nwk,wk);
    assert(INFO==0);

    INFO = blas::ormqrh(n,kr,l,H2,n,tauH2,VR,n,nwk,wk);
    assert(INFO==0);

    // create CDH
    blas::gemmh(l,kr,l,1.0,UR,m,VR,n,M,ll);

    // create R1DH and decompose
    blas::gemmh(m-l,kr,l,1.0,UR+l,m,VR,n,R1DH,m-l);

    INFO = blas::geqrf(m-l, l, R1DH, tauR1, nwk, wk);
    assert(INFO==0);

    T* const MR1 = M+l;

    for (unsigned j=0; j<l; ++j) {
      unsigned i;
      for (i=0; i<=j; ++i) MR1[i+j*ll] = R1DH[i+j*(m-l)];
      for (; i<l; ++i) MR1[i+j*ll] = (T)0.0;
    }

    // create R2CH and decompose
    blas::gemmh(n-l,kr,l,1.0,VR+l,n,UR,m,R2CH,n-l);

    INFO = blas::geqrf(n-l, l, R2CH, tauR2, nwk, wk);
    assert(INFO==0);

    T*const MR2 = M+ll*l;

    for (unsigned j=0; j<l; ++j) {
      unsigned i;
      for (i=0; i<j; ++i) MR2[i+j*ll] = (T)0.0;
      for (; i<l; ++i) MR2[i+j*ll] = R2CH[i*(n-l)+j];
      for (; i<ll; ++i) MR2[i+j*ll] = (T)0.0;
    }

    double* const SR = new double[ll];

    // SVD of M
    INFO = blas::gesvd(ll, ll, M, SR, FT, ll, nwk, wk);
    assert(INFO==0);

    unsigned rank2=ll;
    while (rank2>0 && SR[rank2-1]<=eps2*S[0]) --rank2;

    khat = rank2+kt;
    data = new T[khat*(m+n)];
    T* const dataR = data+kt*m;

    // copy R to data und scale the columns with S
    for (unsigned j=0; j<kt; ++j) {
      unsigned i;
      for (i=0; i<mk; ++i) data[i+j*m] = S[j] * R[i+j*mk];
      for (; i<m; ++i) data[i+j*m] = (T)0.0;
    }
    for (unsigned j=0; j<rank2; ++j) {
      unsigned i;
      for (i=0; i<ll; ++i) dataR[i+j*m] = SR[j] * M[i+j*ll];
      for (; i<m; ++i) dataR[i+j*m] = (T)0.0;
    }

    INFO = blas::ormqr(m, kt, mk, U, tauH1U, data, nwk, wk);
    assert(INFO==0);

    INFO = blas::ormqr(m-l, rank2, mll, R1DH, m-l, tauR1, dataR+l, m, nwk, wk);
    assert(INFO==0);

    INFO = blas::ormqr(m, rank2, l, H1, tauH1, dataR, nwk, wk);
    assert(INFO==0);

    T* const dataV = data + khat*m;
    T* const dataVR = dataV + kt*n;

    // copy V to dataV
    for (unsigned j=0; j<kt; ++j) {
      unsigned i;
      for (i=0; i<nk; ++i) dataV[i+j*n] = VT[i*mnk+j];
      for (; i<n; ++i) dataV[i+j*n] = (T)0.0;
    }
    for (unsigned j=0; j<rank2; ++j) {
      unsigned i;
      for (i=0; i<ll; ++i) dataVR[i+j*n] = FT[i*ll+j];
      for (; i<n; ++i) dataVR[i+j*n] = (T)0.0;
    }

    INFO = blas::ormqr(n, kt, nk, V, tauH2V, dataV, nwk, wk);
    assert(INFO==0);

    INFO = blas::ormqr(n-l, rank2, nll, R2CH, n-l, tauR2, dataVR+l, n, nwk, wk);
    assert(INFO==0);

    INFO = blas::ormqr(n, rank2, l, H2, tauH2, dataVR, nwk, wk);
    assert(INFO==0);
    delete [] tmp2;
    delete [] SR;
    /*{
      T* test = new T[m*n];

      blas::gemmh(m,khat,n,(T)1.0,data,m,data+m*khat,n,test,m);
      double nrm(0.0), error(0.0);
      nrm+=blas::nrm2(m*n,test);
      T* Utmp = new T[m*k];
      T* Vtmp = new T[n*k];
      for(unsigned j=0;j<k;j++){
      unsigned i;
      for(i=0;i<=j;i++) Utmp[i+j*m]=U[i+j*m];
      for(;i<m;i++) Utmp[i+j*m]=(T)0.0;
      for(i=0;i<=j;i++) Vtmp[i+j*n]=V[i+j*n];
      for(;i<n;i++) Vtmp[i+j*n]=(T)0.0;
      }
      blas::ormqr(m,k,mk,U,tauH1U,Utmp,nwk,wk);
      blas::ormqr(n,k,nk,V,tauH2V,Vtmp,nwk,wk);
      blas::gemmha(m,k,n,(T)-1.0,Utmp,m,Vtmp,n,test,m);
      error+=blas::nrm2(m*n,test);
      if(error/nrm>1e-2) std::cout<<"Fehler: "<<error/nrm<<std::endl;

      delete [] Utmp;
      delete [] Vtmp;
      delete [] test;
      }*/
  }
  delete [] S;
  delete [] tmp;
}

// U is mxk and V is nxk
template<class T>
void createLowRankMatHouseholderNB(const double eps, const unsigned kgoal,
                                   const unsigned m,
                                   const unsigned k, const unsigned n,
                                   T* const U, T* const V,
                                   contLowLevel<T>* haarInfo,
                                   const T* const X, const unsigned ldX,
                                   const T* const Y_, const unsigned ldY_,
                                   unsigned& khat, T*& data)
{
  if(haarInfo->col==1){
    createLowRankMatHouseholderNBSingle_(eps, kgoal, m, k, n, U, V,
					 haarInfo, X, Y_, khat, data);
  }else{ 
    createLowRankMatHouseholderNBMulti_(eps, kgoal, m, k, n, U, V, 
					haarInfo, X, ldX, Y_, ldY_, khat, data);   
  }
}

// decompose RURV^H
/*template<class T>
void decomposeSubBlockNB(const double eps, const unsigned kgoal,
                         const double Norm2, const unsigned m, const unsigned n,
                         const unsigned l, const unsigned k,
                         const unsigned mlk, const unsigned nlk,
                         T* const H1U, T* const H2V,
                         unsigned& Pm, unsigned& Pn,
                         T*&E, T*&P, T*&FT, const unsigned nwk, T* wk)
{
  // create M and SVD
  T* const M = new T[mlk*nlk];
  assert(M!=NULL);
  double* const SingVal = new double[MIN(mlk,nlk)];
  assert(SingVal!=NULL);
  blas::utrmmh(mlk,k,nlk,&H1U[l],m,&H2V[l],n,M);
  E = new T[SQR(mlk)+SQR(nlk)];
  FT = E+SQR(mlk);
  int INFO;
  INFO = blas::gesvd(mlk, nlk, M, SingVal, E, mlk, FT, nlk, nwk, wk);
  assert(INFO==0);
  unsigned rank=MIN(MIN(mlk,nlk),kgoal);
  while(rank>0 && SingVal[rank-1]<eps*Norm2) --rank;
  const unsigned mlkr = mlk - rank;
  const unsigned nlkr = nlk - rank;
  T* const tmp=new T[(mlk+nlk)*l+MAX(mlk,nlk)*l+MIN(nlkr,l)+MIN(mlkr,l)];
  assert(tmp!=NULL);
  T* const Chat = tmp;                			//nlk*l
  T* const CDtemp = Chat+nlk*l;                         //MAX(mlk,nlk)*l
  T* const tauChat = CDtemp+MAX(mlk,nlk)*l;             //MIN(nlkr,l)
  T* const Dhat = tauChat+MIN(nlkr,l);                  //mlk*l
  T* const tauDhat = Dhat+mlk*l;                        //MIN(mlkr,l)

  // create Chat
  blas::utrgemmh(nlk,k,l,&H2V[l],n,H1U,m,CDtemp,nlk);
  blas::gemm(nlk,nlk,l,1.0,FT,nlk,CDtemp,nlk,Chat,nlk);
  blas::geqrf(nlkr,l,&Chat[rank],nlk,tauChat,nwk,wk);

  // create Dhat
  blas::utrgemmh(mlk,k,l,&H1U[l],m,H2V,n,CDtemp,mlk);
  blas::gemhm(mlk,mlk,l,1.0,E,mlk,CDtemp,mlk,Dhat,mlk);
  blas::geqrf(mlkr,l,&Dhat[rank],mlk,tauDhat,nwk,wk);

  // create P
  Pm = l+rank+MIN(mlkr,l);
  Pn = l+rank+MIN(nlkr,l);
  P = new T[Pm*Pn];

  // create CH*D
  blas::gemmh(l,k,l,1.0,H1U,m,H2V,n,P,Pm);

  // copy Dhat
  for(unsigned j=0;j<l;j++){
    for(unsigned i=0;i<rank+MIN(mlkr,l);i++){
      if(i<=rank+MIN(mlkr,j))
	P[i+j*Pm+l] = Dhat[i+j*mlk];
      else
	P[i+j*Pm+l] = 0.0;
    }
  }
  // copy Chat
  for(unsigned j=0;j<l;j++){
    for(unsigned i=0;i<rank+MIN(nlkr,l);i++){
      if(i<=rank+MIN(nlkr,j))
	P[j+(i+l)*Pm] = Chat[i+j*nlk];
      else
	P[j+(i+l)*Pm] = 0.0;
    }
  }
  // copy Sigma
  for(unsigned j=0;j<Pn-l;j++){
    for(unsigned i=0;i<Pm-l;i++){
      if(i==j && i<rank)
	P[i+j*Pm+l*(Pm+1)] = (T) SingVal[j];
      else
	P[i+j*Pm+l*(Pm+1)] = 0.0;
    }
  }

  // multiply reflector onto E
  blas::morqr(mlk,mlkr,MIN(mlkr,l),&Dhat[rank],mlk,tauDhat,&E[(mlk)*rank],mlk,
	      nwk,wk);
  // multiply reflector onto FT
  blas::ormqrh(nlkr,nlk,MIN(nlkr,l),&Chat[rank],nlk,tauChat,&FT[rank],nlk,
	       nwk,wk);

  delete [] SingVal;
  delete [] M;
  delete [] tmp;
}

// U is mxk and V is nxk
template<class T>
void createLowRankMatHouseholderNB2(const double eps, const unsigned kgoal,
				    const unsigned m,
				    const unsigned k, const unsigned n,
				    const T* const U, const T* const V,
				    const unsigned l,
				    const T* const X, const unsigned ldX,
				    const T* const Y_, const unsigned ldY_,
				    unsigned& khat, T*& data)
{
  //assert(k>=l);
  if(l<m and l<n and k>=l){
    unsigned nwk = 2*MAX(5*MAX(n,m),3*MIN(n,m)+MAX(n,m));
    int INFO;
    const unsigned mlk = MIN(m-l,k);
    const unsigned nlk = MIN(n-l,k);
    T* tmp = new T[3*SQR(k)+(l+k)*(2+m+n)+MAX(m,n)*(l+MAX(mlk,nlk))+nwk];
    T* wk = tmp;			//nwk
    T* UHU = wk+nwk;			//SQR(k)
    T* VHV = UHU+SQR(k);		//SQR(k)
    T* UHUVHV = VHV+SQR(k);		//SQR(k)
    T* tauH1 = UHUVHV+SQR(k);		//l
    T* H1 = tauH1+l;			//ml
    T* tauH2 = H1+m*l;			//l
    T* H2 = tauH2+l;			//nl
    T* H1U = H2+n*l;			//mk
    T* H2V = H1U+m*k;			//nk
    T* tauH1U = H2V+n*k;		//k
    T* tauH2V = tauH1U+k;		//k
    T* UVhattemp = tauH2V+k;		//MAX(m,n)*(l+MAX(mlk,nlk))

    // norm of UVH
    blas::gemhm(m,k,k,(T)1.0,U,m,U,m,UHU,k);
    blas::gemhm(n,k,k,(T)1.0,V,n,V,n,VHV,k);
    blas::gemm(k,k,k,(T)1.0,UHU,k,VHV,k,UHUVHV,k);

    const double nrm2 = sqrt(vectorIter(5,k,UHUVHV));
    // qr-decomposition
    // create H1
    for(unsigned i=0;i<l;i++)	blas::copy(m,&Y_[ldY_*i],&H1[m*i]);
    INFO = blas::geqrf(m,l, H1,tauH1,nwk,wk);
    assert(INFO==0);
    // create H2
    for(unsigned i=0;i<l;i++) blas::copy(n,&X[ldX*i],&H2[n*i]);
    INFO = blas::geqrf(n,l,H2,tauH2,nwk,wk);
    assert(INFO==0);

    // multiply H1U
    blas::copy(m*k,U,H1U);
    INFO = blas::ormqrh(m,k,l,H1,tauH1,H1U,nwk,wk);
    assert(INFO==0);

    // multiply H2V
    blas::copy(n*k,V,H2V);
    INFO = blas::ormqrh(n,k,l,H2,tauH2,H2V,nwk,wk);
    assert(INFO==0);

    // qr-decomposition H1U
    INFO = blas::geqrf(m-l,k,&H1U[l],m,tauH1U,nwk,wk);
    assert(INFO==0);

    // qr-decomposition H2V
    INFO = blas::geqrf(n-l,k,&H2V[l],n,tauH2V,nwk,wk);
    assert(INFO==0);

    T *E, *FT, *P;
    unsigned Pn, Pm;
    decomposeSubBlockNB(eps, kgoal, nrm2, m, n, l, k, mlk, nlk,
			H1U, H2V, Pm, Pn, E, P, FT, nwk, wk);

    if(l>0){
      // reduce rank by svd
      const unsigned minP = MIN(Pn,Pm);
      double* const PS = new double[minP];
      T* const PVT = new T[minP*Pn];

      // SVD of P
      INFO = blas::gesvd(Pm, Pn, P, PS, PVT, minP, nwk, wk);
      assert(INFO==0);

      khat = minP;
      while (khat>0 && (PS[khat-1]<=1e-14*nrm2 || PS[khat-1]<1e-32)) --khat;
      for(unsigned j=0;j<khat;j++)
	for(unsigned i=0;i<Pm;i++) P[i+j*Pm]*=(T)PS[j];

      if(khat>0){
	data = new T[(n+m)*khat];

	// create U hat
	blas::setzero(m*(l+mlk),UVhattemp);
	for(unsigned i=0;i<l;i++) UVhattemp[i+i*m]=1.0;
	for(unsigned i=0;i<mlk;i++)
	  for(unsigned j=0;j<mlk;j++)
	    UVhattemp[i+j*m+(m+1)*l]=E[i+j*mlk];
	INFO = blas::ormqr(m-l,mlk,mlk,&H1U[l],m,tauH1U,&UVhattemp[(m+1)*l],m,nwk,wk);
	assert(INFO==0);
	INFO = blas::ormqr(m,l+mlk,l,H1,m,tauH1,UVhattemp,m,nwk,wk);
	assert(INFO==0);
	blas::gemm(m,Pm,khat,1.0,UVhattemp,m,P,Pm,data,m);

	// create V hat
	blas::setzero(n*(l+nlk),UVhattemp);
	for(unsigned i=0;i<l;i++) UVhattemp[i+i*n]=1.0;
	for(unsigned i=0;i<nlk;i++)
	  for(unsigned j=0;j<nlk;j++)
	    UVhattemp[i+j*n+(n+1)*l]=FT[j+i*nlk];
	INFO = blas::ormqr(n-l,nlk,nlk,&H2V[l],n,tauH2V,&UVhattemp[(n+1)*l],n,nwk,wk);
	assert(INFO==0);
	INFO = blas::ormqr(n,l+nlk,l,H2,n,tauH2,UVhattemp,n,nwk,wk);
	assert(INFO==0);
	blas::gemmh(n,Pn,khat,1.0,UVhattemp,n,PVT,minP,data+m*khat,n);
      }else{
	data=NULL;
      }

      delete [] tmp;
      delete [] P;
      delete [] E;
      delete [] PS;
      delete [] PVT;
    }else{
      khat = MIN(Pn,Pm);
      data = new T[(n+m)*khat];

      // create U hat
      T* Unew = data;
      if(Pn<Pm){
	blas::setzero(m*(l+mlk),UVhattemp);
	for(unsigned i=0;i<l;i++)
	  UVhattemp[i+i*m]=1.0;
	for(unsigned i=0;i<mlk;i++)
	  for(unsigned j=0;j<mlk;j++)
	    UVhattemp[i+j*m+(m+1)*l]=E[i+j*mlk];
	INFO = blas::ormqr(m-l,mlk,mlk,&H1U[l],m,tauH1U,
			   &UVhattemp[(m+1)*l],m,nwk,wk);
	assert(INFO==0);
	INFO = blas::ormqr(m,l+mlk,l,H1,m,tauH1,UVhattemp,m,nwk,wk);
	assert(INFO==0);
	blas::gemm(m,Pm,khat,1.0,UVhattemp,m,P,Pm,data,m);
      }else{
	blas::setzero(m*khat,Unew);
	for(unsigned i=0;i<l;i++) Unew[i+i*m]=1.0;
	for(unsigned i=0;i<mlk;i++)
	  for(unsigned j=0;j<khat-l;j++)
	    Unew[i+j*m+(m+1)*l]=E[i+j*mlk];
	INFO = blas::ormqr(m-l,khat-l,mlk,&H1U[l],m,tauH1U,&Unew[(m+1)*l],m,nwk,wk);
	assert(INFO==0);
	INFO = blas::ormqr(m,khat,l,H1,m,tauH1,Unew,m,nwk,wk);
	assert(INFO==0);
      }

      // create V hat
      T* Vnew = data+m*khat;
      if(Pn<Pm){
	blas::setzero(n*khat,Vnew);
	for(unsigned i=0;i<l;i++) Vnew[i+i*n]=1.0;
	for(unsigned i=0;i<nlk;i++)
	  for(unsigned j=0;j<khat-l;j++)
	    Vnew[i+j*n+(n+1)*l]=FT[j+i*nlk];
	INFO = blas::ormqr(n-l,khat-l,nlk,&H2V[l],n,tauH2V,&Vnew[(n+1)*l],n,nwk,wk);
	assert(INFO==0);
	INFO = blas::ormqr(n,khat,l,H2,n,tauH2,Vnew,n,nwk,wk);
	assert(INFO==0);
      }else{
	blas::setzero(n*(l+nlk),UVhattemp);
	for(unsigned i=0;i<l;i++) UVhattemp[i+i*n]=1.0;
	for(unsigned i=0;i<nlk;i++)
	  for(unsigned j=0;j<nlk;j++)
	    UVhattemp[i+j*n+(n+1)*l]=FT[j+i*nlk];
	INFO = blas::ormqr(n-l,nlk,nlk,&H2V[l],n,tauH2V,
			   &UVhattemp[(n+1)*l],n,nwk,wk);
	assert(INFO==0);
	INFO = blas::ormqr(n,l+nlk,l,H2,n,tauH2,UVhattemp,n,nwk,wk);
	assert(INFO==0);
	blas::gemmh(n,Pn,Pm,1.0,UVhattemp,n,P,Pm,Vnew,n);
      }
      delete [] tmp;
      delete [] P;
      delete [] E;
    }
  }else{
    khat = k;
    data=new T[k*(m+n)];
    blas::copy(m*k,U,data);
    blas::copy(n*k,V,data+m*k);
  }
}*/

/*// U is mxk and V is nxk
void createLowRankMatNB_OLD(const double eps, const unsigned m, const unsigned k,
                            const unsigned n,
                            const double* const U, const double* V,
                            const unsigned l, const double* const X,
                            const double* const Y_,
                            unsigned& khat, double*& data)
{
	// norm of UVH
	double* UHU = new double[SQR(k)];
	double* VHV = new double[SQR(k)];
	double* UHUVHV = new double[SQR(k)];
	blas::gemhm(m,k,k,1.0,U,m,U,m,UHU,k);
	blas::gemhm(n,k,k,1.0,V,n,V,n,VHV,k);
	blas::gemm(k,k,k,1.0,UHU,k,VHV,k,UHUVHV,k);
	const double nrm2 = sqrt(vectorIter(5,k,UHUVHV));
	// create (XH*X)^-1
	double* invX = new double[l*(l+1)/2];
	blas::setzero(l*(l+1)/2, invX);
	blas::symhm(n, l, X, 1.0, invX);
	lapack::geinv_herm(l, invX);

	//create (Y_H*Y_)^-1
  double* invY_ = new double[l*(l+1)/2];
  blas::setzero(l*(l+1)/2, invY_);
  blas::symhm(m, l, Y_, 1.0, invY_);
  lapack::geinv_herm(l, invY_);

  // create orthoU
  double* orthoU = new double[m*k];
  double* orthoUtemp1 = new double[l*k];
  double* orthoUtemp2 = new double[l*k];
  blas::setzero(l*k,orthoUtemp2);
  blas::copy(m*k, U, orthoU);
  blas::gemhm(m, l, k, 1.0, Y_, m, U, m, orthoUtemp1, l);
  blas::sygemma(l, k, invY_, orthoUtemp1, 1.0, orthoUtemp2);
  blas::gemma(m, l, k, -1.0, Y_, m, orthoUtemp2, l, orthoU, m);
  delete [] orthoUtemp1;
  delete [] orthoUtemp2;

  // create orthoV
  double* orthoV = new double[n*k];
  double* orthoVtemp1 = new double[l*k];
  double* orthoVtemp2 = new double[l*k];
  blas::setzero(l*k,orthoVtemp2);
  blas::copy(n*k, V, orthoV);

  blas::gemhm(n, l, k, 1.0, X, n, V, n, orthoVtemp1, l);
  blas::sygemma(l, k, invX, orthoVtemp1, 1.0, orthoVtemp2);
  blas::gemma(n, l, k, -1.0, X, n, orthoVtemp2, l, orthoV, n);
  delete [] orthoVtemp1;
  delete [] orthoVtemp2;

  // apporximate orthoU and orthoVH
  int INFO;
  unsigned LWORK = 10*MAX(n,m);
  double* WORK = new double[LWORK];
  double* tauU = new double[k];
  double* tauV = new double[k];
  double* tmp1 = new double[m*k];
  double* tmp2 = new double[n*k];
  blas::copy(m*k,orthoU,tmp1);
  blas::copy(n*k,orthoV,tmp2);
  INFO = blas::geqrf(m,k,tmp1,tauU,LWORK,WORK);
  assert(INFO==0);
  INFO = blas::geqrf(n,k,tmp2,tauV,LWORK,WORK);
  assert(INFO==0);
  double* R = new double[SQR(k)];
  blas::utrmmh(k,k,k,tmp1,m,tmp2,n,R);
  double* SingVal = new double[k];
  double* VH = new double[m*k];
  INFO = blas::gesvd(k,k,R,SingVal,VH,k,LWORK,WORK);
  assert(INFO==0);

  unsigned rank=k;
  while(rank>0 && SingVal[rank-1]<eps*nrm2) --rank;
  for(unsigned j=0;j<rank;j++){
  unsigned i(0);
  for(;i<k;i++)
  orthoU[i+j*m] = SingVal[j]*R[i+j*k];
  blas::setzero(m-k,&orthoU[k+j*m]);
  }
  for(unsigned j=0;j<rank;j++){
  unsigned i(0);
  for(;i<k;i++)
  orthoV[i+j*n] = VH[i*k+j];
  blas::setzero(n-k,&orthoV[k+j*n]);
  }
  INFO = blas::ormqr(m, rank, k, tmp1, tauU, orthoU, LWORK, WORK);
  assert(INFO==0);

  INFO = blas::ormqr(n, rank, k, tmp2, tauV, orthoV, LWORK, WORK);
  assert(INFO==0);
  delete [] WORK;
  delete [] tauU;
  delete [] tauV;
  delete [] tmp1;
  delete [] tmp2;
  delete [] R;
  delete [] SingVal;
  delete [] VH;

  // create Ytemp
  double* Ytemp = new double[k*l];
  blas::gemhm(n,k,l,1.0,V,n,X,n,Ytemp,k);

  // create Uhat
  data = new double[(m+n)*(2*l+rank)];
  blas::gemm(m,k,l,1.0,U,m,Ytemp,k,data,m);
  blas::setzero(m*l,&data[l*m]);
  blas::gesymma(m,l,Y_,invY_,1.0,&data[l*m]);
  blas::copy(rank*m,orthoU,&data[2*l*m]);

  // create X_
  double* X_ = new double[n*l];
  double* X_temp = new double[k*l];
  blas::gemhm(m,k,l,1.0,U,m,Y_,m,X_temp,k);
  blas::gemm(n,k,l,1.0,V,n,X_temp,k,X_,n);

  // create (XH*X)^-1*XH*X_
  double* XHX_ = new double[l*l];
  blas::gemhm(n,l,l,1.0,X,n,X_,n,XHX_,l);
  double* invXHXtXHX_ = new double[l*l];
  blas::setzero(l*l,invXHXtXHX_);
  blas::sygemma(l,l,invX,XHX_,1.0,invXHXtXHX_);

  // create Vhat
  blas::setzero(n*l,&data[m*(2*l+rank)]);
  blas::gesymma(n,l,X,invX,1.0,&data[m*(2*l+rank)]);
  blas::copy(n*l,X_,&data[m*(2*l+rank)+n*l]);
  blas::gemma(n,l,l,-1.0,X,n,invXHXtXHX_,l,&data[m*(2*l+rank)+n*l],n);
  blas::copy(rank*n,orthoV,&data[m*(2*l+rank)+2*l*n]);

  khat = 2*l+rank;
  delete [] X_;
  delete [] X_temp;
  delete [] Ytemp;
  delete [] XHX_;
  delete [] invXHXtXHX_;

  delete [] orthoU;
  delete [] orthoV;
  delete [] invX;
  delete [] invY_;
  }*/

// A=U1 V1H, U1 mxk1 V1 nxk1 and B=data, data mxk2xn
template<class T>
void addLowRankNB(const double eps, const unsigned kgoal, const unsigned m,
                  const unsigned k1, const unsigned k2,
                  const unsigned n,
                  const T* const U1, const unsigned ldU1,
                  const T* const V1, const unsigned ldV1,
                  contLowLevel<T>* haarInfo, 
		  const T* const X, const unsigned ldX,
                  const T* const Y_, const unsigned ldY_,
                  unsigned& khat, T*& data)
{
  const unsigned ksum = k1+k2;
  T* Unew = new T[m*ksum];
  T* Vnew = new T[n*ksum];
  for (unsigned i=0; i<k1; i++)
    blas::copy(m,&U1[i*ldU1],&Unew[i*m]);
  blas::copy(m*k2,data,&Unew[m*k1]);
  for (unsigned i=0; i<k1; i++)
    blas::copy(n,&V1[i*ldV1],&Vnew[i*n]);
  blas::copy(n*k2,&data[m*k2],&Vnew[n*k1]);
  delete [] data;
  data = NULL;
  createLowRankMatHouseholderNB(eps, kgoal, m, ksum, n, Unew, Vnew, haarInfo,
                                X, ldX, Y_, ldY_, khat, data);
  delete [] Unew;
  delete [] Vnew;
}

// U is mxk and V is nxk
template<class T>
void createLowRankMatHouseholderNB_OLD_(const T eps, const unsigned m,
                                        const unsigned k, const unsigned n,
                                        const T* const U, const T* const V,
                                        const unsigned l,
                                        const T* const X, const T* const Y_,
                                        unsigned& khat, T*& data)
{
  // qr-decomposition
  unsigned nwk = 10*MAX(n,m);
  T* wk = new T[nwk];
  int INFO;
  // create H1
  T* tauH1 = new T[l];
  T* H1 = new T[m*l];
  blas::copy(m*l,Y_,H1);
  INFO = blas::geqrf(m,l,H1,tauH1,nwk,wk);
  assert(INFO==0);
  // create H2
  T* tauH2 = new T[l];
  T* H2 = new T[n*l];
  blas::copy(n*l,X,H2);
  INFO = blas::geqrf(n,l,H2,tauH2,nwk,wk);
  assert(INFO==0);

  // multiply H1U
  T* H1U = new T[m*k];
  blas::copy(m*k,U,H1U);
  INFO = blas::ormqrh(m,k,l,H1,tauH1,H1U,nwk,wk);
  assert(INFO==0);

  // multiply H2V
  T* H2V = new T[n*k];
  blas::copy(n*k,V,H2V);
  INFO = blas::ormqrh(n,k,l,H2,tauH2,H2V,nwk,wk);
  assert(INFO==0);

  // qr-decomposition H1U
  T* tauH1U = new T[k];
  INFO = blas::geqrf(m-l,k,&H1U[l],m,tauH1U,nwk,wk);
  assert(INFO==0);

  // qr-decomposition H2V
  T* tauH2V = new T[k];
  INFO = blas::geqrf(n-l,k,&H2V[l],n,tauH2V,nwk,wk);
  assert(INFO==0);

  // create M and SVD
  const unsigned mlk = MIN(m-l,k);
  const unsigned nlk = MIN(n-l,k);
  T* M = new T[mlk*nlk];
  blas::utrmmh(mlk,k,nlk,&H1U[l],m,&H2V[l],n,M);
  T* SingVal = new T[MIN(mlk,nlk)];
  T* FT = new T[SQR(nlk)];
  blas::setzero(SQR(nlk),FT);
  T* E = new T[SQR(mlk)];
  blas::setzero(SQR(mlk),E);
  INFO = blas::gesvd(mlk, nlk, M, SingVal, E, mlk, FT, nlk, nwk, wk);
  assert(INFO==0);
  unsigned rank=MIN(mlk,nlk);
  while (rank>0 && SingVal[rank-1]<eps*SingVal[0]) --rank;
  std::cout<<"Rank: "<<rank<<std::endl;

  const unsigned mlkr = mlk - rank;
  const unsigned nlkr = nlk - rank;

  // create Chat
  T* Chat = new T[nlk*l];
  T* CDtemp = new T[MAX(mlk,nlk)*l];
  blas::utrgemmh(nlk,k,l,&H2V[l],n,H1U,m,CDtemp,nlk);
  blas::gemm(nlk,nlk,l,1.0,FT,nlk,CDtemp,nlk,Chat,nlk);
  T* tauChat = new T[MIN(nlkr,l)];
  blas::geqrf(nlkr,l,&Chat[rank],nlk,tauChat,nwk,wk);

  // create Dhat
  T* Dhat = new T[mlk*l];
  blas::utrgemmh(mlk,k,l,&H1U[l],m,H2V,n,CDtemp,mlk);
  blas::gemhm(mlk,mlk,l,1.0,E,mlk,CDtemp,mlk,Dhat,mlk);
  T* tauDhat = new T[MIN(mlkr,l)];
  blas::geqrf(mlkr,l,&Dhat[rank],mlk,tauDhat,nwk,wk);
  delete [] CDtemp;

  // create P
  const unsigned Pm = l+rank+MIN(mlkr,l);
  const unsigned Pn = l+rank+MIN(nlkr,l);
  T* P = new T[Pm*Pn];

  // create CH*D
  blas::gemmh(l,k,l,1.0,H1U,m,H2V,n,P,Pm);
  // copy Dhat
  for (unsigned j=0; j<l; j++)
    for (unsigned i=0; i<rank+MIN(mlkr,l); i++)
      if (i<=rank+MIN(mlkr,j))
        P[i+j*Pm+l] = Dhat[i+j*mlk];
      else
        P[i+j*Pm+l] = 0.0;
  // copy Chat
  for (unsigned j=0; j<l; j++)
    for (unsigned i=0; i<rank+MIN(nlkr,l); i++)
      if (i<=rank+MIN(nlkr,j))
        P[j+(i+l)*Pm] = Chat[i+j*nlk];
      else
        P[j+(i+l)*Pm] = 0.0;
  // copy Sigma
  for (unsigned j=0; j<Pn-l; j++)
    for (unsigned i=0; i<Pm-l; i++)
      if (i==j && i<rank)
        P[i+j*Pm+l*(Pm+1)] = SingVal[j];
      else
        P[i+j*Pm+l*(Pm+1)] = 0.0;

  // svd of P
  T eps2 = 1e-14;
  T* SingVal2 = new T[MIN(Pm,Pn)];
  T* QH = new T[MIN(Pm,Pn)*Pn];
  INFO = blas::gesvd(Pm,Pn,P,SingVal2,QH,MIN(Pm,Pn),nwk,wk);
  assert(INFO==0);
  unsigned rank2=MIN(Pm,Pn);
  while (rank2>0 && SingVal2[rank-1]<eps2*SingVal2[0]) --rank2;
  std::cout<<"Rank2: "<<rank2<<std::endl;
  std::cout<<"Pm: "<<Pm<<" Pn: "<<Pn<<std::endl;
  for (unsigned i=0; i<Pm; i++)
    for (unsigned j=0; j<rank2; j++)
      P[i+j*Pm] *= SingVal2[j];

  // multiply reflector onto E
  blas::morqr(mlk,mlkr,MIN(mlkr,l),&Dhat[rank],mlk,tauDhat,&E[(mlk)*rank],mlk,
              nwk,wk);
  // multiply reflector onto FT
  blas::ormqrh(nlkr,nlk,MIN(nlkr,l),&Chat[rank],nlk,tauChat,&FT[rank],nlk,
               nwk,wk);

  // create V hat
  T* UVhattemp = new T[MAX(m,n)*(l+MAX(mlk,nlk))];
  blas::setzero(n*(l+nlk),UVhattemp);
  for (unsigned i=0; i<l; i++)
    UVhattemp[i+i*n]=1.0;
  for (unsigned i=0; i<nlk; i++)
    for (unsigned j=0; j<nlk; j++)
      UVhattemp[i+j*n+(n+1)*l]=FT[j+i*nlk];
  INFO = blas::ormqr(n-l,nlk,nlk,&H2V[l],n,tauH2V,&UVhattemp[(n+1)*l],n,nwk,wk);
  assert(INFO==0);
  INFO = blas::ormqr(n,l+nlk,l,H2,n,tauH2,UVhattemp,n,nwk,wk);
  assert(INFO==0);

  data = new T[(n+m)*rank2];
  blas::gemmh(n,Pn,rank2,1.0,UVhattemp,n,QH,MIN(Pm,Pn),&data[m*rank2],n);

  // create U hat
  blas::setzero(m*(l+mlk),UVhattemp);
  for (unsigned i=0; i<l; i++)
    UVhattemp[i+i*m]=1.0;
  for (unsigned i=0; i<mlk; i++)
    for (unsigned j=0; j<mlk; j++)
      UVhattemp[i+j*m+(m+1)*l]=E[i+j*mlk];
  INFO = blas::ormqr(m-l,mlk,mlk,&H1U[l],m,tauH1U,&UVhattemp[(m+1)*l],m,nwk,wk);
  assert(INFO==0);
  INFO = blas::ormqr(m,l+mlk,l,H1,m,tauH1,UVhattemp,m,nwk,wk);
  assert(INFO==0);
  blas::gemm(m,Pm,rank2,1.0,UVhattemp,m,P,Pm,data,m);

  khat = rank2;

  delete [] SingVal2;
  delete [] QH;
  delete [] UVhattemp;
  delete [] P;
  delete [] tauChat;
  delete [] tauDhat;
  delete [] Chat;
  delete [] Dhat;
  delete [] wk;
  delete [] tauH1;
  delete [] tauH2;
  delete [] H1;
  delete [] H2;
  delete [] H1U;
  delete [] H2V;
  delete [] tauH1U;
  delete [] tauH2V;
  delete [] M;
  delete [] SingVal;
  delete [] E;
  delete [] FT;
}

#endif
