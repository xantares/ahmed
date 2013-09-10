/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// solve L X = B for X
// forward substitution, B is destroyed
template<class T> static
void LtHGeM_solve_(blcluster* blL, mblock<T>** L,
                   unsigned p, T* B, unsigned ldB)
{
  if (blL->isleaf()) L[blL->getidx()]->ltr_solve(p, B, ldB);
  else {
    unsigned ns = blL->getnrs();
    T* Bi = B;
    for (unsigned i=0; i<ns; ++i) {
      T* Bk = B;
      for (unsigned k=0; k<i; ++k) {
        blcluster *sonL = blL->getson(i, k);
        mltaGeHGeM((T) -1.0, sonL, L, p, Bk, ldB, Bi, ldB);
        Bk += sonL->getn2();
      }
      blcluster *son = blL->getson(i, i);
      LtHGeM_solve_(son, L, p, Bi, ldB);
      Bi += son->getn1();
    }
  }
}


// solve L^H X = B for X
// forward substitution, B is destroyed
template<class T> static
void LtHhGeM_solve_(blcluster* blL, mblock<T>** L,
                    unsigned p, T* B, unsigned ldB)
{
  if (blL->isleaf()) L[blL->getidx()]->ltrh_solve(p, B, ldB);
  else {
    unsigned ns = blL->getncs();
    T* Bi = B + blL->getn2();
    for (unsigned i=ns; i>0; --i) {
      blcluster *son = blL->getson(i-1, i-1);
      T* Bk = Bi;
      Bi -= son->getn2();
      for (unsigned k=i; k<ns; ++k) {
	blcluster *sonL = blL->getson(k, i-1);
	mltaGeHhGeM((T) -1.0, sonL, L, p, Bk, ldB, Bi, ldB);
	Bk += sonL->getn1();
      }
      LtHhGeM_solve_(son, L, p, Bi, ldB);
    }
  }
}


// solve U X = B for X
// backward substitution, B is destroyed
template<class T> static
void UtHGeM_solve_(blcluster* blU, mblock<T>** U,
                   unsigned p, T* B, unsigned ldB)
{
  if (blU->isleaf()) U[blU->getidx()]->utr_solve(p, B, ldB);
  else {
    unsigned ns = blU->getnrs();
    T* Bi = B + blU->getn1();
    for (unsigned i=ns; i>0; --i) {
      blcluster *son = blU->getson(i-1, i-1);
      T* Bk = Bi;
      Bi -= son->getn1();
      for (unsigned k=i; k<ns; ++k) {
	blcluster *sonU = blU->getson(i-1, k);
	mltaGeHGeM((T) -1.0, sonU, U, p, Bk, ldB, Bi, ldB);
	Bk += sonU->getn2();
      }
      UtHGeM_solve_(son, U, p, Bi, ldB);
    }
  }
}


// solve U^H X = B for X
// forward substitution, B contains solution on exit
template<class T> static
void UtHhGeM_solve_(blcluster* blU, mblock<T>** U,
                    unsigned p, T* B, unsigned ldB)
{
  if (blU->isleaf()) U[blU->getidx()]->utrh_solve(p, B, ldB);
  else {
    unsigned ns = blU->getncs();
    T* Bi = B;
    for (unsigned i=0; i<ns; ++i) {
      T* Bk = B;
      for (unsigned k=0; k<i; ++k) {
        blcluster *sonU = blU->getson(k, i);
        mltaGeHhGeM((T) -1.0, sonU, U, p, Bk, ldB, Bi, ldB);
        Bk += sonU->getn1();
      }
      blcluster *son = blU->getson(i, i);
      UtHhGeM_solve_(son, U, p, Bi, ldB);
      Bi += son->getn2();
    }
  }
}

// solve L X = B for X, B is destroyed
template<class T> static
void LtHGeH_solve_(blcluster* blL, mblock<T>** L, blcluster* blB, mblock<T>** B,
                 mblock<T>** X, double eps, unsigned rankmax,
		 contBasis<T>* haar=NULL)
{
  if (blB->isleaf()) {    // B is a leaf
    T* dataB = blB->data(B);
    const unsigned mB = blB->getn1();

    if (blB->isLrM(B)) {
      const unsigned rankB = blB->rank(B);
      LtHGeM_solve_(blL, L, rankB, dataB, mB);
      X[blB->getidx()]->cpyLrM(rankB, dataB, dataB+rankB*mB);
    } else {
      X[blB->getidx()]->setGeM();
      const unsigned nB = blB->getn2();
      LtHGeM_solve_(blL, L, nB, dataB, mB);
      blas::copy(mB*nB, dataB, blB->data(X));
    }
  } else {    // then L cannot be a leaf
    unsigned n1 = blB->getnrs(), n2 = blB->getncs();
    contBasis<T>* haarSon = NULL;
    for (unsigned i=0; i<n1; ++i) {
      blcluster *son = blL->getson(i, i);
      for (unsigned j=0; j<n2; ++j) {
	blcluster* sonB = blB->getson(i, j);
	if(haar) haarSon = haar->son(i,j);
	for (unsigned k=0; k<i; ++k) {
	  blcluster *sonL = blL->getson(i, k);
	  blcluster *sonX = blB->getson(k, j);
	  mltaGeHGeH((T) -1.0, sonL, L, sonX, B, sonB, B, eps, rankmax,
		  haarSon);
	}
	LtHGeH_solve_(son, L, sonB, B, X, eps, rankmax, haarSon);
	delete haarSon;
      }
    }
  }
}

// solve U X = B for X and store X in B, U upper triangular H-Matrix
template<class T> static
void UtHGeH_solve_(blcluster* blU, mblock<T>** U,
                 blcluster* blB, mblock<T>** B, double eps, unsigned rankmax)
{
  if (blB->isleaf()){
    T* dataB = blB->data(B);
    unsigned k, mB = blB->getn1();
    
    if (blB->isLrM(B))
      k = blB->rank(B);
    else
      k = blB->getn2();
    
    UtHGeM_solve(blU, U, k, dataB, mB);
  }else{   // then U cannot be a leaf    
    unsigned ni = blB->getnrs(), nj = blB->getncs();
    for (unsigned i=ni; i>0;){
      --i;
      blcluster *son = blU->getson(i, i);
      for (unsigned j=0; j<nj; ++j){
	blcluster* sonB = blB->getson(i, j);
	for (unsigned k=i+1; k<ni; ++k){
	  blcluster *sonU = blU->getson(k, i);
	  blcluster *sonX = blB->getson(k, j);
	  mltaGeHGeH((T) -1.0, sonU, U, sonX, B, sonB, B, eps, rankmax);
	}
	UtHGeH_solve_(son, U, sonB, B, eps, rankmax);
      }
    }
  }
}


// solve U^H X = B for X and store X in B, U upper triangular H-Matrix
template<class T> static
void UtHhGeH_solve_(blcluster* blU, mblock<T>** U,
                  blcluster* blB, mblock<T>** B, double eps, unsigned rankmax,
                  contBasis<T>* haar=NULL)
{
  if (blB->isleaf()) {
    T* dataB = blB->data(B);
    unsigned k, mB = blB->getn1();

    if (blB->isLrM(B))
      k = blB->rank(B);
    else
      k = blB->getn2();

    UtHhGeM_solve(blU, U, k, dataB, mB);
  } else {   // then U cannot be a leaf
    unsigned ni = blB->getnrs(), nj = blB->getncs();
    contBasis<T>* haarSon = NULL;
    for (unsigned i=0; i<ni; ++i) {
      blcluster *son = blU->getson(i, i);
      for (unsigned j=0; j<nj; ++j) {
	blcluster* sonB = blB->getson(i, j);
	if (haar) haarSon = haar->son(i,j);
	for (unsigned k=0; k<i; ++k) {
	  blcluster *sonU = blU->getson(k, i);
	  blcluster *sonX = blB->getson(k, j);
	  mltaGeHhGeH((T) -1.0, sonU, U, sonX, B, sonB, B, eps, rankmax,
		      haarSon);
	}
	UtHhGeH_solve_(son, U, sonB, B, eps, rankmax, haarSon);
	delete haarSon;
      }
    }
  }
}

// solve X U = B for X (B H-matrix, U UtH-matrix), B is destroyed
template<class T> static
void GeHUtH_solve_(blcluster* blU, mblock<T>** U, blcluster* blB,
                 mblock<T>** B, mblock<T>** X, double eps, unsigned rankmax,
		 contBasis<T>* haar=NULL)
{
  if (blB->isleaf())    // solve U^h X^H = B^H
    {
      const unsigned mB = blB->getn1(), nB = blB->getn2();
      T* dataB = blB->data(B);

      if (blB->isLrM(B)) {
	const unsigned rankB = blB->rank(B);
	UtHhGeM_solve(blU, U, rankB, dataB + rankB*mB, nB);
	X[blB->getidx()]->cpyLrM(rankB, dataB, dataB+rankB*mB);
      }
      else
	{
	  X[blB->getidx()]->setGeM();

	  if (blU->isleaf())
	    U[blU->getidx()]->utr_solve_left(mB, dataB, mB, blB->data(X), mB);
	  else
	    {
	      T* tmp = new T[nB*mB];
	      assert(tmp!=NULL);

	      blas::transpose(mB, nB, dataB, tmp);
	      UtHhGeM_solve(blU, U, mB, tmp, nB);
	      blas::transpose(nB, mB, tmp, blB->data(X));
	      delete [] tmp;
	    }
	}
    }
  else    // then U cannot be a leaf
    {
      unsigned ni = blB->getnrs(), nj = blB->getncs();
      contBasis<T>* haarSon = NULL;
      for (unsigned j=0; j<nj; ++j) {
	blcluster *son = blU->getson(j, j);
	for (unsigned i=0; i<ni; ++i) {
	  blcluster* sonB = blB->getson(i, j);
	  if(haar) haarSon = haar->son(i,j);
	  for (unsigned k=0; k<j; ++k) {
	    blcluster *sonX = blB->getson(i, k);
	    blcluster *sonU = blU->getson(k, j);
	    mltaGeHGeH((T) -1.0, sonX, X, sonU, U, sonB, B, eps, rankmax,
		    haarSon);
	  }
	  GeHUtH_solve_(son, U, sonB, B, X, eps, rankmax, haarSon);
	  delete haarSon;
	}
      }
    }
}

// solve U^H D X = B (U as returned by SPTRF, B dense), X is stored in B
template<class T> static
void solve_utrhDdns(unsigned mB, unsigned nB, T* U, int* piv,
		    T* B, unsigned ldB)
{
  int k = 0, kc = 1;
  while (k<(int)mB) {
    if (piv[k]>0) { // 1x1 diagonal block
      int kp = piv[k]-1;
      if (kp!=k) blas::swap(nB, B+k, ldB, B+kp, ldB);
      
      unsigned idx1 = kc;
      for (unsigned i=k+1; i<mB; ++i) {
	unsigned idx2 = 0;
	for (unsigned j=0; j<nB; ++j) {
	  B[idx2+i] -= U[idx1] * B[idx2+k];
	  idx2 += ldB;
	}
	idx1 += i+1;
      }

      // D^-1 X
      blas::scal(nB, (T)1.0/U[kc-k-1], B+k, ldB);

      ++k;
      kc += k+2;
    } else { // 2x2 diagonal block
      int kp = -piv[k]-1;
      if (kp!=k+1) blas::swap(nB, B+k+1, ldB, B+kp, ldB);

      unsigned idx1 = kc+k+2;
      for (unsigned i=k+2; i<mB; ++i) {
	unsigned idx2 = 0;
	for (unsigned j=0; j<nB; ++j) {
	  B[idx2+i] -= (U[idx1]*B[idx2+k] + U[idx1+1]*B[idx2+k+1]);
	  idx2 += ldB;
	}
	idx1 += i+1;
      }

      // D^-1 X
      T d12 = -U[kc], d11 = U[kc+1], d22 = U[kc-k-1];
      T dinv = d11*d22 - d12*d12;
      unsigned idx = k;
      for (unsigned j=0; j<nB; ++j) {
	T tmp = (d12*B[idx] + d22*B[idx+1]) / dinv;
	B[idx] = (d11*B[idx] + d12*B[idx+1]) / dinv;
	B[idx+1] = tmp;
	idx += ldB;
      }

      kc += 2*k+7;
      k += 2;
    }
  }
}

// solve U X = B (U as returned by SPTRF, B dense), X is stored in B
template<class T> static
void solve_utrdns(unsigned mB, unsigned nB, T* U, int* piv,
		  T* B, unsigned ldB)
{
  int k = mB-1;
  int kc = mB*(mB+3)/2-1;

  while (k>=0) {
    if (piv[k]>0) { // 1x1 diagonal block
      
      unsigned idx1 = kc;
      for (unsigned i=k+1; i<mB; ++i) {
	unsigned idx2 = 0;
	for (unsigned j=0; j<nB; ++j) {
	  B[idx2+k] -= U[idx1] * B[idx2+i];
	  idx2 += ldB;
	}
	idx1 += i+1;
      }

      int kp = piv[k]-1;
      if (kp!=k) blas::swap(nB, B+k, ldB, B+kp, ldB);

      --k;
    } else { // 2x2 diagonal block

      unsigned idx1 = kc;
      for (unsigned i=k+1; i<mB; ++i) {
	unsigned idx2 = 0;
	for (unsigned j=0; j<nB; ++j) {
	  B[idx2+k-1] -= U[idx1-1] * B[idx2+i];
	  B[idx2+k] -= U[idx1] * B[idx2+i];
	  idx2 += ldB;
	}
	idx1 += i+1;
      }

      int kp = -piv[k]-1;
      if (kp!=k) blas::swap(nB, B+k, ldB, B+kp, ldB);
      
      kc -= k+2;
      k -= 2;
    }
    kc -= k+3;
  }
}


// solve U X = B (U H-matrix (blocks as returned by SPTRF), B dense)
// backward substitution; B contains solution on exit
template<class T> static
void UtHGeM_solve_(blcluster* blU, mblock<T>** U, int* piv,
		   unsigned p, T* B, unsigned ldB)
{
  if (blU->isleaf())
    solve_utrdns(blU->getn1(), p, blU->data(U), piv+blU->getb1(), B, ldB); 
  else {
    unsigned ns = blU->getnrs();
    T* Bi = B + blU->getn1();
    for (unsigned i=ns; i>0; --i) {
      blcluster *son = blU->getson(i-1, i-1);
      T* Bk = Bi;
      Bi -= son->getn1();
      for (unsigned k=i; k<ns; ++k) {
        blcluster *sonU = blU->getson(i-1, k);
        mltaGeHGeM((T) -1.0, sonU, U, p, Bk, ldB, Bi, ldB);
        Bk += sonU->getn2();
      }
      UtHGeM_solve_(son, U, piv, p, Bi, ldB);
    }
  }
}

// solve U^H D X = B (U H-matrix (blocks as returned by SPTRF), B dense)
// forward substitution; B contains solution on exit
template<class T> static
void UtHhDdns_solve_(blcluster* blU, mblock<T>** U, int* piv,
		     unsigned p, T* B, unsigned ldB)
{
  if (blU->isleaf())
    solve_utrhDdns(blU->getn1(), p, blU->data(U), piv+blU->getb1(), B, ldB);
  else {
    unsigned ns = blU->getncs();
    T* Bi = B;
    for (unsigned i=0; i<ns; ++i) {
      T* Bk = B;
      for (unsigned k=0; k<i; ++k) {
        blcluster *sonU = blU->getson(k, i);
	blcluster *sonD = blU->getson(k, k);
        mltaGeHhDiHGeM((T) -1.0, U, sonU, sonD, piv, p, Bk, ldB, Bi, ldB);
        Bk += sonU->getn1();
      }
      blcluster *son = blU->getson(i, i);
      UtHhDdns_solve_(son, U, piv, p, Bi, ldB);
      Bi += son->getn2();
    }
  }
}

// solve U^H D X = B (U H-matrix (blocks as returned by SPTRF), B H-matrix)
template<class T> static
void UtHhDH_solve_(mblock<T>** U, blcluster* blB, blcluster* blU, int* piv,
		   double eps, unsigned rankmax)
{
  if (blB->isleaf()) {
    T* dataB = blB->data(U);
    unsigned mB = blB->getn1(), k;

    if (blB->isLrM(U)) k = blB->rank(U);
    else k = blB->getn2();

    UtHhDdns_solve(blU, U, piv, k, dataB, mB);

  } else { // then U cannot be a leaf
    unsigned n1 = blB->getnrs(), n2 = blB->getncs();
    for (unsigned i=0; i<n1; ++i) {
      blcluster* son = blU->getson(i, i);
      for (unsigned j=0; j<n2; ++j) {
	blcluster* sonB = blB->getson(i, j);
	for (unsigned k=0; k<i; ++k) {
	  blcluster* sonU = blU->getson(k, i);
	  blcluster* sonX = blB->getson(k, j);
	  blcluster* sonD = blU->getson(k, k);
	  mltaGeHhDiHGeH((T) -1.0, U, sonU, sonD, piv, sonX, sonB, eps, rankmax);
	}
	UtHhDH_solve(U, sonB, son, piv, eps, rankmax);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen

void LtHGeM_solve(blcluster* blL, mblock<double>** L, unsigned p, double* B,
                  unsigned ldB)
{
  LtHGeM_solve_(blL, L, p, B, ldB);
}
void LtHGeM_solve(blcluster* blL, mblock<float>** L, unsigned p, float* B,
                  unsigned ldB)
{
  LtHGeM_solve_(blL, L, p, B, ldB);
}
void LtHGeM_solve(blcluster* blL, mblock<dcomp>** L, unsigned p, dcomp* B,
                  unsigned ldB)
{
  LtHGeM_solve_(blL, L, p, B, ldB);
}
void LtHGeM_solve(blcluster* blL, mblock<scomp>** L, unsigned p, scomp* B,
                  unsigned ldB)
{
  LtHGeM_solve_(blL, L, p, B, ldB);
}

void LtHhGeM_solve(blcluster* blL, mblock<double>** L, unsigned p,
                   double* B, unsigned ldB)
{
  LtHhGeM_solve_(blL, L, p, B, ldB);
}
void LtHhGeM_solve(blcluster* blL, mblock<float>** L, unsigned p, float* B,
                   unsigned ldB)
{
  LtHhGeM_solve_(blL, L, p, B, ldB);
}
void LtHhGeM_solve(blcluster* blL, mblock<dcomp>** L, unsigned p, dcomp* B,
                   unsigned ldB)
{
  LtHhGeM_solve_(blL, L, p, B, ldB);
}
void LtHhGeM_solve(blcluster* blL, mblock<scomp>** L, unsigned p, scomp* B,
                   unsigned ldB)
{
  LtHhGeM_solve_(blL, L, p, B, ldB);
}

void UtHGeM_solve(blcluster* blL, mblock<double>** L, unsigned p, double* B,
                  unsigned ldB)
{
  UtHGeM_solve_(blL, L, p, B, ldB);
}
void UtHGeM_solve(blcluster* blL, mblock<float>** L, unsigned p, float* B,
                  unsigned ldB)
{
  UtHGeM_solve_(blL, L, p, B, ldB);
}
void UtHGeM_solve(blcluster* blL, mblock<dcomp>** L, unsigned p, dcomp* B,
                  unsigned ldB)
{
  UtHGeM_solve_(blL, L, p, B, ldB);
}
void UtHGeM_solve(blcluster* blL, mblock<scomp>** L, unsigned p, scomp* B,
                  unsigned ldB)
{
  UtHGeM_solve_(blL, L, p, B, ldB);
}

void LtHGeH_solve(blcluster* blL, mblock<double>** L, blcluster* blB,
                mblock<double>** B, mblock<double>** X,
                double eps, unsigned rankmax, contBasis<double>* haar)
{
  LtHGeH_solve_(blL, L, blB, B, X, eps, rankmax, haar);
}
void LtHGeH_solve(blcluster* blL, mblock<float>** L, blcluster* blB,
                mblock<float>** B, mblock<float>** X,
                double eps, unsigned rankmax, contBasis<float>* haar)
{
  LtHGeH_solve_(blL, L, blB, B, X, eps, rankmax, haar);
}
void LtHGeH_solve(blcluster* blL, mblock<dcomp>** L, blcluster* blB,
                mblock<dcomp>** B, mblock<dcomp>** X,
                double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  LtHGeH_solve_(blL, L, blB, B, X, eps, rankmax, haar);
}
void LtHGeH_solve(blcluster* blL, mblock<scomp>** L, blcluster* blB,
                mblock<scomp>** B, mblock<scomp>** X,
                double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  LtHGeH_solve_(blL, L, blB, B, X, eps, rankmax, haar);
}

void UtHGeH_solve(blcluster* blU, mblock<double>** U, blcluster* blB,
                mblock<double>** B, double eps, unsigned rankmax)
{
  UtHGeH_solve_(blU, U, blB, B, eps, rankmax);
}

void UtHGeH_solve(blcluster* blU, mblock<float>** U, blcluster* blB,
                mblock<float>** B, double eps, unsigned rankmax)
{
  UtHGeH_solve_(blU, U, blB, B, eps, rankmax);
}

void UtHGeH_solve(blcluster* blU, mblock<dcomp>** U, blcluster* blB,
                mblock<dcomp>** B, double eps, unsigned rankmax)
{
  UtHGeH_solve_(blU, U, blB, B, eps, rankmax);
}

void UtHGeH_solve(blcluster* blU, mblock<scomp>** U, blcluster* blB,
                mblock<scomp>** B, double eps, unsigned rankmax)
{
  UtHGeH_solve_(blU, U, blB, B, eps, rankmax);
}

void GeHUtH_solve(blcluster* blU, mblock<double>** U, blcluster* blB,
                mblock<double>** B, mblock<double>** X,
                double eps, unsigned rankmax, contBasis<double>* haar)
{
  GeHUtH_solve_(blU, U, blB, B, X, eps, rankmax, haar);
}
void GeHUtH_solve(blcluster* blU, mblock<float>** U, blcluster* blB,
                mblock<float>** B, mblock<float>** X,
                double eps, unsigned rankmax, contBasis<float>* haar)
{
  GeHUtH_solve_(blU, U, blB, B, X, eps, rankmax, haar);
}
void GeHUtH_solve(blcluster* blU, mblock<dcomp>** U, blcluster* blB,
                mblock<dcomp>** B, mblock<dcomp>** X,
                double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  GeHUtH_solve_(blU, U, blB, B, X, eps, rankmax, haar);
}
void GeHUtH_solve(blcluster* blU, mblock<scomp>** U, blcluster* blB,
                mblock<scomp>** B, mblock<scomp>** X,
                double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  GeHUtH_solve_(blU, U, blB, B, X, eps, rankmax, haar);
}

void UtHhGeH_solve(blcluster* blU, mblock<double>** U, blcluster* blB,
                 mblock<double>** B, double eps, unsigned rankmax,
                 contBasis<double>* haar)
{
  UtHhGeH_solve_(blU, U, blB, B, eps, rankmax,haar);
}
void UtHhGeH_solve(blcluster* blU, mblock<float>** U, blcluster* blB,
                 mblock<float>** B, double eps, unsigned rankmax,
                 contBasis<float>* haar)
{
  UtHhGeH_solve_(blU, U, blB, B, eps, rankmax, haar);
}
void UtHhGeH_solve(blcluster* blU, mblock<dcomp>** U, blcluster* blB,
                 mblock<dcomp>** B, double eps, unsigned rankmax,
                 contBasis<dcomp>* haar)
{
  UtHhGeH_solve_(blU, U, blB, B, eps, rankmax, haar);
}
void UtHhGeH_solve(blcluster* blU, mblock<scomp>** U, blcluster* blB,
                 mblock<scomp>** B, double eps, unsigned rankmax,
                 contBasis<scomp>* haar)
{
  UtHhGeH_solve_(blU, U, blB, B, eps, rankmax, haar);
}

void UtHhGeM_solve(blcluster* blU, mblock<double>** U,
                   unsigned p, double* B, unsigned ldB)
{
  UtHhGeM_solve_(blU, U, p, B, ldB);
}
void UtHhGeM_solve(blcluster* blU, mblock<float>** U,
                   unsigned p, float* B, unsigned ldB)
{
  UtHhGeM_solve_(blU, U, p, B, ldB);
}
void UtHhGeM_solve(blcluster* blU, mblock<dcomp>** U,
                   unsigned p, dcomp* B, unsigned ldB)
{
  UtHhGeM_solve_(blU, U, p, B, ldB);
}
void UtHhGeM_solve(blcluster* blU, mblock<scomp>** U,
                   unsigned p, scomp* B, unsigned ldB)
{
  UtHhGeM_solve_(blU, U, p, B, ldB);
}


void UtHhDH_solve(mblock<double>** U, blcluster* blB, blcluster* blU,
		  int* piv, double eps, unsigned rankmax)
{
  UtHhDH_solve_(U, blB, blU, piv, eps, rankmax);
}

void UtHhDH_solve(mblock<float>** U, blcluster* blB, blcluster* blU,
		  int* piv, double eps, unsigned rankmax)
{
  UtHhDH_solve_(U, blB, blU, piv, eps, rankmax);
}

void UtHhDH_solve(mblock<dcomp>** U, blcluster* blB, blcluster* blU,
		  int* piv, double eps, unsigned rankmax)
{
  UtHhDH_solve_(U, blB, blU, piv, eps, rankmax);
}

void UtHhDH_solve(mblock<scomp>** U, blcluster* blB, blcluster* blU,
		  int* piv, double eps, unsigned rankmax)
{
  UtHhDH_solve_(U, blB, blU, piv, eps, rankmax);
}

void UtHGeM_solve(blcluster* blU, mblock<double>** U, int* piv,
		  unsigned p, double* B, unsigned ldB)
{
  UtHGeM_solve_(blU, U, piv, p, B, ldB);
}

void UtHGeM_solve(blcluster* blU, mblock<float>** U, int* piv,
		  unsigned p, float* B, unsigned ldB)
{
  UtHGeM_solve_(blU, U, piv, p, B, ldB);
}

void UtHGeM_solve(blcluster* blU, mblock<dcomp>** U, int* piv,
		  unsigned p, dcomp* B, unsigned ldB)
{
  UtHGeM_solve_(blU, U, piv, p, B, ldB);
}
void UtHGeM_solve(blcluster* blU, mblock<scomp>** U, int* piv,
		  unsigned p, scomp* B, unsigned ldB)
{
  UtHGeM_solve_(blU, U, piv, p, B, ldB);
}

void UtHhDdns_solve(blcluster* blU, mblock<double>** U, int* piv,
		    unsigned p, double* B, unsigned ldB)
{
  UtHhDdns_solve_(blU, U, piv, p, B, ldB);
}
void UtHhDdns_solve(blcluster* blU, mblock<float>** U, int* piv,
		    unsigned p, float* B, unsigned ldB)
{
  UtHhDdns_solve_(blU, U, piv, p, B, ldB);
}
void UtHhDdns_solve(blcluster* blU, mblock<dcomp>** U, int* piv,
		    unsigned p, dcomp* B, unsigned ldB)
{
  UtHhDdns_solve_(blU, U, piv, p, B, ldB);
}
void UtHhDdns_solve(blcluster* blU, mblock<scomp>** U, int* piv,
		    unsigned p, scomp* B, unsigned ldB)
{
  UtHhDdns_solve_(blU, U, piv, p, B, ldB);
}
