/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A^H D B  (A,C H-matrices, B low-rank, D diagonal)
// idea: A^H D B = (A^H D U) V^H
template<class T> static
void mltaGeHhDiHLrM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
		       mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
		       blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->isleaf() && blB->isLrM(B) &&
         blC->isnleaf());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[(nA+mB)*rankB], * const tmp1 = tmp + mB*rankB;
    assert(tmp!=NULL);

    mltDiHGeM(rankB, D, blD, piv, dataB, tmp, mB);

    blas::setzero(nA*rankB, tmp1);
    if (mltaGeHhGeM(d, blA, A, rankB, tmp, mB, tmp1, nA))
      addGeHLrM(blC, C, eps, rankmax, rankB, tmp1, nA, dataB+mB*rankB, nB);

    delete [] tmp;
  }
}


// C += d A^H D B  (A H-matrix, B low-rank, C mblock, D diagonal)
// idea: A^H D B = (A^H D U) V^H
template<class T> static
void mltaGeHhDiHLrM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			 mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			 mblock<T>* mblC, double eps, unsigned rankmax)
{
  assert(blA->getn1()==blB->getn1() && blB->isleaf() && blB->isLrM(B));

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[(nA+mB)*rankB], * const tmp1 = tmp + mB*rankB;
    assert(tmp!=NULL);

    mltDiHGeM(rankB, D, blD, piv, dataB, tmp, mB);

    blas::setzero(nA*rankB, tmp1);
    if (mltaGeHhGeM(d, blA, A, rankB, tmp, mB, tmp1, nA))
      mblC->addLrM(rankB, tmp1, nA, dataB+mB*rankB, nB, eps, rankmax);
    delete [] tmp;
  }
}


// C += d A^H D B  (A,C H-matrices, B dense, D diagonal)
template<class T> static
void mltaGeHhDiHGeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
		       mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
		       blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->isleaf() && blB->isGeM(B) &&
         blC->isnleaf());

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[nB*(nA+mB)], * const tmp1 = tmp + nB*mB;
  assert(tmp!=NULL);

  mltDiHGeM(nB, D, blD, piv, dataB, tmp, mB);

  blas::setzero(nA*nB, tmp1);
  if (mltaGeHhGeM(d, blA, A, nB, tmp, mB, tmp1, nA))
    addGeHGeM(blC, C, tmp1, nA, eps, rankmax);

  delete [] tmp;
}


// C += d A^H D B  (A,C H-matrices, B dense, C mblock, D diagonal)
template<class T> static
void mltaGeHhDiHGeM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			 mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			 mblock<T>* mblC, double eps, unsigned rankmax)
{
  assert(blA->getn1()==blB->getn1() && blB->isleaf() && blB->isGeM(B));

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[nB*(nA+mB)], * const tmp1 = tmp + nB*mB;
  assert(tmp!=NULL);

  mltDiHGeM(nB, D, blD, piv, dataB, tmp, mB);

  blas::setzero(nA*nB, tmp1);
  if (mltaGeHhGeM(d, blA, A, nB, tmp, mB, tmp1, nA)) 
    mblC->addGeM(tmp1, nA, eps, rankmax);

  delete [] tmp;
}


// C += d A^H D B  (A low-rank, B,C H-matrices, D diagonal)
// idea: A^H D B = V (B^H D U)^H
template<class T> static
void mltaLrMhDiHgeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
		       mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
		       blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->isleaf() && blA->isLrM(A) &&
         blC->isnleaf());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned nB = blB->getn2();

    T* const tmp = new T[(mA+nB)*rankA], * const tmp1 = tmp + mA*rankA;
    assert(tmp!=NULL);

    mltDiHGeM(rankA, D, blD, piv, dataA, tmp, mA);

    blas::setzero(nB*rankA, tmp1);
    if (mltaGeHhGeM(d, blB, B, rankA, tmp, mA, tmp1, nB)) 
      addGeHLrM(blC, C, eps, rankmax, rankA, dataA+mA*rankA, nA, tmp1, nB);

    delete [] tmp;
  }
}


// C += d A^H D B  (A low-rank, B H-matrix, C mblock, D diagonal)
// idea: A^H D B = V (B^H D U)^H
template<class T> static
void mltaLrMhDiHGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			 mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			 mblock<T>* mblC, double eps, unsigned rankmax)
{
  assert(blA->getn1()==blB->getn1() && blA->isleaf() && blA->isLrM(A));

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned nB = blB->getn2();

    T* const tmp = new T[(mA+nB)*rankA], * const tmp1 = tmp + mA*rankA;
    assert(tmp!=NULL);

    mltDiHGeM(rankA, D, blD, piv, dataA, tmp, mA);

    blas::setzero(nB*rankA, tmp1);
    if (mltaGeHhGeM(d, blB, B, rankA, tmp, mA, tmp1, nB)) 
      mblC->addLrM(rankA, dataA+mA*rankA, nA, tmp1, nB, eps, rankmax);

    delete [] tmp;
  }
}


// C += d A^H D B  (A dense, B H-matrix, D diagonal)
// idea: A^H D B = (B^H D A)^H
template<class T> static
inline void mltaGeMhDiHGeH_toGeH_(T d, blcluster* blA, mblock<T>** A,
			      blcluster* blD, mblock<T>** D, int* piv,
			      blcluster* blB, mblock<T>** B, blcluster* blC,
			      mblock<T>** C, double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->isleaf() && blA->isGeM(A) &&
         blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn2();

  T *const tmp = new T[nA*(MAX(mA,nB)+nB)], * const tmp1 = tmp+nA*MAX(mA,nB);
  assert(tmp!=NULL);

  mltDiHGeM(nA, D, blD, piv, dataA, tmp, mA);

  blas::setzero(nB*nA, tmp1);
  if (mltaGeHhGeM(d, blB, B, nA, tmp, mA, tmp1, nB)) {
    blas::transpose(nB, nA, tmp1, tmp);
    addGeHGeM(blC, C, tmp, nA, eps, rankmax);
  }

  delete [] tmp;
}


// C += d A^H D B  (A dense, B H-matrix, C mblock, D diagonal)
// idea: A^H D B = (B^H D A)^H
template<class T> static
void mltaGeMhDiHGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			 mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			 mblock<T>* mblC, double eps, unsigned rankmax)
{
  assert(blA->getn1()==blB->getn1() && blA->isleaf() && blA->isGeM(A));

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[nA*(MAX(mA,nB)+nB)], * const tmp1 = tmp+MAX(mA,nB)*nA;
  assert(tmp!=NULL);

  mltDiHGeM(nA, D, blD, piv, dataA, tmp, mA);

  blas::setzero(nB*nA, tmp1);
  if (mltaGeHhGeM(d, blB, B, nA, tmp, mA, tmp1, nB)) {
    blas::transpose(nB, nA, tmp1, tmp);
    mblC->addGeM(tmp, nA, eps, rankmax);
  }

  delete [] tmp;
}


// C += d A^H D B  (A,B H-matrices, C mblock, D diagonal)
template<class T> static
void mltaGeHhDiHGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blD, 
		       mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
		       mblock<T>* mblC, double eps, unsigned rankmax)
{
  if (blB->isleaf() && blB->isLrM(B))
    mltaGeHhDiHLrM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
  else if (blA->isleaf() && blA->isLrM(A))
    mltaLrMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
  else if (blB->isleaf() && blB->isGeM(B))
    mltaGeHhDiHGeM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
  else if (blA->isleaf() && blA->isGeM(A))
    mltaGeMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
  else { // A and B are no leaves, hence the clusters can be subdivided
    assert(blA->getnrs()==blB->getnrs());
    unsigned ns1 = blA->getncs(), nsp = blA->getnrs(), ns2 = blB->getncs();

    // unify first block row
    const unsigned m0 = blA->getson(0,0)->getn2();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0,0)->getn2();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2, *bl3;

    for (unsigned k=0; k<nsp; ++k) {
      bl1 = blA->getson(k, 0);
      bl2 = blB->getson(k, 0);
      bl3 = blD->getson(k, k);
      mltaGeHhDiHGeH_toMbl_(d, bl1, A, bl3, D, piv, bl2, B, R, eps, rankmax);
    }

    for (unsigned j=1; j<ns2; ++j) {
      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);

      for (unsigned k=0; k<nsp; ++k) {
        bl1 = blA->getson(k, 0);
        bl2 = blB->getson(k, j);
	bl3 = blD->getson(k, k);
        mltaGeHhDiHGeH_toMbl_(d, bl1, A, bl3, D, piv, bl2, B, &R2, eps, rankmax);
      }

      mblock<T>* R1 = R;
      n += nj;
      R = new mblock<T>(m0, n);
      R->unify_cols(eps, rankmax, *R1, R2);

      delete R1;
    }

    for (unsigned i=1; i<ns1; ++i) {
      const unsigned mi = blA->getson(0, i)->getn2();
      n = n0;
      mblock<T>* Rz = new mblock<T>(mi, n0);

      for (unsigned k=0; k<nsp; ++k) {
        bl1 = blA->getson(k, i);
        bl2 = blB->getson(k, 0);
	bl3 = blD->getson(k, k);
        mltaGeHhDiHGeH_toMbl_(d, bl1, A, bl3, D, piv, bl2, B, Rz, eps, rankmax);
      }

      for (unsigned j=1; j<ns2; ++j) {
        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);

        for (unsigned k=0; k<nsp; ++k) {
          bl1 = blA->getson(k, i);
          bl2 = blB->getson(k, j);
	  bl3 = blD->getson(k, k);
          mltaGeHhDiHGeH_toMbl_(d, bl1, A, bl3, D, piv, bl2, B, &R2, eps, rankmax);
        }

        mblock<T>* R1 = Rz;
        n += nj;
        Rz = new mblock<T>(mi, n);
        Rz->unify_cols(eps, rankmax, *R1, R2);

        delete R1;
      }

      mblock<T>* R1 = R;
      m += mi;
      R = new mblock<T>(m, n);
      R->unify_rows(eps, rankmax, *R1, *Rz);

      delete R1;
      delete Rz;
    }

    mblC->addMbl(eps, rankmax, R);
    delete R;
  }
}


// C += d A^H D B  (A,B,C H-matrices, D diagonal)
template<class T> static
void mltaGeHhDiHGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
		mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
		blcluster* blC, mblock<T>** C, double eps, unsigned rankmax)
{
  if (blC->isleaf())                 // C is a leaf
    mltaGeHhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, C[blC->getidx()],
		      eps, rankmax);
  else {                             // C is not a leaf
    if (blA->isnleaf() && blB->isnleaf()) {
      unsigned ns1 = blA->getncs(), nsp = blA->getnrs(), ns2 = blB->getncs();
      for (unsigned i=0; i<ns1; ++i) {
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
          for (unsigned k=0; k<nsp; ++k) {
            blcluster* bl1 = blA->getson(k, i);
            blcluster* bl2 = blB->getson(k, j);
            blcluster* bl3 = blD->getson(k, k);
            mltaGeHhDiHGeH_(d, bl1, A, bl3, D, piv, bl2, B, bl, C, eps, rankmax);
          }
        }
      }
    }
    else if (blB->isleaf() && blB->isLrM(B))
      mltaGeHhDiHLrM_toGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMhDiHgeH_toGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHhDiHGeM_toGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    else {                           // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMhDiHGeH_toGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    }
  }
}

// multiplies diagonal matrix U (as generated by DSPTRF) 
// by dense matrix B and stores result in C
// C = U B, B is mB x nB
template<class T> static
void mltDiMGeM_(unsigned mB, unsigned nB, T* U, int* const piv,
	       T* B, T* C, unsigned ldB)
{
  unsigned k=0, kc=0;
  while (k<mB) {
    if (piv[k]>0) { // 1x1 diagonal block
      T d = U[kc];      
      unsigned idx = k;
      for (unsigned i=0; i<nB; ++i) {
	C[idx] = d * B[idx];
	idx += ldB;
      }
      ++k;
    } else { // 2x2 diagonal block
      // multiply by upper triang. diagonal matrix
      T d11 = U[kc], d22 = U[kc+k+1], d12 = U[kc+k];
      unsigned idx = k;
      for (unsigned j=0; j<nB; ++j) {
	C[idx] = d11 * B[idx] + d12 * B[idx+1];
	C[idx+1] = d12 * B[idx] + d22 * B[idx+1];
	idx += ldB;
      } 
      kc += k+2;
      k += 2;
    }
    kc += k+1;
  }
}

// multiplies diagonal H matrix A (as returned by DSPTRF)
// by dense matrix B and stores result in C
// C = A B, nB -> cols of B
template<class T> static
void mltDiHGeM_(unsigned nB, mblock<T>** A, blcluster* blD,
		   int* const piv, T* B, T* C, unsigned ldB)
{
  if (blD->isleaf())
    mltDiMGeM_(blD->getn1(), nB, blD->data(A), piv+blD->getb1(), B, C, ldB);
  else {
    unsigned ns = blD->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = blD->getson(i, i);
      mltDiHGeM_(nB, A, son, piv, B, C, ldB);
      unsigned n1 = son->getn1();
      C += n1; B += n1;
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHhDiHGeH(double d, blcluster* blA, mblock<double>** A, blcluster* blD,
	       mblock<double>** D, int* piv, blcluster* blB, mblock<double>** B,
	       blcluster* blC, mblock<double>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}

void mltaGeHhDiHGeH(float d, blcluster* blA, mblock<float>** A, blcluster* blD,
	       mblock<float>** D, int* piv, blcluster* blB, mblock<float>** B,
	       blcluster* blC, mblock<float>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}

void mltaGeHhDiHGeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blD,
	       mblock<scomp>** D, int* piv, blcluster* blB, mblock<scomp>** B,
	       blcluster* blC, mblock<scomp>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}

void mltaGeHhDiHGeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blD,
	       mblock<dcomp>** D, int* piv, blcluster* blB, mblock<dcomp>** B,
	       blcluster* blC, mblock<dcomp>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}



void mltaGeHhDiHGeH(double d, mblock<double>** A, blcluster* blA, blcluster* blD,
	       int* piv, blcluster* blB, blcluster* blC,
	       double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

void mltaGeHhDiHGeH(float d, mblock<float>** A, blcluster* blA, blcluster* blD,
	       int* piv, blcluster* blB, blcluster* blC,
	       double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

void mltaGeHhDiHGeH(scomp d, mblock<scomp>** A, blcluster* blA, blcluster* blD,
	       int* piv, blcluster* blB, blcluster* blC,
	       double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

void mltaGeHhDiHGeH(dcomp d, mblock<dcomp>** A, blcluster* blA, blcluster* blD,
	       int* piv, blcluster* blB, blcluster* blC,
	       double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}




////////////////////////////////////
//not required to be extern atm
void mltaGeHhDiHGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
			  blcluster* blD, mblock<double>** D, int* piv,
			  blcluster* blB, mblock<double>** B, mblock<double>* mblC,
			  double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
		      blcluster* blD, mblock<float>** D, int* piv,
		      blcluster* blB, mblock<float>** B, mblock<float>* mblC,
		      double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
		      blcluster* blD, mblock<dcomp>** D, int* piv,
		      blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
		      double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
		      blcluster* blD, mblock<scomp>** D, int* piv,
		      blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
		      double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}





void mltaGeHhDiHLrM_toMbl(double d, blcluster* blA, mblock<double>** A,
			blcluster* blD, mblock<double>** D, int* piv,
			blcluster* blB, mblock<double>** B, mblock<double>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHLrM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHLrM_toMbl(float d, blcluster* blA, mblock<float>** A,
			blcluster* blD, mblock<float>** D, int* piv,
			blcluster* blB, mblock<float>** B, mblock<float>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHLrM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHLrM_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
			blcluster* blD, mblock<dcomp>** D, int* piv,
			blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHLrM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHLrM_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
			blcluster* blD, mblock<scomp>** D, int* piv,
			blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHLrM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeM_toMbl(double d, blcluster* blA, mblock<double>** A,
			blcluster* blD, mblock<double>** D, int* piv,
			blcluster* blB, mblock<double>** B, mblock<double>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHGeM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeM_toMbl(float d, blcluster* blA, mblock<float>** A,
			blcluster* blD, mblock<float>** D, int* piv,
			blcluster* blB, mblock<float>** B, mblock<float>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHGeM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeM_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
			blcluster* blD, mblock<dcomp>** D, int* piv,
			blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHGeM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeHhDiHGeM_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
			blcluster* blD, mblock<scomp>** D, int* piv,
			blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeHhDiHGeM_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeMhDiHGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
			blcluster* blD, mblock<double>**D, int* piv,
			blcluster* blB, mblock<double>**B, mblock<double>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeMhDiHGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
			blcluster* blD, mblock<float>**D, int* piv,
			blcluster* blB, mblock<float>**B, mblock<float>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeMhDiHGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
			blcluster* blD, mblock<dcomp>**D, int* piv,
			blcluster* blB, mblock<dcomp>**B, mblock<dcomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaGeMhDiHGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
			blcluster* blD, mblock<scomp>**D, int* piv,
			blcluster* blB, mblock<scomp>**B, mblock<scomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaGeMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaLrMhDiHGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
			blcluster* blD, mblock<double>** D, int* piv,
			blcluster* blB, mblock<double>** B, mblock<double>* mblC,
			double eps, unsigned rankmax)
{
  mltaLrMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaLrMhDiHGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
			blcluster* blD, mblock<float>** D, int* piv,
			blcluster* blB, mblock<float>** B, mblock<float>* mblC,
			double eps, unsigned rankmax)
{
  mltaLrMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaLrMhDiHGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
			blcluster* blD, mblock<dcomp>** D, int* piv,
			blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaLrMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltaLrMhDiHGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
			blcluster* blD, mblock<scomp>** D, int* piv,
			blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
			double eps, unsigned rankmax)
{
  mltaLrMhDiHGeH_toMbl_(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
}

void mltDiHGeM(unsigned nB, mblock<double>** A, blcluster* blD,
		  int* const piv, double* B, double* C, unsigned ldB)
{
  mltDiHGeM_(nB, A, blD, piv, B, C, ldB);
}

void mltDiHGeM(unsigned nB, mblock<float>** A, blcluster* blD,
		  int* const piv, float* B, float* C, unsigned ldB)
{
  mltDiHGeM_(nB, A, blD, piv, B, C, ldB);
}

void mltDiHGeM(unsigned nB, mblock<scomp>** A, blcluster* blD,
		  int* const piv, scomp* B, scomp* C, unsigned ldB)
{
  mltDiHGeM_(nB, A, blD, piv, B, C, ldB);
}

void mltDiHGeM(unsigned nB, mblock<dcomp>** A, blcluster* blD,
		  int* const piv, dcomp* B, dcomp* C, unsigned ldB)
{
  mltDiHGeM_(nB, A, blD, piv, B, C, ldB);
}
