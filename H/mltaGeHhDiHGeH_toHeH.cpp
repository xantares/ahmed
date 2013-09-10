/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A^H D B  (A H-matrices, B low-rank, C symm. H-matrix)
// idea: A^H D B = (A^H D U) V^H
template<class T> static
void mltaGeHhDiHLrM_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			  mblock<T>** D, int* piv, blcluster* blB,
			  mblock<T>** B, blcluster* blC, mblock<T>** C,
			  double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->isleaf() && blB->isLrM(B) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1(), nB = blB->getn2();

    T* const tmp = new T[(mB+nB)*rankB], * const tmp1 = tmp + mB*rankB;
    assert(tmp!=NULL);

    mltDiHGeM(rankB, D, blD, piv, dataB, tmp, mB);

    blas::setzero(nB*rankB, tmp1);
    if (mltaGeHhGeM(d, blA, A, rankB, tmp, mB, tmp1, nB))
      addHeHLrM(blC, C, eps, rankmax, rankB, tmp1, nB, dataB+mB*rankB, nB);

    delete [] tmp;
  }
}


// C += d A^H D B  (A H-matrices, B dense, C symm. H-matrix, D diagonal)
template<class T> static
void mltaGeHhDiHGeM_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			  mblock<T>** D, int* piv, blcluster* blB, mblock<T>**B,
			  blcluster* blC, mblock<T>** C,
			  double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->isleaf() && blB->isGeM(B) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1(), nB = blB->getn2();

  T* const tmp = new T[nB*(nB+mB)], * const tmp1 = tmp + mB*nB;
  assert(tmp!=NULL);

  mltDiHGeM(nB, D, blD, piv, dataB, tmp, mB);

  blas::setzero(nB*nB, tmp1);
  if (mltaGeHhGeM(d, blA, A, nB, tmp, mB, tmp1, nB))
    addHeHGeM(blC, C, tmp1, nB, eps, rankmax);

  delete [] tmp;
}


// C += d A^H D B  (A low-rank, B H-matrices, C symm. H-matrix, D diagonal)
// idea: A^H D B = V (B^H D U)^H
template<class T> static
void mltaLrMhDiHGeH_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			  mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			  blcluster* blC, mblock<T>** C,
			  double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->isleaf() && blA->isLrM(A) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1(), nA = blA->getn2();

    T* const tmp = new T[(mA+nA)*rankA], * const tmp1 = tmp + mA*rankA;
    assert(tmp!=NULL);

    mltDiHGeM(rankA, D, blD, piv, dataA, tmp, mA);

    blas::setzero(nA*rankA, tmp1);
    if (mltaGeHhGeM(d, blB, B, rankA, tmp, mA, tmp1, nA))
      addHeHLrM(blC, C, eps, rankmax, rankA, dataA+mA*rankA, nA, tmp1, nA);

    delete [] tmp;
  }
}


// C += d A^H D B  (A dense, B H-matrix, C symm. H-matrix, D diagonal)
// idea: A^H D B = B^H D A
template<class T> static
void mltaGeMhDiHGeH_toHeH_(T d, blcluster* blA, mblock<T>**A, blcluster* blD,
			  mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			  blcluster* blC, mblock<T>** C,
			  double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->isleaf() && blA->isGeM(A) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1(), nA = blA->getn2();

  T* const tmp = new T[nA*(nA+mA)], * const tmp1 = tmp + mA*nA;
  assert(tmp!=NULL);

  mltDiHGeM(nA, D, blD, piv, dataA, tmp, mA);

  blas::setzero(nA*nA, tmp1);
  if (mltaGeHhGeM(d, blB, B, nA, tmp, mA, tmp1, nA)) 
    addHeHGeM(blC, C, tmp1, nA, eps, rankmax);

  delete [] tmp;
}


// C += d A^H D B  (A,B H-matrices, C symm. H-matrix, D diagonal)
template<class T> static
void mltaGeHhDiHGeH_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blD,
			mblock<T>** D, int* piv, blcluster* blB, mblock<T>** B,
			blcluster* blC, mblock<T>** C,
			double eps, unsigned rankmax)
{
  if (blC->isleaf()) {
    mblock<T>* const mblC = C[blC->getidx()];             // C is a leaf

    // C is a diag block. Hence, one of the blocks blA, blB must be a leaf
    if (blB->isleaf() && blB->isLrM(B))
      mltaGeHhDiHLrM_toMbl(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMhDiHGeH_toMbl(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHhDiHGeM_toMbl(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
    else {                // blA must be the leaf and it must be dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMhDiHGeH_toMbl(d, blA, A, blD, D, piv, blB, B, mblC, eps, rankmax);
    }
  } else {                           // C is not a leaf
    if (blA->isnleaf() && blB->isnleaf()) {
      unsigned ns1 = blA->getncs(), nsp = blA->getnrs(), ns2 = blB->getncs();
      for (unsigned i=0; i<ns1; ++i) {
        blcluster* bl = blC->getson(i, i);

        for (unsigned k=0; k<nsp; ++k) {
          blcluster* bl1 = blA->getson(k, i);
          blcluster* bl2 = blB->getson(k, i);
          blcluster* bl3 = blD->getson(k, k);
          mltaGeHhDiHGeH_toHeH_(d, bl1, A, bl3, D, piv, bl2, B, bl, C,
			     eps, rankmax);
        }

        for (unsigned j=i+1; j<ns2; ++j) {
          bl = blC->getson(i, j);

          for (unsigned k=0; k<nsp; ++k) {
            blcluster* bl1 = blA->getson(k, i);
            blcluster* bl2 = blB->getson(k, j);
	    blcluster* bl3 = blD->getson(k, k);
            mltaGeHhDiHGeH(d, bl1, A, bl3, D, piv, bl2, B, bl, C, eps, rankmax);
          }
        }
      }
    }
    else if (blB->isleaf() && blB->isLrM(B))
      mltaGeHhDiHLrM_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMhDiHGeH_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHhDiHGeM_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    else {                           // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMhDiHGeH_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
    }
  }
}



///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHhDiHGeH_toHeH(double d, blcluster* blA, mblock<double>** A,
		       blcluster* blD, mblock<double>** D, int* piv,
		       blcluster* blB, mblock<double>** B, blcluster* blC,
		       mblock<double>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}

void mltaGeHhDiHGeH_toHeH(float d, blcluster* blA, mblock<float>** A,
		       blcluster* blD, mblock<float>** D, int* piv,
		       blcluster* blB, mblock<float>** B, blcluster* blC,
		       mblock<float>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}

void mltaGeHhDiHGeH_toHeH(dcomp d, blcluster* blA, mblock<dcomp>** A,
		       blcluster* blD, mblock<dcomp>** D, int* piv,
		       blcluster* blB, mblock<dcomp>** B, blcluster* blC,
		       mblock<dcomp>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}

void mltaGeHhDiHGeH_toHeH(scomp d, blcluster* blA, mblock<scomp>** A,
		       blcluster* blD, mblock<scomp>** D, int* piv,
		       blcluster* blB, mblock<scomp>** B, blcluster* blC,
		       mblock<scomp>** C, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, D, piv, blB, B, blC, C, eps, rankmax);
}


void mltaGeHhDiHGeH_toHeH(double d, mblock<double>** A, blcluster* blA, 
		       blcluster* blD, int* piv, blcluster* blB,
		       blcluster* blC, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

void mltaGeHhDiHGeH_toHeH(float d, mblock<float>** A, blcluster* blA, 
		       blcluster* blD, int* piv, blcluster* blB,
		       blcluster* blC, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

void mltaGeHhDiHGeH_toHeH(dcomp d, mblock<dcomp>** A, blcluster* blA,
		       blcluster* blD, int* piv, blcluster* blB,
		       blcluster* blC, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

void mltaGeHhDiHGeH_toHeH(scomp d, mblock<scomp>** A, blcluster* blA, 
		       blcluster* blD, int* piv, blcluster* blB,
		       blcluster* blC, double eps, unsigned rankmax)
{
  mltaGeHhDiHGeH_toHeH_(d, blA, A, blD, A, piv, blB, A, blC, A, eps, rankmax);
}

