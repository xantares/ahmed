/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B^H  (A H-matrices, B low-rank, C symm. H-matrix)
// idea: A B^H = (A V) U^H
template<class T> static
void mltaGeHLrMh_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                         mblock<T>** B, blcluster* blC, mblock<T>** C,
                         double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blB->isleaf() && blB->isLrM(B) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1(), nB = blB->getn2();

    T* const tmp = new T[mB*rankB];
    assert(tmp!=NULL);

    blas::setzero(mB*rankB, tmp);
    if (mltaGeHGeM(d, blA, A, rankB, dataB + rankB*mB, nB, tmp, mB))
      addHeHLrM(blC, C, eps, rankmax, rankB, tmp, mB, dataB, mB);

    delete [] tmp;
  }
}


// C += d A B^H  (A H-matrices, B dense, C symm. H-matrix)
template<class T> static
void mltaGeHGeMh_toHeH_(T d, blcluster* blA, mblock<T>** A,
                         blcluster* blB, mblock<T>** B, blcluster* blC,
                         mblock<T>** C, double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blB->isleaf() && blB->isGeM(B) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1(), nB = blB->getn2();

  T* const tmp = new T[mB*(mB+nB)];
  assert(tmp!=NULL);

  blas::transpose(mB, nB, dataB, tmp);
  T* const tmp1 = tmp + nB*mB;
  blas::setzero(mB*mB, tmp1);
  if (mltaGeHGeM(d, blA, A, mB, tmp, nB, tmp1, mB))
    addHeHGeM(blC, C, tmp1, mB, eps, rankmax);

  delete [] tmp;
}


// C += d A B^H  (A low-rank, B H-matrices, C symm. H-matrix)
// idea: A B^H = U (B V)^H
template<class T> static
void mltaLrMGeHh_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                         mblock<T>** B, blcluster* blC, mblock<T>** C,
                         double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blA->isleaf() && blA->isLrM(A) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1(), nA = blA->getn2();

    T* const tmp = new T[mA*rankA];
    assert(tmp!=NULL);

    blas::setzero(mA*rankA, tmp);
    if (mltaGeHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, mA))
      addHeHLrM(blC, C, eps, rankmax, rankA, dataA, mA, tmp, mA);

    delete [] tmp;
  }
}


// C += d A B^H  (A dense, B H-matrix, C symm. H-matrix)
// idea: A B^H = B A^H
template<class T> static
void mltaGeMGeHh_toHeH_(T d, blcluster* blA, mblock<T>** A,
                         blcluster* blB, mblock<T>** B, blcluster* blC,
                         mblock<T>** C, double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blA->isleaf() && blA->isGeM(A) &&
         blC->isnleaf() && blC->getn1()==blC->getn2());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1(), nA = blA->getn2();

  T* const tmp = new T[mA*(mA+nA)];
  assert(tmp!=NULL);
  blas::transpose(mA, nA, dataA, tmp);

  T* const tmp1 = tmp + nA*mA;
  blas::setzero(mA*mA, tmp1);
  if (mltaGeHGeM(d, blB, B, mA, tmp, nA, tmp1, mA)) 
    addHeHGeM(blC, C, tmp1, mA, eps, rankmax);

  delete [] tmp;
}


// C += d A B^H  (A,B H-matrices, C symm. H-matrix)
template<class T> static
void mltaGeHGeHh_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, blcluster* blC, mblock<T>** C,
                       double eps, unsigned rankmax)
{
  if (blC->isleaf()) {
    mblock<T>* const mblC = C[blC->getidx()];             // C is a leaf

    // C is a diag block. Hence, one of the blocks blA, blB must be a leaf
    if (blB->isleaf() && blB->isLrM(B))
      mltaGeHLrMh_toMbl(d, blA, A, blB, B, mblC, eps, rankmax);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMGeHh_toMbl(d, blA, A, blB, B, mblC, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHGeMh_toMbl(d, blA, A, blB, B, mblC, eps, rankmax);
    else {                // blA must be the leaf and it must be dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMGeHh_toMbl(d, blA, A, blB, B, mblC, eps, rankmax);
    }
  } else {                            // C is not a leaf
    if (blA->isnleaf() && blB->isnleaf()) {
      unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getnrs();
      for (unsigned i=0; i<ns1; ++i) {
        blcluster* bl = blC->getson(i, i);

        for (unsigned k=0; k<nsp; ++k) {
          blcluster* bl1 = blA->getson(i, k);
          blcluster* bl2 = blB->getson(i, k);
          mltaGeHGeHh_toHeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax);
        }

        for (unsigned j=i+1; j<ns2; ++j) {
          bl = blC->getson(i, j);

          for (unsigned k=0; k<nsp; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(j, k);
            mltaGeHGeHh(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }
        }
      }
    }
    else if (blB->isleaf() && blB->isLrM(B))
      mltaGeHLrMh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMGeHh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHGeMh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else {                                       // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMGeHh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    }
  }
}



///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHGeHh_toHeH(double d, blcluster* blA, mblock<double>** A,
                      blcluster* blB, mblock<double>** B, blcluster* blC,
                      mblock<double>** C, double eps, unsigned rankmax)
{
  mltaGeHGeHh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHGeHh_toHeH(float d, blcluster* blA, mblock<float>** A,
                      blcluster* blB, mblock<float>** B, blcluster* blC,
                      mblock<float>** C, double eps, unsigned rankmax)
{
  mltaGeHGeHh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHGeHh_toHeH(dcomp d, blcluster* blA, mblock<dcomp>** A,
                      blcluster* blB, mblock<dcomp>** B, blcluster* blC,
                      mblock<dcomp>** C, double eps, unsigned rankmax)
{
  mltaGeHGeHh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHGeHh_toHeH(scomp d, blcluster* blA, mblock<scomp>** A,
		       blcluster* blB, mblock<scomp>** B, blcluster* blC,
		       mblock<scomp>** C, double eps, unsigned rankmax)
{
  mltaGeHGeHh_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}
