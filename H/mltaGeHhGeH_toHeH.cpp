/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A^H B  (A H-matrices, B low-rank, C symm. H-matrix)
template<class T> static
void mltaGeHhLrM_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                         mblock<T>** B, blcluster* blC, mblock<T>** C,
                         double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blC->getn1()==blC->getn2() &&
         blB->isleaf() && blB->isLrM(B) && blC->isnleaf());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1(), nB = blB->getn2();

    T* const tmp = new T[nB*rankB];
    assert(tmp!=NULL);

    blas::setzero(nB*rankB, tmp);
    if (mltaGeHhGeM(d, blA, A, rankB, dataB, mB, tmp, nB))
      addHeHLrM(blC, C, eps, rankmax, rankB, tmp, nB, dataB+rankB*mB, nB,
		 haar);
    delete [] tmp;
  }
}


// C += d A^H B  (A H-matrices, B dense, C symm. H-matrix)
template<class T> static
void mltaGeHhGeM_toHeH_(T d, blcluster* blA, mblock<T>** A,
                         blcluster* blB, mblock<T>** B, blcluster* blC,
                         mblock<T>** C, double eps, unsigned rankmax,
                         contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blC->getn1()==blC->getn2() &&
         blB->isleaf() && blB->isGeM(B) && blC->isnleaf());

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1(), nB = blB->getn2();

  T* const tmp = new T[nB*nB];
  assert(tmp!=NULL);

  blas::setzero(nB*nB, tmp);
  if (mltaGeHhGeM(d, blA, A, nB, dataB, mB, tmp, nB)) 
    addHeHGeM(blC, C, tmp, nB, eps, rankmax, haar);
  delete [] tmp;
}


// C += d A^H B  (A low-rank, B H-matrix, C symm. H-matrix)
// idea: A^H B = V (B^H U)^H
template<class T> static
void mltaLrMhGeH_toHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                         mblock<T>** B, blcluster* blC, mblock<T>** C,
                         double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blC->getn1()==blC->getn2() &&
         blA->isleaf() && blA->isLrM(A) && blC->isnleaf());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1(), nA = blA->getn2();

    T* const tmp = new T[nA*rankA];
    assert(tmp!=NULL);

    blas::setzero(rankA*nA, tmp);
    if (mltaGeHhGeM(d, blB, B, rankA, dataA, mA, tmp, nA)) 
      addHeHLrM(blC, C, eps, rankmax, rankA, dataA+rankA*mA, nA, tmp, nA,
		 haar);
    delete [] tmp;
  }
}


// C += d A^H B  (A dense, B H-matrix, C symm. H-matrix)
// idea: A^H B = B^H A
template<class T> static
void mltaGeMhGeH_toHeH_(T d, blcluster* blA, mblock<T>** A,
                         blcluster* blB, mblock<T>** B,
                         blcluster* blC, mblock<T>** C,
                         double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blC->getn1()==blC->getn2() &&
         blA->isleaf() && blA->isGeM(A) && blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1(), nA = blA->getn2();

  T* const tmp = new T[nA*nA];
  assert(tmp!=NULL);

  blas::setzero(nA*nA, tmp);
  if (mltaGeHhGeM(d, blB, B, nA, dataA, mA, tmp, nA))
    addHeHGeM(blC, C, tmp, nA, eps, rankmax, haar);
  delete [] tmp;
}


// C += d A^H B  (A,B H-matrices, C symm. H-matrix)
template<class T> static
void mltaGeHhGeH_toHeH_(T d, blcluster* blA, mblock<T>** A,
                       blcluster* blB, mblock<T>** B, blcluster* blC,
                       mblock<T>** C, double eps, unsigned rankmax,
                       contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blC->getn2());

  if (blC->isleaf()) {
    mblock<T>* const mblC = C[blC->getidx()]; // C is a leaf -> dense symm

    // one of the blocks blA, blB must be a leaf
    if (blB->isleaf() && blB->isLrM(B)){
      mltaGeHhLrM_toMbl(d, blA, A, blB, B, mblC, eps, rankmax, haar);
    } else if (blA->isleaf() && blA->isLrM(A)){
      mltaLrMhGeH_toMbl(d, blA, A, blB, B, mblC, eps, rankmax, haar);
    } else if (blB->isleaf() && blB->isGeM(B)){
      mltaGeHhGeM_toMbl(d, blA, A, blB, B, mblC, eps, rankmax, haar);
    } else {                // blA must be the leaf and it must be dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMhGeH_toMbl(d, blA, A, blB, B, mblC, eps, rankmax, haar);
    }
  } else {                            // C is not a leaf
    if (blA->isnleaf() && blB->isnleaf()) {
      unsigned ns1 = blA->getncs(), nsp = blA->getnrs(), ns2 = blB->getncs();
      contBasis<T>* haarSon = NULL;
      for (unsigned i=0; i<ns1; ++i) {
	blcluster* bl = blC->getson(i, i);
	if (haar) haarSon = haar->son(i,i);
	for (unsigned k=0; k<nsp; ++k) {
	  blcluster* bl1 = blA->getson(k, i);
	  blcluster* bl2 = blB->getson(k, i);
	  mltaGeHhGeH_toHeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
	}
	delete haarSon;

	for (unsigned j=i+1; j<ns2; ++j) {
	  bl = blC->getson(i, j);
	  if (haar) haarSon = haar->son(i,j);

	  for (unsigned k=0; k<nsp; ++k) {
	    blcluster* bl1 = blA->getson(k, i);
	    blcluster* bl2 = blB->getson(k, j);
	    mltaGeHhGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
	  }
	  delete haarSon;
	}
      }
    }
    else if (blB->isleaf() && blB->isLrM(B)){
      mltaGeHhLrM_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    } else if (blA->isleaf() && blA->isLrM(A)){
      mltaLrMhGeH_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    } else if (blB->isleaf() && blB->isGeM(B)){
      mltaGeHhGeM_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    } else {                                       // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMhGeH_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHhGeH_toHeH(double d, blcluster* blA, mblock<double>** A,
                      blcluster* blB, mblock<double>** B, blcluster* blC,
                      mblock<double>** C, double eps, unsigned rankmax,
                      contBasis<double>* haar)
{
  mltaGeHhGeH_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHhGeH_toHeH(float d, blcluster* blA, mblock<float>** A,
                      blcluster* blB, mblock<float>** B, blcluster* blC,
                      mblock<float>** C, double eps, unsigned rankmax,
                      contBasis<float>* haar)
{
  mltaGeHhGeH_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHhGeH_toHeH(dcomp d, blcluster* blA, mblock<dcomp>** A,
                      blcluster* blB, mblock<dcomp>** B, blcluster* blC,
                      mblock<dcomp>** C, double eps, unsigned rankmax,
                      contBasis<dcomp>* haar)
{
  mltaGeHhGeH_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax,  haar);
}

void mltaGeHhGeH_toHeH(scomp d, blcluster* blA, mblock<scomp>** A,
                      blcluster* blB, mblock<scomp>** B, blcluster* blC,
                      mblock<scomp>** C, double eps, unsigned rankmax,
                      contBasis<scomp>* haar)
{
  mltaGeHhGeH_toHeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}
