/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A^H B  (A,C H-matrices, B low-rank)
template<class T> static
void mltaGeHhLrM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
			mblock<T>** B, blcluster* blC, mblock<T>** C,
			double eps, unsigned rankmax, contBasis<T>* haar=NULL)
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
    
    T* const tmp = new T[nA*rankB];
    assert(tmp!=NULL);

    blas::setzero(nA*rankB, tmp);
    if (mltaGeHhGeM(d, blA, A, rankB, dataB, mB, tmp, nA))
      addGeHLrM(blC, C, eps, rankmax, rankB, tmp, nA, dataB+rankB*mB, nB,
		haar);
    
    delete [] tmp;
  }
}


// C += d A^H B  (A H-matrix, B low-rank, C mblock)
template<class T> static
void mltaGeHhLrM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax, 
                        contBasis<T>* haar=NULL)
{
  assert(blA->getn1()==blB->getn1() && blB->isleaf() && blB->isLrM(B));
  
  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[nA*rankB];
    assert(tmp!=NULL);

    blas::setzero(nA*rankB, tmp);

    if (mltaGeHhGeM(d, blA, A, rankB, dataB, mB, tmp, nA)) {

      contLowLevel<T>* haarInfo = NULL;
      unsigned cols = 0;
      T* X = NULL;
      if (haar!=NULL) {
	haarInfo = new contLowLevel<T>(haar);
	cols = haar->getcols();
	X = haar->getX();
	haar->initX();
      }
      mblC->addLrM(rankB, tmp, nA, dataB+rankB*mB, nB, eps, rankmax, haarInfo,
		   X, blB->getn2(), X+cols*blB->getn2(), blA->getn2());
      delete haarInfo;
    }
    
    delete [] tmp;
  }
}


// C += d A^H B  (A,C H-matrices, B dense)
template<class T> static
void mltaGeHhGeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
			mblock<T>** B, blcluster* blC, mblock<T>** C,
			double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
	 blB->getn2()==blC->getn2() && blB->isleaf() && blB->isGeM(B) &&
	 blC->isnleaf());
  
  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned nA = blA->getn2();
  
  T* const tmp = new T[nA*nB];
  assert(tmp!=NULL);

  blas::setzero(nA*nB, tmp);
  if (mltaGeHhGeM(d, blA, A, nB, dataB, mB, tmp, nA))
    addGeHGeM(blC, C, tmp, nA, eps, rankmax, haar);
  
  delete [] tmp;
}


// C += d A^H B  (A,C H-matrices, B dense, C mblock)
template<class T> static
void mltaGeHhGeM_toMbl_(T d, blcluster* blA, mblock<T>** A,
                        blcluster* blB, mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->getn1()==blB->getn1() && blB->isleaf() && blB->isGeM(B));
  
  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned nA = blA->getn2();
  
  T* const tmp = new T[nA*nB];
  assert(tmp!=NULL);
  
  blas::setzero(nA*nB, tmp);
  if (mltaGeHhGeM(d, blA, A, nB, dataB, mB, tmp, nA)) {
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar!=NULL) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    mblC->addGeM(tmp, nA, eps, rankmax, haarInfo, X, blB->getn2(),
		 X+cols*blB->getn2(), blA->getn2());
    delete haarInfo;
  }
  
  delete [] tmp;
}


// C += d A^H B  (A low-rank, B,C H-matrices)
// idea: A^H B = V (B^H U)^H
template<class T> static
void mltaLrMhGeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
			mblock<T>** B, blcluster* blC, mblock<T>** C,
			double eps, unsigned rankmax, contBasis<T>* haar=NULL)
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
    
    T* const tmp = new T[nB*rankA];
    assert(tmp!=NULL);
    
    blas::setzero(nB*rankA, tmp);
    if (mltaGeHhGeM(d, blB, B, rankA, dataA, mA, tmp, nB))
      addGeHLrM(blC, C, eps, rankmax, rankA, dataA+rankA*mA, nA, tmp, nB, haar);
    
    delete [] tmp;
  }
}


// C += d A^H B  (A low-rank, B H-matrix, C mblock<double>)
// idea: A^H B = V (B^H U)^H
template<class T> static
void mltaLrMhGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->getn1()==blB->getn1() && blA->isleaf() && blA->isLrM(A));
	
  
  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned nB = blB->getn2();
    
    T* const tmp = new T[nB*rankA];
    assert(tmp!=NULL);
    
    blas::setzero(nB*rankA, tmp);
    if (mltaGeHhGeM(d, blB, B, rankA, dataA, mA, tmp, nB)) {
      contLowLevel<T>* haarInfo = NULL;
      unsigned cols = 0;
      T* X = NULL;
      if (haar!=NULL) {
	haarInfo = new contLowLevel<T>(haar);
	cols = haar->getcols();
	X = haar->getX();
	haar->initX();
      }
      mblC->addLrM(rankA, dataA+rankA*mA, nA, tmp, nB, eps, rankmax, haarInfo,
		   X, blB->getn2(), X+cols*blB->getn2(), blA->getn2());
      delete haarInfo;
    }
		
    delete [] tmp;
  }
}


// C += d A^H B  (A dense, B H-matrix)
// idea: A^H B = (B^H A)^H
template<class T> static
void mltaGeMhGeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
			mblock<T>** B, blcluster* blC, mblock<T>** C,
			double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn2() && blA->getn1()==blB->getn1() &&
	 blB->getn2()==blC->getn2() && blA->isleaf() && blA->isGeM(A) &&
	 blC->isnleaf());
  
  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn2();
  
  T* const tmp = new T[2*nA*nB];
  assert(tmp!=NULL);
  blas::setzero(nB*nA, tmp);
  if (mltaGeHhGeM(d, blB, B, nA, dataA, mA, tmp, nB)) {
  
    T* const tmp1 = tmp + nA*nB;
    blas::transpose(nB, nA, tmp, tmp1);
    addGeHGeM(blC, C, tmp1, nA, eps, rankmax, haar);
  }
  
  delete [] tmp;
}


// C += d A^H B  (A dense, B H-matrix, C mblock<double>)
// idea: A^H B = (B^H A)^H
template<class T> static
void mltaGeMhGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->getn1()==blB->getn1() && blA->isleaf() && blA->isGeM(A));
  
  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn2();
  
  
  T* const tmp = new T[2*nA*nB];
  assert(tmp!=NULL);
  blas::setzero(nB*nA, tmp);
  if (mltaGeHhGeM(d, blB, B, nA, dataA, mA, tmp, nB)) {  
    T* const tmp1 = tmp + nA*nB;
    blas::transpose(nB, nA, tmp, tmp1);
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar!=NULL) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    mblC->addGeM(tmp1, nA, eps, rankmax, haarInfo, X, blB->getn2(), 
		 X+cols*blB->getn2(), blA->getn2());
    delete haarInfo;
  }
  
  delete [] tmp;
}


// C += d A^H B  (A,B H-matrices, C mblock)
template<class T> static
void mltaGeHhGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
			mblock<T>** B, mblock<T>* mblC,
			double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  if (blB->isleaf() && blB->isLrM(B)){
    mltaGeHhLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  } else if (blA->isleaf() && blA->isLrM(A)){
    mltaLrMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  } else if (blB->isleaf() && blB->isGeM(B)){
    mltaGeHhGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  } else if (blA->isleaf() && blA->isGeM(A)){
    mltaGeMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  } else {         // A and B are no leaves, hence the clusters can be subdivided
    assert(blA->getnrs()==blB->getnrs());
    unsigned ns1 = blA->getncs(), nsp = blA->getnrs(), ns2 = blB->getncs();

    // unify first block row
    const unsigned m0 = blA->getson(0, 0)->getn2();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0, 0)->getn2();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2;

    contBasis<T>* haarSon = NULL;
    if(haar) haarSon = haar->son(0,0);

    for (unsigned k=0; k<nsp; ++k) {
      bl1 = blA->getson(k, 0);
      bl2 = blB->getson(k, 0);
      mltaGeHhGeH_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax, haarSon);
    }
    delete haarSon;

    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    if(haar!=NULL){
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
    }

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);	
      if(haar) haarSon = haar->son(0,j);
      for (unsigned k=0; k<nsp; ++k) {
	bl1 = blA->getson(k, 0);
	bl2 = blB->getson(k, j);
	mltaGeHhGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax, haarSon);
      }
      delete haarSon;

      T* X=NULL;
      if(haar!=NULL){
	X=haar->getX();
	haar->initX();
      }

      mblock<T>* R1 = R;
      n += nj;
      R = new mblock<T>(m0, n);
      R->unify_cols(eps, rankmax, *R1, R2, haarInfo, X, blB->getn2(), // hier nachbessern
		    X+cols*blB->getn2(), blA->getn2());


      delete R1;
    }
    for (unsigned i=1; i<ns1; ++i) {

      const unsigned mi = blA->getson(0, i)->getn2();
      n = n0;
      mblock<T>* Rz = new mblock<T>(mi, n0);
      if(haar) haarSon = haar->son(i,0);
      for (unsigned k=0; k<nsp; ++k) {
	bl1 = blA->getson(k, i);
	bl2 = blB->getson(k, 0);
	mltaGeHhGeH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);
      }
      delete haarSon;

      for (unsigned j=1; j<ns2; ++j) {

	const unsigned nj = blB->getson(0, j)->getn2();
	mblock<T> R2(mi, nj);
	for (unsigned k=0; k<nsp; ++k) {
	  bl1 = blA->getson(k, i);
	  bl2 = blB->getson(k, j);
	  if(haar) haarSon = haar->son(i,j);
	  mltaGeHhGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax, haarSon);
	  delete haarSon;
	}
	T* X = NULL;
	if(haar){
	  X = haar->getX();
	  haar->initX();
	}

	mblock<T>* R1 = Rz;
	n += nj;
	Rz = new mblock<T>(mi, n);
	Rz->unify_cols(eps, rankmax, *R1, R2, haarInfo, X, blB->getn2(),
		       X+cols*blB->getn2()+bl1->getb2()-blA->getb2(), 
		       blA->getn2());
	delete R1;
      }

      mblock<T>* R1 = R;
      m += mi;
      R = new mblock<T>(m, n);
      T* X = NULL;
      if(haar){
	X = haar->getX();
	haar->initX();
      }
      R->unify_rows(eps, rankmax, *R1, *Rz, haarInfo, X, blB->getn2(), 
		    X+cols*blB->getn2(), blA->getn2());
      delete R1;
      delete Rz;
    }
    T* X = NULL;
    if(haar){
      X = haar->getX();
      haar->initX();
    }
    mblC->addMbl(eps, rankmax, R, haarInfo, X, blB->getn2(), 
		 X+cols*blB->getn2(), blA->getn2());
    delete haarInfo;
    delete R;
  }
}


// C += d A B  (A,B,C H-matrices)
template<class T> static
void mltaGeHhGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		  mblock<T>** B, blcluster* blC, mblock<T>** C,
		  double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  if (blC->isleaf()){                // C is a leaf
    mltaGeHhGeH_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax, haar);
  } else {                            // C is not a leaf
    if (!blA->isleaf() && !blB->isleaf()) {
      unsigned ns1 = blA->getncs(), nsp = blA->getnrs(), ns2 = blB->getncs();
      contBasis<T>* haarSon = NULL;
      for (unsigned i=0; i<ns1; ++i)
	for (unsigned j=0; j<ns2; ++j) {
	  blcluster* bl = blC->getson(i, j);
	  if (haar) haarSon = haar->son(i,j);
	  for (unsigned k=0; k<nsp; ++k) {
	    blcluster* bl1 = blA->getson(k, i);
	    blcluster* bl2 = blB->getson(k, j);
	    mltaGeHhGeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
	  }
	  delete haarSon;
	}
    }
    else if (blB->isleaf() && blB->isLrM(B)){
      mltaGeHhLrM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    } else if (blA->isleaf() && blA->isLrM(A)){
      mltaLrMhGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    } else if (blB->isleaf() && blB->isGeM(B)){
      mltaGeHhGeM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    } else {                                       // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMhGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHhGeH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
		 mblock<double>** B, blcluster* blC, mblock<double>** C,
		 double eps, unsigned rankmax, contBasis<double>* haar)
{
  mltaGeHhGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHhGeH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
		 mblock<float>** B, blcluster* blC, mblock<float>** C,
		 double eps, unsigned rankmax, contBasis<float>* haar)
{
  mltaGeHhGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHhGeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
		 mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
		 double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  mltaGeHhGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHhGeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
		 mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
		 double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  mltaGeHhGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}


void mltaGeHhGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
		       blcluster* blB, mblock<double>** B, mblock<double>* mblC,
		       double eps, unsigned rankmax)
{
  mltaGeHhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHhGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
		       blcluster* blB, mblock<float>** B, mblock<float>* mblC,
		       double eps, unsigned rankmax)
{
  mltaGeHhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHhGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
		       blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
		       double eps, unsigned rankmax)
{
  mltaGeHhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHhGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
		       blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
		       double eps, unsigned rankmax)
{
  mltaGeHhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHhLrM_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B, mblock<double>* mblC,
                       double eps, unsigned rankmax, contBasis<double>* haar)
{
  mltaGeHhLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHhLrM_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B, mblock<float>* mblC,
                       double eps, unsigned rankmax, contBasis<float>* haar)
{
  mltaGeHhLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHhLrM_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
                       double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  mltaGeHhLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHhLrM_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
                       double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  mltaGeHhLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}


void mltaGeHhGeM_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax,
                       contBasis<double>* haar)
{
  mltaGeHhGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHhGeM_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B,
                       mblock<float>* mblC, double eps, unsigned rankmax, 
                       contBasis<float>* haar)
{
  mltaGeHhGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHhGeM_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B,
                       mblock<dcomp>* mblC, double eps, unsigned rankmax,
                       contBasis<dcomp>* haar)
{
  mltaGeHhGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, 
                     haar);
}

void mltaGeHhGeM_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B,
                       mblock<scomp>* mblC, double eps, unsigned rankmax, 
                       contBasis<scomp>* haar)
{
  mltaGeHhGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}


void mltaGeMhGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax, 
                       contBasis<double>* haar)
{
  mltaGeMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeMhGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B,
                       mblock<float>* mblC, double eps, unsigned rankmax, 
                       contBasis<float>* haar)
{
  mltaGeMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeMhGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B,
                       mblock<dcomp>* mblC, double eps, unsigned rankmax, 
                       contBasis<dcomp>* haar)
{
  mltaGeMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeMhGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B,
                       mblock<scomp>* mblC, double eps, unsigned rankmax, 
                       contBasis<scomp>* haar)
{
  mltaGeMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}


void mltaLrMhGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax, 
                       contBasis<double>* haar)
{
  mltaLrMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaLrMhGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B,
                       mblock<float>* mblC, double eps, unsigned rankmax, 
                       contBasis<float>* haar)
{
  mltaLrMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaLrMhGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B,
                       mblock<dcomp>* mblC, double eps, unsigned rankmax, 
                       contBasis<dcomp>* haar)
{
  mltaLrMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaLrMhGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B,
                       mblock<scomp>* mblC, double eps, unsigned rankmax, 
                       contBasis<scomp>* haar)
{
  mltaLrMhGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}
