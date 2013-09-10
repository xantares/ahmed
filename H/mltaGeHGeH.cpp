/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B  (A,C H-matrices, B low-rank mblock)
// idea: A B = (A U) V^H
template<class T> static
void mltaGeHLrM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax, contBasis<T>* haar)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->isleaf() && blB->isLrM(B) &&
         blC->isnleaf());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[rankB*mA];
    assert(tmp!=NULL);
    blas::setzero(rankB*mA, tmp);
    if (mltaGeHGeM(d, blA, A, rankB, dataB, mB, tmp, mA)) 
      addGeHLrM(blC, C, eps, rankmax, rankB, tmp, mA, dataB+rankB*mB, nB, haar);
    delete [] tmp;
  }
}


// C += d A B  (A H-matrix, B low-rank mblock, C mblock)
// idea: A B = (A U) V^H
template<class T> static
void mltaGeHLrM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blB->isLrM(B) && blA->getn2()==blB->getn1());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* dataB = blB->data(B);
    const unsigned nB = blB->getn2();
    const unsigned mB = blB->getn1();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[rankB*mA];
    assert(tmp!=NULL);
    blas::setzero(rankB*mA, tmp);
    if (mltaGeHGeM(d, blA, A, rankB, dataB, mB, tmp, mA)){
      contLowLevel<T>* haarInfo = NULL;
      unsigned cols = 0;
      T* X = NULL;
      if (haar!=NULL) {
	haarInfo = new contLowLevel<T>(haar);	
	cols = haar->getcols();
	X = haar->getX();
	haar->initX();
      }

      mblC->addLrM(rankB, tmp, mA, dataB+mB*rankB, nB, eps, rankmax,
		   haarInfo, X, blB->getn2(), 
		   X+cols*blB->getn2(), blA->getn1());
      delete haarInfo;
    }
    delete [] tmp;
  }
}


// C += d A B  (A low-rank mblock, B,C H-matrices)
// idea: A B = U (B^H V)^H
template<class T> static
void mltaLrMGeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->isleaf() && blA->isLrM(A) &&
         blC->isnleaf());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned nB = blB->getn2();

    T* const tmp = new T[rankA*nB];
    assert(tmp!=NULL);
    blas::setzero(rankA*nB, tmp);
    if (mltaGeHhGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nB))
      addGeHLrM(blC, C, eps, rankmax, rankA, dataA, mA, tmp, nB, haar);
    delete [] tmp;
  }
}


// C += d A B  (A low-rank mblock, B H-matrix, C mblock)
// idea: A B = U (B^H V)^H
template<class T> static
void mltaLrMGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->isLrM(A) && blA->getn2()==blB->getn1());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned nB = blB->getn2();

    T* const dataA = blA->data(A);
    T* const tmp = new T[rankA*nB];
    assert(tmp!=NULL);
    blas::setzero(rankA*nB, tmp);
    if (mltaGeHhGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nB)){
      contLowLevel<T>* haarInfo = NULL;
      unsigned cols = 0;
      T* X = NULL;
      if (haar!=NULL) {
	haarInfo = new contLowLevel<T>(haar);
	cols = haar->getcols();
	X = haar->getX();
	haar->initX();
      }
      mblC->addLrM(rankA, dataA, mA, tmp, nB, eps, rankmax, haarInfo,
		   X, blB->getn2(), X+cols*blB->getn2(), blA->getn1());
      delete haarInfo;
    }
    delete [] tmp;
  }
}


// C += d A B  (A,C H-matrices, B dense mblock)
template<class T> static
void mltaGeHGeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->isleaf() && blB->isGeM(B) &&
         blC->isnleaf());

  T* dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mA = blA->getn1();

  T* const tmp = new T[mA*nB];
  assert(tmp!=NULL);
  blas::setzero(mA*nB, tmp);
  if (mltaGeHGeM(d, blA, A, nB, dataB, mB, tmp, mA))
    addGeHGeM(blC, C, tmp, mA, eps, rankmax, haar);
  delete [] tmp;
}


// C += d A B  (A H-matrix, B dense mblock, C mblock)
template<class T> static
void mltaGeHGeM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blB->isGeM(B) && blA->getn2()==blB->getn1());

  const unsigned mA = blA->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mB = blB->getn1();
  T* dataB = blB->data(B);

  T* const tmp = new T[mA*nB];
  assert(tmp!=NULL);
  blas::setzero(mA*nB, tmp);
  if (mltaGeHGeM(d, blA, A, nB, dataB, mB, tmp, mA)){    
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar!=NULL) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }    
    mblC->addGeM(tmp, mA, eps, rankmax, haarInfo, X, blB->getn2(),
		 X+cols*blB->getn2(), blA->getn1());
    delete haarInfo;
  }
  delete [] tmp;
}


// C += d A B  (A dense mblock, B,C H-matrices)
// idea: A B = (B^H A^H)^H
template<class T> static
void mltaGeMGeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, blcluster* blC, mblock<T>** C,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->isleaf() && blA->isGeM(A) &&
         blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mA*(MAX(nA,nB)+nB)];
  assert(tmp!=NULL);
  blas::transpose(mA, nA, dataA, tmp);
  T* const tmp1 = tmp + mA*MAX(nA,nB);
  blas::setzero(nB*mA, tmp1);
  if (mltaGeHhGeM(d, blB, B, mA, tmp, nA, tmp1, nB)) {
    blas::transpose(nB, mA, tmp1, tmp);
    addGeHGeM(blC, C, tmp, mA, eps, rankmax, haar);
  }
  delete [] tmp;
}


// C += d A B  (A dense mblock, B H-matrix, C mblock)
// idea: A B = (B^H A^H)^H
template<class T> static
void mltaGeMGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->isGeM(A) && blA->getn2()==blB->getn1());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mA*(MAX(nA,nB)+nB)];
  assert(tmp!=NULL);
  blas::transpose(mA, nA, dataA, tmp);
  T* const tmp1 = tmp + mA*MAX(nA,nB);
  blas::setzero(nB*mA, tmp1);
  if (mltaGeHhGeM(d, blB, B, mA, tmp, nA, tmp1, nB)) {
    blas::transpose(nB, mA, tmp1, tmp);
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar!=NULL) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    mblC->addGeM(tmp, mA, eps, rankmax, haarInfo, X, blB->getn2(),
		 X+cols*blB->getn2(), blA->getn1());
    delete haarInfo;
  }
  delete [] tmp;
}


// C += d A B  (A,B H-matrices, C mblock)
template<class T> static
void mltaGeHGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, mblock<T>* mblC,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  if (blB->isleaf() && blB->isLrM(B))
    mltaGeHLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  else if (blA->isleaf() && blA->isLrM(A))
    mltaLrMGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  else if (blB->isleaf() && blB->isGeM(B))
    mltaGeHGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  else if (blA->isleaf() && blA->isGeM(A))
    mltaGeMGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  else {        // A and B are no leaves, hence the clusters can be subdivided
    assert(blA->getncs()==blB->getnrs());
    unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getncs();

    // unify first block row
    const unsigned m0 = blA->getson(0, 0)->getn1();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0, 0)->getn2();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2;

    contBasis<T>* haarSon = NULL;
    if(haar) haarSon = haar->son(0,0);

    for (unsigned k=0; k<nsp; ++k) {
      bl1 = blA->getson(0, k);
      bl2 = blB->getson(k, 0);
      mltaGeHGeH_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax, haarSon);
    }
    delete haarSon;

    unsigned cols =0;
    contLowLevel<T>* haarInfo = NULL;
    if(haar!=NULL){
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
    }

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);
      if(haar) haarSon = haar->son(0,j);
      for (unsigned k=0; k<nsp; ++k) {
        bl1 = blA->getson(0, k);
        bl2 = blB->getson(k, j);
        mltaGeHGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax, haarSon);
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
      R->unify_cols(eps, rankmax, *R1, R2, haarInfo, X, blB->getn2(),
		    X+cols*blB->getn2(), blA->getn2());
      delete R1;
    }

    for (unsigned i=1; i<ns1; ++i) {

      const unsigned mi = blA->getson(i, 0)->getn1();
      n = n0;
      mblock<T>* Rz = new mblock<T>(mi, n0);
      if(haar) haarSon = haar->son(i,0);
      for (unsigned k=0; k<nsp; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(k, 0);
        mltaGeHGeH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);
      }
      delete haarSon;

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);
	if(haar) haarSon = haar->son(i,j);
        for (unsigned k=0; k<nsp; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(k, j);
          mltaGeHGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax,
			    haarSon);
        }
	delete haarSon;

	T* X = NULL;
	if(haar!=NULL){
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
      if(haar!=NULL){
	X = haar->getX();
	haar->initX();
      }
      R->unify_rows(eps, rankmax, *R1, *Rz, haarInfo, X, blB->getn2(),
		    X+cols*blB->getn2(), blA->getn2());
      delete R1;
      delete Rz;
    }
    T* X = NULL;
    if(haar!=NULL){
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
void mltaGeHGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		 mblock<T>** B, blcluster* blC, mblock<T>** C,
		 double eps, unsigned rankmax, contBasis<T>* haar)
{
  if (blC->isleaf()) {              // C is a leaf
    mltaGeHGeH_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax,
		      haar);
  }
  else {                           // C is not a leaf
    if (blA->isnleaf() && blB->isnleaf()) {
      unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getncs();
      contBasis<T>* haarSon = NULL;
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
	  if(haar) haarSon = haar->son(i,j);
          for (unsigned k=0; k<nsp; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(k, j);
            mltaGeHGeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
          }
	  delete haarSon;
        }
    }
    else if (blB->isleaf() && blB->isLrM(B))
      mltaGeHLrM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHGeM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    else {                                        // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHGeH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
		mblock<double>** B, blcluster* blC, mblock<double>** C,
		double eps, unsigned rankmax, contBasis<double>* haar)
{
  mltaGeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHGeH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
		mblock<float>** B, blcluster* blC, mblock<float>** C,
		double eps, unsigned rankmax, contBasis<float>* haar)
{
  mltaGeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHGeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
		mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
		double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  mltaGeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHGeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
		mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
		double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  mltaGeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}


void mltaGeHGeH_toMbl(double d, blcluster* blA, mblock<double>** A,
		      blcluster* blB, mblock<double>** B, mblock<double>* mblC,
		      double eps, unsigned rankmax, contBasis<double>* haar)
{
  mltaGeHGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHGeH_toMbl(float d, blcluster* blA, mblock<float>** A,
		      blcluster* blB, mblock<float>** B, mblock<float>* mblC,
		      double eps, unsigned rankmax, contBasis<float>* haar)
{
  mltaGeHGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHGeH_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
		      blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
		      double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  mltaGeHGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}

void mltaGeHGeH_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
		      blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
		      double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  mltaGeHGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
}
