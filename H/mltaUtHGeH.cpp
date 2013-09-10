/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B  (A UtH-matrix, B low-rank mblock<T>, C H-matrix)
// idea: A B = (A U) V^T
template<class T> static
void mltaUtHLrM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, blcluster* blC, mblock<T>** C,
                       double eps, unsigned rankmax, contBasis<T>* haar)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isLrM(B) && blC->isnleaf());

  const unsigned rankB = blB->rank(B);

  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[mA*rankB];
    assert(tmp!=NULL);

    blas::setzero(mA*rankB, tmp);
    mltaUtHGeM(d, blA, A, rankB, dataB, mB, tmp, mA);

    addGeHLrM(blC, C, eps, rankmax, rankB, tmp, mA, dataB+rankB*mB, nB, haar);
    delete [] tmp;
  }
}


// C += d A B  (A UtH-matrix, B low-rank mblock<T>, C mblock<T>)
// idea: A B = (A U) V^T
template<class T> static
void mltaUtHLrM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, mblock<T>* mblC,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->getn2()==blB->getn1() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isLrM(B));

  const unsigned rankB = blB->rank(B);

  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[mA*rankB];
    assert(tmp!=NULL);

    blas::setzero(mA*rankB, tmp);
    mltaUtHGeM(d, blA, A, rankB, dataB, mB, tmp, mA);
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar!=NULL) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    
    mblC->addLrM(rankB, tmp, mA, dataB+rankB*mB, nB, eps, rankmax, haarInfo,
		 X, blB->getn2(), X+cols*blB->getn2(), blA->getn1());
    delete haarInfo;
    delete [] tmp;
  }
}


// C += d A B  (A UtH-matrix, B dense mblock<T>, C H-matrix)
template<class T> static
void mltaUtHGeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, blcluster* blC, mblock<T>** C,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isGeM(B) && blC->isnleaf());

  T* const dataB = blB->data(B);
  const unsigned mA = blA->getn1();
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mA*nB];
  assert(tmp!=NULL);

  blas::setzero(mA*nB, tmp);
  mltaUtHGeM(d, blA, A, nB, dataB, mB, tmp, mA);

  addGeHGeM(blC, C, tmp, mA, eps, rankmax, haar);

  delete [] tmp;
}


// C += d A B  (A UtH-matrix, B dense mblock<T>, C mblock<T>)
template<class T> static
void mltaUtHGeM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, mblock<T>* mblC,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->getn2()==blB->getn1() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isGeM(B));

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mA = blA->getn1();

  T* const tmp = new T[mA*nB];
  assert(tmp!=NULL);

  blas::setzero(mA*nB, tmp);
  mltaUtHGeM(d, blA, A, nB, dataB, mB, tmp, mA);
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
  delete [] tmp;
}


// C += d A B  (A UtH-matrix, B H-matrices, C mblock<T>)
template<class T> static
void mltaUtHGeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  if (blB->isleaf()) {
    if (blB->isLrM(B))
      mltaUtHLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
    else
      mltaUtHGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  }
  else {         // A and B are no leaves, hence the clusters can be subdivided

    assert(blA->getncs()==blB->getnrs());
    unsigned ns1 = blA->getnrs(), ns2 = blB->getncs();

    // unify first block row
    const unsigned m0 = blA->getson(0, 0)->getn1();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0, 0)->getn2();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2;
    
    contBasis<T>* haarSon = NULL;
    if(haar) haarSon = haar->son(0,0);

    bl1 = blA->getson(0, 0);
    bl2 = blB->getson(0, 0);
    mltaUtHGeH_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax, haarSon);

    for (unsigned k=1; k<ns1; ++k) {
      bl1 = blA->getson(0, k);
      bl2 = blB->getson(k, 0);
      mltaGeHGeH_toMbl(d, bl1, A, bl2, B, R, eps, rankmax, haarSon);
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
      bl1 = blA->getson(0, 0);
      bl2 = blB->getson(0, j);
      mltaUtHGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax, haarSon);

      for (unsigned k=1; k<ns1; ++k) {
        bl1 = blA->getson(0, k);
        bl2 = blB->getson(k, j);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax, haarSon);
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
      bl1 = blA->getson(i, i);
      bl2 = blB->getson(i, 0);
      mltaUtHGeH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);

      for (unsigned k=i+1; k<ns1; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(k, 0);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);
      }
      delete haarSon;

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);
	if(haar) haarSon = haar->son(i,j);
        bl1 = blA->getson(i, i);
        bl2 = blB->getson(i, j);
        mltaUtHGeH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax,
			  haarSon);

        for (unsigned k=i+1; k<ns1; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(k, j);
          mltaGeHGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax,
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


// C += d A B  (A UtH-matrix, B,C H-matrices)
template<class T> static
void mltaUtHGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		 mblock<T>** B, blcluster* blC, mblock<T>** C,
		 double eps, unsigned rankmax, contBasis<T>* haar)
{
  if (blC->isleaf())
    mltaUtHGeH_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax, haar);
  else {                            // C is not a leaf
    if (blB->isnleaf()) {           // then A has sons
      unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getncs();
      contBasis<T>* haarSon = NULL;
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
          blcluster* bl1 = blA->getson(i, i);
          blcluster* bl2 = blB->getson(i, j);
	  if(haar) haarSon = haar->son(i,j);
          mltaUtHGeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
          for (unsigned k=i+1; k<nsp; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(k, j);
            mltaGeHGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
          }
	  delete haarSon;
        }
    }
    else {  // B is a leaf
      if (blB->isLrM(B))  // B is lwr
        mltaUtHLrM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
      else // B is dense
        mltaUtHGeM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    }
  }
}


// ----------------------------------------------------------------------------
// Instanzen

void mltaUtHGeH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
		mblock<double>** B, blcluster* blC, mblock<double>** C,
		double eps, unsigned rankmax, contBasis<double>* haar)
{
  mltaUtHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaUtHGeH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
		mblock<float>** B, blcluster* blC, mblock<float>** C,
		double eps, unsigned rankmax, contBasis<float>* haar)
{
  mltaUtHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaUtHGeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
		mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
		double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  mltaUtHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaUtHGeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
		mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
		double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  mltaUtHGeH_(d, blA, A, blB, B, blC, C,eps, rankmax, haar);
}



