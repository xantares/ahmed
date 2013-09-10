/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B  (B LtH-matrix, A low-rank mblock<T>, C H-matrix)
// idea: A B = U (B^H V)^H
template<class T> static
void mltaLrMLtH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, blcluster* blC, mblock<T>** C,
                       double eps, unsigned rankmax, contBasis<T>* haar)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isLrM(A) && blC->isnleaf());

  const unsigned rankA = blA->rank(A);

  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[nA*rankA];
    assert(tmp!=NULL);

    blas::setzero(nA*rankA, tmp);
    mltaLtHhGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nA);

    addGeHLrM(blC, C, eps, rankmax, rankA, dataA, mA, tmp, nA, haar);
    delete [] tmp;
  }
}


// C += d A B  (A low-rank mblock<T>, B LtH-matrix, C mblock<T>)
// idea: A B = U (B^H V)^H
template<class T> static
void mltaLrMLtH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, mblock<T>* mblC,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->getn2()==blB->getn1() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isLrM(A));

  const unsigned rankA = blA->rank(A);

  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[nA*rankA];
    assert(tmp!=NULL);

    blas::setzero(nA*rankA, tmp);
    mltaLtHhGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nA);
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar!=NULL) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }

    mblC->addLrM(rankA, dataA, mA, tmp, nA, eps, rankmax,
		 haarInfo, X, blB->getn2(), 
		 X+cols*blB->getn2(), blA->getn1());
    delete haarInfo;
    delete [] tmp;
  }
}


// C += d A B  (B LtH-matrix, A dense mblock<T>, C H-matrix)
// idea: A B = (B^H A^H)^H
template<class T> static
void mltaGeMLtH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, blcluster* blC, mblock<T>** C,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn2()==blC->getn2() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isGeM(A) && blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[2*mA*nA];
  assert(tmp!=NULL);

  blas::transpose(mA, nA, dataA, tmp);
  blas::setzero(nA*mA, tmp+mA*nA);
  mltaLtHhGeM(d, blB, B, mA, tmp, nA, tmp+mA*nA, nA);
  blas::transpose(nA, mA, tmp+mA*nA, tmp);

  addGeHGeM(blC, C, tmp, mA, eps, rankmax, haar);
  delete [] tmp;
}


// C += d A B  (A dense mblock<T>, B LtH-matrix, C mblock<T>)
// idea: A B = (B^H A^H)^H
template<class T> static
void mltaGeMLtH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		       mblock<T>** B, mblock<T>* mblC,
		       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(blA->isleaf() && blA->isGeM(A) && blA->getn2()==blB->getn1());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned nB = blB->getn1();
  
  T* const tmp = new T[mA*(MAX(nA,nB)+nB)];
  assert(tmp!=NULL);
  
  blas::transpose(mA, nA, dataA, tmp);
  T* const tmp1 = tmp+MAX(nA,nB)*mA;
  blas::setzero(nB*mA, tmp1);
  mltaLtHhGeM(d, blB, B, mA, tmp, nA, tmp1, nB);
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
  delete [] tmp;
}


// C += d A B  (A H-matrix, B LtH-matrix, C mblock<T>)
template<class T> static
void mltaGeHLtH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  if (blA->isleaf()) {
    if (blA->isLrM(A))
      mltaLrMLtH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
    else
      mltaGeMLtH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax, haar);
  }
  else {         // A and B are no leaves, hence the clusters can be subdivided

    assert(blA->getncs()==blB->getnrs());
    unsigned ns1 = blA->getnrs(), ns2 = blA->getncs();

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
    mltaGeHLtH_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax, haarSon);

    for (unsigned k=1; k<ns2; ++k) {
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

      const unsigned nj = blB->getson(j, 0)->getn2();
      mblock<T> R2(m0, nj);
      if(haar) haarSon = haar->son(0,j);
      bl1 = blA->getson(0, j);
      bl2 = blB->getson(j, j);
      mltaGeHLtH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax, haarSon);

      for (unsigned k=j+1; k<ns2; ++k) {
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
      bl1 = blA->getson(i, 0);
      bl2 = blB->getson(0, 0);
      mltaGeHLtH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);

      for (unsigned k=1; k<ns2; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(k, 0);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);
      }
      delete haarSon;

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);
	if(haar) haarSon = haar->son(i,j);
        bl1 = blA->getson(i, j);
        bl2 = blB->getson(j, j);
        mltaGeHLtH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);

        for (unsigned k=i+1; k<ns2; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(k, j);
          mltaGeHGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax, haarSon);
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


// C += d A B  (A,C H-matrices, B LtH-matrix)
template<class T> static
void mltaGeHLtH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
		 mblock<T>** B, blcluster* blC, mblock<T>** C,
		 double eps, unsigned rankmax, contBasis<T>* haar)
{
  if (blC->isleaf())
    mltaGeHLtH_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax,
		      haar);
  else {                            // C is not a leaf
    if (blA->isnleaf()) {           // then B has sons
      unsigned ns1 = blA->getnrs(), ns2 = blA->getncs();
      contBasis<T>* haarSon = NULL;
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
          blcluster* bl1 = blA->getson(i, j);
          blcluster* bl2 = blB->getson(j, j);
	  if(haar) haarSon = haar->son(i,j);
          mltaGeHLtH_(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
          for (unsigned k=j+1; k<ns2; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(k, j);
            mltaGeHGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax, haarSon);
          }
	  delete haarSon;
        }
    }
    else {  // A is a leaf
      if (blA->isLrM(A))  // A is LrM
        mltaLrMLtH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
      else // A is dense
        mltaGeMLtH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
    }
  }
}


// ----------------------------------------------------------------------------
// Instanzen

void mltaGeHLtH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
		mblock<double>** B, blcluster* blC, mblock<double>** C,
		double eps, unsigned rankmax, contBasis<double>* haar)
{
  mltaGeHLtH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHLtH(float  d, blcluster* blA, mblock<float>** A, blcluster* blB,
		mblock<float>** B, blcluster* blC, mblock<float>** C,
		double eps, unsigned rankmax, contBasis<float>* haar)
{
  mltaGeHLtH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHLtH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
		mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
		double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  mltaGeHLtH_(d, blA, A, blB, B, blC, C, eps, rankmax, haar);
}

void mltaGeHLtH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
		mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
		double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  mltaGeHLtH_(d, blA, A, blB, B, blC, C,eps, rankmax, haar);
}



