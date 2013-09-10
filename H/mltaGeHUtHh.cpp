/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"
#include "blas.h"

// C += d A B^H  (B UtH-matrix, A low-rank mblock<T>, C H-matrix)
// idea: A B^H = U (B V)^H
template<class T> static
void mltalwrUtHh_to_H_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
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
    mltaUtHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nA);

    addGeHLrM(blC, C, eps, rankmax, rankA, dataA, mA, tmp, nA);
    delete [] tmp;
  }
}


// C += d A B^H  (B UtH-matrix, A low-rank mblock<T>, C mblock<T>)
// idea: A B^H = U (B V)^H
template<class T> static
void mltalwrUtHh_to_mbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                          mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn2() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isLrM(A));

  const unsigned rankA = blA->rank(A);

  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[nA*rankA];
    assert(tmp!=NULL);

    blas::setzero(nA*rankA, tmp);
    mltaUtHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nA);

    mblC->addLrM(rankA, dataA, mA, tmp, nA, eps, rankmax);
    delete [] tmp;
  }
}


// C += d A B^H  (B UtH-matrix, A dense mblock<T>, C H-matrix)
// idea: A B^H = (B A^H)^H
template<class T> static
void mltadnsUtHh_to_H_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isGeM(A) && blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[2*mA*nA];
  assert(tmp!=NULL);

  blas::transpose(mA, nA, dataA, tmp);
  blas::setzero(nA*mA, tmp+mA*nA);
  mltaUtHGeM(d, blB, B, mA, tmp, nA, tmp+mA*nA, nA);
  blas::transpose(nA, mA, tmp+mA*nA, tmp);

  addGeHGeM(blC, C, tmp, mA, eps, rankmax);
  delete [] tmp;
}


// C += d A B^H  (B UtH-matrix, A dense mblock<T>, C mblock<T>)
// idea: A B^H = (B A^H)^H
template<class T> static
void mltadnsUtHh_to_mbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                          mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->isleaf() && blA->isGeM(A) && blA->getn2()==blB->getn1());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned mB = blB->getn1();

  T* const tmp = new T[mA*(MAX(nA,mB)+mB)];
  assert(tmp!=NULL);

  blas::transpose(mA, nA, dataA, tmp);
  T* const tmp1 = tmp+MAX(nA,mB)*mA;
  blas::setzero(mB*mA, tmp1);
  mltaUtHGeM(d, blB, B, mA, tmp, nA, tmp1, mB);
  blas::transpose(mB, mA, tmp1, tmp);
  mblC->addGeM(tmp, mA, eps, rankmax);

  delete [] tmp;
}


// C += d A B^H  (A H-matrix, B UtH-matrices, C mblock<T>)
template<class T> static
void mltaGeHUtHh_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax)
{
  if (blA->isleaf())    //  B should be upper triangular
    if (blA->isLrM(A))
      mltalwrUtHh_to_mbl_(d, blA, A, blB, B, mblC, eps, rankmax);
    else
      mltadnsUtHh_to_mbl_(d, blA, A, blB, B, mblC, eps, rankmax);
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

    bl1 = blA->getson(0, 0);
    bl2 = blB->getson(0, 0);
    mltaGeHUtHh_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax);

    for (unsigned k=1; k<ns2; ++k) {
      bl1 = blA->getson(0, k);
      bl2 = blB->getson(0, k);
      mltaGeHGeHh_toMbl(d, bl1, A, bl2, B, R, eps, rankmax);
    }

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);

      bl1 = blA->getson(0, j);
      bl2 = blB->getson(j, j);
      mltaGeHUtHh_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

      for (unsigned k=j+1; k<ns2; ++k) {
        bl1 = blA->getson(0, k);
        bl2 = blB->getson(j, k);
        mltaGeHGeHh_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
      }

      mblock<T>* R1 = R;
      n += nj;
      R = new mblock<T>(m0, n);
      R->unify_cols(eps, rankmax, *R1, R2);

      delete R1;
    }

    for (unsigned i=1; i<ns1; ++i) {

      const unsigned mi = blA->getson(i, 0)->getn1();
      n = n0;
      mblock<T>* Rz = new mblock<T>(mi, n0);

      bl1 = blA->getson(i, 0);
      bl2 = blB->getson(0, 0);
      mltaGeHUtHh_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax);

      for (unsigned k=1; k<ns2; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(0, k);
        mltaGeHGeHh_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax);
      }

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);

        bl1 = blA->getson(i, j);
        bl2 = blB->getson(j, j);
        mltaGeHUtHh_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax);

        for (unsigned k=j+1; k<ns2; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(j, k);
          mltaGeHGeHh_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax);
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


// C += d A B^H  (B UtH-matrix, A,C H-matrices)
template<class T> static
void mltaGeHUtHh_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                 mblock<T>** B, blcluster* blC, mblock<T>** C,
                 double eps, unsigned rankmax)
{
  if (blC->isleaf())
    mltaGeHUtHh_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax);
  else {                            // C is not a leaf
    if (blA->isnleaf()) {           // then B has sons
      unsigned ns1 = blA->getnrs(), ns2 = blA->getncs();
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
          blcluster* bl1 = blA->getson(i, j);
          blcluster* bl2 = blB->getson(j, j);
          mltaGeHUtHh_(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          for (unsigned k=j+1; k<ns2; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(j, k);
            mltaGeHGeHh(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }
        }
    }else {                 // A is a leaf
      if (blA->isLrM(A))    // A is lwr
        mltalwrUtHh_to_H_(d, blA, A, blB, B, blC, C, eps, rankmax);
      else                  // A is dense
        mltadnsUtHh_to_H_(d, blA, A, blB, B, blC, C, eps, rankmax);
    }
  }
}


// ----------------------------------------------------------------------------
// Instanzen

void mltaGeHUtHh(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
                mblock<double>** B, blcluster* blC, mblock<double>** C,
                double eps, unsigned rankmax)
{
  mltaGeHUtHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHUtHh(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
                mblock<float>** B, blcluster* blC, mblock<float>** C,
                double eps, unsigned rankmax)
{
  mltaGeHUtHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHUtHh(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
                mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
                double eps, unsigned rankmax)
{
  mltaGeHUtHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHUtHh(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
                mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
                double eps, unsigned rankmax)
{
  mltaGeHUtHh_(d, blA, A, blB, B, blC, C,eps, rankmax);
}





