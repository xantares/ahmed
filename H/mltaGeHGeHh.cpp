/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B^H  (A,C H-matrices, B low-rank)
// idea: A B^H = (A V) U^H
template<class T> static
void mltaGeHLrMh_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                      mblock<T>** B, blcluster* blC, mblock<T>** C,
                      double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blB->isleaf() && blB->isLrM(B) &&
         blC->isnleaf());

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[mA*rankB];
    assert(tmp!=NULL);

    blas::setzero(mA*rankB, tmp);
    if (mltaGeHGeM(d, blA, A, rankB, dataB+rankB*mB, nB, tmp, mA))
      addGeHLrM(blC, C, eps, rankmax, rankB, tmp, mA, dataB, mB);

    delete [] tmp;
  }
}


// C += d A B^H  (A H-matrix, B low-rank, C mblock)
// idea: A B^H = (A V) U^H
template<class T> static
void mltaGeHLrMh_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn2() && blB->isleaf() && blB->isLrM(B));

  const unsigned rankB = blB->rank(B);
  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[mA*rankB];
    assert(tmp!=NULL);

    blas::setzero(mA*rankB, tmp);
    if (mltaGeHGeM(d, blA, A, rankB, dataB+rankB*mB, nB, tmp, mA)) 
      mblC->addLrM(rankB, tmp, mA, dataB, mB, eps, rankmax);

    delete [] tmp;
  }
}


// C += d A B^H  (A,C H-matrices, B dense)
template<class T> static
void mltaGeHGeMh_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                      mblock<T>** B, blcluster* blC, mblock<T>** C,
                      double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blB->isleaf() && blB->isGeM(B) &&
         blC->isnleaf());

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mA = blA->getn1();

  T* const tmp = new T[mB*(mA+nB)];
  assert(tmp!=NULL);
  blas::transpose(mB, nB, dataB, tmp);

  T* const tmp1 = tmp + nB*mB;
  blas::setzero(mA*mB, tmp1);
  if (mltaGeHGeM(d, blA, A, mB, tmp, nB, tmp1, mA))
    addGeHGeM(blC, C, tmp1, mA, eps, rankmax);

  delete [] tmp;
}


// C += d A B^H  (A,C H-matrices, B dense, C mblock)
template<class T> static
void mltaGeHGeMh_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn2() && blB->isleaf() && blB->isGeM(B));

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mA = blA->getn1();

  T* const tmp = new T[mB*(mA+nB)];
  assert(tmp!=NULL);
  blas::transpose(mB, nB, dataB, tmp);

  T* const tmp1 = tmp + mB*nB;
  blas::setzero(mA*mB, tmp1);
  if (mltaGeHGeM(d, blA, A, mB, tmp, nB, tmp1, mA)) 
    mblC->addGeM(tmp1, mA, eps, rankmax);

  delete [] tmp;
}


// C += d A B^H  (A low-rank, B,C H-matrices)
// idea: A B^H = U (B V)^H
template<class T> static
void mltaLrMGeHh_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                      mblock<T>** B, blcluster* blC, mblock<T>** C,
                      double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blA->isleaf() && blA->isLrM(A) &&
         blC->isnleaf());

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned mB = blB->getn1();

    T* const tmp = new T[mB*rankA];
    assert(tmp!=NULL);

    blas::setzero(mB*rankA, tmp);
    if (mltaGeHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, mB)) 
      addGeHLrM(blC, C, eps, rankmax, rankA, dataA, mA, tmp, mB);

    delete [] tmp;
  }
}


// C += d A B^H  (A low-rank, B H-matrix, C mblock)
// idea: A B^t = U (B V)^t
template<class T> static
void mltaLrMGeHh_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn2() && blA->isleaf() && blA->isLrM(A));

  const unsigned rankA = blA->rank(A);
  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();
    const unsigned mB = blB->getn1();

    T* const tmp = new T[mB*rankA];
    assert(tmp!=NULL);

    blas::setzero(mB*rankA, tmp);
    if (mltaGeHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, mB)) 
      mblC->addLrM(rankA, dataA, mA, tmp, mB, eps, rankmax);

    delete [] tmp;
  }
}


// C += d A B^H  (A dense, B H-matrix)
// idea: A B^H = (B A^H)^H
template<class T> static
inline void mltaGeMGeHh_toGeH_(T d, blcluster* blA, mblock<T>** A,
                             blcluster* blB, mblock<T>** B,
                             blcluster* blC, mblock<T>** C,
                             double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn2() &&
         blB->getn1()==blC->getn2() && blA->isleaf() && blA->isGeM(A) &&
         blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned mB = blB->getn1();

  T* const tmp = new T[mA*(MAX(nA,mB)+mB)];
  assert(tmp!=NULL);

  blas::transpose(mA, nA, dataA, tmp);
  T* const tmp1 = tmp + MAX(nA,mB)*mA;
  blas::setzero(mB*mA, tmp1);
  if (mltaGeHGeM(d, blB, B, mA, tmp, nA, tmp1, mB)) {
    blas::transpose(mB, mA, tmp1, tmp);
    addGeHGeM(blC, C, tmp, mA, eps, rankmax);
  }

  delete [] tmp;
}


// C += d A B^H  (A dense, B H-matrix, C mblock)
// idea: A B^H = (B A^H)^H
template<class T> static
void mltaGeMGeHh_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn2() && blA->isleaf() && blA->isGeM(A));

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();
  const unsigned mB = blB->getn1();

  T* const tmp = new T[mA*(MAX(nA,mB)+mB)];
  assert(tmp!=NULL);

  blas::transpose(mA, nA, dataA, tmp);
  T* const tmp1 = tmp + MAX(nA,mB)*mA;
  blas::setzero(mB*mA, tmp1);
  if (mltaGeHGeM(d, blB, B, mA, tmp, nA, tmp1, mB)) {
    blas::transpose(mB, mA, tmp1, tmp);
    mblC->addGeM(tmp, mA, eps, rankmax);
  }

  delete [] tmp;
}


// C += d A B^H  (A,B H-matrices, C mblock)
template<class T> static
void mltaGeHGeHh_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                      mblock<T>** B, mblock<T>* mblC,
                      double eps, unsigned rankmax)
{
  if (blB->isleaf() && blB->isLrM(B))
    mltaGeHLrMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blA->isleaf() && blA->isLrM(A))
    mltaLrMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blB->isleaf() && blB->isGeM(B))
    mltaGeHGeMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blA->isleaf() && blA->isGeM(A))
    mltaGeMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else {         // A and B are no leaves, hence the clusters can be subdivided
    assert(blA->getncs()==blB->getncs());
    unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getnrs();

    // unify first block row
    const unsigned m0 = blA->getson(0, 0)->getn1();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0, 0)->getn1();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2;

    for (unsigned k=0; k<nsp; ++k) {
      bl1 = blA->getson(0, k);
      bl2 = blB->getson(0, k);
      mltaGeHGeHh_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax);
    }

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(j, 0)->getn1();
      mblock<T> R2(m0, nj);

      for (unsigned k=0; k<nsp; ++k) {
        bl1 = blA->getson(0, k);
        bl2 = blB->getson(j, k);
        mltaGeHGeHh_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);
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

      for (unsigned k=0; k<nsp; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(0, k);
        mltaGeHGeHh_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax);
      }

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(j, 0)->getn1();
        mblock<T> R2(mi, nj);

        for (unsigned k=0; k<nsp; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(j, k);
          mltaGeHGeHh_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);
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


// C += d A B^H  (A,B,C H-matrices)
template<class T> static
void mltaGeHGeHh_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
               mblock<T>** B, blcluster* blC, mblock<T>** C,
               double eps, unsigned rankmax)
{
  if (blC->isleaf())                // C is a leaf
    mltaGeHGeHh_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax);
  else {                            // C is not a leaf
    if (blA->isnleaf() && blB->isnleaf()) {
      unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getnrs();
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
          for (unsigned k=0; k<nsp; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(j, k);
            mltaGeHGeHh_(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }
        }
    }
    else if (blB->isleaf() && blB->isLrM(B))
      mltaGeHLrMh_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMGeHh_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHGeMh_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else {                                       // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMGeHh_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHGeHh(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
              mblock<double>** B, blcluster* blC, mblock<double>** C,
              double eps, unsigned rankmax)
{
  mltaGeHGeHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHGeHh(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
              mblock<float>** B, blcluster* blC, mblock<float>** C,
              double eps, unsigned rankmax)
{
  mltaGeHGeHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHGeHh(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
              mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
              double eps, unsigned rankmax)
{
  mltaGeHGeHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHGeHh(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
              mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
              double eps, unsigned rankmax)
{
  mltaGeHGeHh_(d, blA, A, blB, B, blC, C, eps, rankmax);
}


void mltaGeHGeHh_toMbl(double d, blcluster* blA, mblock<double>** A,
                     blcluster* blB, mblock<double>** B, mblock<double>* mblC,
                     double eps, unsigned rankmax)
{
  mltaGeHGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHGeHh_toMbl(float d, blcluster* blA, mblock<float>** A,
                     blcluster* blB, mblock<float>** B, mblock<float>* mblC,
                     double eps, unsigned rankmax)
{
  mltaGeHGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHGeHh_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
                     mblock<dcomp>** B, mblock<dcomp>* mblC,
                     double eps, unsigned rankmax)
{
  mltaGeHGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHGeHh_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                     blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
                     double eps, unsigned rankmax)
{
  mltaGeHGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}


void mltaGeHGeMh_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax)
{
  mltaGeHGeMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHGeMh_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B, mblock<float>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeHGeMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHGeMh_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeHGeMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHGeMh_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeHGeMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}


void mltaGeMGeHh_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax)
{
  mltaGeMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeMGeHh_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B, mblock<float>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeMGeHh_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeMGeHh_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}


void mltaLrMGeHh_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax)
{
  mltaLrMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaLrMGeHh_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B, mblock<float>* mblC,
                       double eps, unsigned rankmax)
{
  mltaLrMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaLrMGeHh_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaLrMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaLrMGeHh_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaLrMGeHh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}


void mltaGeHLrMh_toMbl(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B,
                       mblock<double>* mblC, double eps, unsigned rankmax)
{
  mltaGeHLrMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHLrMh_toMbl(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B, mblock<float>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeHLrMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHLrMh_toMbl(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B, mblock<dcomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeHLrMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}

void mltaGeHLrMh_toMbl(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B, mblock<scomp>* mblC,
                       double eps, unsigned rankmax)
{
  mltaGeHLrMh_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
}
