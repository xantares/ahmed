/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"


// C += d A B  (A symm. H-matrix, C H-matrix, B low-rank)
// idea: A B = (A U) V^H
template<class T> static
void mltaHeHLrM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isLrM(B) && blC->isnleaf());

  const unsigned rankB = blB->rank(B);

  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();

    T* const tmp = new T[rankB*mB];
    assert(tmp!=NULL);
    blas::setzero(rankB*mB, tmp);
    mltaHeHGeM(d, blA, A, rankB, dataB, mB, tmp, mB);

    addGeHLrM(blC, C, eps, rankmax, rankB, tmp, mB, dataB+rankB*mB, nB);
    delete [] tmp;
  }
}


// C += d A B  (A symm. H-matrix, B low-rank, C mblock)
// idea: A B = (A U) V^H
template<class T> static
void mltaHeHLrM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                          mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isLrM(B));

  const unsigned rankB = blB->rank(B);

  if (rankB>0) {
    T* const dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();

    T* const tmp = new T[rankB*mB];
    assert(tmp!=NULL);

    blas::setzero(rankB*mB, tmp);
    mltaHeHGeM(d, blA, A, rankB, dataB, mB, tmp, mB);
    mblC->addLrM(rankB, tmp, mB, dataB+rankB*mB, nB, eps, rankmax);

    delete [] tmp;
  }
}


// C += d A B  (A symm. H-matrix, C H-matrix, B dense)
template<class T> static
void mltaHeHGeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isGeM(B) && blC->isnleaf());

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mB*nB];
  assert(tmp!=NULL);

  blas::setzero(mB*nB, tmp);
  mltaHeHGeM(d, blA, A, nB, dataB, mB, tmp, mB);
  addGeHGeM(blC, C, tmp, mB, eps, rankmax);

  delete [] tmp;
}


// C += d A B  (A symm. H-matrix, B dense, C mblock)
template<class T> static
void mltaHeHGeM_toMbl_(T d, blcluster* blA, mblock<T>** A,
                          blcluster* blB, mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isGeM(B));

  T* const dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mB*nB];
  assert(tmp!=NULL);

  blas::setzero(mB*nB, tmp);
  mltaHeHGeM(d, blA, A, nB, dataB, mB, tmp, mB);
  mblC->addGeM(tmp, mB, eps, rankmax);

  delete [] tmp;
}


// C += d A B  (B,C H-matrices, A symm. dense)
// idea: A B = (B^H A)^H
template<class T> static
void mltaHeMGeH_toGeH_(T d, blcluster* blA, mblock<T>** A,
                        blcluster* blB, mblock<T>** B,
                        blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blA->isleaf() && blA->isGeM(A) && blC->isnleaf());

  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mB*(MAX(mB,nB)+nB)];
  assert(tmp!=NULL);

  A[blA->getidx()]->convHeM_toGeM(tmp, mB);
  T* const tmp1 = tmp + mB*MAX(mB,nB);
  blas::setzero(nB*mB, tmp1);
  if (mltaGeHhGeM(d, blB, B, mB, tmp, mB, tmp1, nB)) {
    blas::transpose(nB, mB, tmp1, tmp);
    addGeHGeM(blC, C, tmp, mB, eps, rankmax);
  }
  delete [] tmp;
}


// C += d A B  (B H-matrix, A symm. dense, C mblock)
// idea: A B = (B^t A)^t
template<class T> static
void mltaHeMGeH_toMbl_(T d, blcluster* blA, mblock<T>** A,
                          blcluster* blB, mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->isGeM(A) && blA->getn2()==blB->getn1() &&
         blA->getn1()==blA->getn2());

  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();

  T* const tmp = new T[mB*(MAX(mB,nB)+nB)];
  assert(tmp!=NULL);

  A[blA->getidx()]->convHeM_toGeM(tmp, mB);
  T* const tmp1 = tmp + mB*MAX(mB,nB);
  blas::setzero(nB*mB, tmp1);
  if (mltaGeHhGeM(d, blB, B, mB, tmp, mB, tmp1, nB)) {
    blas::transpose(nB, mB, tmp1, tmp);
    mblC->addGeM(tmp, mB, eps, rankmax);
  }

  delete [] tmp;
}


// C += d A B  (A symm. H-matrix, B H-matrices, C mblock)
template<class T> static
inline void mltaHeHGeH_toMbl_(T d, blcluster* blA, mblock<T>** A,
                               blcluster* blB, mblock<T>** B,
                               mblock<T>* mblC, double eps, unsigned rankmax)
{
  assert(blA->getn1()==blA->getn2());

  if (blB->isleaf() && blB->isLrM(B))
    mltaHeHLrM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blB->isleaf() && blB->isGeM(B))
    mltaHeHGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blA->isleaf() && blA->isGeM(A))
    mltaHeMGeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
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

    bl1 = blA->getson(0, 0);
    bl2 = blB->getson(0, 0);
    mltaHeHGeH_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax);

    for (unsigned k=1; k<nsp; ++k) {
      bl1 = blA->getson(0, k);
      bl2 = blB->getson(k, 0);
      mltaGeHGeH_toMbl(d, bl1, A, bl2, B, R, eps, rankmax);
    }

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);

      bl1 = blA->getson(0, 0);
      bl2 = blB->getson(0, j);
      mltaHeHGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

      for (unsigned k=1; k<nsp; ++k) {
        bl1 = blA->getson(0, k);
        bl2 = blB->getson(k, j);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
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

      for (unsigned k=0; k<i; ++k) {
        bl1 = blA->getson(k, i);
        bl2 = blB->getson(k, 0);
        mltaGeHhGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax);
      }

      bl1 = blA->getson(i, i);
      bl2 = blB->getson(i, 0);
      mltaHeHGeH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax);

      for (unsigned k=i+1; k<nsp; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(k, 0);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax);
      }

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);

        for (unsigned k=0; k<i; ++k) {
          bl1 = blA->getson(k, i);
          bl2 = blB->getson(k, j);
          mltaGeHhGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
        }

        bl1 = blA->getson(i, i);
        bl2 = blB->getson(i, j);
        mltaHeHGeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

        for (unsigned k=i+1; k<nsp; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(k, j);
          mltaGeHGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
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


// C += d A B  (A symm. H-matrix, B,C H-matrices)
template<class T> static
void mltaHeHGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                 mblock<T>** B, blcluster* blC, mblock<T>** C,
                 double eps, unsigned rankmax)
{
  assert(blA->getn1()==blA->getn2());

  if (blC->isleaf())                // C is a leaf
    mltaHeHGeH_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax);
  else {                            // C is not a leaf
    if (!blA->isleaf() && !blB->isleaf()) {
      unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getncs();
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster *bl1, *bl2, *bl = blC->getson(i, j);

          for (unsigned k=0; k<i; ++k) {
            bl1 = blA->getson(k, i);
            bl2 = blB->getson(k, j);
            mltaGeHhGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }

          bl1 = blA->getson(i, i);
          bl2 = blB->getson(i, j);
          mltaHeHGeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax);

          for (unsigned k=i+1; k<nsp; ++k) {
            bl1 = blA->getson(i, k);
            bl2 = blB->getson(k, j);
            mltaGeHGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }
        }
    }
    else if (blB->isleaf() && blB->isLrM(B))
      mltaHeHLrM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaHeHGeM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else {                                       // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaHeMGeH_toGeH(d, blA, A, blB, B, blC, C, eps, rankmax);
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaHeHGeH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
                mblock<double>** B, blcluster* blC, mblock<double>** C,
                double eps, unsigned rankmax)
{
  mltaHeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaHeHGeH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
                mblock<float>** B, blcluster* blC, mblock<float>** C,
                double eps, unsigned rankmax)
{
  mltaHeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaHeHGeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
                mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
                double eps, unsigned rankmax)
{
  mltaHeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaHeHGeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
                mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
                double eps, unsigned rankmax)
{
  mltaHeHGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}




void mltaHeMGeH_toGeH(double d, blcluster* blA, mblock<double>** A,
                       blcluster* blB, mblock<double>** B, blcluster* blC,
                       mblock<double>** C, double eps, unsigned rankmax)
{
  mltaHeMGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaHeMGeH_toGeH(float d, blcluster* blA, mblock<float>** A,
                       blcluster* blB, mblock<float>** B, blcluster* blC,
                       mblock<float>** C, double eps, unsigned rankmax)
{
  mltaHeMGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaHeMGeH_toGeH(dcomp d, blcluster* blA, mblock<dcomp>** A,
                       blcluster* blB, mblock<dcomp>** B, blcluster* blC,
                       mblock<dcomp>** C, double eps, unsigned rankmax)
{
  mltaHeMGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaHeMGeH_toGeH(scomp d, blcluster* blA, mblock<scomp>** A,
                       blcluster* blB, mblock<scomp>** B, blcluster* blC,
                       mblock<scomp>** C, double eps, unsigned rankmax)
{
  mltaHeMGeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}
