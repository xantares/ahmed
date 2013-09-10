/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B  (A LtH-matrix, B low-rank mblock<T>, C H-matrix)
// idea: A B = (A U) V^t
template<class T> static
void mltaLtHlwr_to_H_(T d, blcluster* blA, mblock<T>** A,
                       blcluster* blB, mblock<T>** B, blcluster* blC,
                       mblock<T>** C, double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isLrM(B) && blC->isnleaf());

  const unsigned rankB = blB->rank(B);

  if (rankB>0) {
    T* dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[mA*rankB];
    assert(tmp!=NULL);

    blas::setzero(mA*rankB, tmp);
    mltaLtHGeM(d, blA, A, rankB, dataB, mB, tmp, mA);

    addGeHLrM(blC, C, eps, rankmax, rankB, tmp, mA, dataB, nB);
    delete [] tmp;
  }
}


// C += d A B  (A LtH-matrix, B low-rank mblock<T>, C mblock<T>)
// idea: A B = (A U) V^T
template<class T> static
void mltaLtHlwr_to_mbl_(T d, blcluster* blA, mblock<T>** A,
                         blcluster* blB, mblock<T>** B,
                         mblock<T>* mblC, double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blB->isleaf() && blB->isLrM(B) &&
         blA->getn1()==blA->getn2());

  const unsigned rankB = blB->rank(B);

  if (rankB>0) {
    T* dataB = blB->data(B);
    const unsigned mB = blB->getn1();
    const unsigned nB = blB->getn2();
    const unsigned mA = blA->getn1();

    T* const tmp = new T[mA*rankB];
    assert(tmp!=NULL);

    blas::setzero(mA*rankB, tmp);
    mltaLtHGeM(d, blA, A, rankB, dataB, mB, tmp, mA);

    mblC->addLrM(rankB, tmp, mA, dataB, nB, eps, rankmax);

    delete [] tmp;
  }
}


// C += d A B  (A LtH-matrix, B dense mblock<T>, C H-matrix)
template<class T> static
void mltaLtHGeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, blcluster* blC, mblock<T>** C,
                       double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blA->getn1()==blA->getn2() &&
         blB->isleaf() && blB->isGeM(B) && blC->isnleaf());


  T* dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mA = blA->getn1();

  T* const tmp = new T[mA*nB];
  assert(tmp!=NULL);

  blas::setzero(mA*nB, tmp);
  mltaLtHGeM(d, blA, A, nB, dataB, mB, tmp, mA);

  addGeHGeM(blC, C, tmp, mA, eps, rankmax);
  delete [] tmp;
}


// C += d A B  (A LtH-matrix, B dense mblock<T>, C mblock<T>)
template<class T> static
void mltaLtHGeM_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                         mblock<T>** B, mblock<T>* mblC,
                         double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blB->isleaf() && blB->isGeM(B) &&
         blA->getn1()==blA->getn2());

  T* dataB = blB->data(B);
  const unsigned mB = blB->getn1();
  const unsigned nB = blB->getn2();
  const unsigned mA = blA->getn1();

  T* const tmp = new T[mA*nB];
  assert(tmp!=NULL);

  blas::setzero(mA*nB, tmp);
  mltaLtHGeM(d, blA, A, nB, dataB, mB, tmp, mA);

  mblC->addGeM(tmp, mA, eps, rankmax);

  delete [] tmp;
}

// C += d A B  (A LtH-matrix, B H-matrix, C mblock<T>)
template<class T> static
void mltaLtHH_to_mbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                       mblock<T>** B, mblock<T>* mblC,
                       double eps, unsigned rankmax)
{
  if (blB->isleaf()) {
    if (blB->isLrM(B))
      mltaLtHlwr_to_mbl_(d, blA, A, blB, B, mblC, eps, rankmax);
    else
      mltaLtHGeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  }
  else {        // A and B are no leaves, hence the clusters can be subdivided

    assert(blA->getncs()==blB->getnrs());
    unsigned ns1 = blA->getnrs(), ns2 = blB->getncs();

    // unify first block row
    const unsigned m0 = blA->getson(0, 0)->getn1();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0, 0)->getn2();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2;

    bl1 = blA->getson(0, 0);
    bl2 = blB->getson(0, 0);
    mltaLtHH_to_mbl_(d, bl1, A, bl2, B, R, eps, rankmax);

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);

      bl1 = blA->getson(0, 0);
      bl2 = blB->getson(0, j);
      mltaLtHH_to_mbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

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
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(k, 0);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax);
      }

      bl1 = blA->getson(i, i);
      bl2 = blB->getson(i, 0);
      mltaLtHH_to_mbl_(d, bl1, A, bl2, B, Rz, eps, rankmax);

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);

        for (unsigned k=0; k<i; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(k, j);
          mltaGeHGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
        }

        bl1 = blA->getson(i, i);
        bl2 = blB->getson(i, j);
        mltaLtHH_to_mbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

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




// C += d A B  (A LtH-matrix, B,C H-matrices)
template<class T> static
void mltaLtHH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                mblock<T>** B, blcluster* blC, mblock<T>** C,
                double eps, unsigned rankmax)
{
  if (blC->isleaf())
    mltaLtHH_to_mbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax);
  else {                                               // C is not a leaf
    if (blB->isnleaf()) { // then A has sons
      unsigned ns1 = blA->getnrs(), ns2 = blB->getncs();
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* bl = blC->getson(i, j);
          for (unsigned k=0; k<i; ++k) {
            blcluster* bl1 = blA->getson(i, k);
            blcluster* bl2 = blB->getson(k, j);
            mltaGeHGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }
          blcluster* bl1 = blA->getson(i, i);
          blcluster* bl2 = blB->getson(i, j);
          mltaLtHH_(d, bl1, A, bl2, B, bl, C, eps, rankmax);
        }
    }
    else {
      if (blB->isLrM(B))
        mltaLtHlwr_to_H_(d, blA, A, blB, B, blC, C, eps, rankmax);
      else
        mltaLtHGeM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    }
  }
}


// ----------------------------------------------------------------------------
// Instanzen

void mltaLtHH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
               mblock<double>** B, blcluster* blC, mblock<double>** C,
               double eps, unsigned rankmax)
{
  mltaLtHH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaLtHH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
               mblock<float>** B, blcluster* blC, mblock<float>** C,
               double eps, unsigned rankmax)
{
  mltaLtHH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaLtHH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
               mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
               double eps, unsigned rankmax)
{
  mltaLtHH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaLtHH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
               mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
               double eps, unsigned rankmax)
{
  mltaLtHH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}




