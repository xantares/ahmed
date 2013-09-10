/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d A B  (A low-rank, C H-matrices, B symm. H-matrix)
// idea: A B = U (B V)^H
template<class T> static
void mltaLrMHeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isLrM(A) && blC->isnleaf());

  const unsigned rankA = blA->rank(A);

  if (rankA>0) {
    T* const dataA = blA->data(A);
    const unsigned mA = blA->getn1();
    const unsigned nA = blA->getn2();

    T* const tmp = new T[nA*rankA];
    assert(tmp!=NULL);
    blas::setzero(nA*rankA, tmp);
    mltaHeHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nA);
    addGeHLrM(blC, C, eps, rankmax, rankA, dataA, mA, tmp, nA);
    delete [] tmp;
  }
}


// C += d A B  (B symm. H-matrix, A low-rank, C mblock)
// idea: A B = U (B V)^H
template<class T> static
void mltaLrMHeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                          mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
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
    mltaHeHGeM(d, blB, B, rankA, dataA+rankA*mA, nA, tmp, nA);
    mblC->addLrM(rankA, dataA, mA, tmp, nA, eps, rankmax);

    delete [] tmp;
  }
}


// C += d A B  (A,C H-matrices, B dense Herm.)
template<class T> static
void mltaGeHHeM_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->getn1()==blB->getn2() &&
         blB->isleaf() && blB->isGeM(B) && blC->isnleaf());

  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[(mA+nA)*nA];
  T* const tmp1 = tmp + nA*nA;
  assert(tmp!=NULL);

  B[blB->getidx()]->convHeM_toGeM(tmp, nA);
  blas::setzero(mA*nA, tmp1);
  if (mltaGeHGeM(d, blA, A, nA, tmp, nA, tmp1, mA))
    addGeHGeM(blC, C, tmp1, mA, eps, rankmax);
  delete [] tmp;
}


// C += d A B  (A H-matrix, B dense symm., C mblock)
template<class T> static
void mltaGeHHeM_toMbl_(T d, blcluster* blA, mblock<T>** A,
                          blcluster* blB, mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blB->getn1()==blB->getn2() &&
         blB->isleaf() && blB->isGeM(B));

  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[(mA+nA)*nA];
  T* const tmp1 = tmp + nA*nA;
  assert(tmp!=NULL);

  B[blB->getidx()]->convHeM_toGeM(tmp, nA);
  blas::setzero(mA*nA, tmp1);
  if (mltaGeHGeM(d, blA, A, nA, tmp, nA, tmp1, mA))
    mblC->addGeM(tmp1, mA, eps, rankmax);
  delete [] tmp;
}


// C += d A B  (B symm. H-matrix, C H-matrix, A dense)
// idea: A B = (B A^H)^H
template<class T> static
void mltaGeMHeH_toGeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, blcluster* blC, mblock<T>** C,
                        double eps, unsigned rankmax)
{
  assert(blC->getn1()==blA->getn1() && blA->getn2()==blB->getn1() &&
         blB->getn2()==blC->getn2() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isGeM(A) && blC->isnleaf());

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[2*mA*nA];
  assert(tmp!=NULL);
  blas::transpose(mA, nA, dataA, tmp);

  T* const tmp1 = tmp + nA*mA;
  blas::setzero(nA*mA, tmp1);
  mltaHeHGeM(d, blB, B, mA, tmp, nA, tmp1, nA);
  blas::transpose(nA, mA, tmp1, tmp);
  addGeHGeM(blC, C, tmp, mA, eps, rankmax);
  delete [] tmp;
}


// C += d A B  (B symm. H-matrix, A dense, C block)
// idea: A B = (B A^H)^H
template<class T> static
void mltaGeMHeH_toMbl_(T d, blcluster* blA, mblock<T>** A,
                          blcluster* blB, mblock<T>** B, mblock<T>* mblC,
                          double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blB->getn1()==blB->getn2() &&
         blA->isleaf() && blA->isGeM(A));

  T* const dataA = blA->data(A);
  const unsigned mA = blA->getn1();
  const unsigned nA = blA->getn2();

  T* const tmp = new T[2*mA*nA];
  assert(tmp!=NULL);
  blas::transpose(mA, nA, dataA, tmp);

  T* const tmp1 = tmp + nA*mA;
  blas::setzero(nA*mA, tmp1);
  mltaHeHGeM(d, blB, B, mA, tmp, nA, tmp1, nA);
  blas::transpose(nA, mA, tmp1, tmp);
  mblC->addGeM(tmp, mA, eps, rankmax);
  delete [] tmp;
}


// C += d A B  (A H-matrix, B herm. H-Matrix, C mblock)
template<class T> static
void mltaGeHHeH_toMbl_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                        mblock<T>** B, mblock<T>* mblC,
                        double eps, unsigned rankmax)
{
  assert(blA->getn2()==blB->getn1() && blB->getn1()==blB->getn2());

  if (blA->isleaf() && blA->isLrM(A))
    mltaLrMHeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blB->isleaf() && blB->isGeM(B))
    mltaGeHHeM_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else if (blA->isleaf() && blA->isGeM(A))
    mltaGeMHeH_toMbl_(d, blA, A, blB, B, mblC, eps, rankmax);
  else {        // A and B are no leaves, hence the clusters can be subdivided
    assert(blA->getncs()==blB->getnrs());
    unsigned k, ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getncs();

    // unify first block row
    const unsigned m0 = blA->getson(0, 0)->getn1();
    unsigned m = m0;
    const unsigned n0 = blB->getson(0, 0)->getn2();
    unsigned n = n0;
    mblock<T>* R = new mblock<T>(m0, n0);
    blcluster *bl1, *bl2;

    bl1 = blA->getson(0, 0);
    bl2 = blB->getson(0, 0);
    mltaGeHHeH_toMbl_(d, bl1, A, bl2, B, R, eps, rankmax);

    for (k=1; k<nsp; ++k) {
      bl1 = blA->getson(0, k);
      bl2 = blB->getson(k, 0);
      mltaGeHGeH_toMbl(d, bl1, A, bl2, B, R, eps, rankmax);
    }

    for (unsigned j=1; j<ns2; ++j) {

      const unsigned nj = blB->getson(0, j)->getn2();
      mblock<T> R2(m0, nj);

      for (k=0; k<j; ++k) {
        bl1 = blA->getson(0, k);
        bl2 = blB->getson(k, j);
        mltaGeHGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
      }

      bl1 = blA->getson(0, j);
      bl2 = blB->getson(j, j);
      mltaGeHHeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

      for (k=j+1; k<nsp; ++k) {
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
      mltaGeHHeH_toMbl_(d, bl1, A, bl2, B, Rz, eps, rankmax);

      for (k=1; k<nsp; ++k) {
        bl1 = blA->getson(i, k);
        bl2 = blB->getson(0, k);
        mltaGeHGeHh_toMbl(d, bl1, A, bl2, B, Rz, eps, rankmax);
      }

      for (unsigned j=1; j<ns2; ++j) {

        const unsigned nj = blB->getson(0, j)->getn2();
        mblock<T> R2(mi, nj);

        for (k=0; k<j; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(k, j);
          mltaGeHGeH_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
        }

        bl1 = blA->getson(i, j);
        bl2 = blB->getson(j, j);
        mltaGeHHeH_toMbl_(d, bl1, A, bl2, B, &R2, eps, rankmax);

        for (k=j+1; k<nsp; ++k) {
          bl1 = blA->getson(i, k);
          bl2 = blB->getson(j, k);
          mltaGeHGeHh_toMbl(d, bl1, A, bl2, B, &R2, eps, rankmax);
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


// C += d A B  (A,C H-matrices, B symm. H-matrix)
template<class T> static
void mltaGeHHeH_(T d, blcluster* blA, mblock<T>** A, blcluster* blB,
                 mblock<T>** B, blcluster* blC, mblock<T>** C,
                 double eps, unsigned rankmax)
{
  if (blC->isleaf())               // C is a leaf
    mltaGeHHeH_toMbl_(d, blA, A, blB, B, C[blC->getidx()], eps, rankmax);
  else {                            // C is not a leaf
    if (!blA->isleaf() && !blB->isleaf()) {
      unsigned ns1 = blA->getnrs(), nsp = blA->getncs(), ns2 = blB->getncs();
      unsigned k;
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster *bl1, *bl2, *bl = blC->getson(i, j);

          for (k=0; k<j; ++k) {
            bl1 = blA->getson(i, k);
            bl2 = blB->getson(k, j);
            mltaGeHGeH(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }

          bl1 = blA->getson(i, j);
          bl2 = blB->getson(j, j);
          mltaGeHHeH_(d, bl1, A, bl2, B, bl, C, eps, rankmax);

          for (k=j+1; k<nsp; ++k) {
            bl1 = blA->getson(i, k);
            bl2 = blB->getson(j, k);
            mltaGeHGeHh(d, bl1, A, bl2, B, bl, C, eps, rankmax);
          }
        }
    }
    else if (blA->isleaf() && blA->isLrM(A))
      mltaLrMHeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else if (blB->isleaf() && blB->isGeM(B))
      mltaGeHHeM_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    else {                                       // A must be a leaf and dense
      assert(blA->isleaf() && blA->isGeM(A));
      mltaGeMHeH_toGeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Instanzen
//

void mltaGeHHeH(double d, blcluster* blA, mblock<double>** A, blcluster* blB,
                mblock<double>** B, blcluster* blC, mblock<double>** C,
                double eps, unsigned rankmax)
{
  mltaGeHHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHHeH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
                mblock<float>** B, blcluster* blC, mblock<float>** C,
                double eps, unsigned rankmax)
{
  mltaGeHHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHHeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
                mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
                double eps, unsigned rankmax)
{
  mltaGeHHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

void mltaGeHHeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
                mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
                double eps, unsigned rankmax)
{
  mltaGeHHeH_(d, blA, A, blB, B, blC, C, eps, rankmax);
}

