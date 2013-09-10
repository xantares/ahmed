/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// type 'L' LtH, 'H' H(default), 'U' UtH
template<class T> static
unsigned Hmax_rank_(blcluster* bl, mblock<T>** A, char type)
{
  unsigned k = 0;

  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    if (A[idx]->isLrM() && k<A[idx]->rank()) k = A[idx]->rank();
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      if (type!='U' || bl->isndbl())
        for (unsigned j=0; j<i; ++j) {
          blcluster* son = bl->getson(i, j);
          if (son) {
            unsigned r = Hmax_rank_(son, A, type);
            if (k<r) k = r;
          }
        }

      unsigned r = Hmax_rank_(bl->getson(i, i), A, type);
      if (k<r) k = r;

      if (type!='L' || bl->isndbl())
        for (unsigned j=i+1; j<ns2; ++j) {
          blcluster* son = bl->getson(i, j);
          if (son) {
            unsigned r = Hmax_rank_(son, A, type);
            if (k<r) k = r;
          }
        }
    }
  }
  return k;
}

// 'U' for upper, 'L' for lower triangular matrix, 'H' is default
template<class T> static
unsigned long sizeH_(blcluster* bl, mblock<T>** A, char type)
{
  if (bl->isleaf()) return A[bl->getidx()]->size();
  else {
    unsigned long size = 0;
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();

    if (type=='U' && bl->isdbl()) {
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=i; j<ns2; ++j)
          size += sizeH_(bl->getson(i, j), A, type);
    } else if (type=='L' && bl->isdbl()) {
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<=i; ++j)
          size += sizeH_(bl->getson(i, j), A, type);
    } else {
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j) {
          blcluster* son = bl->getson(i, j);
          if (son) size += sizeH_(son, A, type);
        }
    }
    return size;
  }
}


template<class T> static
unsigned long sizeH_(unsigned n, mblock<T>** A)
{
  unsigned long sum = 0;
  for (unsigned i=0; i<n; ++i)
    if (A[i]!=NULL) sum += A[i]->size();

  return sum;
}


void copymbl_(mblock<double>* A, mblock<float>* B)
{
  unsigned n1 = B->n1 = A->n1, n2 = B->n2 = A->n2;

  if (A->isLrM()) {
    const unsigned k = A->rank();
    B->setrank(k);
    if (k>0) blas::copy(k*(n1+n2), A->data, B->data);
  } else { // mbl is dense
    assert(A->isGeM());
    if (A->isHeM()) {
      B->setHeM();
      blas::copy(n1*(n1+1)/2, A->data, B->data);
    }
    else if (A->isLtM()) {
      B->setLtM();
      blas::copy(n1*(n1+1)/2, A->data, B->data);
    }
    else if (A->isUtM()) {
      B->setUtM();
      blas::copy(n1*(n1+1)/2, A->data, B->data);
    }
    else {
      B->setGeM();
      blas::copy(n1*n2, A->data, B->data);
    }
  }
}

void copymbl_(mblock<dcomp>* A, mblock<scomp>* B)
{
  unsigned n1 = B->n1 = A->n1, n2 = B->n2 = A->n2;

  if (A->isLrM()) {
    const unsigned k = A->rank();
    B->setrank(k);
    if (k>0) blas::copy(k*(n1+n2), A->data, B->data);
  } else { // mbl is dense
    assert(A->isGeM());
    if (A->isHeM()) {
      B->setHeM();
      blas::copy(n1*(n1+1)/2, A->data, B->data);
    }
    else if (A->isLtM()) {
      B->setLtM();
      blas::copy(n1*(n1+1)/2, A->data, B->data);
    }
    else if (A->isUtM()) {
      B->setUtM();
      blas::copy(n1*(n1+1)/2, A->data, B->data);
    }
    else {
      B->setGeM();
      blas::copy(n1*n2, A->data, B->data);
    }
  }
}

// copies A to B
// B is an array of nleaves(bl) mblock<float>*
void copyH(blcluster* bl, mblock<double>** A, mblock<float>** B)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    delete B[idx];
    B[idx] = new mblock<float>(bl->getn1(), bl->getn2());
    assert(B[idx]!=NULL);
    copymbl_(A[idx], B[idx]);
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) copyH(son, A, B);
      }
  }
}

// copies A to B
// B is an array of n mblock<float>*
void copyH(unsigned n, mblock<double>** A, mblock<float>** B)
{
  for (unsigned i=0; i<n; ++i) {
    delete B[i];
    B[i] = new mblock<float>(A[i]->getn1(), A[i]->getn2());
    assert(B[i]!=NULL);
    copymbl_(A[i], B[i]);
  }
}

// copies A to B
// B is an array of nleaves(bl) mblock<scomp>*
void copyH(blcluster* bl, mblock<dcomp>** A, mblock<scomp>** B)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    delete B[idx];
    B[idx] = new mblock<scomp>(bl->getn1(), bl->getn2());
    assert(B[idx]!=NULL);
    copymbl_(A[idx], B[idx]);
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) copyH(son, A, B);
      }
  }
}

// copies A to B
// B is an array of nleaves(bl) mblock<double>*
template<class T> static
void copyH_(blcluster* bl, mblock<T>** A, mblock<T>** B)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    delete B[idx];
    B[idx] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(B[idx]!=NULL);
    B[idx]->copy(*A[idx]);
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) copyH_(son, A, B);
      }
  }
}

template<class T> static void freembls_recursive_(blcluster* bl, mblock<T>** &A)
{
  if (bl->isleaf()) {
    delete A[bl->getidx()];
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) freembls_recursive_(bl->getson(i, j), A);
  }
}

template<class T> static void initGeH_0_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    A[bl->getidx()] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(A[bl->getidx()]!=NULL);
    if (bl->isdbl()) A[bl->getidx()]->init0_GeM(bl->getn1(), bl->getn2());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) initGeH_0_(bl->getson(i, j), A);
  }
}

template<class T> static void initHeH_0_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    A[bl->getidx()] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(A[bl->getidx()]!=NULL);
    A[bl->getidx()]->init0_HeM(bl->getn1());
  } else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      initHeH_0_(bl->getson(i, i), A);
      for (unsigned j=i+1; j<ns; ++j) initGeH_0_(bl->getson(i, j), A);
    }
  }
}

template<class T> static void initLtH_0_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    A[bl->getidx()] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(A[bl->getidx()]!=NULL);
    if (bl->isdbl()) A[bl->getidx()]->init0_LtM(bl->getn1());
  } else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      for (unsigned j=0; j<i; ++j) initGeH_0_(bl->getson(i, j), A);
      initLtH_0_(bl->getson(i, i), A);
    }
  }
}

template<class T> static void initUtH_0_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    A[bl->getidx()] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(A[bl->getidx()]!=NULL);
    if (bl->isdbl()) A[bl->getidx()]->init0_UtM(bl->getn1());
  } else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      initUtH_0_(bl->getson(i, i), A);
      for (unsigned j=i+1; j<ns; ++j) initGeH_0_(bl->getson(i, j), A);
    }
  }
}

template<class T> static void initGeH_Id_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    A[bl->getidx()] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(A[bl->getidx()]!=NULL);
    if (bl->isdbl()) A[bl->getidx()]->initId_GeM(bl->getn1());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j)
        initGeH_Id_(bl->getson(i, j), A);
  }
}

template<class T> static void initHeH_Id_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) {
    A[bl->getidx()] = new mblock<T>(bl->getn1(), bl->getn2());
    assert(A[bl->getidx()]!=NULL);
    A[bl->getidx()]->initId_HeM(bl->getn1());
  } else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      initHeH_Id_(bl->getson(i, i), A);
      for (unsigned j=i+1; j<ns; ++j) initGeH_Id_(bl->getson(i, j), A);
    }
  }
}

template<class T>
void setGeHzero_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) A[bl->getidx()]->freedata();
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j) setGeHzero_(bl->getson(i, j), A);
  }
}

template<class T>
void setHeHzero_(blcluster* bl, mblock<T>** A)
{
  if (bl->isleaf()) A[bl->getidx()]->freedata();
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      setHeHzero_(bl->getson(i, i), A);
      for (unsigned j=i+1; j<ns2; ++j) setGeHzero_(bl->getson(i, j), A);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Instanzen

unsigned Hmax_rank(blcluster* bl, mblock<double>** A, char type)
{
  return Hmax_rank_(bl, A, type);
}
unsigned Hmax_rank(blcluster* bl, mblock<float>** A, char type)
{
  return Hmax_rank_(bl, A, type);
}
unsigned Hmax_rank(blcluster* bl, mblock<dcomp>** A, char type)
{
  return Hmax_rank_(bl, A, type);
}
unsigned Hmax_rank(blcluster* bl, mblock<scomp>** A, char type)
{
  return Hmax_rank_(bl, A, type);
}

unsigned long sizeH(blcluster* bl, mblock<double>** A, char type)
{
  return sizeH_(bl, A, type);
}
unsigned long sizeH(blcluster* bl, mblock<float>** A, char type)
{
  return sizeH_(bl, A, type);
}
unsigned long sizeH(blcluster* bl, mblock<dcomp>** A, char type)
{
  return sizeH_(bl, A, type);
}
unsigned long sizeH(blcluster* bl, mblock<scomp>** A, char type)
{
  return sizeH_(bl, A, type);
}

unsigned long sizeH(unsigned n, mblock<double>** A)
{
  return sizeH_(n, A);
}
unsigned long sizeH(unsigned n, mblock<float>** A)
{
  return sizeH_(n, A);
}
unsigned long sizeH(unsigned n, mblock<dcomp>** A)
{
  return sizeH_(n, A);
}
unsigned long sizeH(unsigned n, mblock<scomp>** A)
{
  return sizeH_(n, A);
}

void copyH(blcluster* bl, mblock<double>** A, mblock<double>** B)
{
  copyH_(bl, A, B);
}
void copyH(blcluster* bl, mblock<float>** A, mblock<float>** B)
{
  copyH_(bl, A, B);
}
void copyH(blcluster* bl, mblock<dcomp>** A, mblock<dcomp>** B)
{
  copyH_(bl, A, B);
}
void copyH(blcluster* bl, mblock<scomp>** A, mblock<scomp>** B)
{
  copyH_(bl, A, B);
}

void freembls_recursive(blcluster* bl, mblock<double>** &A)
{
  freembls_recursive_(bl, A);
}
void freembls_recursive(blcluster* bl, mblock<float>** &A)
{
  freembls_recursive_(bl, A);
}
void freembls_recursive(blcluster* bl, mblock<dcomp>** &A)
{
  freembls_recursive_(bl, A);
}
void freembls_recursive(blcluster* bl, mblock<scomp>** &A)
{
  freembls_recursive_(bl, A);
}

void initGeH_0(blcluster* bl, mblock<double>** &A)
{
  allocmbls(bl, A);
  initGeH_0_(bl, A);
}
void initGeH_0(blcluster* bl, mblock<float>** &A)
{
  allocmbls(bl, A);
  initGeH_0_(bl, A);
}
void initGeH_0(blcluster* bl, mblock<dcomp>** &A)
{
  allocmbls(bl, A);
  initGeH_0_(bl, A);
}
void initGeH_0(blcluster* bl, mblock<scomp>** &A)
{
  allocmbls(bl, A);
  initGeH_0_(bl, A);
}

void initGeH_0_withoutAlloc(blcluster* bl, mblock<double>** &A)
{
  initGeH_0_(bl, A);
}
void initGeH_0_withoutAlloc(blcluster* bl, mblock<float>** &A)
{
  initGeH_0_(bl, A);
}
void initGeH_0_withoutAlloc(blcluster* bl, mblock<dcomp>** &A)
{
  initGeH_0_(bl, A);
}
void initGeH_0_withoutAlloc(blcluster* bl, mblock<scomp>** &A)
{
  initGeH_0_(bl, A);
}

void initHeH_0_withoutAlloc(blcluster* bl, mblock<double>** &A)
{
  initHeH_0_(bl, A);
}
void initHeH_0_withoutAlloc(blcluster* bl, mblock<float>** &A)
{
  initHeH_0_(bl, A);
}
void initHeH_0_withoutAlloc(blcluster* bl, mblock<dcomp>** &A)
{
  initHeH_0_(bl, A);
}
void initHeH_0_withoutAlloc(blcluster* bl, mblock<scomp>** &A)
{
  initHeH_0_(bl, A);
}

void initHeH_0(blcluster* bl, mblock<double>** &A)
{
  allocmbls(bl, A);
  initHeH_0_(bl, A);
}
void initHeH_0(blcluster* bl, mblock<float>** &A)
{
  allocmbls(bl, A);
  initHeH_0_(bl, A);
}
void initHeH_0(blcluster* bl, mblock<dcomp>** &A)
{
  allocmbls(bl, A);
  initHeH_0_(bl, A);
}
void initHeH_0(blcluster* bl, mblock<scomp>** &A)
{
  allocmbls(bl, A);
  initHeH_0_(bl, A);
}

void initLtH_0(blcluster* bl, mblock<double>** &A)
{
  allocmbls(bl, A);
  initLtH_0_(bl, A);
}
void initLtH_0(blcluster* bl, mblock<float>** &A)
{
  allocmbls(bl, A);
  initLtH_0_(bl, A);
}
void initLtH_0(blcluster* bl, mblock<dcomp>** &A)
{
  allocmbls(bl, A);
  initLtH_0_(bl, A);
}
void initLtH_0(blcluster* bl, mblock<scomp>** &A)
{
  allocmbls(bl, A);
  initLtH_0_(bl, A);
}

void initLtH_0_withoutAlloc(blcluster* bl, mblock<double>** &A)
{
  initLtH_0_(bl, A);
}
void initLtH_0_withoutAlloc(blcluster* bl, mblock<float>** &A)
{
  initLtH_0_(bl, A);
}
void initLtH_0_withoutAlloc(blcluster* bl, mblock<dcomp>** &A)
{
  initLtH_0_(bl, A);
}
void initLtH_0_withoutAlloc(blcluster* bl, mblock<scomp>** &A)
{
  initLtH_0_(bl, A);
}

void initUtH_0(blcluster* bl, mblock<double>** &A)
{
  allocmbls(bl, A);
  initUtH_0_(bl, A);
}
void initUtH_0(blcluster* bl, mblock<float>** &A)
{
  allocmbls(bl, A);
  initUtH_0_(bl, A);
}
void initUtH_0(blcluster* bl, mblock<dcomp>** &A)
{
  allocmbls(bl, A);
  initUtH_0_(bl, A);
}
void initUtH_0(blcluster* bl, mblock<scomp>** &A)
{
  allocmbls(bl, A);
  initUtH_0_(bl, A);
}

void initUtH_0_withoutAlloc(blcluster* bl, mblock<double>** &A)
{
  initUtH_0_(bl, A);
}
void initUtH_0_withoutAlloc(blcluster* bl, mblock<float>** &A)
{
  initUtH_0_(bl, A);
}
void initUtH_0_withoutAlloc(blcluster* bl, mblock<dcomp>** &A)
{
  initUtH_0_(bl, A);
}
void initUtH_0_withoutAlloc(blcluster* bl, mblock<scomp>** &A)
{
  initUtH_0_(bl, A);
}

void initGeH_Id(blcluster* bl, mblock<double>** &A)
{
  allocmbls(bl, A);
  initGeH_Id_(bl, A);
}
void initGeH_Id(blcluster* bl, mblock<float>** &A)
{
  allocmbls(bl, A);
  initGeH_Id_(bl, A);
}
void initGeH_Id(blcluster* bl, mblock<dcomp>** &A)
{
  allocmbls(bl, A);
  initGeH_Id_(bl, A);
}
void initGeH_Id(blcluster* bl, mblock<scomp>** &A)
{
  allocmbls(bl, A);
  initGeH_Id_(bl, A);
}

void initHeH_Id(blcluster* bl, mblock<double>** &A)
{
  allocmbls(bl, A);
  initHeH_Id_(bl, A);
}
void initHeH_Id(blcluster* bl, mblock<float>** &A)
{
  allocmbls(bl, A);
  initHeH_Id_(bl, A);
}
void initHeH_Id(blcluster* bl, mblock<dcomp>** &A)
{
  allocmbls(bl, A);
  initHeH_Id_(bl, A);
}
void initHeH_Id(blcluster* bl, mblock<scomp>** &A)
{
  allocmbls(bl, A);
  initHeH_Id_(bl, A);
}

void setGeHzero(blcluster* bl, mblock<double>** A)
{
  setGeHzero_(bl, A);
}
void setGeHzero(blcluster* bl, mblock<float>** A)
{
  setGeHzero_(bl, A);
}
void setGeHzero(blcluster* bl, mblock<dcomp>** A)
{
  setGeHzero_(bl, A);
}
void setGeHzero(blcluster* bl, mblock<scomp>** A)
{
  setGeHzero_(bl, A);
}

void setHeHzero(blcluster* bl, mblock<double>** A)
{
  setHeHzero_(bl, A);
}
void setHeHzero(blcluster* bl, mblock<float>** A)
{
  setHeHzero_(bl, A);
}
void setHeHzero(blcluster* bl, mblock<dcomp>** A)
{
  setHeHzero_(bl, A);
}
void setHeHzero(blcluster* bl, mblock<scomp>** A)
{
  setHeHzero_(bl, A);
}
