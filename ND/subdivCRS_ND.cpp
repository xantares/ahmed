/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "parallel.h"

template<class T> static
void getCRSBlock1_(const CRSblock<T>* const A, CRSblock<T>*& B,
                   const unsigned k, const unsigned l)
{
  B = new CRSblock<T>[4];
  B[0].iA = new unsigned[k+1];
  B[1].iA = new unsigned[k+1];
  B[2].iA = new unsigned[(A->size-k-l)+1];
  B[3].iA = new unsigned[(A->size-k-l)+1];
  
  B[0].size = B[1].size = k;
  B[2].size = B[3].size = A->size-k-l;
  
  //DryRun und Belegen der iAs
  //Block 0 und 1
  unsigned i, co1, co2; //Counter
  B[0].iA[0] = B[1].iA[0] = 0;
  for (i = 0; i < k; i++) {
    co1 = co2 = 0;
    for (unsigned j = A->iA[i]; j < A->iA[i+1]; j++) {
      if (A->jA[j] < k) ++co1;
      else ++co2;
    }
    B[0].iA[i+1] = B[0].iA[i] + co1;
    B[1].iA[i+1] = B[1].iA[i] + co2;
  }
  //Block 2 und 3
  B[2].iA[0] = B[3].iA[0] = 0;
  for (i = k+l; i < A->size; i++) {
    co1 = co2 = 0;
    for (unsigned j = A->iA[i]; j < A->iA[i+1]; j++) {
      if (A->jA[j] < k) ++co1;
      else if (A->jA[j] >= k+l) ++co2;
    }
    B[2].iA[i-(k+l)+1] = B[2].iA[(i-(k+l))] + co1;
    B[3].iA[i-(k+l)+1] = B[3].iA[(i-(k+l))] + co2;
  }



  //restliche Felder anlegen
  B[0].A = new T[B[0].iA[k]];
  B[0].jA = new unsigned[B[0].iA[k]];

  B[1].A = new T[B[1].iA[k]];
  B[1].jA = new unsigned[B[1].iA[k]];

  B[2].A = new T[B[2].iA[A->size-(k+l)]];
  B[2].jA = new unsigned[B[2].iA[A->size-(k+l)]];

  B[3].A = new T[B[3].iA[A->size-(k+l)]];
  B[3].jA = new unsigned[B[3].iA[A->size-(k+l)]];
  //Felder beschreiben
  //Block 0 und 1
  co1 = co2 = 0;
  for (i = 0; i < A->iA[k]; i++) {
    if (A->jA[i] < k) {
      B[0].A[co1] = A->A[i];
      B[0].jA[co1++] = A->jA[i];
    }
    else {
      B[1].A[co2] = A->A[i];
      B[1].jA[co2++] = A->jA[i]-(k+l);
    }
  }
  //Block 2 und 3
  co1 = co2 = 0;
  for (i = A->iA[k+l]; i < A->iA[A->size];i++) {
    if (A->jA[i] < k) {
      B[2].A[co1] = A->A[i];
      B[2].jA[co1++] = A->jA[i];
    }
    else if (A->jA[i] < k+l);
    else {
      B[3].A[co2] = A->A[i];
      B[3].jA[co2++] = A->jA[i]-(k+l);
    }
  }
}

template<class T> static
void getCRSBlock2_(const CRSblock<T>* const A, CRSblock<T>*& C,
                   const unsigned k, const unsigned l)
{
  C = new CRSblock<T>[3];
  C[0].iA = new unsigned[l+1];
  C[1].iA = new unsigned[l+1];
  C[2].iA = new unsigned[(A->size-k-l)+1];

  C[0].size = C[1].size = l;
  C[2].size = A->size-k-l;

  //DryRun und Belegen der iAs
  //Block 0 und 1
  unsigned i, co1, co2; //Counter
  C[0].iA[0] = C[1].iA[0] = 0;
  for (i = k; i < k+l; i++) {
    co1 = co2 = 0;
    for (unsigned j = A->iA[i]; j < A->iA[i+1]; j++) {
      if (A->jA[j] >= k+l) ++co2;
      else ++co1;
    }
    C[0].iA[i-k+1] = C[0].iA[i-k] + co1;
    C[1].iA[i-k+1] = C[1].iA[i-k] + co2;
  }
  //Block 2
  C[2].iA[0] = 0;
  for (i = k+l; i < A->size; i++) {
    co1 = 0;
    for (unsigned j = A->iA[i]; j < A->iA[i+1]; j++) {
      if (A->jA[j] < k);
      else if (A->jA[j] < k+l) ++co1;
    }
    C[2].iA[i-(k+l)+1] = C[2].iA[(i-(k+l))] + co1;
  }

  //restliche Felder anlegen
  C[0].A = new T[C[0].iA[l]];
  C[0].jA = new unsigned[C[0].iA[l]];

  C[1].A = new T[C[1].iA[l]];
  C[1].jA = new unsigned[C[1].iA[l]];

  C[2].A = new T[C[2].iA[A->size-(k+l)]];
  C[2].jA = new unsigned[C[2].iA[A->size-(k+l)]];
  //Felder beschreiben
  //Block 0 und 1
  co1 = co2 = 0;
  for (i = A->iA[k]; i < A->iA[k+l]; i++) {
    if (A->jA[i] >= k+l) {
      C[1].A[co1] = A->A[i];
      C[1].jA[co1++] = A->jA[i]-(k+l);
    }
    else {
      C[0].A[co2] = A->A[i];
      C[0].jA[co2++] = A->jA[i]-k;
    }
  }
  //Block 2
  co1 = 0;
  for (i = A->iA[k+l]; i < A->iA[A->size];i++) {
    if (A->jA[i] < k);
    else if (A->jA[i]<k+l) {
      C[2].A[co1] = A->A[i];
      C[2].jA[co1++] = A->jA[i]-k;
    }
  }
}

//für die waagerechten
template<class T> static
void getCRS1Block1_(const CRSblock<T>* const A, CRSblock<T>*& D,
                    const unsigned k, const unsigned l)
{
  D = new CRSblock<T>[2];
  D[0].iA = new unsigned[A->size+1];
  D[1].iA = new unsigned[A->size+1];

  D[0].size = D[1].size = A->size;

  unsigned i, co1, co2; //Counter
  D[0].iA[0] = D[1].iA[0] = 0;
  for (i = 0; i < A->size; i++) {
    co1 = co2 = 0;
    for (unsigned j = A->iA[i]; j < A->iA[i+1]; j++) {
      if (A->jA[j] < k) ++co1;
      else if (A->jA[j] < k+l);
      else ++co2;
    }
    D[0].iA[i+1] = D[0].iA[i] + co1;
    D[1].iA[i+1] = D[1].iA[i] + co2;
  }

  D[0].A = new T[D[0].iA[A->size]];
  D[0].jA = new unsigned[D[0].iA[A->size]];
  D[1].A = new T[D[1].iA[A->size]];
  D[1].jA = new unsigned[D[1].iA[A->size]];

  co1 = co2 = 0;
  for (i = 0; i < A->iA[A->size]; i++) {
    if (A->jA[i] < k) {
      D[0].A[co1] = A->A[i];
      D[0].jA[co1++] = A->jA[i];
    }
    else if (A->jA[i] < k+l);
    else {
      D[1].A[co2] = A->A[i];
      D[1].jA[co2++] = A->jA[i]-(k+l);
    }
  }
}

template<class T> static
void getCRS1Block2_(const CRSblock<T>* const A, CRSblock<T>*& E,
                    const unsigned k, const unsigned l)
{
  E = new CRSblock<T>[1];

  E[0].iA = new unsigned[A->size+1];
  E[0].size = A->size;
  unsigned i, co1; //Counter
  E[0].iA[0] = 0;
  for (i = 0; i < A->size; i++) {
    co1 = 0;
    for (unsigned j = A->iA[i]; j < A->iA[i+1]; j++) {
      if (A->jA[j] < k);
      else if (A->jA[j] < k+l) ++co1;

    }
    E[0].iA[i+1] = E[0].iA[i] + co1;
  }

  E[0].A = new T[E[0].iA[A->size]];
  E[0].jA = new unsigned[E[0].iA[A->size]];

  co1 = 0;
  for (i = 0; i < A->iA[A->size]; i++) {
    if (A->jA[i] < k);
    else if (A->jA[i] < k+l) {
      E[0].A[co1] = A->A[i];
      E[0].jA[co1++] = A->jA[i]-k;
    }
  }
}

//für die senkrechten
template<class T> static
void getCRS2Block1_(const CRSblock<T>* const A, CRSblock<T>*& G,
                    const unsigned k, const unsigned l)
{
  G = new CRSblock<T>[2];
  G[0].iA = new unsigned[k+1];
  G[1].iA = new unsigned[(A->size-k-l)+1];

  G[0].size = k;
  G[1].size = A->size-k-l;

  unsigned i, co = 0; //Counter

  G[0].iA[0] = G[1].iA[0] = 0;
  for (i = 0; i < k; i++) {
    co = A->iA[i+1]-A->iA[i];
    G[0].iA[i+1] = G[0].iA[i] + co;
  }

  for (i = k+l; i < A->size; i++) {
    co = A->iA[i+1]-A->iA[i];
    G[1].iA[i-k-l+1] = G[1].iA[i-k-l] + co;
  }


  G[0].A = new T[G[0].iA[k]];
  G[0].jA = new unsigned[G[0].iA[k]];
  G[1].A = new T[G[1].iA[A->size-(k+l)]];
  G[1].jA = new unsigned[G[1].iA[A->size-(k+l)]];

  co = 0;
  for (i = 0; i < A->iA[k]; i++) {
    G[0].A[co] = A->A[i];
    G[0].jA[co++] = A->jA[i];
  }
  co = 0;
  for (i = A->iA[k+l]; i < A->iA[A->size]; i++) {
    G[1].A[co] = A->A[i];
    G[1].jA[co++] = A->jA[i];
  }
}

template<class T> static
void getCRS2Block2_(const CRSblock<T>* const A, CRSblock<T>*& F,
                    const unsigned k, const unsigned l)
{
  F = new CRSblock<T>[1];
  F[0].iA = new unsigned[l+1];

  F[0].size = l;

  unsigned i;
  unsigned co = 0; //Counter

  F[0].iA[0] = 0;
  for (i = k; i < k+l; i++) {
    co = A->iA[i+1]-A->iA[i];
    F[0].iA[i-k+1] = F[0].iA[i-k] + co;
  }

  F[0].A = new T[F[0].iA[l]];
  F[0].jA = new unsigned[F[0].iA[l]];

  co = 0;
  for (i = A->iA[k]; i < A->iA[k+l]; i++) {
    F[0].A[co] = A->A[i];
    F[0].jA[co++] = A->jA[i];
  }

}
/*****************************************************************************
 *****************************************************************************/

template<class T> static
void subdivideCRS1_ND_(const unsigned begp, const unsigned p,
                       const blcluster* const rootH, CRSblock<T>* H)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p!=1 && !rootH->isleaf()) {
    if (begp<=rank && rank<begp+p/2) {
      getCRS1Block1_(H, H->sons, rootH->getson(0,0)->getn1(),
                     rootH->getson(1,0)->getn1());
      subdivideCRS1_ND_(begp, p/2, rootH->getson(0,0), &H->sons[0]);
    }else {
      getCRS1Block2_(H, H->sons, rootH->getson(0,0)->getn1(),
                     rootH->getson(1,0)->getn1());
      subdivideCRS1_ND_(begp+p/2, p/2, rootH->getson(1,0), &H->sons[0]);
    }
  }
}

template<class T> static
void subdivideCRS2_ND_(const unsigned begp, const unsigned p,
                       const blcluster* const rootH, CRSblock<T>* H)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p!=1 && !rootH->isleaf()) {
    if (begp<=rank && rank<begp+p/2) {
      getCRS2Block1_(H, H->sons, rootH->getson(0,0)->getn1(),
                     rootH->getson(1,0)->getn1());
      subdivideCRS2_ND_(begp, p/2, rootH->getson(0,0), &H->sons[0]);
    }
    else {
      getCRS2Block2_(H, H->sons, rootH->getson(0,0)->getn1(),
                     rootH->getson(1,0)->getn1());
      subdivideCRS2_ND_(begp+p/2, p/2, rootH->getson(1,0), &H->sons[0]);
    }
  }
}

template<class T> static
void subdivideCRS_ND_(const unsigned begp, const unsigned p,
                      const blcluster* const rootH, CRSblock<T>* H)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p!=1 && !rootH->isleaf()) {
    unsigned nsons = rootH->getncs();
    if (begp<=rank && rank<begp+p/2){
      
      getCRSBlock1_(H, H->sons, rootH->getson(0,0)->getn1(),
                    rootH->getson(1,1)->getn1());
      
      subdivideCRS_ND_(begp, p/2, rootH->getson(0,0), &H->sons[0]);
      if (nsons==3) {      
	subdivideCRS1_ND_(begp, p/2, rootH->getson(0,2), &H->sons[2]);
      	subdivideCRS2_ND_(begp, p/2, rootH->getson(0,2), &H->sons[1]);
      }
    } else {
      getCRSBlock2_(H, H->sons, rootH->getson(0,0)->getn1(),
                    rootH->getson(1,1)->getn1());
      
      subdivideCRS_ND_(begp+p/2, p/2, rootH->getson(1,1), &H->sons[0]);
      if (nsons==3) {      
	subdivideCRS1_ND_(begp+p/2, p/2, rootH->getson(1,2), &H->sons[2]);
      	subdivideCRS2_ND_(begp+p/2, p/2, rootH->getson(1,2), &H->sons[1]);
      }
    }
  }
}


void subdivideCRS_ND(const unsigned p, const blcluster* const rootH,
                     CRSblock<float>* H)
{
  subdivideCRS_ND_(0, p, rootH, H);
}

void subdivideCRS_ND(const unsigned p, const blcluster* const rootH,
                     CRSblock<double>* H)
{
  subdivideCRS_ND_(0, p, rootH, H);
}

void subdivideCRS_ND(const unsigned p, const blcluster* const rootH,
                     CRSblock<scomp>* H)
{
  subdivideCRS_ND_(0, p, rootH, H);
}

void subdivideCRS_ND(const unsigned p, const blcluster* const rootH,
                     CRSblock<dcomp>* H)
{
  subdivideCRS_ND_(0, p, rootH, H);
}
