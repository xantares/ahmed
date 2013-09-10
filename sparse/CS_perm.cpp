/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "sparse.h"

// permute CRS matrix according to op_perm(orig->perm) and po_perm(perm->orig)
template<class T>
void CS_perm_(unsigned n, T* A, unsigned* jA, unsigned* iA,
              unsigned* po_perm, unsigned* op_perm,
              T*& B, unsigned*& jB, unsigned*& iB)
{
  unsigned nnz = iA[n];
  B = new T[nnz];
  jB = new unsigned[nnz];
  iB = new unsigned[n+1];
  T* pB = B;
  unsigned *pjB = jB, *piB = iB;

  nnz = 0;
  for (unsigned i=0; i<n; ++i) {
    litem<T>* list = NULL;
    unsigned oi = op_perm[i];
    for (unsigned k=iA[oi]; k<iA[oi+1]; ++k)
      add_item(list, po_perm[jA[k]], A[k]);

    pack_list(list, pB, pjB);
    *piB++ = nnz;
    nnz += iA[oi+1]-iA[oi];
  }

  *piB = nnz;
}

void CS_perm(unsigned n, double* A, unsigned* jA, unsigned* iA,
             unsigned* po_perm, unsigned* op_perm,
             double*& B, unsigned*& jB, unsigned*& iB)
{
  CS_perm_(n, A, jA, iA, po_perm, op_perm, B, jB, iB);
}

void CS_perm(unsigned n, dcomp* A, unsigned* jA, unsigned* iA,
             unsigned* po_perm, unsigned* op_perm,
             dcomp*& B, unsigned*& jB, unsigned*& iB)
{
  CS_perm_(n, A, jA, iA, po_perm, op_perm, B, jB, iB);
}
