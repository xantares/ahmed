/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "basmod.h"

// generate adjacency matrix from a CRS matrix
// INPUT:  the CRS matrix given by the arrays iA, jA
// OUTPUT: the symmetrized sparsity pattern returned in CRS format (upper part)
void genAdjM(unsigned n, unsigned* &iA, unsigned* &jA)
{
  unsigned i;
  // count entries of each row
  unsigned* iAn = new unsigned[n+1];
  for (i=0; i<=n; ++i) iAn[i] = 0;

  // go through all strictly lower triangular entries (i,j) and check
  // whether (j,i) exists in the upper triangular part

  // set n pointers to the beginning of each row
  unsigned* co = new unsigned[n];
  for (i=0; i<n; ++i) co[i] = iA[i];

  for (i=0; i<n; ++i)
    for (unsigned k=iA[i]; k<iA[i+1]; ++k) {
      const unsigned j = jA[k];
      if (i<j) ++iAn[i+1]; // upper triangular entries count
      else {               // lower triangular only if there is no counter part
        const unsigned k1 = iA[j], k2 = iA[j+1];
        if (i<jA[k1] || i>jA[k2-1]) ++iAn[j+1]; // i is out of bounds
        else {             // go through all uninspected entries in the jth row
          while (co[j]<k2 && i>jA[co[j]]) ++co[j];
          if (co[j]==k2 || i<jA[co[j]]) ++iAn[j+1];
        }
      }
    }

  // construct array iAn by summing up the contents of iAn
  // con is a set of pointer refering to iAn
  unsigned* con = new unsigned[n];
  co[0] = con[0] = 0;
  for (i=1; i<n; ++i) {
    co[i] = iA[i];
    con[i] = iAn[i];
    iAn[i+1] += iAn[i];
  }

  unsigned *jAn = new unsigned[iAn[n]];
  for (i=1; i<n; ++i)
    for (unsigned k=iA[i]; k<iA[i+1]; ++k) {
      const unsigned j = jA[k];
      // copy all transposed lower triangular entries and all upper
      // triangular elements up to that position
      if (j<i) {
        while (co[j]<iA[j+1] && i>jA[co[j]]) {
          if (jA[co[j]]>j) jAn[con[j]++] = jA[co[j]];
          ++co[j];
        }

        if (co[j]==iA[j+1] || i<=jA[co[j]]) {
          jAn[con[j]++] = i;
          ++co[i];
          if (i==jA[co[j]]) ++co[j];
        }
      }
    }

  // finish rows
  for (i=0; i<n; ++i)
    for (unsigned k=co[i]; k<iA[i+1]; ++k)
      if (i<jA[k]) jAn[con[i]++] = jA[k];

  swap(jA, jAn);
  swap(iA, iAn);

  delete [] jAn;
  delete [] con;
  delete [] co;
  delete [] iAn;

}

void genFullAdjMat (unsigned n, unsigned* &iA, unsigned* &jA)
{
  unsigned i;
  // count entries of each column
  unsigned* cnt = new unsigned[n];
  for (i=0; i<n; ++i) cnt[i] = 0;

  for (i=0; i<n; ++i) {
    unsigned j = iA[i];
    const unsigned idx = iA[i+1];
    while (j < idx) {
      cnt[jA[j]]++;
      j++;
    }
  }

  // summing up entries
  for (i=2; i<n; ++i) cnt[i] += cnt[i-1];

  unsigned* iAn = new unsigned[n+1]; // VALGRIND meldet hier Fehler
  iAn[0] = 0;
  for (i=1; i<=n; ++i) iAn[i] = iA[i] + cnt[i-1];

  unsigned *jAn = new unsigned [iAn[n]];
  for (unsigned k=0; k<n; k++) cnt[k] = iAn[k];

  for (i=0; i<n; ++i) {
    unsigned j = iA[i];
    const unsigned idx = iA[i+1];
    while (j < idx) {
      jAn[cnt[i]++] = jA[j];
      jAn[cnt[jA[j]]++] = i;
      j++;
    }
  }

  swap(jA, jAn);
  swap(iA, iAn);

  delete [] jAn;
  delete [] iAn;
  delete [] cnt;
}

// generate a weighted adjacency matrix from a CRS matrix
// INPUT:  the CRS matrix given by the arrays iA, jA, A
// OUTPUT: the symmetrized matrix returned in CRS format (upper part)
void genAdjM(unsigned n, unsigned* &iA, unsigned* &jA, unsigned* &A)
{
  unsigned i;
  // count entries of each row
  unsigned* iAn = new unsigned[n+1];
  for (i=0; i<=n; ++i) iAn[i] = 0;

  // go through all strictly lower triangular entries (i,j) and check
  // whether (j,i) exists in the upper triangular part

  // set n pointers to the beginning of each row
  unsigned* co = new unsigned[n];
  for (i=0; i<n; ++i) co[i] = iA[i];

  for (i=0; i<n; ++i) {
    for (unsigned k=iA[i]; k<iA[i+1]; ++k) {
      const unsigned j = jA[k];
      if (i<j) ++iAn[i+1]; // upper triangular entries count
      else {               // lower triangular only if there is no counter part
        const unsigned k1 = iA[j], k2 = iA[j+1];
        if (i<jA[k1] || i>jA[k2-1]) ++iAn[j+1]; // i is out of bounds
        else {             // go through all uninspected entries in the jth row
          while (co[j]<k2 && i>jA[co[j]]) ++co[j];
          if (co[j]==k2 || i<jA[co[j]]) ++iAn[j+1];
        }
      }
    }
  }

  // construct array iAn by summing up the contents of iAn
  // con is a set of pointer refering to iAn
  unsigned* con = new unsigned[n];
  co[0] = con[0] = 0;
  for (i=1; i<n; ++i) {
    co[i] = iA[i];
    con[i] = iAn[i];
    iAn[i+1] += iAn[i];
  }

  unsigned *jAn = new unsigned[iAn[n]];
  unsigned *An = new unsigned[iAn[n]];

  for (i=1; i<n; ++i)
    for (unsigned k=iA[i]; k<iA[i+1]; ++k) {
      const unsigned j = jA[k];
      // copy all transposed lower triangular entries and all upper
      // triangular elements up to that position
      if (j<i) {
        while (co[j]<iA[j+1] && i>jA[co[j]]) {
          if (jA[co[j]]>j) {
            An[con[j]] = A[co[j]]; // entry of jth row is copied
            jAn[con[j]++] = jA[co[j]];
          }
          ++co[j];
        }

        if (co[j]==iA[j+1] || i<=jA[co[j]]) {
          if (jA[co[j]] == i) An[con[j]] = (A[co[i]] + A[co[j]]) / 2;
          else An[con[j]] = A[co[i]]; // entry of ith row is copied in jth row
          jAn[con[j]++] = i;
          ++co[i];
          if (i==jA[co[j]]) ++co[j];
        }
      }
    }

  // finish rows
  for (i=0; i<n; ++i)
    for (unsigned k=co[i]; k<iA[i+1]; ++k)
      if (i<jA[k]) {
        An[con[i]] = A[k];
        jAn[con[i]++] = jA[k];
      }

  delete [] con;
  delete [] co;
  swap(jA, jAn);
  swap(iA, iAn);
  swap(A, An);
  delete [] iAn;
  delete [] jAn;
  delete [] An;
}

void genFullAdjMat(unsigned n, unsigned* &iA, unsigned* &jA, unsigned* &A)
{
  genAdjM (n, iA, jA, A);

  unsigned i;
  // count entries of each row
  unsigned* cnt = new unsigned[n+1];
  for (i=0; i<=n; ++i) cnt[i] = 0;

  for (i=0; i<n; ++i) {
    unsigned j = iA[i];
    const unsigned idx = iA[i+1];
    while (j < idx) {
      cnt[jA[j]]++;
      j++;
    }
  }

  // summing up entries
  for (i=1; i<n; ++i) cnt[i+1] += cnt[i];

  unsigned* iAn = new unsigned[n+1];
  iAn[0] = 0;
  for (i=1; i<=n; ++i) iAn[i] = iA[i] + cnt[i-1];

  unsigned *jAn = new unsigned[iAn[n]];
  unsigned *An = new unsigned[iAn[n]];
  for (unsigned k=0; k<=n; k++) cnt[k] = iAn[k];

  for (i=0; i<n; ++i) {
    unsigned j = iA[i];
    const unsigned idx = iA[i+1];
    while (j < idx) {
      An[cnt[i]] = A[j];
      An[cnt[jA[j]]] = A[j];
      jAn[cnt[i]++] = jA[j];
      jAn[cnt[jA[j]]++] = i;
      j++;
    }
  }

  swap(jA, jAn);
  swap(iA, iAn);
  swap(A, An);
  delete [] cnt;
  delete [] iAn;
  delete [] jAn;
  delete [] An;
}


