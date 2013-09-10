/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "basmod.h"

unsigned getRemoteNode (const unsigned* const perm, const unsigned* const iperm,
		unsigned beg, unsigned end, unsigned s, unsigned nrows,
		unsigned *row_ptr, unsigned *col_idx, unsigned &diam)
{
  unsigned remotly_node (s);
  bool *visited(new bool[nrows]);

  for (unsigned k=0; k<nrows; k++) visited[k] = false;
  visited[s] = true;

  // stores row indices in the actual component
  unsigned *nbrsf = new unsigned[2*nrows];
  unsigned *new_nbrs = nbrsf+nrows, *old_nbrs = nbrsf;
  unsigned nnbs(0), nonbs(1);
  old_nbrs[0] = s;

  unsigned cnt(1);

  while (cnt<nrows && nonbs>0) {
    nnbs = 0;
    for (unsigned j(0); j<nonbs; ++j) {
      unsigned r (perm[old_nbrs[j]]), idx (row_ptr[r+1]);
      for (unsigned k(row_ptr[r]); k<idx; ++k) {
        unsigned c(iperm[col_idx[k]]);
        // if this node isn't visited yet =>
        if (!visited[c] && beg <= c && c < end) {
          visited[c] = true;
          new_nbrs[nnbs++] = c;
          remotly_node = c;
        }
      }
    }
    swap(new_nbrs, old_nbrs);
    nonbs = nnbs;
    diam++;
    cnt += nnbs;
  }

  delete [] nbrsf;
  delete [] visited;

  return remotly_node;
}

unsigned getRemotlyNode (unsigned s, unsigned n, unsigned *iA, unsigned *jA, unsigned &diam)
{
  unsigned remotly_node (s);
  bool *visited(new bool[n]);

  for (unsigned k=0; k<n; k++) visited[k] = false;
  visited[s] = true;

  // stores row indices in the actual component
  unsigned *nbrsf = new unsigned[2*n];
  unsigned *new_nbrs = nbrsf+n, *old_nbrs = nbrsf;
  unsigned nnbs(0), nonbs(1);
  old_nbrs[0] = s;

  unsigned cnt(1);

  while (cnt<n && nonbs>0) {
    nnbs = 0;
    for (unsigned j(0); j<nonbs; ++j) {
      unsigned r (old_nbrs[j]), idx (iA[r+1]);
      for (unsigned k(iA[r]); k<idx; ++k) {
        unsigned c(jA[k]);
        // if this node isn't visited yet =>
        if (!visited[c]) {
          visited[c] = true;
          new_nbrs[nnbs++] = c;
          remotly_node = c;
        }
      }
    }
    swap(new_nbrs, old_nbrs);
    nonbs = nnbs;
    diam++;
    cnt += nnbs;
  }

  delete [] nbrsf;
  delete [] visited;

  return remotly_node;
}
