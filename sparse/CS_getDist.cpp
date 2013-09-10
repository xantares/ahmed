/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "basmod.h"
#include <limits>

bool getDist (unsigned *membership, unsigned *nbrsf, unsigned size_t1, unsigned &max_dist, 
	unsigned nrows, unsigned *row_ptr, unsigned *col_idx)
{
  bool nfound (true);
  unsigned *new_nbrs = nbrsf+nrows, *old_nbrs = nbrsf;
  unsigned nnbs(0), nonbs(size_t1);

  unsigned cnt(size_t1), dist(0);

  while (cnt<nrows && nonbs>0 && nfound && dist < max_dist) {
    nnbs = 0;
    for (unsigned j(0); j<nonbs && nfound; ++j) {
      unsigned r (old_nbrs[j]), idx (row_ptr[r+1]);
      for (unsigned k(row_ptr[r]); k<idx && nfound; ++k) {
        unsigned c(col_idx[k]);
        // node c is member of cluster t2
        if (membership[c] == 2) nfound = false;
        // if this node isn't member of t1 or t2 yet =>
        if (membership[c] == 0) {
          membership[c] = 3; // make node c member of the neighbourhood of t1
          new_nbrs[nnbs++] = c;
        }
      }
    }
    swap(new_nbrs, old_nbrs);
    nonbs = nnbs;
    dist++;
    cnt += nnbs;
  }

  if (!nfound) max_dist = dist;
  return !nfound;
}

unsigned getDist (unsigned s, unsigned t, unsigned nrows, unsigned *row_ptr, unsigned *col_idx)
{
  bool nfound (true);
  unsigned *visited(new unsigned[nrows]); 

  for (unsigned k=0; k<nrows; k++) visited[k] = 0;
  visited[s] = 1;

  // stores row indices in the actual component
  unsigned *nbrsf = new unsigned[2*nrows];
  unsigned *new_nbrs = nbrsf+nrows, *old_nbrs = nbrsf;
  unsigned nnbs(0), nonbs(1);
  old_nbrs[0] = s;

  unsigned cnt(1), dist(0);

  while (cnt<nrows && 0<nonbs && nfound) {
    nnbs = 0;
    for (unsigned j(0); j<nonbs; ++j) {
      unsigned r (old_nbrs[j]), idx (row_ptr[r+1]);
      for (unsigned k(row_ptr[r]); k<idx && nfound; ++k) {
        unsigned c(col_idx[k]);
        // if this node isn't visited yet =>
        if (visited[c] == 0) {
          visited[c] = 1;
          new_nbrs[nnbs++] = c;
        }
        if (c == t) nfound = false;
      }
    }
    swap(new_nbrs, old_nbrs);
    nonbs = nnbs;
    dist++;
    cnt += nnbs;
  }
  	
  delete [] nbrsf;
  delete [] visited;

  if (!nfound) return dist;
  else return std::numeric_limits<unsigned>::max();
}

void reverseBFS (unsigned *membership, unsigned *nbrsf, unsigned size_t1,  
	unsigned nrows, unsigned *row_ptr, unsigned *col_idx)
{
  unsigned *new_nbrs = nbrsf+nrows, *old_nbrs = nbrsf;
  unsigned nnbs(0), nonbs(size_t1);

  unsigned cnt(size_t1);

  while (cnt<nrows && nonbs>0) {
    nnbs = 0;
    for (unsigned j(0); j<nonbs; ++j) {
      unsigned r (old_nbrs[j]), idx (row_ptr[r+1]);
      for (unsigned k(row_ptr[r]); k<idx; ++k) {
        unsigned c(col_idx[k]);
        // node c belongs to cluster U(t1)
        if (membership[c] == 3) { // c belongs to the neighbourhood of t1
          membership[c] = 0; // cancel the membership 
          new_nbrs[nnbs++] = c;
        }
      }
    }
    swap(new_nbrs, old_nbrs);
    nonbs = nnbs;
    cnt += nnbs;
  }
}