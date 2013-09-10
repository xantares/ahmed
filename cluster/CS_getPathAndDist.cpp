/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "CS_getPathAndDist.h"
#include "basmod.h"

unsigned getPathAndDist (unsigned s, unsigned t, int* tree,
   unsigned nrows, unsigned *row_ptr, unsigned *col_idx)
{
    bool nfound (true);

    // stores row indices in the actual component
    unsigned *nbrsf = new unsigned[2*nrows];
    unsigned  *old_nbrs = nbrsf, *new_nbrs = old_nbrs+nrows;
    unsigned nnbs = 0, nonbs = 1;
    old_nbrs[0] = s;

    unsigned cnt = 1, dist = 0;

    while (cnt<nrows && nonbs>0 && nfound) {
      nnbs = 0;
      for (unsigned j=0; j<nonbs && nfound; ++j) {
	unsigned r = old_nbrs[j], idx = row_ptr[r+1];
	for (unsigned k=row_ptr[r]; k<idx && nfound; ++k) {
	  unsigned c = col_idx[k];
	  // if this node isn't visited yet =>
	  if (tree[c] == -1) {
	    tree[c] = r;
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
    return dist;
}
