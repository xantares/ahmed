/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <limits>
#include "basmod.h"

// checks whether a given CRS matrix is connected and
// returns an approximation of the diameter
bool isConnected(unsigned n, unsigned *iA, unsigned *jA, unsigned &diam)
{
  unsigned *nbrsf = new unsigned[2*n];
  unsigned *new_nbrs = nbrsf+n, *old_nbrs = nbrsf;
  unsigned nnbs, nonbs;

  // store whether an index belongs to the component of i0
  bool *comp = new bool[n];
  for (unsigned i=1; i<n; i++) comp[i] = false;
  comp[0] = true;
  old_nbrs[0] = 0;
  nonbs = 1;

  diam = 0;
  unsigned cnt = 1;
  while (cnt<n && nonbs>0) {
    nnbs = 0;
    // go through all active neighbors
    for (unsigned j=0; j<nonbs; ++j) {
      const unsigned jp = old_nbrs[j], nbeg = iA[jp], nend = iA[jp+1];
      for (unsigned k=nbeg; k<nend; ++k) {
        const unsigned idx = jA[k];
        // has this index been checked ?
        if (!comp[idx]) {
          comp[idx] = true;
          new_nbrs[nnbs++] = idx;
        }
      }
    }
    ++diam;
    cnt += nnbs;
    swap(new_nbrs, old_nbrs);
    nonbs = nnbs;
  }

  delete [] comp;
  delete [] nbrsf;

  if (cnt == n) return true;
  else {
    diam = std::numeric_limits<unsigned>::max();
    return false;
  }
}
