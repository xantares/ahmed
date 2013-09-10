/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "AdjMatrix.h"
#include "quicksort.h"
//#include "specialSort.h"
//#include "binSearch.h" // for getSubMat

int AdjMat::isElement(unsigned r, unsigned c) const
{
  if (r>=n) return -1;
  const unsigned end = iA[r+1];
  for (unsigned k=iA[r]; k<end; k++) if (jA[k]==c) return k;
  return -1;
}

unsigned AdjMat::getElement(unsigned r, unsigned c) const
{
  const unsigned end = iA[r+1];
  for (unsigned k=iA[r]; k<end; k++) if (jA[k] == c) return A[k];
  return std::numeric_limits<unsigned>::max();
}

void AdjMat::getDiam(unsigned &diam) const
{
  unsigned node0 = 0, tmp_diam = 0;
  diam = 0;
  unsigned node1 = getRemotlyNode(node0, n, iA, jA, tmp_diam);
  node0 = getRemotlyNode(node1, n, iA, jA, diam);
  if (tmp_diam > diam) diam = tmp_diam;
  
  tmp_diam = 0;
  node1 = getRemotlyNode(node0, n, iA, jA, tmp_diam);
  if (tmp_diam > diam) diam = tmp_diam;
  
  tmp_diam = 0;
  node0 = getRemotlyNode(node1, n, iA, jA, tmp_diam);
  if (tmp_diam > diam) diam = tmp_diam;
}

void AdjMat::getDiam(unsigned &diam, const unsigned* const perm,
		     const unsigned* const iperm,
		     unsigned beg, unsigned end) const
{
  unsigned node0 = 0, tmp_diam = 0;
  diam = 0;
  unsigned node1 = getRemoteNode(perm, iperm, beg, end, node0, n, iA, jA,
				 tmp_diam);
  node0 = getRemoteNode(perm, iperm, beg, end, node1, n, iA, jA, diam);
  if (tmp_diam > diam) diam = tmp_diam;

  tmp_diam = 0;
  node1 = getRemoteNode(perm, iperm, beg, end, node0, n, iA, jA, tmp_diam);
  if (tmp_diam > diam) diam = tmp_diam;

  tmp_diam = 0;
  node0 = getRemoteNode(perm, iperm, beg, end, node1, n, iA, jA, tmp_diam);
  if (tmp_diam > diam) diam = tmp_diam;
}

/*
AdjMat* AdjMat::getMat(unsigned beg, unsigned end,
                       const unsigned* const op_perm,
                       const unsigned* const po_perm,
                       unsigned &nedges, unsigned* &edges) const
{
  const unsigned s = end-beg;
  unsigned nnz = 0;
  unsigned *iA_blk = new unsigned[s+1], *jA_blk, *A_blk;
  iA_blk[0] = 0;

  for (unsigned k=0; k<s; k++) {
    const unsigned r = op_perm[beg+k], idx = iA[r+1];
    for (unsigned j=iA[r]; j<idx; ++j) {
      const unsigned pc = po_perm[jA[j]];
      if (beg <= pc && pc < end) nnz++;
    }
    iA_blk[k+1] = nnz;
  }

  // determine the number of components ncomp
  // and the nodes of the components in the field comp
  unsigned *comp = new unsigned[s], ncomp;
  getComponents(beg, end, op_perm, po_perm, ncomp, comp);

  if (ncomp>1) {
    edges = new unsigned[2*ncomp];

    // determine nodes to insert edges between in order to make the components
    // connected
    unsigned *comp_nodes = new unsigned[ncomp];
    for (unsigned k=0; k<ncomp; ++k) comp_nodes[k] = 0;
    for (unsigned k=0; k<s; ++k)
      if (comp_nodes[comp[k]-1] == 0) comp_nodes[comp[k]-1] = op_perm[beg+k];

    delete [] comp;

    // count the edges for every node
    unsigned *cnt = new unsigned[s];
    for (unsigned k=0; k<s; ++k) cnt[k] = iA_blk[k+1]-iA_blk[k];
    // count the edges necessary for connectivity
    for (unsigned k=0; k+1<ncomp; ++k) {
      // linking component k with component k+1
      ++cnt[po_perm[comp_nodes[k]]-beg];
      ++cnt[po_perm[comp_nodes[k+1]]-beg];
      edges[2*nedges] = po_perm[comp_nodes[k]];
      edges[2*nedges+1] = po_perm[comp_nodes[k+1]];
      nedges++;
    }

    // rebuild iA_blk
    for (unsigned k=1; k<=s; ++k) iA_blk[k] = iA_blk[k-1] + cnt[k-1];

    jA_blk = new unsigned[iA_blk[s]];
    A_blk = new unsigned[iA_blk[s]];
    nnz=0;

    for (unsigned k=0; k<s; ++k) cnt[k] = iA_blk[k];
    // insert "normal" edges
    for (unsigned k=0; k<s; ++k) {
      unsigned r = op_perm[beg+k], idx = iA[r+1];
      for (unsigned j=iA[r]; j<idx; ++j) {
        unsigned pc = po_perm[jA[j]];
        if (beg <= pc && pc < end) {
          jA_blk[cnt[k]] = pc-beg;
          A_blk[cnt[k]] = 1;
          cnt[k]++;
          nnz++;
        }
      }
    }
    // insert additional edges
    for (unsigned k=0; k+1<ncomp; ++k) {
      jA_blk[cnt[po_perm[comp_nodes[k]]-beg]] = po_perm[comp_nodes[k+1]]-beg;
      jA_blk[cnt[po_perm[comp_nodes[k+1]]-beg]] = po_perm[comp_nodes[k]]-beg;
      ++cnt[po_perm[comp_nodes[k]]-beg];
      ++cnt[po_perm[comp_nodes[k+1]]-beg];
    }

    delete [] cnt;
    delete [] comp_nodes;
  } else {
    jA_blk = new unsigned[iA_blk[s]];
    A_blk = new unsigned[iA_blk[s]];
    nnz=0;

    for (unsigned k=0; k<s; k++) {
      unsigned r = op_perm[beg+k], idx = iA[r+1];
      for (unsigned j=iA[r]; j<idx; ++j) {
        unsigned pc = po_perm[jA[j]];
        if (beg <= pc && pc < end) {
          jA_blk[nnz] = pc - beg;
          A_blk[nnz] = 1;
          ++nnz;
        }
      }
    }
    delete [] comp;
  }

  for (unsigned k=0; k<s; ++k) quickSort(jA_blk, iA_blk[k], iA_blk[k+1]);

  return new AdjMat(s, iA_blk, jA_blk, A_blk);
}
*/

AdjMat* AdjMat::getMat(unsigned beg, unsigned end,
		       const unsigned* const op_perm,
		       const unsigned* const po_perm) const
{
  const unsigned nsize = end-beg; // size of new matrix
  unsigned i;
  
  unsigned *pos = new unsigned[nsize];
  for (i=0; i<nsize; i++) pos[i] = 0;
	
  for (i=beg; i<end; i++) {
    const unsigned r = op_perm[i], idx = iA[r+1]; // row r in original matrix
    for (unsigned j=iA[r]; j<idx; j++) {
      const unsigned c = po_perm[jA[j]]; // column c in permuted matrix
      if (beg<=c && c<end) ++pos[i-beg];
    }
  }
  unsigned *iAn = new unsigned[nsize+1];
  iAn[0] = 0;
  for (i=0; i<nsize; i++) {
    iAn[i+1] = iAn[i] + pos[i];
    pos[i] = iAn[i];
  }
  
  unsigned *jAn = new unsigned[iAn[nsize]];
  for (i=beg; i<end; i++) {
    const unsigned r = op_perm[i], idx = iA[r+1]; // row r in original matrix
    for (unsigned j=iA[r]; j<idx; j++) {
      const unsigned c = po_perm[jA[j]]; // column in permuted matrix
      if (beg<=c && c<end) jAn[pos[i-beg]++] = c-beg;
    }		
  }
  
  delete [] pos;

  for (i=0; i<nsize; ++i) quickSort(jAn, iAn[i], iAn[i+1]);
  return new AdjMat(nsize, iAn, jAn);
}


/*
AdjMat* AdjMat::getMat(unsigned beg, unsigned end, 
		       unsigned* op_perm, unsigned* po_perm,
		       unsigned &nadj_mats, AdjMat** &adj_mats) const
{
  unsigned s = end-beg, k;
  unsigned *comp = new unsigned[s];
  unsigned *reordering = new unsigned[s];
  
  // determine the number of components nadj_mats
  // and the nodes of the components in the field comp
  getComponents(beg, end, op_perm, po_perm, nadj_mats, comp);
  // init perm
  for (k=0; k<s; k++) reordering[k] = k;
  quickSort (comp, 0, s, reordering);
  // apply local reordering to op_perm and po_perm
  unsigned *t_op_perm = new unsigned [s];
  for (unsigned k=0; k<s; ++k) 
    t_op_perm[k] = op_perm[beg + reordering[k]];
  
  for (unsigned k=0; k<s; ++k) op_perm[beg+k] = t_op_perm[k];
  for (unsigned k=beg; k<end; ++k) po_perm[op_perm[k]] = k;
  delete [] t_op_perm;
  delete [] reordering;

  adj_mats = new AdjMat*[nadj_mats];
  	
  if (nadj_mats == 1) {
    adj_mats[0] = getMat (beg, end, op_perm, po_perm);
    delete [] comp;
    return adj_mats[0];
  }
  
  for (k=0; k<nadj_mats; k++)
    adj_mats[k] = getSubMat (beg, end, op_perm, po_perm, k+1, comp);
  
  delete [] comp;
  
  return adj_mats[0];
}
*/

 /*
AdjMat* AdjMat::getSubMat(unsigned beg, unsigned end, 
			  const unsigned* const op_perm,
			  const unsigned* const po_perm,
			  unsigned k, const unsigned* const comp) const
{
   unsigned s=0, l, size=end-beg, cnt=0;
   // count elements in k-th component
   for (l=0; l<size; l++) if (comp[l] == k) s++;
   
   unsigned *local_comp_idx = new unsigned [s];
   // write indices of k-th component to local_comp_idx
   for (l=0; l<size; l++) 
      if (comp[l] == k) // comp respects the permutation
         local_comp_idx[cnt++] = beg+l;
   
   unsigned i, c; // row and col idx in permuted matrix
   unsigned j, idx; // pointer in jA
   unsigned r; // row idx in original matrix
	
   unsigned *iAn = new unsigned [s+1];
   iAn[0] = 0;

   unsigned *pos = new unsigned [s+1];
   for (i=0; i<=s; i++) pos[i] = 0;
	
   for (i=0; i<s; i++) {
      r = op_perm[local_comp_idx[i]];
      idx = iA[r+1];
      for (j=iA[r]; j<idx; j++) {
         c = po_perm[jA[j]];
         if (beg<=c && c<end) ++pos[i];
      }
   }
   
   for (i=0; i<s; i++) iAn[i+1] = iAn[i] + pos[i];
   for (i=0; i<s; i++) pos[i] = iAn[i];

   unsigned *jAn = new unsigned[iAn[s]];
   for (i=0; i<s; i++) {
      r = op_perm[local_comp_idx[i]];
      idx = iA[r+1];
      for (j=iA[r]; j<idx; j++) {
         c = po_perm[jA[j]];
         if (beg <= c && c < end) {
         	unsigned p; // position of c inside local_comp_idx
         	if (binSearch (0, s, local_comp_idx, c, p)) jAn[pos[i]++] = p;
         }
      }		
   }

   delete [] pos;
   delete [] local_comp_idx;
   for (i=0; i<s; ++i) quickSort(jAn, iAn[i], iAn[i+1]);
   return new AdjMat(s, iAn, jAn, NULL); 	
}
 */

  /*
void AdjMat::addEdges(std::list<unsigned>* edges)
{
  unsigned nedges = edges->size() / 2;
  // *** make two copies of the list of edges ***
  unsigned *edges0s = new unsigned [2*nedges]; // start points of first copy
  unsigned *edges0e = new unsigned [2*nedges]; // end points of first copy
  unsigned k=0, j=0, i=0, r=0, c=0;

  std::list<unsigned>::iterator it = edges->begin();
  while (it != edges->end()) {
    edges0s[k] = (*it);
    edges0e[k+1] = (*it);
    ++it;
    ++k;
    edges0s[k] = (*it);
    edges0e[k-1] = (*it);
    ++it;
    k++;
  }
  // done

  // *** sorting edges with respect to the first component
  // and the second component ***
  specialSort(edges0s, 0, 2*nedges, edges0e);

  unsigned *iAn = new unsigned[n+1], *cnt = new unsigned[n];
  for (k=0; k<n; ++k) cnt[k] = iA[k+1]-iA[k];

  for (k=0; k<2*nedges; ++k) ++cnt[edges0s[k]];

  iAn[0] = 0;
  for (k=0; k<n; ++k) iAn[k+1] = iAn[k]+cnt[k];

  unsigned  *jAn = new unsigned[iAn[n]], *An = new unsigned[iAn[n]];

  for (r=0; r<n; ++r) {
    // consider j-th edge, which is aquivalent to row edges0s[j]
    if (r<edges0s[j]) {
      // transfer all rows until new edge
      // i ... position in fields jAn, An
      unsigned idx = iA[edges0s[j]];

      for (c=iA[r]; c<idx; ++c) {
        jAn[i] = jA[c];
        An[i] = A[c];
        i++;
      }
      // skip some rows
      r = edges0s[j]-1;
    }
    while (r == edges0s[j]) {
      unsigned idx = iA[r+1], idxn = iAn[r+1];
      // search for position of the new entry edges0e{j] in row r
      c = iA[r];
      while (c<idx && r == edges0s[j] && i<idxn) {
        if (jA[c] < edges0e[j]) {
          jAn[i] = jA[c];
          An[i] = A[c];
          c++;
          i++;
        } else {
          jAn[i] = edges0e[j];
          An[i] = 1; // new edge weight
          j++;
          i++;
        }
      }
      // new entries at the end of the row
      while (c==idx && r == edges0s[j] && i<idxn) {
        jAn[i] = edges0e[j];
        An[i] = 1; // new edge weight
        i++;
        j++;
      }
      // transmit remaining row
      while (c<idx  && i<idxn) {
        jAn[i] = jA[c];
        An[i] = A[c];
        i++;
        c++;
      }
    }
    if (j == 2*nedges && r<n) {
      // all edges inserted - transmit remaining matrix entries
      unsigned idx(iA[n]);
      for (c = iA[r+1]; c<idx; ++c) {
        jAn[i] = jA[c];
        An[i] = A[c];
        i++;
      }
      r = n;
    }
  }

  swap(A, An);
  swap(jA, jAn);
  swap (iA, iAn);

  delete [] An;
  delete [] jAn;
  delete [] iAn;
  delete [] cnt;

  delete [] edges0s;
  delete [] edges0e;
}

*/

void AdjMat::getComponents(unsigned beg, unsigned end,
                           const unsigned* const op_perm,
			   const unsigned* const po_perm,
                           unsigned &ncomp, unsigned* comp) const
{
  ncomp = 0; // number of components
  const unsigned s = end-beg;
  // init: all nodes belong to the null-component
  for (unsigned k=0; k<s; k++) comp[k] = ncomp;

  // counter, start idx, number of new / old neighbours
  unsigned cnt=0, start=0, nnnbrs=0, nonbrs=0;
  // fields for storing neighbour indices
  unsigned *nbrsf = new unsigned[2*s];
  unsigned *new_nbrs = nbrsf+s, *old_nbrs = nbrsf;

  while (cnt < s) {
    // search for first node not belonging to a component
    while (start<s && comp[start]!=0) start++;
    ncomp++; // determine next component
    comp[start] = ncomp;
    cnt++;

    old_nbrs[0] = start;
    nonbrs = 1;

    // determine the connected component ncomp using breadth first search
    while (cnt<s && nonbrs>0) {
      // go through all active neighbours
      for (unsigned j=0; j<nonbrs; ++j) {
        // transform to original index -> original row
        const unsigned r = op_perm[beg+old_nbrs[j]];
        // iterate over original row
        const unsigned idx1 = iA[r+1];
        for (unsigned k=iA[r]; k<idx1; k++) {
          // get column c and permuted column pc
          const unsigned pc = po_perm[jA[k]];
          // if permuted column pc is in the intervall (block)
          if (beg <= pc && pc < end) {
            // determine the component index
            const unsigned comp_idx = pc-beg;
            // has this index been checked?
            if (!comp[comp_idx]) {
              comp[comp_idx] = ncomp;
              new_nbrs[nnnbrs++] = comp_idx;
            }
          }
        }
      }
      swap(new_nbrs, old_nbrs);
      nonbrs = nnnbrs;
      cnt += nnnbrs;
      nnnbrs = 0;
    }
  }

  delete [] nbrsf;
}

