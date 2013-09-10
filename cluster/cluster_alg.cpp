/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <limits>
#include "blas.h"
#include "cluster_alg.h"
#include "Separator.h"
#include "AdjMatrix.h"
#include "metis.h"

// METIS function
#if METIS_VERSION != 5
extern "C" void METIS_NodeComputeSeparator(const unsigned*, unsigned*,
					   unsigned*, unsigned*, unsigned*,
					   const int*, const unsigned*,
					   unsigned*);
#endif

void cluster_alg::subdivide(unsigned bmin)
{
  // if cluster is a connected component then diam is lower than
  // std::numeric_limits<unsigned>::max()
  if (diam!=std::numeric_limits<unsigned>::max() && size()>bmin) {

    // for METIS    
    const unsigned nrow = _local_adj_mat->getRows();
    unsigned *xadj = _local_adj_mat->getRowIdxArray();
    unsigned *adjncy = _local_adj_mat->getColIdxArray();
    const unsigned sepsize = 0, nnz = xadj[nrow];
    const unsigned maxn = MAX(nnz, nrow+1);
    unsigned *wgt = new unsigned[maxn];
    for (unsigned k=0; k<maxn; k++) wgt[k] = 1;
    unsigned *part = new unsigned[nrow+1];

    // subdivide the index set into three parts employing METIS
#if METIS_VERSION == 5
    METIS_ComputeVertexSeparator((idx_t*)&nrow, (idx_t*)xadj, (idx_t*)adjncy,
				 (idx_t*)wgt, NULL, (idx_t*)&sepsize,
				 (idx_t*)part);
#else
    int options[8];
    options[0] = 0; // use the default values
    options[7] = -1; // initialize the random number generator
    METIS_NodeComputeSeparator(&nrow, xadj, adjncy, wgt, wgt, options,
			       &sepsize, part);
#endif

    delete [] wgt;

    // create and init local permutations
    _l_op_perm = new unsigned[nend-nbeg];
    for (unsigned i=0; i<nend-nbeg; ++i) _l_op_perm[i] = i;
    _l_po_perm = new unsigned[nend-nbeg];

    unsigned isep1, isep2;
    updatePerm(part, isep1, isep2);
    delete [] part;

    // next recursion step
    if ((isep1>=bmin) && (isep2>=bmin+isep1)) {

      // construct adj matrices for [0, isep1), [isep1,isep2), [isep2, nend)
      AdjMat *l_adj0 = _local_adj_mat->getMat(0, isep1, _l_op_perm, _l_po_perm);
      AdjMat *l_adj1 = _local_adj_mat->getMat(isep1, isep2, _l_op_perm,
					      _l_po_perm);
      AdjMat *l_adj2 = _local_adj_mat->getMat(isep2, nend-nbeg, _l_op_perm,
					      _l_po_perm);

      delete _local_adj_mat;
      _local_adj_mat = NULL;
	
      delete [] _l_op_perm;
      delete [] _l_po_perm;
      _l_op_perm = _l_po_perm = NULL;

      isep1 += nbeg;
      isep2 += nbeg;
      
      // if isp2==nend then the separator is empty
      if (isep2==nend) std::cout << "empty separator." << std::endl;
      nsons = (isep2==nend) ? 2 : 3;
      sons = (cluster**) new ClusterAlg*[nsons];
      
      // constructing child nodes for index cluster tree
      sons[0] = (cluster*) new cluster_alg(this, nbeg, isep1, op_perm,
					   po_perm, depth+1, adj_mat, l_adj0);
      ((ClusterAlg*) sons[0])->setVirtualDepth(depth+1);

      sons[1] = (cluster*) new cluster_alg(this, isep1, isep2, op_perm,
					   po_perm, depth+1, adj_mat, l_adj1);
      ((ClusterAlg*) sons[1])->setVirtualDepth(depth+1);

      if (isep2!=nend) {
	sons[2] = (cluster*) new Separator(this, isep2, nend, op_perm,
					   po_perm, depth+1, adj_mat, l_adj2);
	((ClusterAlg*) sons[2])->setVirtualDepth(depth+1);
      }

      // continue recursion
      for (unsigned k=0; k<nsons; k++)
	((ClusterAlg*) sons[k])->subdivide(bmin);

    } else {
      delete _local_adj_mat;
      _local_adj_mat = NULL;

      delete[] _l_op_perm;
      delete[] _l_po_perm;
      _l_op_perm = _l_po_perm = NULL;
    } // end if next recursion step
  } // end if ( connected && size () > bmin )
}

void cluster_alg::updatePerm(unsigned* reordering, unsigned &isep0,
			     unsigned &isep1)
{
  unsigned beg = 0, end = nend-nbeg;
  while (beg<end) {
    if (reordering[beg]>=1) {
      --end;
      while (beg<end && reordering[end]>=1) --end;

      swap(reordering[beg], reordering[end]);
      // local permutation
      swap(_l_op_perm[beg], _l_op_perm[end]);
      // global permutation
      swap(op_perm[nbeg+beg], op_perm[nbeg+end]);
    }
    ++beg;
  }

  isep0 = (beg>end) ? beg-1 : end;

  beg = isep0, end = nend - nbeg;
  while (beg<end) {
    if (reordering[beg] == 2) {
      --end;
      while (beg<end && reordering[end]==2) --end;

      swap(reordering[beg], reordering[end]);
      // local permutation
      swap(_l_op_perm[beg], _l_op_perm[end]);
      // global permutation
      swap(op_perm[nbeg+beg], op_perm[nbeg+end]);
    }
    ++beg;
  }

  isep1 = (beg>end) ? beg-1 : end;

  // invert local and global permutation
  for (unsigned k=nbeg; k<nend; ++k) {
    po_perm[op_perm[k]] = k;
    _l_po_perm[_l_op_perm[k-nbeg]] = k-nbeg;
  }
}

bool cluster_alg::isadm(double eta2, cluster* cl, bl_info& info)
{
  // case: diagonal element
  if (cl == this) {
    info.is_adm = info.is_sep = false;
    return info.is_adm;
  }

  // this is a cluster by definition, check whether cl is a cluster
  // then nested dissection reordering makes a zero block
  if (!((ClusterAlg*) cl)->isSeparator()) {
    info.is_sep = info.is_adm = true;
    return info.is_adm;
  }

  // cl is a separator, check adm for separators
  return cl->isadm(eta2, this, info);
}

void cluster_alg::computeNeighbor()
{
  for (unsigned k=0; k<nsons; k++) ((ClusterAlg*)sons[k])->computeNeighbor();
}

void cluster_alg::refine(unsigned max_level)
{
  if (depth < max_level) {
    if (nsons == 0) {
      nsons = 1;
      sons = (cluster**) new ClusterAlg*[nsons];
      sons[0] = new cluster_alg(this, nbeg, nend, op_perm, po_perm, depth+1,
				adj_mat, NULL);
      ((ClusterAlg*)sons[0])->setVirtualDepth(vdepth);
      ((ClusterAlg*)sons[0])->setDiam(diam);
    }
    for (unsigned k=0; k<nsons; ++k)
      ((ClusterAlg*)sons[k])->refine(max_level);
  }
}

void cluster_alg::createClusterTree(unsigned bmin, unsigned* aop_perm,
				    unsigned* apo_perm)
{
  op_perm = aop_perm;
  po_perm = apo_perm;
  // create top local problem
  unsigned n = adj_mat->getRows(), nnz = adj_mat->getNNZ();
  unsigned *iAn = new unsigned[n+1];
  for (unsigned i=0; i<=n; ++i) iAn[i] = adj_mat->getRowIdxArray()[i];
  unsigned *jAn = new unsigned[nnz];
  for (unsigned i=0; i<nnz; ++i) jAn[i] = adj_mat->getColIdxArray()[i];
  
  _local_adj_mat = new AdjMat(n, iAn, jAn);
    
  // create cluster tree
  double time = cputime(0.0); // debug
  subdivide(bmin);
  time = cputime(time); // debug
  std::cout << "subdivide: " << time << "s." << std::endl; // debug

  // refine
  unsigned max_depth = getMaxDepth();
  refine(max_depth);

  // compute boundaries
  computeNeighbor();

  // transform neighbors (list -> sorted array)
  neighborsToArray();
}
