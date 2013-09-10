/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "Separator.h"
#include "AdjMatrix.h"
#include "metis.h"

#if METIS_VERSION != 5
extern "C" void METIS_PartGraphRecursive(const unsigned*, unsigned*, unsigned*,
 					 unsigned*, const unsigned*,
					 const unsigned*, const unsigned*,
 					 const unsigned*, const unsigned*,
 					 unsigned*, unsigned*);
#endif

void Separator::subdivide(unsigned bmin)
{
  if (size()>bmin) {
    const unsigned s = adj_mat->getRows(); // global matrix size
    if (s < size()*pow2(depth+2)) {
      
      // if cluster is connected then diam is lower than
      // std::limits<unsigned>::max()
      if (diam != std::numeric_limits<unsigned>::max()) {
	const unsigned nrow = _local_adj_mat->getRows();
	unsigned *reorder = new unsigned[nrow];
	unsigned *xadj = _local_adj_mat->getRowIdxArray();
	unsigned *adjncy = _local_adj_mat->getColIdxArray();
	const unsigned nparts = 2;
	unsigned edgecut;

#if METIS_VERSION == 5
	const unsigned ncon = 1;
	METIS_PartGraphRecursive((idx_t*)&nrow, (idx_t*)&ncon, (idx_t*)xadj,
				 (idx_t*)adjncy, NULL, NULL, NULL,
	 			 (idx_t*)&nparts, NULL, NULL,NULL,
				 (idx_t*)&edgecut, (idx_t*)reorder);
#else
	const unsigned wgtflag = 0, numflag = 0, options = 0;
	METIS_PartGraphRecursive(&nrow, xadj, adjncy, NULL, NULL, &wgtflag,
				 &numflag, &nparts, &options, &edgecut,
				 reorder);
#endif

	_l_op_perm = new unsigned[nend-nbeg];
	_l_po_perm = new unsigned[nend-nbeg];
	// init local permutations
	for (unsigned i=0; i<nend-nbeg; ++i) _l_op_perm[i] = i;
	
	unsigned isep = updatePerm(reorder) + nbeg;
	delete [] reorder;
	
	// invert local and global permutation
	for (unsigned k=nbeg; k<nend; ++k) {
	  po_perm[op_perm[k]] = k;
	  _l_po_perm[_l_op_perm[k-nbeg]] = k-nbeg;
	}
		
	// do next recursion step if size is equal the appropriate level
	if ((isep>=bmin+nbeg) && (nend>=bmin+isep)) {

	  AdjMat *l_adj0 = _local_adj_mat->getMat(0, isep-nbeg, _l_op_perm,
						  _l_po_perm);
	  AdjMat *l_adj1 = _local_adj_mat->getMat(isep-nbeg, nend-nbeg,
						  _l_op_perm, _l_po_perm);

	  delete _local_adj_mat;
	  _local_adj_mat = NULL;
	
	  delete [] _l_op_perm;
	  delete [] _l_po_perm;
	  _l_op_perm = _l_po_perm = NULL;

	  nsons = 2;
	  sons = (cluster**) new ClusterAlg*[nsons];
	  // constructing child nodes for index cluster tree
	  sons[0] = new Separator(this, nbeg, isep, op_perm, po_perm,
				  depth+1, adj_mat, l_adj0);
	  sons[1] = new Separator(this, isep, nend, op_perm, po_perm,
				  depth+1, adj_mat, l_adj1);
	  
	  // recursive subdivision
	  ((ClusterAlg*)sons[0])->subdivide(bmin);
	  ((ClusterAlg*)sons[1])->subdivide(bmin);	  
	} else {
	  delete _local_adj_mat;
	  _local_adj_mat = NULL;

	  delete [] _l_op_perm;
	  delete [] _l_po_perm;
	  _l_op_perm = _l_po_perm = NULL;
	}
      }
    } else {
      nsons = 1;
      sons = (cluster**) new ClusterAlg*[nsons];
      
      _l_op_perm = new unsigned[nend-nbeg];
      _l_po_perm = new unsigned[nend-nbeg];
      // init local permutations
      for (unsigned i=0; i<nend-nbeg; ++i) _l_op_perm[i] = _l_po_perm[i] = i;
      AdjMat *l_adj0 = _local_adj_mat->getMat(0, nend-nbeg, _l_op_perm,
					      _l_po_perm);
      delete _local_adj_mat;
      _local_adj_mat = NULL;

      delete [] _l_op_perm;
      delete [] _l_po_perm;
      _l_op_perm = _l_po_perm = NULL;

      sons[0] = new Separator(this, nbeg, nend, op_perm, po_perm, depth+1,
			      adj_mat, l_adj0);

      ((ClusterAlg*)sons[0])->setVirtualDepth(vdepth);
      ((ClusterAlg*)sons[0])->subdivide(bmin);
      ((ClusterAlg*)sons[0])->setDiam(diam);
    }
  }
}


unsigned Separator::updatePerm(unsigned* reordering)
{
  unsigned beg = 0, end = nend-nbeg;
  while (beg<end) {
    if (reordering[beg]==1) {
      --end;
      while (beg<end && reordering[end]==1) --end;
      // local permutation
      swap(_l_op_perm[beg], _l_op_perm[end]);
      // global permutation
      swap(op_perm[nbeg+beg], op_perm[nbeg+end]);
    }
    ++beg;
  }
  return ((beg>end) ? beg-1 : end);
}

void Separator::refine(unsigned max_level)
{
  if (depth<max_level) {
    if (nsons==0) {
      nsons = 1;
      sons = (cluster**) new ClusterAlg*[nsons];
      sons[0] = new Separator(this, nbeg, nend, op_perm, po_perm, depth+1,
			      adj_mat, NULL);
      ((ClusterAlg*)sons[0])->setVirtualDepth(vdepth);
      ((ClusterAlg*)sons[0])->setDiam(diam);
    }

    for (unsigned k=0; k<nsons; ++k)
      ((ClusterAlg*)sons[k])->refine(max_level);
  }
}

void Separator::computeNeighbor()
{
  if (!((ClusterAlg*)parent)->isSeparator()) {
    // case: separator arises from a cluster
    if (parent->getns() == 3) {
      const unsigned s = nend-nbeg;
      ClusterAlg* son0 = (ClusterAlg*)parent->getson(0);
      son0->addNeighborElement(Neighbor(this, s, s));
      addNeighborElement(Neighbor(son0, s, s));
      ClusterAlg* son1 = (ClusterAlg*)parent->getson(1);
      son1->addNeighborElement(Neighbor(this, s, s));
      addNeighborElement(Neighbor(son1, s, s));
    }
  } else {
    // case: separator arises from a separator by subdivision
    if (parent->getns() == 2) {
      // compute the neighbour size
      unsigned size_edge_cut = 0;
      if (parent->getson(0)->getnbeg() == nbeg) {
        unsigned nnodes0 = 0, nnodes1 = 0;
	ClusterAlg* son1 = (ClusterAlg*)parent->getson(1);
        size_edge_cut = isNeighbour(son1, nnodes0, nnodes1);
        if (size_edge_cut > 0) {
          son1->addNeighborElement(Neighbor(this, size_edge_cut, nnodes1));
          addNeighborElement(Neighbor(son1, size_edge_cut, nnodes0));
        }
      } else {
        unsigned nnodes0 = 0, nnodes1 = 0;
	ClusterAlg* son0 = (ClusterAlg*)parent->getson(0);
        size_edge_cut = isNeighbour(son0, nnodes0, nnodes1);
        if (size_edge_cut > 0) {
          son0->addNeighborElement(Neighbor(this, size_edge_cut, nnodes0));
          addNeighborElement(Neighbor(son0, size_edge_cut, nnodes1));
        }
      }
    }
  }

  // get children of the neighbors of father and test neighborhood
  const std::list<Neighbor> previous_neighbours = parent->getNeighbor();
  std::list<Neighbor>::const_iterator it = previous_neighbours.begin();
  while (it != previous_neighbours.end()) {
    const unsigned ns = (*it).ngbr->getns();
    for (unsigned k=0; k<ns; ++k) {
      unsigned nnodes0 = 0, nnodes1 = 0;
      ClusterAlg* pk = (ClusterAlg*) (it->ngbr)->getson(k);
      const unsigned size_edge_cut = isNeighbour(pk, nnodes0, nnodes1);
      if (size_edge_cut > 0) {
        addNeighborElement(Neighbor(pk, size_edge_cut, nnodes0));
        pk->addNeighborElement(Neighbor(this, size_edge_cut, nnodes1));
      }
    }
    ++it;
  }

  for (unsigned k=0; k<nsons; ++k) ((Separator*)sons[k])->computeNeighbor();
}

bool Separator::isadm(double eta2, cluster* cl, bl_info& info)
{
  if (cl == this) { // case: diagonal element
    info.is_adm = info.is_sep = false;
    return info.is_adm;
  }

  if (hasNeighbour((ClusterAlg*)cl)) {
    info.is_adm = info.is_sep = false;
    return info.is_adm;
  } else {
    info.is_sep = true;

#ifdef NO_DIST
    info.is_adm = true;
    return info.is_adm;
#else

    unsigned dist = 0;
    unsigned tmp_diam = max(diam, ((ClusterAlg*)(cl))->getDiam());

#if defined(VIRTUAL_DEPTH_WEIGHT) || defined(NEIGHBOUR_SIZE_WEIGHT) || defined(REFINEMENT)
    dist = getDist((ClusterAlg*)cl);
#endif

    if (tmp_diam < eta2 * dist) info.is_adm = true;
    else info.is_adm = false;

    return info.is_adm;
#endif
  }
}
