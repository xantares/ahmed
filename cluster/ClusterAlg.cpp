/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blas.h"
#include "basmod.h" // fuer funktion cputime
#include "ClusterAlg.h"
#include "../cluster/MGraph.h"
#include "AdjMatrix.h"
#include <fstream>


ClusterAlg::ClusterAlg(unsigned n, unsigned* iA, unsigned* jA) :
  cluster(0, n), parent(NULL), op_perm(NULL), po_perm(NULL), _l_op_perm(NULL), _l_po_perm(NULL), diam(0), depth(0)
{
  unsigned nnz = iA[n];
  
  // create adjacency matrix
  unsigned *row_ptr = new unsigned[n+1];
  for (unsigned k=0; k<=n; ++k) row_ptr[k] = iA[k];
  unsigned *col_idx  = new unsigned[nnz];
  for (unsigned k=0; k<nnz; ++k) col_idx[k] = jA[k];
  
  adj_mat = new AdjMat(n, row_ptr, col_idx);
  adj_mat->makeSymmetric();

  _local_adj_mat = NULL;
}

ClusterAlg::ClusterAlg(ClusterAlg *father, unsigned beg, unsigned end,
                       unsigned* op_p, unsigned* po_p,
                       unsigned d, AdjMat* global_mat, AdjMat* local_mat)
  : cluster(beg, end), parent(father), op_perm(op_p), po_perm(po_p),
    _l_op_perm(NULL), _l_po_perm(NULL), diam(0), depth(d), adj_mat(global_mat), _local_adj_mat(local_mat), neighbors(NULL), nneighbors(0), vdepth(depth)
{
  // compute an approximation of diameter
  if (_local_adj_mat) {
    if (_local_adj_mat->getNNZ() > 1) {
      _local_adj_mat->getDiam(diam);
    } else diam = 0;
  }
}


ClusterAlg::~ClusterAlg()
{
  if (parent == NULL) delete adj_mat;
  delete _local_adj_mat;
  delete [] neighbors;
}

unsigned ClusterAlg::isNeighbour(ClusterAlg* const other, unsigned &nnodes0,
				 unsigned &nnodes1) const
{
  nnodes0 = nnodes1 = 0;
  unsigned nbrs_size = 0;
  if (this == other) return nbrs_size;
  
  unsigned beg (other->getnbeg()), end (other->getnend());
#ifdef NEIGHBOUR_SIZE_WEIGHT
  bool *is_ngbr_0 (new bool [nend-nbeg]);
  bool *is_ngbr_1 (new bool [end-beg]);
  for (unsigned k(0); k<nend-nbeg; ++k) is_ngbr_0[k] = false;
  for (unsigned k(0); k<end-beg; ++k) is_ngbr_1[k] = false;
#endif

  for (unsigned k=beg; k<end; ++k) {
    unsigned orow (op_perm[k]);
    unsigned idx0(adj_mat->getBeginRow(orow)), idx1(adj_mat->getBeginRow(orow+1));
    for (unsigned j=idx0; j<idx1; ++j) {
      unsigned pc (po_perm[adj_mat->getCol(j)]);
      if (nbeg <= pc && pc < nend) {
        nbrs_size++;
#ifdef NEIGHBOUR_SIZE_WEIGHT
        is_ngbr_0[pc-nbeg] = true;
        is_ngbr_1[k-beg] = true;
#endif
#if defined(VIRTUAL_DEPTH_WEIGHT) || defined(UNWEIGHTED)
        nnodes0 = nnodes1 = 1;
        return nbrs_size;
#endif
      }
    }
  }
  
#ifdef NEIGHBOUR_SIZE_WEIGHT
  for (unsigned k(0); k<nend-nbeg; ++k) if (is_ngbr_0[k]) nnodes0++;
  for (unsigned k(0); k<end-beg; ++k) if (is_ngbr_1[k]) nnodes1++;
  delete [] is_ngbr_0;
  delete [] is_ngbr_1;
#endif

  return nbrs_size;
}

unsigned ClusterAlg::hasNeighbour(ClusterAlg* const other) const
{
  if (nneighbors == 0) return 0;
  
  // binary search into neighbors array
  bool nfound = true, interval_changed = true;
  unsigned key = other->getnbeg();
  unsigned beg = 0, end = nneighbors;
  unsigned idx((beg + end) / 2);
  if (key < (neighbors[0].ngbr)->getnbeg())
    return 0;
  if ((neighbors[nneighbors - 1].ngbr)->getnbeg() < key)
    return 0;
  
  while (beg < end && nfound && interval_changed) {
    if (key == (neighbors[idx].ngbr)->getnbeg())
      nfound = false;
    if (key < (neighbors[idx].ngbr)->getnbeg()) {
      end = idx;
      idx = (beg + end) / 2;
      if (idx == end)
	interval_changed = false;
    }
    if (key > (neighbors[idx].ngbr)->getnbeg()) {
      beg = idx;
      idx = (beg + end) / 2;
      if (idx == beg)
	interval_changed = false;
    }
  }
  
  if (nfound)
    return 0;
  else
    return neighbors[idx].nnodes;
}

void ClusterAlg::neighborsToArray()
{
  nneighbors = list_boundary.size();
  neighbors = new Neighbor[nneighbors];

  std::list<Neighbor>::const_iterator it(list_boundary.begin());
  unsigned idx = 0;
  while (it != list_boundary.end()) {
    neighbors[idx].ngbr = it->ngbr;
    neighbors[idx].s_edge_cut = it->s_edge_cut;
    neighbors[idx].nnodes = it->nnodes;
    idx++;
    it++;
  }
  list_boundary.clear();
  sortNeighbors(neighbors, 0, nneighbors);
  for (unsigned k=0; k<nsons; ++k) ((ClusterAlg*)getson(k))->neighborsToArray();
}

#ifndef NO_DIST
#ifdef REFINEMENT
unsigned ClusterAlg::getDist(ClusterAlg* other)
{
  // construct graph
  MGraph<ClusterAlg*, unsigned> *graph (constructGraph ( other ));

  // compute path(s)
  std::list<ClusterAlg*> path0, path1, path2;
  unsigned dist0 (graph->getPathsAndDist ( this, other, path0, path1, path2 ));
  unsigned dist1(std::numeric_limits<unsigned>::max()), dist2(dist1), tmp_dist;
  delete graph;

  // *** first level refinement
  // compute refined graph for path0
  unsigned max_depth (depth), refine_level(3);
  getMaxDepth(max_depth);
  if (max_depth - depth < refine_level) refine_level = max_depth - depth;
  if (refine_level != 0) {
    refine_level += depth;
    graph = constructRefinedGraph (refine_level, path0);

    // choose a new src and target node in the new graph
    // starting from old src and target
    ClusterAlg*src (*path0.begin()), *target (*(--path0.end()));
    getNode (refine_level, src);
    other->getNode (refine_level, target);

    // compute distance in refined graph
    std::list<ClusterAlg*> path00, path01, path02;
    tmp_dist = graph->getPathsAndDist ( src, target, path00, path01, path02 );

    if (tmp_dist < dist1) dist1 = tmp_dist;

    // tidy up
    if (graph) delete graph;

    // *** second level refinement
    // *** refinement path p00
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3) {
      refine_level += depth;
      graph = constructRefinedGraph(refine_level, path00);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path00.begin();
      target = *(--path00.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path000, path001, path002;
      tmp_dist = graph->getPathsAndDist ( src, target, path000, path001, path002 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }
    // *** refinement path p01
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3 && path01.size() != 0) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path01);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path01.begin();
      target = *(--path01.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path010, path011, path012;
      tmp_dist = graph->getPathsAndDist ( src, target, path010, path011, path012 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }

    // *** refinement path p02
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3 && path02.size() != 0) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path02);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path02.begin();
      target = *(--path02.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path020, path021, path022;
      tmp_dist = graph->getPathsAndDist ( src, target, path020, path021, path022 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }
  }

  // *** first level refinement
  // compute refined graph for path1
  refine_level=3;
  if (max_depth - depth < refine_level) refine_level = max_depth - depth;
  if (refine_level != 0 && path1.size() > 0) {
    refine_level += depth;
    graph = constructRefinedGraph (refine_level, path1);

    // choose a new src and target node in the new graph
    // starting from old src and target
    ClusterAlg *src (*path1.begin()), *target (*(--path1.end()));
    other->getNode (refine_level, target);
    getNode (refine_level, src);

    // std::cout << "\033[32m |p1| " << dist << " \033[0m" << std::flush;
    // compute distance in refined graph
    std::list<ClusterAlg*> path10, path11, path12;
    tmp_dist = graph->getPathsAndDist ( src, target, path10, path11, path12 );

    if (tmp_dist < dist1) dist1 = tmp_dist;

    // tidy up
    if (graph) delete graph;

    // *** second level refinement
    // *** refinement path p10
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path10);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path10.begin();
      target = *(--path10.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path100, path101, path102;
      tmp_dist = graph->getPathsAndDist ( src, target, path100, path101, path102 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }
    // *** refinement path p11
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3 && path11.size() != 0) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path11);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path11.begin();
      target = *(--path11.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path110, path111, path112;
      tmp_dist = graph->getPathsAndDist ( src, target, path110, path111, path112 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }

    // *** refinement path p12
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3 && path12.size() != 0) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path12);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path12.begin();
      target = *(--path12.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path120, path121, path122;
      tmp_dist = graph->getPathsAndDist ( src, target, path120, path121, path122 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }
  }
  // *** first level refinement
  // compute refined graph for path2
  refine_level=3;
  if (max_depth - depth < refine_level) refine_level = max_depth - depth;
  if (refine_level != 0 && path2.size() > 0) {
    refine_level += depth;
    graph = constructRefinedGraph (refine_level, path2);

    // choose a new src and target node in the new graph
    // starting from old src and target
    ClusterAlg *src (*path2.begin()), *target (*(--path2.end()));
    other->getNode (refine_level, target);
    getNode (refine_level, src);

    // compute distance in refined graph
    std::list<ClusterAlg*> path20, path21, path22;
    tmp_dist = graph->getPathsAndDist ( src, target, path20, path21, path22 );

    if (tmp_dist < dist1) dist1 = tmp_dist;

    // tidy up
    if (graph) delete graph;

    // *** second level refinement
    // *** refinement path p20
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path20);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path20.begin();
      target = *(--path20.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path200, path201, path202;
      tmp_dist = graph->getPathsAndDist ( src, target, path200, path201, path202 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }
    // *** refinement path p21
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3 && path21.size() != 0) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path21);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path21.begin();
      target = *(--path21.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path210, path211, path212;
      tmp_dist = graph->getPathsAndDist ( src, target, path210, path211, path212 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }

    // *** refinement path p22
    refine_level=6;
    if (max_depth - depth < refine_level) refine_level = max_depth - depth;
    if (refine_level != 0 && refine_level > 3 && path22.size() != 0) {
      refine_level += depth;
      graph = constructRefinedGraph (refine_level, path22);

      // choose a new src and target node in the new graph
      // starting from old src and target
      src = *path22.begin();
      target = *(--path22.end());
      getNode (refine_level, target);
      other->getNode (refine_level, src);

      // compute distance in refined graph
      std::list<ClusterAlg*> path220, path221, path222;
      tmp_dist = graph->getPathsAndDist ( src, target, path220, path221, path222 );

      if (tmp_dist < dist2) dist2 = tmp_dist;

      // tidy up
      if (graph) delete graph;
    }
  }

  if (dist2 < std::numeric_limits<unsigned>::max())
    return dist2;
  if (dist1 < std::numeric_limits<unsigned>::max())
    return dist1;

  return dist0;
}
#endif

#endif

#ifndef NO_DIST
#ifndef REFINEMENT
unsigned ClusterAlg::getDist(ClusterAlg* other)
{
  // construct graph
  MGraph<ClusterAlg*, unsigned> *graph (constructGraph ( other ));

  // calc dist
  if (graph) {
    unsigned dist;
    if (graph->isGraphConnected ()) {
      dist = graph->getDist( this, other );
    } else {
      delete graph;
      graph = constructGraph (other, 1);

      if (graph) {
        if (graph->isGraphConnected ()) dist = graph->getDist( this, other );
        else {
          delete graph;
          graph = constructGraph (other, 2);
          if (graph) {
            if (graph->isGraphConnected ()) dist = graph->getDist( this, other );
            else {
              delete graph;

              graph = constructGraph (other, 3);

              if (graph) {
                if (graph->isGraphConnected ()) dist = graph->getDist( this, other );
                else {
                  delete graph;
                  graph = constructGraph (other, 4);
                  if (graph) {
                    if (graph->isGraphConnected ()) dist = graph->getDist( this, other );
                    else {
                      delete graph;
                      graph = NULL;
                      dist = std::numeric_limits<unsigned>::max();
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // tidy up
    if (graph) delete graph;
    return dist;
  }
  return 0;
}
#endif
#endif

MGraph<ClusterAlg*, unsigned>* ClusterAlg::constructGraph(ClusterAlg* other, unsigned up)
{
  std::list<ClusterAlg*> list_graph_nodes;

  ClusterAlg *p0 = this, *p1 = other;
  // search for roots of subtrees in the cluster tree
  findRoots (p0, p1);

  if (p0 == NULL || p1 == NULL) return NULL;
  if (p0->getDepth() < up || p1->getDepth() < up) return NULL;
  for (unsigned l=0; l<up; ++l) {
    p0 = (ClusterAlg*)p0->getParent();
    p1 = (ClusterAlg*)p1->getParent();
  }

  // fill list_graph_nodes with the appropriate graph nodes
  if (p0 == p1) p1->getNodes (depth, list_graph_nodes);
  else {
    p0->getNodes (depth, list_graph_nodes);
    p1->getNodes (depth, list_graph_nodes);
  }

  // create and fill vertex array for graph (adjacency matrix)
  unsigned row_size = list_graph_nodes.size(), idx = 0;
  ClusterAlg **graph_nodes = new ClusterAlg* [row_size];

  std::list<ClusterAlg*>::iterator node_it = list_graph_nodes.begin ();
  while (node_it != list_graph_nodes.end()) {
    graph_nodes[idx++] = *node_it;
    node_it++;
  }

  // count nnz-elements, fill liA, ljA and lA
  std::list<unsigned> liA, ljA, lA;
  liA.push_back (0);
  unsigned col = 0, col_cnt;
  std::list<ClusterAlg*>::iterator it0(list_graph_nodes.begin ()), it1;
  while (it0 != list_graph_nodes.end()) {
    it1 = list_graph_nodes.begin ();
    col_cnt = 0;
    while (it1 != list_graph_nodes.end()) {
      if (it0 != it1) {
        unsigned nbr_size ((*it0)->hasNeighbour(*it1));
#ifdef VIRTUAL_DEPTH_WEIGHT
        if (nbr_size > 0) {
          ljA.push_back (col_cnt);
          col++;
          // determine weights
          if ((*it0)->isSeparator() && (*it1)->isSeparator()) {
            unsigned diff0 ((*it0)->getDepth()-(*it0)->getVirtualDepth()+1);
            unsigned diff1 ((*it1)->getDepth()-(*it1)->getVirtualDepth()+1);
            if (diff0 < diff1) lA.push_back (diff0);
            else lA.push_back(diff1);
          } else {
            if ((*it0)->isSeparator() && !(*it1)->isSeparator())
              lA.push_back ((*it0)->getDepth()-(*it0)->getVirtualDepth()+1);
            else lA.push_back ((*it1)->getDepth()-(*it1)->getVirtualDepth()+1);
          }
        }
#endif

#ifdef NEIGHBOUR_SIZE_WEIGHT
        if (nbr_size) {
          ljA.push_back (col_cnt);
          col++;
          // determine weights
          if ((*it0)->isSeparator() && (*it1)->isSeparator()) {
            unsigned mnbr_size (nbr_size > (*it1)->hasNeighbour(*it0) ? nbr_size : (*it1)->hasNeighbour(*it0));
            unsigned min_diam ((*it0)->getDiam() < (*it1)->getDiam() ? (*it0)->getDiam() : (*it1)->getDiam());
            lA.push_back (ceil((double)min_diam / mnbr_size));
          } else {
            if ((*it0)->isSeparator() && !(*it1)->isSeparator()) {
              if (ceil((double)(*it0)->getDiam()/nbr_size) == 1)
                lA.push_back (ceil(sqrt((double)(*it0)->getDiam())));
              else lA.push_back (ceil((double)(*it0)->getDiam()/nbr_size));
            } else {
              if (ceil((double)(*it1)->getDiam()/nbr_size) == 1)
                lA.push_back (ceil(sqrt((double)(*it1)->getDiam())));
              else lA.push_back (ceil((double)(*it1)->getDiam()/nbr_size));
            }
          }
        }
#endif

#ifdef UNWEIGHTED
        if (nbr_size > 0) {
          ljA.push_back (col_cnt);
          lA.push_back (1);
          col++;
        }
#endif
      }
      it1++;
      col_cnt++;
    }
    liA.push_back (col);
    it0++;
  }

  // create AdjMat in CRS-format: iA, jA, A
  unsigned *iA = new unsigned[row_size+1];
  unsigned col_ptr_size = ljA.size();
  unsigned *jA = new unsigned[col_ptr_size];
  unsigned  *A = new unsigned[col_ptr_size];

  // fill iA
  std::list<unsigned>::const_iterator it2 = liA.begin ();
  idx = 0;
  while (it2 != liA.end()) {
    iA[idx++] = *it2;
    it2++;
  }
  iA[row_size] = col_ptr_size;

  // fill jA and A
  it2 = ljA.begin ();
  std::list<unsigned>::const_iterator it3 (lA.begin());
  idx = 0;
  while (it2 != ljA.end()) {
    jA[idx] = *it2;
    A[idx] = *it3;
    it2++;
    it3++;
    idx++;
  }

  AdjMat *tmp_adj(new AdjMat (row_size, iA, jA, A));

  MGraph <ClusterAlg*, unsigned> *graph
  (new MGraph<ClusterAlg*, unsigned> (graph_nodes, tmp_adj));

  delete [] graph_nodes;
  return graph;
}

void ClusterAlg::findRoots(ClusterAlg* &p0, ClusterAlg* &p1) const
{
  bool n_found (true);
  while ((p0 != NULL || p1 != NULL) && p0 != p1 && n_found) {
    if (p0->hasNeighbour (p1)) n_found = false;
    else {
      p0 = (ClusterAlg*)(p0->getParent());
      p1 = (ClusterAlg*)(p1->getParent());
    }
  }
}

void ClusterAlg::getNodes(unsigned level, std::list<ClusterAlg*> &nodes )
{
  if (level == depth) nodes.push_back (this);
  else
    for (unsigned k=0; k<nsons; ++k)
      ((ClusterAlg*)sons[k])->getNodes (level, nodes);
}

unsigned partitionNeighbors(Neighbor* neighbors, unsigned beg, unsigned end)
{
  unsigned i = beg+1, j = end-1;
  unsigned m = (neighbors[beg].ngbr)->getnbeg();

  for (;;) {
    while ((i<end) && ((neighbors[i].ngbr)->getnbeg() < m)) i++;
    while ((j>beg) && !((neighbors[j].ngbr)->getnbeg() < m)) j--;

    if (i >= j) break;
    swap(neighbors[i], neighbors[j]);
  }

  swap(neighbors[beg], neighbors[j]);
  return j;
}

void sortNeighbors(Neighbor* neighbors, unsigned beg, unsigned end)
{
  if (beg < end) {
    unsigned p = partitionNeighbors (neighbors, beg, end);
    sortNeighbors(neighbors, beg, p);
    sortNeighbors(neighbors, p+1, end);
  }
}
