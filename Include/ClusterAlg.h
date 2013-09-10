/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file  ClusterAlg.h
    \brief Include file for the class ClusterAlg.
    \author M. Bebendorf and T. Fischer
    \date   5/4/2007
*/
#ifndef CLUSTERALG_H
#define CLUSTERALG_H

#include "cluster.h"
#include <list>

#define NO_DIST   // isadm does not compute distance - it is based only on neighbors
#ifdef NO_DIST
#define UNWEIGHTED
#endif

#define VIRTUAL_DEPTH_WEIGHT
//#define NEIGHBOUR_SIZE_WEIGHT

//#define REFINEMENT
#ifdef REFINEMENT
#define UNWEIGHTED
#endif

class ClusterAlg;
class AdjMat;
template <class NodeDataType, class EdgeWeightType> class MGraph;

struct Neighbor {
  Neighbor(ClusterAlg* neighbour, unsigned sec, unsigned nn)
      : ngbr(neighbour), s_edge_cut(sec), nnodes(nn) {};

  Neighbor() : ngbr(NULL), s_edge_cut(0), nnodes(0) {};

  ClusterAlg* ngbr;
  unsigned s_edge_cut;
  unsigned nnodes;
};

/** \brief Base class for storing clusters of degrees of freedom without
 geometric information.
*/

class ClusterAlg : public cluster
{
public:
  /**
   * Constructor creates the root of the cluster tree
   * @param n
   * @param iA
   * @param jA
   * @return
   */
  ClusterAlg(unsigned n, unsigned* iA, unsigned* jA);
  /*!
    \brief Constructor
    \param father parent node in cluster tree
    \param beg beginning index of the cluster (see class DOFPoint)
    \param end beginning index of the next cluster (see class DOFPoint)
    \param op_perm
    \param po_perm
    \param depth distance between this node and the root
    \param global_mat reference to the global adjacency matrix of the matrix graph in crs format
    \param local_mat pointer to the local adjacency matrix of the matrix graph in crs format
    format
  */
  ClusterAlg(ClusterAlg* father, unsigned beg, unsigned end,
             unsigned* op_perm, unsigned* po_perm, unsigned depth,
             AdjMat* global_mat, AdjMat* local_mat);

  /** \brief Destructor.
   * Destructor frees all form the objects allocated memory.
   * */
  virtual ~ClusterAlg();

  //
  /**
   * only for compatibility. empty implementation.
   * @param
   * @param
   * @param
   */
  // virtual void subdivide(unsigned int, unsigned*, unsigned*, unsigned int &) {};


  /**
    \brief subdivide this cluster recursively, needed to generate the
    cluster tree
    \param bmin the minimal size of clusters (input)
  */
  virtual void subdivide(unsigned bmin) = 0;

  /**
   * Method creates recursively the cluster tree, i.e. changes the permutation op_perm and po_perm
   * and create child cluster trees. For this task only the adjacency matrix is used.
   * @param op_perm permutation: permutated_row = op_perm[original_row]
   * @param po_perm reverse permutation: original_row = po_perm[permutated_idx]
   * @param bmin threshold value for stopping further refinement
   * @return a cluster tree
   */
  virtual void createClusterTree(unsigned bmin, unsigned* op_perm,
				 unsigned* po_perm) { }

protected:
  /** \brief Method returns the pointer to the parent cluster.
   \returns parent cluster */
  ClusterAlg* getParent() const { return parent; }

  /** hasNeighbour iterates through the boundary list searching for pointer to other
  \param other the ClusterAlg to check with
  \returns the size of neighbour
  */
  unsigned hasNeighbour(ClusterAlg* const other) const;

public:
  /** Method traverses the ClusterAlg-tree and computes the boundary
   clusters of the actual cluster.
   */
  virtual void computeNeighbor() = 0;

  /** \brief method checks if there exists edges in \f$(i, j) \in E(G(A))\f$
  with \f$i\f$ in this cluster and \f$j\f$ in cluster other.
  \param input: other ClusterAlg object, in particular the represented index set
  \param nnodes0 output: number of nodes in this cluster adjacent with any node into the
       other cluster
  \param nnodes1 output: number of nodes in other cluster adjacent with any node into this
       cluster
  \returns the number of edges between the clusters
  */
  unsigned isNeighbour(ClusterAlg* const other, unsigned &nnodes0,
                       unsigned &nnodes1) const;

  /** transforms the list of neighbors into a sorted array of neighbors */
  void neighborsToArray();

  /** \brief Method returns the list of ClusterAlg* that are in the
  neighborhood of the actual ClusterAlg.
  \returns the boundary of the actual cluster
  */
  const std::list<Neighbor>& getNeighbor() const {
    return list_boundary;
  }


  /** Method adds Cluster or Separator to the boundary of the actual cluster.
  */
  void addNeighborElement (Neighbor b) { list_boundary.push_back(b); }

  /**
   * \brief This method checks the admissibility of cluster this with
   * the cluster cl.
   * */
  virtual bool isadm(double eta2, cluster* cl, bl_info& info) = 0;

  /**
   * Method returns the status of this ClusterAlg object.
   * */
  virtual bool isSeparator() const = 0;

  /**
   * The attribute depth takes the position of this cluster into the
   * cluster-tree.
   * \returns the value of depth
   * */
  unsigned getDepth() const { return depth; }

  /**@name Public Member Functions - Graph related methods for
  	admissibility checking

  The following methods base on graph-methods. The graph is generated with
  constructGraph from nodes consisting of cluster-tree-entities and
  edges between clusters if they are neighbored.
  */
  //@{
  /** \brief Method computes an estimation of the distance between the cluster
  represented by this object and the cluster other based on graphs for the
  admissibility check.

  In order to compute the estimation getDist() sets up a graph and
  uses the algorithm of Dijkstra to compute the distance, whereas
  the graph is build by method constructGraph().
  \param other cluster object to whom distance should be computed
  \returns distance in the graph
  */

#ifndef NO_DIST
  unsigned getDist(ClusterAlg* other);
#endif

  /** \brief Method constructGraph constructs the graph $G_D$.
     For this purpose two neighbored predecessors $p_0, p_1 \in T_I$ are searched
     using findRoots(). All clusters rooted at $p_0, p_1$ in the actual level are collected
     in a list by getNodes(). These are the vertices of $G_D$. Two vertices are neighbored iff
     the clusters are neighbored. The weight of an edge is difference between the actual level
     and the virtual level.
     \param other an other cluster for the block $b = this \times t_2$
     \return a simple weighted graph
  */
  MGraph<ClusterAlg*, unsigned>* constructGraph(ClusterAlg* t_2, unsigned up = 0);

#ifdef REFINEMENT
  MGraph<ClusterAlg*, unsigned>* constructRefinedGraph(unsigned refine_level, std::list<ClusterAlg*> path);
  void getNode(unsigned refine_level, ClusterAlg* &node);
#endif

  /** \brief Method findRoots returns first neighbored pair of predecessors in the cluster tree.
     \param $p_0$ input: the starting cluster, output: the predecessor neighbored
        with a predecessor of $p_1$
     \param $p_1$ input: the starting cluster, output: the predecessor neighbored
        with a predecessor of $p_0$
  */
  void findRoots(ClusterAlg* &p_0, ClusterAlg* &p_1) const;

  /** \brief Method traverses recursively the cluster-tree until a given
  	level. If this level is achieved the pointer to this cluster is inserted
  	into the list.
  	\param nodes the list of pointers to the appropriate clusters (output)
  */
  void getNodes(unsigned level, std::list<ClusterAlg*> &nodes);

  /**
   * \brief method returns an estimation of the diameter of the cluster.
   * for the estimation holds:
   * radius of graph <= estimated diam <= diameter of graph
   * the value is computed in the constructor
  */
  unsigned getDiam() const { return diam; }
  void setDiam(unsigned d) { diam = d; }

  unsigned getVirtualDepth() const { return vdepth; }
  void setVirtualDepth(unsigned vd) { vdepth = vd; }

  virtual void refine(unsigned max_level) = 0;

protected: // *** attributes of the class ********* //
  /** pointer to parent */
  ClusterAlg *parent;
  /** global permutation: original <- op_perm <- permutation */
  unsigned* op_perm;
  /** global permutation: permutation <- po_perm <- original */
  unsigned* po_perm;
  /** local permutation: original <- l_op_perm <- permutation */
  unsigned* _l_op_perm;
  /** local permutation: permutation <- l_po_perm <- original */
  unsigned* _l_po_perm;
  /** approximation for diameter of graph */
  unsigned diam;
  /** distance between the root node and this cluster */
  unsigned depth;
  /** The AdjMat stores the edge set of the matrix graph $G = (V,E)$.
   (see class AdjMat) */
  AdjMat* adj_mat;
  /**
   * local adjacency matrix
   */
  AdjMat* _local_adj_mat;

  /**
   * list of clusters that are in neighborhood
   * */
  std::list<Neighbor> list_boundary;
  /**
   * array of neighbors
   */
  Neighbor* neighbors;
  /**
   * number of neighbors
   */
  unsigned nneighbors;

  /** level since the index set is persistent */
  unsigned vdepth;
};

/** helper functions for sorting the neighbors in order to accelerate the
  neighbor search */
unsigned partitionNeighbors (Neighbor* neighbors, unsigned beg, unsigned end);
void sortNeighbors ( Neighbor* neighbor, unsigned beg, unsigned end );

#endif
