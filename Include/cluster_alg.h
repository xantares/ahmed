/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file  cluster_alg.h
    \brief Include file for the class cluster_alg.
    \author M. Bebendorf and T. Fischer
    \date   2/23/2007
*/
#ifndef CLUSTER_NON_GEOMETRIC_H
#define CLUSTER_NON_GEOMETRIC_H

#include "ClusterAlg.h"

class Separator;

/** \brief class for storing clusters of degrees of freedom without
    geometric information. This class stores clusters (no separators)
*/
class cluster_alg : public ClusterAlg
{
public:
  /**
   * Constructor creates the root of the cluster tree
   * @param n
   * @param jA
   * @param iA
   * @return
   */
  cluster_alg(unsigned n, unsigned* iA, unsigned* jA) : ClusterAlg (n, iA, jA) {}

protected:
  /** \brief Constructor
      \param father parent node in cluster tree
      \param beg beginning index of the cluster
      \param end beginning index of the next cluster
      \param op_perm permutation
      \param po_perm permutation
      \param depth distance between this node and the root
      \param global_mat reference to adjacency matrix of the matrix graph in
      crs format
      \param local_mat pointer to the local adjacency matrix of the matrix
      graph in crs format
  */
  cluster_alg(ClusterAlg* father, unsigned beg, unsigned end,
              unsigned* op_perm, unsigned* po_perm, unsigned depth,
              AdjMat* global_mat, AdjMat* local_mat)
      : ClusterAlg(father, beg, end, op_perm, po_perm, depth, global_mat, local_mat) {}

public:
  /** Destructor. */
  virtual ~cluster_alg() {}

  /**
   * Method creates recursively the cluster tree, i.e. changes the permutation
   * op_perm and po_perm and create child cluster trees. For this task only the
   * adjacency matrix is used.
   * @param op_perm permutation: permutated_row = op_perm[original_row]
   * @param po_perm reverse permutation: original_row = po_perm[permutated_idx]
   * @param bmin threshold value for stopping further refinement
   * @return a cluster tree
   */
  virtual void createClusterTree(unsigned bmin, unsigned* op_perm,
				 unsigned* po_perm);

protected:
  /** update perm */
  void updatePerm(unsigned* reordering, unsigned &isep0, unsigned &isep1);

  virtual void refine(unsigned max_level);

  /** Method traverses the ClusterAlg tree and computes the boundary
      clusters of the actual cluster that consists of the sibling-separator
      and parts of other separators from higher levels.
  */
  virtual void computeNeighbor();

  void subdivide(unsigned bmin);

public:
  /** \brief This method checks the admissibility of the cluster represented
      by this object and the cluster cl given as parameter to isadm.

      If both clusters are not separators then isadm returns true and the
      parameter info.is_adm is set true, too.
      (because of nested dissection reordering)

      If cluster cl is a separator, then method isadm tests if cl belongs to
      the nearfield or farfield of the cluster represented by this object.
      First the neighbourhood is checked. If cl and this object are not in
      neighbourhood it is necesarry to compute the distance between them.
      Therefor isadm computes an approximation of the distance between the
      accordingly clusters/nodes in dual graph (see getDist()) and an
      approximation of the diameter of the cluster (see getDiam()).

      \param eta2 parameter to justify the admissibility condition (input)
      \param cl check the admissibility for cluster this and cluster cl (input)
      \param info the two attributes info.is_adm and info.is_sep are set
      (output)
      \returns true if the clusters are admissible else false
  */
  bool isadm(double eta2, cluster* cl, bl_info& info);

  /** Method returns the status of this ClusterAlg object. Instances
      of this class are "normal" Clusters.
      \returns false
  */
  virtual bool isSeparator() const {
    return false;
  }
};

#endif
