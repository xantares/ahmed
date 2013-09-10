/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file  Separator.h
    \brief Include file for the class ClusterAlg.
    \author M. Bebendorf and T. Fischer
    \date   2/23/2007
*/
#ifndef CLUSTER_SEPARATOR_H
#define CLUSTER_SEPARATOR_H

#include "ClusterAlg.h"

class cluster_alg;

/** \brief class for storing clusters of degrees of freedom without
 geometric information.
*/
class Separator : public ClusterAlg
{
public:
  /** brief Constructor builds a initial object for clustering
      \param father pointer to the father node in cluster tree
      \param beg index in op_perm and po_perm which describes the begin of the index set of the Separator
      \param end index in op_perm and po_perm which describes the begin of the index set of the next
      ClusterAlg
      \param op_perm permutation 
      \param po_perm inverse permutation 
      \param depth distance between this node and the root
      \param global_mat reference to adjacency matrix of the matrix graph in crs format
      \param local_mat reference to the local adjacency matrix of the matrix graph in crs format
  */
  Separator(ClusterAlg* father, unsigned beg, unsigned end,
	    unsigned* op_perm, unsigned* po_perm, unsigned depth,
	    AdjMat* global_mat, AdjMat* local_mat)
  : ClusterAlg (father, beg, end, op_perm, po_perm, depth, global_mat,
		local_mat)
  {
    if (parent == NULL)
      std::cerr << "Separator::Separator parent == NULL" << std::endl;
  }
  
  /** Destructor. */
  virtual ~Separator() {}

  /** \brief subdivides this cluster recursively, needed to generate the
      cluster tree
      \param bmin the minimal size of clusters (input)
  */
  virtual void subdivide(unsigned bmin);
  
  virtual void refine(unsigned max_level);
  
  /** Method traverses the ClusterTree and computes the Boundaries.
   */
  virtual void computeNeighbor();
  
  /** \brief This method checks the admissibility of cluster this with
      the cluster cl.
      
      - If one of the clusters this or cl are separators, then isadm is looking
      into the graph, whereas nodes are set of clusters and edges are
      boundaries/separators, for a suitable node (or cluster) (see
      ClusterAlg::getAnotherCluster()) to check the admissible condition.
      Therefor it computes the distance between the accordingly clusters/nodes
      (see getDist()) and the diameter of the accordingly
      clusters/nodes (see getDiam())
      
      \param eta2 parameter to justify the admissibility condition (input)
      \param cl check the admissibility for cluster this and cluster cl (input)
      \param info the two attributes info.is_adm and info.is_sep are set
      (output)
      \returns true if the clusters are admissible else false
  */
  virtual bool isadm(double eta2, cluster* cl, bl_info& info);
  
  /** Method returns the status of this ClusterAlg object. Instances
      of this class are Separators. 
      \returns true
  */
  virtual bool isSeparator() const { return true; }
  
 private:
  /** update perm */
  unsigned updatePerm(unsigned *reordering);
};

#endif
