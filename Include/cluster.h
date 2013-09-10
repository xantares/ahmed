/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file  cluster.h
  \brief Include file for the class cluster
  \author M. Bebendorf
  \date   3/24/2004
*/

#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>
#include <limits>
#include <assert.h>

struct bl_info {
  unsigned is_adm : 1;      // is this block admissible ?
  unsigned is_sep : 1;      // do the clusters have neighb. vertices ?
  // this is only for sparse matrices
};

class blcluster_ref;

//! the basis class for storing clusters of degrees of freedom
class cluster
{
  void getncl_(unsigned& n) {
    if (nsons>0) {
      n += nsons;
      for (unsigned i=0; i<nsons; ++i) sons[i]->getncl_(n);
    }
  }

protected:
  //! beginning and ending index (beginning index of next)
  unsigned nbeg, nend;

  //! number of sons, nsons==0 if this is a leaf
  unsigned nsons;

  //! the array of sons of this cluster in the cluster tree
  cluster* *sons;

public:
  /*!
    \brief Constructor
    \param k beginning index of the cluster
    \param l beginning index of the next cluster
  */
  cluster(unsigned k, unsigned l) : nbeg(k), nend(l), nsons(0), sons(NULL) {}

  virtual ~cluster() {
    for (unsigned i=0; i<nsons; ++i) delete sons[i];
    nsons = 0;
    delete [] sons;
    sons = NULL;
  }

  unsigned getnbeg() const {
    return nbeg;
  }
  void setnbeg (unsigned beg) {
    nbeg = beg;
  }
  unsigned getnend() const {
    return nend;
  }
  void setnend(unsigned end) {
    nend = end;
  }
  unsigned size() const {
    return nend-nbeg;
  }

  bool isnleaf() const {
    return (nsons!=0);
  }
  bool isleaf() const {
    return !isnleaf();
  }

  unsigned getns() const {
    return nsons;
  }
  cluster* getson(unsigned i) const {
    assert(i<nsons);
    return sons[i];
  }
  void setsons(unsigned n, cluster** p) {
    for (unsigned i=0; i<nsons; ++i) delete sons[i];
    delete [] sons;
    if ((nsons=n) > 0) {
      sons = new cluster*[nsons];
      for (unsigned i=0; i<nsons; ++i) sons[i] = p[i];
    } else sons = NULL;
  }

  //! is this cluster admissible with the cluster cl ?
  //! bool& only valid for clusterSpND
  virtual bool isadm(double eta2, cluster* cl, bl_info& info)=0;

  /**
   * creates recursively the cluster tree, i.e. changes the permutation op_perm and po_perm
   * and create child cluster trees
   * @param op_perm permutation: permutated_row = op_perm[original_row]
   * @param po_perm reverse permutation: original_row = po_perm[permutated_idx]
   * @param bmin threshold value for stopping further refinement
   */
  virtual void createClusterTree(unsigned bmin, unsigned* op_perm,
				 unsigned* po_perm) = 0;

  void getMaxDepth_(unsigned &act_depth, unsigned &max_depth) {
    if (nsons>0) {
      act_depth++;
      if (act_depth > max_depth) max_depth = act_depth;
      for (unsigned i=0; i<nsons; ++i)
        sons[i]->getMaxDepth_(act_depth, max_depth);
      act_depth--;
    }
  }
  unsigned getMaxDepth() {
    unsigned act_depth = 0, max_depth = 0;
    getMaxDepth_(act_depth,max_depth);
    return max_depth;
  }

  void getMinDepth_(unsigned &act_depth, unsigned &min_depth) {
    if (nsons>0) {
      act_depth++;
      for (unsigned i=0; i<nsons; i++) sons[i]->getMinDepth_(act_depth, min_depth);
      act_depth--;
    } else {
      if (act_depth<min_depth) min_depth = act_depth;
    }
  }
  unsigned getMinDepth() {
    unsigned act_depth = 0, min_depth = std::numeric_limits<unsigned>::max();
    getMinDepth_(act_depth, min_depth);
    return min_depth;
  }

  unsigned getncl() {
    unsigned n = 0;
    getncl_(n);
    return n;
  }

  // the array stores:
  // nsons, sons[0]->getnend(), ..., sons[nsons-1]->getnend(),
  // sons[0]->writeUnsigned(), ... sons[nsons-1]->writeUnsigned()
  unsigned* writeUnsigned(unsigned* clTree) {
    clTree[0] = nsons;
    ++clTree;
    if (nsons > 0) {
      for (unsigned i=0; i<nsons; ++i) clTree[i] = sons[i]->getnend();
      clTree += nsons;
      for (unsigned i=0; i<nsons; ++i)
        clTree = sons[i]->writeUnsigned(clTree);
    }
    return clTree;
  }

  virtual unsigned* readUnsigned(unsigned*) { return 0; }
};

//! basis class for storing clusters with a reference to diagonal blocks
class cluster_ref : public cluster
{
  friend class blcluster_ref;

  blcluster_ref* dbl; // bl is the diagonal block corresponding to this cluster
  // = NULL if no such block exists
public:
  cluster_ref(unsigned k, unsigned l) : cluster(k, l), dbl(NULL) { }
  blcluster_ref* getdbl() const {
    return dbl;
  }
};

//! basis class of a geometric cluster
class cluster_geo : public cluster
{
 protected:
  //! squared diameter
  double diam2;

public:
  cluster_geo(unsigned k, unsigned l) : cluster(k, l), diam2(0.0) { }

  //! returns the spatial dimension
  virtual unsigned dim() const = 0;

  virtual cluster_geo* clone(unsigned*, unsigned, unsigned) const = 0;
  double getdiam2() { return diam2; }
};

#endif
