/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file   blcluster.h
    \brief  Include file for the class blcluster
    \author M. Bebendorf
    \date   4/24/2004                            */

#ifndef BLCLUSTER_H
#define BLCLUSTER_H

#include "mblock.h"
#include "cluster.h"

//! in this class the block structure of the H-matrix is stored
class blcluster
{
  friend blcluster* subdivide_augm(unsigned, unsigned, cluster*, cluster*,
                                   double, unsigned&);
  friend blcluster* subdivide_sym_augm(unsigned, unsigned, cluster*, double,
                                       unsigned&);
protected:

  //! dimensions of this block
  unsigned n1, n2;

  //! position in global matrix (inidices of first row/column)
  unsigned b1, b2;

  //! number of row and column sons
  unsigned ns1, ns2;

  //! the sons of this block cluster
  blcluster** sons;

  //! the index of this block, valid only if this block is a leaf
  unsigned idx;

  //! properties of this block
  bl_info info;

  bl_info& getinfo() {
    return info;
  } 

  void subdivide_(cluster*, cluster*, double, unsigned&, unsigned, unsigned);
  void subdivide_sym_(cluster*, double, unsigned&, unsigned, unsigned);

public:
  //! empty constructor
  blcluster() : ns1(0), ns2(0), sons(NULL) { }
  blcluster(unsigned i1, unsigned i2, unsigned m1, unsigned m2)
      : n1(m1), n2(m2), b1(i1), b2(i2), ns1(0), ns2(0), sons(NULL) {
    info.is_sep = false;
  }

  //! another constructor
  blcluster(cluster* cl1, cluster* cl2) : ns1(0), ns2(0), sons(NULL) {
    n1 = cl1->size();
    n2 = cl2->size();
    b1 = cl1->getnbeg();
    b2 = cl2->getnbeg();
    info.is_sep = false;
  }

  //! copy constructor
  blcluster(blcluster* bl) : n1(bl->n1), n2(bl->n2), b1(bl->b1), b2(bl->b2),
      ns1(bl->ns1), ns2(bl->ns2), idx(bl->idx), info(bl->info) {
    unsigned ns = getns();
    if (ns) {
      sons = new blcluster*[ns];
      for (unsigned i=0; i<ns; ++i) {
        if (bl->sons[i]) sons[i] = new blcluster(bl->sons[i]);
        else sons[i] = NULL;
      }
    } else sons = NULL;
  }

  //! destructor
  virtual ~blcluster() {
    unsigned ns = getns();
    for (unsigned i=0; i<ns; ++i) {
      delete sons[i];
      sons[i] = NULL;
    }
    delete [] sons;
    ns1 = ns2 = 0;
    sons = NULL;
  }

  //! set the sons
  void setsons(unsigned n1, unsigned n2, blcluster** p) {
    if(sons) delete [] sons;
    unsigned ns = (ns1 = n1) * (ns2 = n2);
    if (ns) {
      sons = new blcluster*[ns];
      assert(sons);
      for (unsigned i=0; i<ns; ++i) sons[i] = p[i];
    } else sons = NULL;
  }

  //! the row index of the first entry of this block within the H-matrix
  unsigned getb1() const {
    return b1;
  }
  //! the column index of the first entry of this block within the H-matrix
  unsigned getb2() const {
    return b2;
  }

  //! the number of rows of this block
  unsigned getn1() const {
    return n1;
  }
  //! the number of columns of this block
  unsigned getn2() const {
    return n2;
  }

  //! is this block admissible ?
  bool isadm() const {
    return info.is_adm;
  }
  bool isnadm() const {
    return !isadm();
  }
  void setadm(bool b) {
    info.is_adm = b;
  }

  //! are the clusters of this separated ?
  bool issep() const {
    return info.is_sep;
  }
  void setsep(bool b) {
    info.is_sep = b;
  }

  //! is this a leaf ?
  bool isnleaf() const {
    return (ns1 && ns2);
  }
  bool isleaf() const {
    return !isnleaf();
  }

  //! number of row sons
  unsigned getnrs() const {
    return ns1;
  }

  //! number of column sons
  unsigned getncs() const {
    return ns2;
  }

  //! number of sons
  unsigned getns() const {
    return ns1*ns2;
  }

  blcluster* getson(unsigned i, unsigned j) const {
    return sons[i*ns2+j];
  }

  //! is this a diagonal block ?
  bool isdbl() const {
    return (b1==b2);
  }
  bool isndbl() const {
    return !isdbl();
  }

  unsigned getidx() const {
    return idx;
  }
  void setidx(unsigned i) {
    idx = i;
  }

  // returns the column-middle of blocks, returned idx belongs to lower part
  unsigned long getmid() const {
    if (isleaf()) return 0;
    else {
      unsigned temp1 = 0, temp2 = 0, mid = 0;
      for (unsigned i=0; i<getncs(); i++) {
        temp2 += getson(0,i)->getn2();
        if (std::abs(temp1/(double)getn2()-0.5)
            >std::abs(temp2/(double)getn2()-0.5))
          mid += getson(0,i)->getn2();
        temp1 = temp2;
      }
      return mid;
    }
  }

  //! returns the size of the block cluster tree starting from this
  unsigned long size() const {
    if (isleaf()) return sizeof(blcluster);
    else {
      unsigned long size = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i) if (sons[i]) size += sons[i]->size();
      return size;
    }
  }

  //! returns the number of leaves in the block cluster tree starting from this
  unsigned long nleaves() const {
    if (isleaf()) return 1;
    else {
      unsigned long nl = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i)
        if (sons[i]) nl += sons[i]->nleaves();
      return nl;
    }
  }

  //! returns the number of leaves in the upper triangular part
  unsigned long nupleaves() const {
    if (isleaf() && b1<=b2) return 1;
    else {
      unsigned long nl = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i)
        if (sons[i] && sons[i]->b1<=sons[i]->b2) nl += sons[i]->nleaves();
      return nl;
    }
  }

  //! returns the number of leaves in the upper triangular part
  unsigned long nlwleaves() const {
    if (isleaf() && b1>=b2) return 1;
    else {
      unsigned long nl = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i)
        if (sons[i] && sons[i]->b1>=sons[i]->b2) nl += sons[i]->nleaves();
      return nl;
    }
  }

  //! returns the number of leaves on the block diagonal (i.e. b1==b2)
  unsigned long ndbls() const {
    if (isleaf()) return 1;
    else {
      unsigned long ndbl = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i)
        if (sons[i] && sons[i]->isdbl()) ++ndbl;
      return ndbl;
    }
  }

  unsigned long nadmleaves() const {
    if (isleaf()) return isadm();
    else {
      unsigned long nl = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i)
        if (sons[i]) nl += sons[i]->nadmleaves();
      return nl;
    }
  }

  unsigned long nnadmleaves() const {
    if (isleaf()) return isnadm();
    else {
      unsigned long nl = 0;
      unsigned ns = getns();
      for (unsigned i=0; i<ns; ++i)
        if (sons[i]) nl += sons[i]->nnadmleaves();
      return nl;
    }
  }


  //! is the corresponding block in A stored as a low-rank matrix ?
  bool isLrM(mblock<double>** const A) const {
    return A[idx]->isLrM();
  }
  bool isLrM(mblock<float>** const A) const {
    return A[idx]->isLrM();
  }
  bool isLrM(mblock<scomp>** const A) const {
    return A[idx]->isLrM();
  }
  bool isLrM(mblock<dcomp>** const A) const {
    return A[idx]->isLrM();
  }

  //! is the corresponding block in A stored as a dense matrix ?
  bool isGeM(mblock<double>** const A) const {
    return A[idx]->isGeM();
  }
  bool isGeM(mblock<float>** const A) const {
    return A[idx]->isGeM();
  }
  bool isGeM(mblock<scomp>** const A) const {
    return A[idx]->isGeM();
  }
  bool isGeM(mblock<dcomp>** const A) const {
    return A[idx]->isGeM();
  }

  //! is the corresponding block in A stored as a lower triangular matrix ?
  bool isLtM(mblock<double>** const A) const {
    return A[idx]->isLtM();
  }
  bool isLtM(mblock<float>** const A) const {
    return A[idx]->isLtM();
  }
  bool isLtM(mblock<scomp>** const A) const {
    return A[idx]->isLtM();
  }
  bool isLtM(mblock<dcomp>** const A) const {
    return A[idx]->isLtM();
  }

  //! is the corresponding block in A stored as a upper triangular matrix ?
  bool isUtM(mblock<double>** const A) const {
    return A[idx]->isUtM();
  }
  bool isUtM(mblock<float>** const A) const {
    return A[idx]->isUtM();
  }
  bool isUtM(mblock<scomp>** const A) const {
    return A[idx]->isUtM();
  }
  bool isUtM(mblock<dcomp>** const A) const {
    return A[idx]->isUtM();
  }

  //! returns the rank of the corresponding block in A
  unsigned rank(mblock<double>** const A) const {
    assert(isLrM(A));
    return A[idx]->rank();
  }
  unsigned rank(mblock<float>** const A) const {
    assert(isLrM(A));
    return A[idx]->rank();
  }
  unsigned rank(mblock<scomp>** const A) const {
    assert(isLrM(A));
    return A[idx]->rank();
  }
  unsigned rank(mblock<dcomp>** const A) const {
    assert(isLrM(A));
    return A[idx]->rank();
  }

  //! returns the pointer to the data of the corresponding block in A
  double* data(mblock<double>** const A) const {
    assert(isleaf());
    return A[idx]->data;
  }
  float* data(mblock<float>** const A) const {
    assert(isleaf());
    return A[idx]->data;
  }
  scomp* data(mblock<scomp>** const A) const {
    assert(isleaf());
    return A[idx]->data;
  }
  dcomp* data(mblock<dcomp>** const A) const {
    assert(isleaf());
    return A[idx]->data;
  }

  virtual blcluster* clone(cluster* cl1, cluster* cl2) const {
    return new blcluster(cl1, cl2);
  }
  virtual void subdivide(cluster* cl1, cluster* cl2, double eta2, unsigned& nblcks,
                 unsigned maxdepth=0) {
    unsigned lvl = 0;
    nblcks = 0;
    subdivide_(cl1, cl2, eta2, nblcks, lvl, maxdepth);
  }
  virtual void subdivide_sym(cluster* cl, double eta2, unsigned& nblcks,
                     unsigned maxdepth=0) {
    unsigned lvl = 0;
    nblcks = 0;
    subdivide_sym_(cl, eta2, nblcks, lvl, maxdepth);
  }
};


// a block cluster class with reference to the generating cluster tree
class blcluster_ref : public blcluster
{
  //! the row cluster and the column cluster
protected:
  cluster_ref *rcl, *ccl;

public:
  blcluster_ref(unsigned i1, unsigned i2, unsigned m1, unsigned m2)
      : blcluster(i1, i2, m1, m2) { }
  blcluster_ref(cluster_ref* cl1, cluster_ref* cl2)
      : blcluster(cl1, cl2) {
    rcl = cl1;
    ccl = cl2;
    if (isdbl()) cl1->dbl = this;
  }
  //! copy constructor
  blcluster_ref(blcluster_ref* bl) : blcluster(bl), rcl(bl->rcl),
      ccl(bl->ccl) { }
  cluster_ref* getrcl() const {
    return rcl;
  }
  cluster_ref* getccl() const {
    return ccl;
  }
  //! number of row cluster sons
  //  unsigned getnrsons() const { return rcl->getns(); }
  //! return row cluster sons
  cluster* getrson(unsigned i) const {
    return rcl->getson(i);
  }

  //! number of column cluster sons
  //unsigned getncsons() const { return ccl->getns(); }
  //! return column cluster sons
  cluster* getcson(unsigned i) const {
    return ccl->getson(i);
  }

  blcluster* clone(cluster* cl1, cluster* cl2) const {
    return new blcluster_ref((cluster_ref*)cl1, (cluster_ref*)cl2);
  }
};


class blcluster_geo : public blcluster
{
public:
  blcluster_geo(unsigned i1, unsigned i2, unsigned m1, unsigned m2)
      : blcluster(i1, i2, m1, m2) { }
  blcluster_geo(cluster* cl1, cluster* cl2)
      : blcluster(cl1, cl2) { }

  //! copy constructor
  blcluster_geo(blcluster_geo* bl) : blcluster(bl) { }

  blcluster* clone(cluster* cl1, cluster* cl2) const {
    return new blcluster_geo(cl1, cl2);
  }
};



extern void genpart_from_tree(blcluster*, unsigned, blcluster**&);
extern void genpart(cluster*, cluster*, unsigned, double, unsigned&,
                    blcluster**&);

extern blcluster *genblcltree_augm(unsigned, unsigned, cluster*,
                                   cluster*, double, unsigned&);
extern blcluster *gensymblcltree_augm(unsigned, unsigned, cluster*,
                                      double, unsigned&);


#endif  // BLCLUSTER_H

