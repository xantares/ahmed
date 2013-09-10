/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef CLUSTERBBX_H
#define CLUSTERBBX_H

#include "cluster.h"
#include <cmath>
#include "blas.h"
#include "basmod.h"

/*!
  \brief A class for storing clusters of degrees of freedom.
  Subdivision is based on bounding boxes (bbx)
*/
template<class T>
class cluster_bbx : public cluster_geo
{
 protected:
  double cntrdir;           // the center with respect to the main direction
  double *xminmax;
  unsigned maindir;         // main direction
  T* dofs;

  double dist2(cluster_bbx<T>* cl) const
  {
    double d = 0.0;
    const unsigned n = dim();

    for (unsigned i=0; i<n; ++i) {
      const double e = cl->xminmax[i] - xminmax[i+n];
      if (e>0.0) d += e*e;
      else {
        const double e1 = xminmax[i] - cl->xminmax[i+n];
        if (e1>0.0) d += e1*e1;
      }
    }
    return d;
  }

 public:
  virtual bool isadm(double eta, cluster* cl, bl_info& info) = 0;

  void init(unsigned* op_perm)
  {
    const unsigned n = dim();
    unsigned i;

    xminmax = new double[2*n];
    assert(xminmax!=NULL);

    // compute bounding box
    for (i=0; i<n; ++i) {
      T* const v = cluster_bbx<T>::dofs + op_perm[nbeg];
      xminmax[i] = v->getcenter(i) - sqrt(v->getradius2());
      xminmax[i+n] = v->getcenter(i) + sqrt(v->getradius2());
    }

    for (unsigned j=nbeg+1; j<nend; ++j) {
      T* const v = cluster_bbx<T>::dofs + op_perm[j];
      for (i=0; i<n; ++i) {
        const double x = v->getcenter(i);
        const double r = sqrt(v->getradius2());
        if (x-r<xminmax[i]) xminmax[i] = x-r;
        else if (x+r>xminmax[i+n]) xminmax[i+n] = x+r;
      }
    }

    // calculate diam2 and main direction
    diam2 = 0.0;
    unsigned imax = 0;
    double max = 0.0;

    for (i=0; i<n; ++i) {
      const double e = xminmax[i+n]-xminmax[i];
      diam2 += e*e;
      if (e>max) {
        max = e;
        imax = i;
      }
    }

    if (xminmax[n]-xminmax[0]>0.8*max)
      maindir = 0;
    else
      maindir = imax;

    // calculate cntrdir
    cntrdir = 0.5*(xminmax[maindir]+xminmax[maindir+n]);
  }

  virtual void createClusterTree(unsigned bmin, unsigned* op_perm,
				 unsigned* po_perm)
  {
    unsigned nnbeg = nend;

    for (unsigned i=nbeg; i<nnbeg;) {
      T* const v = cluster_bbx<T>::dofs + op_perm[i];

      if (v->getcenter(maindir) < cntrdir) ++i;
      else {                           // interchange dofs[i] and dofs[nnbeg]
        swap(op_perm[i], op_perm[--nnbeg]);
        po_perm[op_perm[i]] = i;
        po_perm[op_perm[nnbeg]] = nnbeg;
      }
    }

    if (nnbeg-nbeg>bmin && nend-nnbeg>bmin) {

      // first new son
      cluster* son[2];
      son[0] = clone(op_perm, nbeg, nnbeg);
      son[0]->createClusterTree(bmin, op_perm, po_perm);

      // second new son
      son[1] = clone(op_perm, nnbeg, nend);
      son[1]->createClusterTree(bmin, op_perm, po_perm);
      cluster::setsons(2, son);
    }
  }

  cluster_bbx(T* dofs_, unsigned k, unsigned l) : cluster_geo(k, l), dofs(dofs_)
  { }
  virtual ~cluster_bbx() { delete [] xminmax; }
};

//! a class for storing clusters of degrees of freedom in 1d
template<class T>
class cluster1d_bbx : public cluster_bbx<T>
{
protected:
  virtual unsigned dim() const { return 1; }
  virtual cluster_bbx<T>* clone(unsigned* op_perm,
				unsigned beg, unsigned end) const
  { return new cluster1d_bbx<T>(cluster_bbx<T>::dofs, op_perm, beg, end); }

public:
  bool isadm(double eta, cluster* cl, bl_info& info)
  {
    cluster1d_bbx<T>* p = (cluster1d_bbx<T>*) cl;
    assert(p!=NULL);  // cl is actually a cluster_bbx*
    info.is_sep = false;

    const double d2 = MIN(cluster_bbx<T>::diam2, p->diam2);
    if (d2<eta*dist2(p)) return (info.is_adm=true);
    /*
    double di2 = 0.0;
    if (0<xminmax[0]) di2 = xminmax[0]*xminmax[0];
    if (xminmax[1]<0) di2 = xminmax[1]*xminmax[1];

    if (diam2/di2<eta2) return true;

    di2 = 0.0;
    if (0<p->xminmax[0]) di2 = p->xminmax[0]*p->xminmax[0];
    if (p->xminmax[1]<0) di2 = p->xminmax[1]*p->xminmax[1];

    if (p->diam2/di2<eta2) return true;
    */

    // weak admissibility
    if (cluster_bbx<T>::xminmax[1]==p->xminmax[0]) return (info.is_adm=true);
    if (p->xminmax[1]==cluster_bbx<T>::xminmax[0]) return (info.is_adm=true);

    return (info.is_adm=false);
  }

  cluster1d_bbx(T* dofs_, unsigned* op_perm, unsigned k, unsigned l) : cluster_bbx<T>(dofs_, k, l) {
    cluster_bbx<T>::init(op_perm);
  }
};

//! a class for storing clusters of degrees of freedom in 2d
template<class T>
class cluster2d_bbx : public cluster_bbx<T>
{
protected:
  virtual unsigned dim() const { return 2; }
  virtual cluster_bbx<T>* clone(unsigned* op_perm,
				unsigned beg, unsigned end) const
  { return new cluster2d_bbx(cluster_bbx<T>::dofs, op_perm, beg, end); }

public:
  bool isadm(double eta, cluster* cl, bl_info& info)
  {
    cluster2d_bbx<T>* p = (cluster2d_bbx<T>*) cl;
    assert(p!=NULL);  // cl is actually a cluster_bbx*
    info.is_sep = false;

    const double d2 = MIN(cluster_bbx<T>::diam2, p->diam2);
    if (d2<eta*this->dist2(p)) return (info.is_adm=true);
    /*
    double di2 = 0.0;
    if (0<xminmax[0]) di2 = xminmax[0]*xminmax[0];
    if (xminmax[1]<0) di2 = xminmax[1]*xminmax[1];

    if (diam2/di2<eta2) return true;

    di2 = 0.0;
    if (0<p->xminmax[0]) di2 = p->xminmax[0]*p->xminmax[0];
    if (p->xminmax[1]<0) di2 = p->xminmax[1]*p->xminmax[1];

    if (p->diam2/di2<eta2) return true;
    */

    // weak admissibility
    if (cluster_bbx<T>::xminmax[2]==p->xminmax[0]) {
      if (cluster_bbx<T>::xminmax[3]==p->xminmax[1])
	return (info.is_adm=true);
      if (p->xminmax[3]==cluster_bbx<T>::xminmax[1])
	return (info.is_adm=true);
    }
    if (p->xminmax[2]==cluster_bbx<T>::xminmax[0]) {
      if (cluster_bbx<T>::xminmax[3]==p->xminmax[1])
	return (info.is_adm=true);
      if (p->xminmax[3]==cluster_bbx<T>::xminmax[1])
	return (info.is_adm=true);
    }
    return (info.is_adm=false);
  }

  cluster2d_bbx(T* dofs_, unsigned* op_perm, unsigned k, unsigned l) : cluster_bbx<T>(dofs_, k, l) {
    cluster_bbx<T>::init(op_perm);
  }
};

 
//! a class for storing clusters of degrees of freedom in 3d
template<class T>
class cluster3d_bbx : public cluster_bbx<T>
{
protected:
  virtual unsigned dim() const { return 3; }
  virtual cluster_bbx<T>* clone(unsigned* op_perm, unsigned beg, unsigned end) const
  { return new cluster3d_bbx(cluster_bbx<T>::dofs,op_perm, beg, end); }

public:
  bool isadm(double eta, cluster* cl, bl_info& info) 
  {
    cluster3d_bbx<T>* p = (cluster3d_bbx<T>*) cl;
    assert(p!=NULL);  // cl is actually a cluster_bbx*
    info.is_sep = false;

    const double d2 = MIN(cluster_bbx<T>::diam2, p->diam2);
    if (d2<eta*this->dist2(p)) return (info.is_adm=true);
    /*
    double di2 = 0.0;
    if (0<xminmax[0]) di2 = xminmax[0]*xminmax[0];
    if (xminmax[1]<0) di2 = xminmax[1]*xminmax[1];

    if (diam2/di2<eta2) return true;

    di2 = 0.0;
    if (0<p->xminmax[0]) di2 = p->xminmax[0]*p->xminmax[0];
    if (p->xminmax[1]<0) di2 = p->xminmax[1]*p->xminmax[1];

    if (p->diam2/di2<eta2) return true;
    */

    // weak admissibility
    if (cluster_bbx<T>::xminmax[3]==p->xminmax[0]) {
      if (cluster_bbx<T>::xminmax[4]==p->xminmax[1]) {
	if (cluster_bbx<T>::xminmax[5]==p->xminmax[2])
	  return (info.is_adm=true);
	if (p->xminmax[5]==cluster_bbx<T>::xminmax[2])
	  return (info.is_adm=true);
      }
      if (p->xminmax[4]==cluster_bbx<T>::xminmax[1]) {
	if (cluster_bbx<T>::xminmax[5]==p->xminmax[2])
	  return (info.is_adm=true);
	if (p->xminmax[5]==cluster_bbx<T>::xminmax[2])
	  return (info.is_adm=true);
      }
    }
    if (p->xminmax[3]==cluster_bbx<T>::xminmax[0]) {
      if (cluster_bbx<T>::xminmax[4]==p->xminmax[1]) {
	if (cluster_bbx<T>::xminmax[5]==p->xminmax[2])
	  return (info.is_adm=true);
	if (p->xminmax[5]==cluster_bbx<T>::xminmax[2])
	  return (info.is_adm=true);
      }
      if (p->xminmax[4]==cluster_bbx<T>::xminmax[1]) {
	if (cluster_bbx<T>::xminmax[5]==p->xminmax[2])
	  return (info.is_adm=true);
	if (p->xminmax[5]==cluster_bbx<T>::xminmax[2])
	  return (info.is_adm=true);
      }
    }
    return (info.is_adm=false);
  }

  cluster3d_bbx(T* dofs_, unsigned* op_perm, unsigned k, unsigned l) : cluster_bbx<T>(dofs_, k, l) {
    cluster_bbx<T>::init(op_perm);
  }
};

#endif
