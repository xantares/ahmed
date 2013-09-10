/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


/*! \file   cluster_pca.h
    \brief  Include file for the class cluster_pca
    \author M. Bebendorf
    \date   7/14/2004
*/

#ifndef CLUSTERPCA_H
#define CLUSTERPCA_H

#include "cluster.h"
#include "blas.h"
#include "basmod.h"

template<class T> class bemcluster;

/*!
  \brief A class for storing clusters of degrees of freedom.
  Subdivision is based on the principal component analysis (PCA)
*/

template<class T>
class cluster_pca : public cluster_geo
{
  template<class TT> friend class bemcluster;
  double* com;  

 protected:
  T* dofs;

 public:
  //Constructor
  cluster_pca(T* dofs_, unsigned k, unsigned l) : cluster_geo(k, l), dofs(dofs_)
  { }

  //Destructor
  virtual ~cluster_pca() { delete [] com; }

  //get i-th component of the center of mass (0 <= i < dim())
  double getcom(int i) { return com[i]; }


  //check admissibility criterion (for this and cl) and write to info
  bool isadm(double eta2, cluster* cl, bl_info& info)
  {
    cluster_pca* p = (cluster_pca*) cl;  // cl is actually a cluster_pca*
    assert(p!=NULL);
    info.is_sep = false;
    
    const double d2 = MIN(diam2, p->diam2);
    
    if (d2<eta2*dist2(p))
      return (info.is_adm=true);
    else
      return (info.is_adm=false);
  }

  //compute distance between this and c
  double dist2(cluster_pca<T>* cl) const
  {
    double d2 = 0.0;
    const unsigned n = cl->dim();

    for (unsigned j=0; j<n; ++j) {
      const double d = com[j] - cl->com[j];
      d2 += d*d;
    }

    return d2;
  }

  //required for subdivision with virtual functions
  //init() is called implicitly by subdivide() with virtual function clone() 
  void init(unsigned* op_perm)
  {
    unsigned i, j;

    // calculate the center of mass
    unsigned n = dim();
    com = new double[n];

    for (j=0; j<n; ++j) com[j] = 0.0;

    for (i=nbeg; i<nend; ++i)
      for (j=0; j<n; ++j)
        com[j] += dofs[op_perm[i]].getcenter(j);

    const double d = 1.0/(nend-nbeg);
    for (j=0; j<n; ++j) com[j] *= d;

    // compute the radius of this cluster

    diam2 = 0.0;

    for (i=nbeg; i<nend; ++i) {

      double d2 = 0.0;

      for (j=0; j<n; ++j) {
        const double d1 = 2*(dofs[op_perm[i]].getcenter(j)-com[j]);
        d2 += d1*d1;
      }

      if (d2>diam2) diam2 = d2;
    }
  }

  
  unsigned subdivide(const unsigned bmin, unsigned* op_perm, unsigned* po_perm)
  {
    const unsigned n = dim();
    const unsigned m = n*(n+1)/2;
    unsigned i;
    double* const c = new double[m];     // components of the covariant matrix
    blas::setzero(m, c);

    double* const d = new double[n];

    for (i=nbeg; i<nend; i++) {

      unsigned j;

      for (j=0; j<n; ++j) d[j] = dofs[op_perm[i]].getcenter(j) - com[j];

      unsigned l = 0;
      for (j=0; j<n; ++j)
        for (unsigned k=j; k<n; ++k)
          c[l++] += d[j] * d[k];
    }

    double* const u = new double[n];
    for (i=0; i<n; ++i) u[i] = 1.0/sqrt((double) n);       // starting vector u

    i = 0;

    do {
      unsigned j;
      for (j=0; j<n; ++j) d[j] = u[j];

      for (j=0; j<n; ++j) {

        unsigned k;
        u[j] = 0.0;

        for (k=0; k<j; ++k) {
          const unsigned l = k*n - k*(k+1)/2 + j;
          u[j] += c[l] * d[k];
        }
        for (; k<n; ++k) {
          const unsigned l = j*n - j*(j+1)/2 + k;
          u[j] += c[l] * d[k];
        }
      }

      double modu = 0.0;
      for (j=0; j<n; ++j) modu += u[j]*u[j];

      modu = 1.0/sqrt(modu);
      for (j=0; j<n; ++j) u[j] *= modu;

    } while (++i<10);

    delete [] d;
    delete [] c;

    // u is the main direction of the cluster
    double mdtcl = 0.0;
    for (i=0; i<n; ++i) mdtcl += u[i] * com[i];

    unsigned nsep = nend;

    for (i=nbeg; i<nsep; ) {

      T* const v = dofs + op_perm[i];
      double skp = 0.0;
      for (unsigned j=0; j<n; ++j) skp += v->getcenter(j) * u[j];

      if (skp<mdtcl) ++i;
      else {                              // interchange dofs[i] and dofs[nsep]
        swap(op_perm[i], op_perm[--nsep]);
        po_perm[op_perm[i]] = i;
        po_perm[op_perm[nsep]] = nsep;
      }
    }

    delete [] u;

    if (nsep == nbeg || nsep == nend) nsep = (nend+nbeg)/2;

    return nsep;
  }

  //recursive subdivison via principal component analysis (PCA) until size()<= bmin
  virtual void createClusterTree(const unsigned bmin, unsigned* op_perm,
				 unsigned* po_perm)
  {
    unsigned nsep = subdivide(bmin, op_perm, po_perm);

    if (nsep>=nbeg+bmin && nend>=nsep+bmin) {
      cluster* son[2];
      son[0] = clone(op_perm, nbeg, nsep);
      son[0]->createClusterTree(bmin, op_perm, po_perm);
      son[1] = clone(op_perm, nsep, nend);
      son[1]->createClusterTree(bmin, op_perm, po_perm);
      setsons(2, son);
    }
  }
};


//! a class for storing clusters of degrees of freedom in 2d
template<class T>
class cluster2d_pca : public cluster_pca<T>
{
 protected:
  virtual unsigned dim() const { return 2; }

  virtual cluster_pca<T>* clone(unsigned* op_perm,
				unsigned beg, unsigned end) const
  { return new cluster2d_pca<T>(cluster_pca<T>::dofs, op_perm, beg, end); }

 public:
  cluster2d_pca(T* dofs_, unsigned* op_perm, unsigned k, unsigned l) : cluster_pca<T>(dofs_, k, l) {
    cluster_pca<T>::init(op_perm);
  }
};

//! a class for storing clusters of degrees of freedom in 3d
template<class T>
class cluster3d_pca : public cluster_pca<T>
{
protected:
  virtual unsigned dim() const { return 3; }

  virtual cluster_pca<T>* clone(unsigned* op_perm,
				unsigned beg, unsigned end) const
  { return new cluster3d_pca<T>(cluster_pca<T>::dofs, op_perm, beg, end); }

public:
 cluster3d_pca(T* dofs_, unsigned* op_perm, unsigned k, unsigned l) : cluster_pca<T>(dofs_, k, l) {
    cluster_pca<T>::init(op_perm);
  }
};

#endif
