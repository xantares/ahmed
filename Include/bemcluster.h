/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef BEMCLUSTER_H
#define BEMCLUSTER_H

#include "cluster_pca.h"
template<class T>
class bemcluster : public cluster3d_pca<T>
{
private:
  // index closest to the centroid (needed for the starting index of ACA)
  unsigned icom;

public:
  // cluster with entries k <= i < l
  bemcluster(T* dofs, unsigned* op_perm, unsigned k, unsigned l) : cluster3d_pca<T>(dofs, op_perm, k, l)
  {
    // compute icom
    icom = this -> nbeg;

    const unsigned n = cluster3d_pca<T>::dim();
    unsigned i, j;

    double smin = 0.0;
    for (j=0; j<n; ++j) {
      const double e = cluster3d_pca<T>::getcom(j) - dofs[this->nbeg].getcenter(j);
      smin += e*e;
    }

    for (i=this->nbeg+1; i<this->nend; ++i) {
      double s = 0.0;
      for (j=0; j<n; ++j) {
        const double e = cluster3d_pca<T>::getcom(j) - dofs[i].getcenter(j);
        s += e*e;
      }
      if (s<smin) {
        icom = i;
        smin = s;
      }
    }
  }

  //! required for recursive construction
  virtual bemcluster* clone(unsigned* op_perm, unsigned beg, unsigned end) const {
    return new bemcluster(cluster_pca<T>::dofs, op_perm, beg, end);
  }

  //! returns the icom index
  unsigned geticom() const { return icom; }
};




#endif
