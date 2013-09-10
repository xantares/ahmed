/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef BEMBLCLUSTER_H
#define BEMBLCLUSTER_H

#include "blcluster.h"
#include "bemcluster.h"

template<class T1, class T2>
class bemblcluster : public blcluster_geo
{
  bemcluster<T1>* cl1;
  bemcluster<T2>* cl2;

public:
  bemblcluster(unsigned i1, unsigned i2, unsigned m1, unsigned m2)
      : blcluster_geo(i1, i2, m1, m2) { }

  bemblcluster(bemcluster<T1>* p1, bemcluster<T2>* p2) : blcluster_geo(p1, p2),
      cl1(p1), cl2(p2) { }

  //! copy constructor
  bemblcluster(bemblcluster<T1,T2>* bl) : blcluster_geo(bl), cl1(bl->cl1), cl2(bl->cl2) { }

  blcluster* clone(cluster* p1, cluster* p2) const {
    return new bemblcluster((bemcluster<T1>*)p1, (bemcluster<T2>*)p2);
  }

  bemcluster<T1>* getcl1() const { return cl1; }
  bemcluster<T2>* getcl2() const { return cl2; }
};

#endif
