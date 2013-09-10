/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef HELPER_H
#define HELPER_H

#include "mblock.h"
#include "blcluster.h"

template<class T>
bool nrmiszero(blcluster* const bl, mblock<T>** const A)
{
  if (bl->isleaf() and bl->data(A)) {
    if (A[bl->getidx()]->nrmF2()==0)
      return true;
    else
      return false;
  } else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; i++) {
      for (unsigned k=0; k<ns; k++) {
        if (bl->getson(i,k))
          if (nrmiszero(bl->getson(i,k),A))
            return true;
      }
    }
  }
  return false;
}


#endif
