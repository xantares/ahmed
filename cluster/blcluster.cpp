/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <assert.h>
#include "blcluster.h"
#include "cluster.h"

// generates a block cluster tree starting from block (cl1, cl2)
// the number of generated leaves is added to nblcks
// n is used for the number of generated blocks when starting from each pair
// -> goes to blcluster.n



/*! \brief generate block cluster tree, returns a pointer to this tree
    \param cl1   pointer to a block cluster tree for the degrees of freedoms of
                 the rows (input)
    \param cl2   pointer to a block cluster tree for the degrees of freedoms of
                 the columns (input)
    \param eta    the cluster parameter eta (input)
    \param nblcks the number of generated blcks (output)
    \param lvl    the depth of the current block
    \param maxdepth maximum depth

    \par Description:
    This function generates a block cluster tree by subdividing each block
    into four subblocks. The tree is terminated with blocks that are either
    admissible or one of the dimensions is smaller than bmin. In this case the
    first son of such a cluster will point to NULL. */


void blcluster::subdivide_(cluster* cl1, cluster* cl2,
                           double eta2, unsigned& nblcks,
                           unsigned lvl, unsigned maxdepth)
{
  if (cl1->isadm(eta2, cl2, info)) {       // cluster pair is admissible
    setsons(0, 0, NULL);
    setidx(nblcks++);
  } else {
    if ((!maxdepth || lvl<maxdepth) && cl1->isnleaf() && cl2->isnleaf()) {
      ns1 = cl1->getns();
      ns2 = cl2->getns();
      sons = new blcluster*[getns()];
      ++lvl;
      for (unsigned i=0; i<ns1; ++i) {
        cluster* cl1s = cl1->getson(i);
        for (unsigned j=0; j<ns2; ++j) {
          cluster* cl2s = cl2->getson(j);
          sons[i*ns2+j] = clone(cl1s, cl2s);
          sons[i*ns2+j]->subdivide_(cl1s, cl2s, eta2, nblcks, lvl, maxdepth);
        }
      }
    } else {
      setsons(0, 0, NULL);
      setidx(nblcks++);
    }
  }
}

// ----------------------------------------------------------------------------
// symmetric block cluster tree

/*! \brief generate a symmetric block cluster tree, returns a pointer to
           this tree

    \param cl   pointer to a block cluster tree for the degrees of freedoms of
                 the rows (input)
    \param eta    the cluster parameter eta (input)
    \param nblcks the number of generated blcks (output)
    \param lvl    the depth of the current block
    \param maxdepth maximum depth

    \par Description:
    This function generates a block cluster tree by subdividing each block
    into four subblocks. The tree is terminated with blocks that are either
    admissible or one of the dimensions is smaller than bmin. In this case the
    first son of such a cluster will point to NULL.
    In contrast to genblcltree this function sets the third son of each
    diagonal block to NULL. 
*/

void blcluster::subdivide_sym_(cluster* cl, double eta2, unsigned& nblcks,
                               unsigned lvl, unsigned maxdepth)
{
  setadm(false);

  if ((!maxdepth || lvl<maxdepth) && cl->isnleaf()) {
    ns1 = ns2 = cl->getns();
    sons = new blcluster*[getns()];
    ++lvl;
    for (unsigned i=0; i<ns1; ++i) {
      unsigned j;
      for (j=0; j<i; ++j) sons[i*ns2+j] = NULL;
      cluster* cl1s = cl->getson(i);
      sons[i*(ns2+1)] = clone(cl1s, cl1s);
      sons[i*(ns2+1)]->subdivide_sym_(cl1s, eta2, nblcks, lvl, maxdepth);
      for (j=i+1; j<ns2; ++j) {
        cluster*cl2s = cl->getson(j);
        sons[i*ns2+j] = clone(cl1s, cl2s);
        sons[i*ns2+j]->subdivide_(cl1s, cl2s, eta2, nblcks, lvl, maxdepth);
      }
    }
  } else {
    setsons(0, 0, NULL);
    setidx(nblcks++);
  }
}


