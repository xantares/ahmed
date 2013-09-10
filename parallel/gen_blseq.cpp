/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"


// ((compare with gen_NortonBlSeq !!!) can this be removed ? -> OpenMP ?
static void gen_BlSequence_(blcluster* bl, blcluster**& BlList)
{
  if (bl->isleaf()) *BlList++ = bl;
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) gen_BlSequence_(son, BlList);
      }
    }
  }
}

static void gen_upBlSequence_(blcluster* bl, blcluster**& BlList,
                              unsigned& nupblcks)
{
  if (bl->isleaf() && bl->getb1()<=bl->getb2()) {
    *BlList++ = bl;
    ++nupblcks;
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son && son->getb1()<=son->getb2())
          gen_upBlSequence_(son, BlList, nupblcks);
      }
    }
  }
}

static void gen_lwBlSequence_(blcluster* bl, blcluster**& BlList,
                              unsigned& nlwblcks)
{
  if (bl->isleaf() && bl->getb1()>=bl->getb2()) {
    *BlList++ = bl;
    ++nlwblcks;
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son && son->getb1()>=son->getb2())
          gen_lwBlSequence_(son, BlList, nlwblcks);
      }
    }
  }
}



void gen_BlSequence(blcluster* bl, blcluster**& BlList)
{
  BlList = new blcluster*[bl->nleaves()];
  assert(BlList!=NULL);

  blcluster** bllist = BlList;
  gen_BlSequence_(bl, bllist);
}


// generate sequence of blocks from the upper triangular part
void gen_upBlSequence(blcluster* bl, blcluster**& BlList, unsigned& nupblcks)
{
  BlList = new blcluster*[bl->nupleaves()];
  assert(BlList!=NULL);

  blcluster** bllist = BlList;
  nupblcks = 0;
  gen_upBlSequence_(bl, bllist, nupblcks);
}


// generate sequence of blocks from the lower triangular part
void gen_lwBlSequence(blcluster* bl, blcluster**& BlList, unsigned& nlwblcks)
{
  BlList = new blcluster*[bl->nlwleaves()];
  assert(BlList!=NULL);

  blcluster** bllist = BlList;
  nlwblcks = 0;
  gen_lwBlSequence_(bl, bllist, nlwblcks);
}

// ----------------------------------------------------------------------------

// Construction of Hilbert block sequence
// see Greg Breinholt, Christoph Schierz
// Algorithm 781: generating Hilbert's space-filling curve by recursion
// ACM Transactions on Mathematical Software (TOMS) archive
// Volume 24, Issue 2 (June 1998)
// Pages: 184 - 189, 1998

static void gen_HilbertBlSeq_(blcluster* bl, unsigned i1, unsigned i2,
                              blcluster**& BlList, unsigned& nbl)
{
  if (bl->isleaf()) {
    bl->setidx(nbl++);
    *BlList++ = bl;
  } else {
    assert(bl->getnrs()==2 && bl->getncs()==2);
    blcluster* son;

    gen_HilbertBlSeq_(bl->getson(i1, i1), i1, 1-i2, BlList, nbl);

    if ((son=bl->getson(i2, 1-i2)))
      gen_HilbertBlSeq_(son, i1, i2, BlList, nbl);

    gen_HilbertBlSeq_(bl->getson(1-i1, 1-i1), i1, i2, BlList, nbl);

    if ((son=bl->getson(1-i2, i2)))
      gen_HilbertBlSeq_(son, 1-i1,  i2, BlList, nbl);
  }
}

void gen_HilbertBlSeq(blcluster* bl, blcluster**& BlList)
{
  BlList = new blcluster*[bl->nleaves()];
  assert(BlList!=NULL);

  blcluster** bllist = BlList;
  unsigned nbl = 0;
  gen_HilbertBlSeq_(bl, 0, 0, bllist, nbl);
}

// The Z- or Norton curve
static void gen_NortonBlSeq_(blcluster* bl, blcluster**& BlList, unsigned& nbl)
{
  if (bl->isleaf()) {
    bl->setidx(nbl++);
    *BlList++ = bl;
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        if (son) gen_NortonBlSeq_(son, BlList, nbl);
      }
    }
  }
}

void gen_NortonBlSeq(blcluster* bl, blcluster**& BlList)
{
  BlList = new blcluster*[bl->nleaves()];
  assert(BlList!=NULL);

  blcluster** bllist = BlList;
  unsigned nbl = 0;
  gen_NortonBlSeq_(bl, bllist, nbl);
}
