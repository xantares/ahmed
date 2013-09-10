/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "basmod.h"

extern void gen_HilbertBlSeq(blcluster*, blcluster**&);
extern void gen_NortonBlSeq(blcluster*, blcluster**&);


// computes maximum cost of the minimal partition with complexity p(n-p)
// see
// Bjarne Olstad and Fredrik Manne
// Efficient Partitioning of Sequences
// IEEE Trans. Comp. 44 (11), 1995
//

static unsigned MaxCostMinPart_(unsigned n, unsigned* cost, unsigned p)
{
  if (p<2) {
    unsigned sum = 0;
    for (unsigned i=0; i<n; ++i) sum += cost[i];
    return sum;
  }

  unsigned j, *g0 = new unsigned[n], *g1 = new unsigned[n];

  g0[n-1] = cost[n-1];
  for (unsigned i=n-1; i>=p; i--) g0[i-1] = g0[i] + cost[i-1];

  for (unsigned k=2; k<=p; k++) {

    g1[n-k] = MAX(cost[n-k], g0[n-k+1]);
    j = n-k;
    unsigned fij = cost[n-k];

    for (unsigned i=n-k-1; i+k>=p; --i) {
      fij += cost[i];

      if (fij <= g0[j+1]) g1[i] = g0[j+1];
      else {
        if (cost[i] >= g0[i+1]) {
          fij = g1[i] = cost[i];
          j = i;
        } else {
          while (fij-cost[j] >= g0[j]) fij -= cost[j--];
          g1[i] = MIN(fij, g0[j]);
          if (g1[i] == g0[j]) fij -= cost[j--];
        }
      }
    }

    swap(g0, g1);
  }

  unsigned max = g0[0];
  delete [] g1;
  delete [] g0;
  return max;
}


// default cost functional for genBlSeqPart
static unsigned cost_default_(blcluster& bl)
{
  const unsigned kavg = 10; // average rank

  return bl.isadm() ? kavg*(bl.getn1()+bl.getn2()) : bl.getn1()*bl.getn2(); 
}

static unsigned cost_symm_default_(blcluster& bl)
{
  const unsigned kavg = 10; // average rank

  if (bl.isdbl()) {
    return (bl.getn1()+1)*bl.getn1()/2;
  } else
    return bl.isadm() ? kavg*(bl.getn1()+bl.getn2()) : bl.getn1()*bl.getn2();
}


void genBlSeqPart(blcluster* bl, unsigned nproc, blcluster**& BlList,
                  unsigned*& part, unsigned (*cost_fnct)(blcluster&))
{
  unsigned i;
  // generate block sequence
  gen_HilbertBlSeq(bl, BlList);
  //gen_NortonBlSeq(bl, BlList);

  // compute cost of each block
  unsigned nbl = bl->nleaves();
  unsigned* cost = new unsigned[nbl];
  if (cost_fnct==NULL) cost_fnct = cost_symm_default_;
  for (i=0; i<nbl; ++i) cost[i] = cost_fnct(*BlList[i]);

  const unsigned max = MaxCostMinPart_(nbl, cost, nproc);

  // now compute partition using the cost of the maximum interval
  part = new unsigned[nproc+1];
  part[0] = 0;

  unsigned idx = 1, sum = cost[0];
  i = 0;

  while (idx<nbl) {
    if (sum+cost[idx] <= max) sum += cost[idx];
    else {
      part[++i] = idx;
      sum = cost[idx];
    }
    ++idx;
  }

  part[nproc] = nbl;

  /*
  std::cout << "genBlSeqPart: ";
  for (unsigned i=0; i<nproc; ++i) {
    unsigned sum = 0;
    for (unsigned j=part[i]; j<part[i+1]; ++j) sum += cost[j];
    std::cout << sum << ' ';
  }
  std::cout << std::endl;
  */

  /*
  unsigned maxi = 0;
  for (unsigned i=0; i<nproc; i++) {
    unsigned sum = 0;
    for (unsigned j=part[i]; j<part[i+1]; ++j) sum += cost[j];
    if (sum>maxi) maxi = sum;
  }
  std::cout << maxi << ' ' << max << std::endl;
  */

  delete [] cost;
}

