/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef MATGEN_SQNTL
#define MATGEN_SQNTL

#include <cmath>
#include "bemblcluster.h"
#include "basmod.h"
#include "apprx.h"
#include "H.h"

// The class MATGEN_T has to provide three functions:
// void cmpblck(unsigned b1, unsigned n1, unsigned b2, unsigned n2, T* data)
//      stores the entries of the block defined
//      by b1, n1, b2, n2 (in permuted ordering) in data
// void cmpblcksym(unsigned b1, unsigned n1, T* data)
//      stores the upper part of symmetric blocks on the diagonal
//      columnwise in data
// T scale(unsigned b1, unsigned n1, unsigned b2, unsigned n2)
//      is the expected size of the entries in this block

template<class T,class T1,class T2, class MATGEN_T> static
void matgenGeH_sqntl_(MATGEN_T& MatGen, bemblcluster<T1,T2>* bl, bool recmpr,
		   double eps, unsigned rankmax, unsigned& i, unsigned nblcks,
                   mblock<T>** A)
{
  if (bl->isleaf()) {
    progressbar(std::cout, "Approximating ... ", i++, nblcks, 20, true);
    apprx_unsym(MatGen, A[bl->getidx()], bl, eps, rankmax);
   } else {

    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    bool all_sons = true;              // indicates whether all sons are leaves
    for (unsigned k=0; k<ns1; ++k) {
      for (unsigned l=0; l<ns2; ++l) {
        bemblcluster<T1,T2>* son = (bemblcluster<T1,T2>*) bl->getson(k, l);
        if (son) {
          matgenGeH_sqntl_(MatGen, son, recmpr, eps, rankmax, i, nblcks, A);
          if (son->isnleaf()) all_sons = false;
        } else all_sons = false;
      }
    }

    if (recmpr && all_sons && !bl->isdbl()) unify_sons(bl, A, eps, rankmax);
  }
}

template<class T,class T1, class MATGEN_T> static
void matgen_sqntl_sym_(MATGEN_T& MatGen, bemblcluster<T1,T1>* bl, bool recmpr,
                       double eps, unsigned rankmax, unsigned& i,
                       unsigned nblcks, mblock<T>** A, const bool& cmplx_sym)
{
  if (bl->isleaf()) {
    progressbar(std::cout, "Approximating ... ", i++, nblcks, 20, true);
    apprx_sym(MatGen, A[bl->getidx()], bl, eps, rankmax, cmplx_sym);
  } else {

    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    bool all_sons = true;              // indicates whether all sons are leaves
    for (unsigned k=0; k<ns1; ++k) {
      for (unsigned l=0; l<ns2; ++l) {
        bemblcluster<T1,T1>* son = (bemblcluster<T1,T1>*) bl->getson(k, l);
        if (son) {
          matgen_sqntl_sym_(MatGen, son, recmpr, eps, rankmax, i, nblcks, 
                            A, cmplx_sym);
          if (son->isnleaf()) all_sons = false;
        } else all_sons = false;
      }
    }

    if (recmpr && all_sons && !bl->isdbl()) unify_sons(bl, A, eps, rankmax);
  }
}


// let root be the root of a block cluster tree and A its mblocks
// the following function computes the H-approximant of the descendant bl

template<class T,class T1,class T2, class MATGEN_T>
void matgenGeH_sqntl(MATGEN_T& MatGen, bemblcluster<T1,T2>* root,
		  bemblcluster<T1,T2>* bl, bool recmpr, double eps,
		  unsigned rankmax, mblock<T>** A)
{
  unsigned nblcks = bl->nleaves();
  unsigned i = 0;

  matgenGeH_sqntl_(MatGen, bl, recmpr, eps, rankmax, i, nblcks, A);

  if (recmpr) {
    i = 0;
    fill_gaps((blcluster*)root, A, root->nleaves(), i);
  }
}


template<class T,class T1, class MATGEN_T>
void matgen_sym_sqntl(MATGEN_T& MatGen, bemblcluster<T1,T1>* root,
		      bemblcluster<T1,T1>* bl, bool recmpr, double eps,
		      unsigned rankmax, mblock<T>** A, 
                      const bool& cmplx_sym = false)
{
  unsigned nblcks = bl->nleaves();
  unsigned i = 0;

  matgen_sqntl_sym_(MatGen, bl, recmpr, eps, rankmax, i, nblcks, A, cmplx_sym);

  if (recmpr) {
    i = 0;
    fill_gaps((blcluster*)root, A, root->nleaves(), i);
  }
}

template<class T,class T1, class MATGEN_T>
void matgenHeH_sqntl(MATGEN_T& MatGen, bemblcluster<T1,T1>* root,
		     bemblcluster<T1,T1>* bl, bool recmpr, double eps,
		     unsigned rankmax, mblock<T>** A)
{
  matgen_sym_sqntl(MatGen, root, bl, recmpr, eps, rankmax, A, false);
}
template<class T,class T1, class MATGEN_T>
void matgenSyH_sqntl(MATGEN_T& MatGen, bemblcluster<T1,T1>* root,
		     bemblcluster<T1,T1>* bl, bool recmpr, double eps,
		     unsigned rankmax, mblock<T>** A)
{
  matgen_sym_sqntl(MatGen, root, bl, recmpr, eps, rankmax, A, true);
}


#endif
