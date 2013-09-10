/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef MATGEN_OMP
#define MATGEN_OMP

#include <cmath>
#include "bemblcluster.h"
#include "basmod.h"
#include "mblock.h"
#include "apprx.h"

#ifdef _OPENMP
#include <omp.h>
#include "bllist.h"
#endif

// The class MATGEN_T has to provide three functions:
// void cmpbl(unsigned b1, unsigned n1, unsigned b2, unsigned n2, T* data)
//      stores the entries of the block defined
//      by b1, n1, b2, n2 (in permuted ordering) in data
// void cmpblsym(unsigned b1, unsigned n1, T* data)
//      stores the upper part of symmetric blocks on the diagonal
//      columnwise in data
// T scale(unsigned b1, unsigned n1, unsigned b2, unsigned n2)
//      is the expected size of the entries in this block

template<class T,class T1,class T2, class MATGEN_T> static
void _thr(MATGEN_T& MatGen, unsigned& co, unsigned nblcks,
	  bemblcluster<T1,T2>* bl, double eps, unsigned rankmax, mblock<T>** A)
{
  char st[100];
  sprintf(st, "Approximating using %d threads ... ", omp_get_num_threads());
  progressbar(std::cout, st, co++, nblcks, 20, true);
  apprx_unsym(MatGen, A[bl->getidx()], bl, eps, rankmax);
}


template<class T,class T1, class MATGEN_T> static
void _thr_sym(MATGEN_T& MatGen, unsigned& co, unsigned nblcks,
	      bemblcluster<T1,T1>* bl, double eps, unsigned rankmax,
	      mblock<T>** A, const bool& cmplx_sym)
{
  char st[100];
  sprintf(st, "Approximating using %d threads ... ", omp_get_num_threads());
  progressbar(std::cout, st, co++, nblcks, 20, true);
  apprx_sym(MatGen, A[bl->getidx()], bl, eps, rankmax, cmplx_sym);
}


template<class T,class T1,class T2, class MATGEN_T>
void matgenGeH_omp(MATGEN_T& MatGen, unsigned nblcks, bemblcluster<T1,T2>* bl,
		   double eps, unsigned rankmax, mblock<T>** A)
{
  // generate the list of blocks
  blcluster** BlList;
  gen_BlSequence(bl, BlList);

  unsigned counter = 0;

#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<(int) nblcks; i++)
    _thr(MatGen, counter, nblcks, (bemblcluster<T1,T2>*)BlList[i], eps,
	 rankmax, A);

  delete [] BlList;
}

template<class T,class T1, class MATGEN_T>
void matgen_sym_omp(MATGEN_T& MatGen, unsigned nblcks, bemblcluster<T1,T1>* bl,
                    double eps, unsigned rankmax, mblock<T>** A, 
		    const bool& cmplx_sym = false)
{
  // generate the list of blocks
  blcluster** BlList;
  gen_BlSequence(bl, BlList);

  unsigned counter = 0;
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<(int) nblcks; i++)
    _thr_sym(MatGen, counter, nblcks, (bemblcluster<T1,T1>*)BlList[i],
             eps, rankmax, A, cmplx_sym);

  delete [] BlList;
}

template<class T,class T1, class MATGEN_T>
void matgenHeH_omp(MATGEN_T& MatGen, unsigned nblcks, bemblcluster<T1,T1>* bl,
		   double eps, unsigned rankmax, mblock<T>** A)
{
  matgen_sym_omp(MatGen, nblcks, bl, eps, rankmax, A, false);
}
template<class T,class T1, class MATGEN_T>
void matgenSyH_omp(MATGEN_T& MatGen, unsigned nblcks, bemblcluster<T1,T1>* bl,
		   double eps, unsigned rankmax, mblock<T>** A)
{
  matgen_sym_omp(MatGen, nblcks, bl, eps, rankmax, A, true);
}

#endif
