/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef MATGEN_MPI
#define MATGEN_MPI

#include "parallel.h"
#include <cmath>
#include "bemblcluster.h"
#include "basmod.h"
#include "mblock.h"
#include "apprx.h"

// The class MATGEN_T has to provide three functions:
// void cmpblck(unsigned b1, unsigned n1, unsigned b2, unsigned n2, T* data)
//      stores the entries of the block defined
//      by b1, n1, b2, n2 (in permuted ordering) in data
// void cmpblcksym(unsigned b1, unsigned n1, T* data)
//      stores the upper part of symmetric blocks on the diagonal
//      columnwise in data
// T scale(unsigned b1, unsigned n1, unsigned b2, unsigned n2)
//      is the expected size of the entries in this block

template<class T, class T1, class T2, class MATGEN_T>
void matgenGeH_mpi(MATGEN_T& MatGen, bemblcluster<T1,T2>** blList,
		   unsigned* seq_part, double eps, unsigned rankmax, 
		   mblock<T>** A)
{
  unsigned rank = COMM_AHMED.Get_rank();

  for (unsigned i=seq_part[rank]; i<seq_part[rank+1]; ++i) {
    bemblcluster<T1,T2>* bl = (bemblcluster<T1,T2>*)blList[i];
    unsigned idx = i - seq_part[rank];
    apprx_unsym(MatGen, A[idx], bl, eps, rankmax);
  }
}

template<class T, class T1, class MATGEN_T>
void matgen_sym_mpi(MATGEN_T& MatGen, bemblcluster<T1,T1>** blList,
		    unsigned* seq_part, double eps, unsigned rankmax,
		    mblock<T>** A, const bool& cmplx_sym = false)
{
  unsigned rank = COMM_AHMED.Get_rank();

  for (unsigned i=seq_part[rank]; i<seq_part[rank+1]; ++i) {
    bemblcluster<T1,T1>* bl = (bemblcluster<T1,T1>*)blList[i];
    unsigned idx = i - seq_part[rank];
    apprx_sym(MatGen, A[idx], bl, eps, rankmax, cmplx_sym);
  }
}

template<class T, class T1, class MATGEN_T>
void matgenHeH_mpi(MATGEN_T& MatGen, bemblcluster<T1,T1>** blList,
		   unsigned* seq_part, double eps, unsigned rankmax,
		   mblock<T>** A)
{
  matgen_sym_mpi(MatGen, blList, seq_part, eps, rankmax, A, false);
}
template<class T, class T1, class MATGEN_T>
void matgenSyH_mpi(MATGEN_T& MatGen, bemblcluster<T1,T1>** blList,
		   unsigned* seq_part, double eps, unsigned rankmax,
		   mblock<T>** A)
{
  matgen_sym_mpi(MatGen, blList, seq_part, eps, rankmax, A, true);
}

#endif
