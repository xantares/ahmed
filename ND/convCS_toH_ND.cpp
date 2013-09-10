/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <mpi.h>
#include "blcluster.h"
#include "basmod.h"
#include "H.h"
#include "parallel.h"

template<class T>
unsigned long sizeH_3_2_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                            mblock<T>** A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  unsigned long size = 0;
  if (p == 1 || blclTree->isleaf()) size += sizeH(blclTree, A, 'H');
  else {
    unsigned cols = blclTree->getncs();
    if (begp<=rank && rank< begp+p/2) {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH0i = blclTree->getson(0,i);
	size += sizeH_3_2_ND_(begp, p/2, rootH0i, A);
	if (blclTree->getnrs()==3) {
	  blcluster* rootH2i = blclTree->getson(2,i);
	  size += sizeH(rootH2i, A, 'H');
	}
      }
    } else {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH1i = blclTree->getson(1,i);
        size += sizeH_3_2_ND_(begp+p/2, p/2, rootH1i, A);
	if (blclTree->getnrs()==3) {
	  blcluster* rootH2i = blclTree->getson(2,i);
	  size += sizeH(rootH2i, A, 'H');
	}
      }
    }
  }
  return size;
}

template<class T>
unsigned long sizeH_2_3_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                            mblock<T>** A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  unsigned long size = 0;
  if (p == 1 || blclTree->isleaf()) size += sizeH(blclTree, A, 'H');
  else {
    unsigned rows = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi0 = blclTree->getson(i,0);
        size += sizeH_2_3_ND_(begp, p/2, rootHi0, A);
	if (blclTree->getncs()==3) {
	  blcluster* rootHi2 = blclTree->getson(i,2);
	  size += sizeH(rootHi2, A, 'H');
	}
      }
    } else {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi1 = blclTree->getson(i,1);
        size += sizeH_2_3_ND_(begp+p/2, p/2, rootHi1, A);
	if (blclTree->getncs()==3) {
	  blcluster* rootHi2 = blclTree->getson(i,2);
	  size += sizeH(rootHi2, A, 'H');
	}
      }
    }
  }
  return size;
}

// 'U' for upper, 'L' for lower triangular matrix, 'S' for symmetric
// 'H' is default
template<class T>
unsigned long sizeH_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                        mblock<T>** A, char type)
{
  unsigned rank = COMM_AHMED.Get_rank();
  unsigned long size = 0;
  if (p == 1 || blclTree->isleaf()) size += sizeH(blclTree, A, type);
  else {
    if (type == 'U') {
      unsigned nsons = blclTree->getncs();
      if (begp<=rank && rank<begp+p/2) {
        blcluster* rootU00 = blclTree->getson(0,0);
        size += sizeH_ND_(begp, p/2, rootU00, A, type);
	if (nsons==3) {
	  blcluster* rootU02 = blclTree->getson(0,2);
	  size += sizeH_3_2_ND_(begp, p/2, rootU02, A);
	  if (rank==begp) {
	    blcluster* rootU22 = blclTree->getson(2,2);
	    size += sizeH(rootU22, A, type);
	  }
	}
      } else {
        blcluster* rootU11 = blclTree->getson(1,1);
	size += sizeH_ND_(begp+p/2, p/2, rootU11, A, type);
	if (nsons==3) {
	  blcluster* rootU12 = blclTree->getson(1,2);
	  size += sizeH_3_2_ND_(begp+p/2, p/2, rootU12, A);
	}
      }
    } else if (type == 'L') {
      unsigned nsons = blclTree->getncs();
      if (begp<=rank && rank<begp+p/2) {
        blcluster* rootL00 = blclTree->getson(0,0);
        size += sizeH_ND_(begp, p/2, rootL00, A, type);
	if (nsons==3) {
	  blcluster* rootL20 = blclTree->getson(2,0);
	  size += sizeH_2_3_ND_(begp, p/2, rootL20, A);
	  if (rank==begp) {
	    blcluster* rootL22 = blclTree->getson(2,2);
	    size += sizeH(rootL22, A, type);
	  }
	}
      } else {
        blcluster* rootL11 = blclTree->getson(1,1);
        size += sizeH_ND_(begp+p/2, p/2, rootL11, A, type);
	if (nsons==3) {
	  blcluster* rootL21 = blclTree->getson(2,1);
	  size += sizeH_2_3_ND_(begp+p/2, p/2, rootL21, A);
	}
      }
    } else if (type == 'S') {
      unsigned nsons = blclTree->getncs();
      if (begp<=rank && rank<begp+p/2) {
        blcluster* rootH00 = blclTree->getson(0, 0);
        size += sizeH_ND_(begp, p/2, rootH00, A, type);
	if (nsons==3) {
	  blcluster* rootH02 = blclTree->getson(0, 2);
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  size += sizeH_3_2_ND_(begp, p/2, rootH02, A);
	  size += sizeH(rootH22, A, type);
	}
      } else {
        blcluster* rootH11 = blclTree->getson(1, 1);
        size += sizeH_ND_(begp+p/2, p/2, rootH11, A, type);
	if (nsons==3) {
	  blcluster* rootH12 = blclTree->getson(1, 2);
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  size += sizeH_3_2_ND_(begp+p/2, p/2, rootH12, A);
	  size += sizeH(rootH22, A, type);
	}
      }
    } else {
      unsigned nsons = blclTree->getncs();
      if (begp<=rank && rank<begp+p/2) {
        blcluster* rootH00 = blclTree->getson(0, 0);
        size += sizeH_ND_(begp, p/2, rootH00, A, type);
	if (nsons==3) {
	  blcluster* rootH02 = blclTree->getson(0, 2);
	  blcluster* rootH20 = blclTree->getson(2, 0);
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  size += sizeH_3_2_ND_(begp, p/2, rootH02, A);
	  size += sizeH_2_3_ND_(begp, p/2, rootH20, A);
	  size += sizeH(rootH22, A, type);
	}
      } else {
        blcluster* rootH11 = blclTree->getson(1, 1);
        size += sizeH_ND_(begp+p/2, p/2, rootH11, A, type);
	if (nsons==3) {
	  blcluster* rootH12 = blclTree->getson(1, 2);
	  blcluster* rootH21 = blclTree->getson(2, 1);
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  size += sizeH_3_2_ND_(begp+p/2, p/2, rootH12, A);
	  size += sizeH_2_3_ND_(begp+p/2, p/2, rootH21, A);
	  size += sizeH(rootH22, A, type);
	}
      }
    }
  }
  return size;
}

// 'U' for upper, 'L' for lower triangular matrix, 'S' for symmetric
// 'H' is default
unsigned long sizeH_ND(unsigned nproc, blcluster* blclTree, mblock<double>** A,
                       char type)
{
  return sizeH_ND_(0, nproc, blclTree, A, type);
}

unsigned long sizeH_ND(unsigned nproc, blcluster* blclTree, mblock<float>** A,
                       char type)
{
  return sizeH_ND_(0, nproc, blclTree, A, type);
}

unsigned long sizeH_ND(unsigned nproc, blcluster* blclTree, mblock<dcomp>** A,
                       char type)
{
  return sizeH_ND_(0, nproc, blclTree, A, type);
}

unsigned long sizeH_ND(unsigned nproc, blcluster* blclTree, mblock<scomp>** A,
                       char type)
{
  return sizeH_ND_(0, nproc, blclTree, A, type);
}


/* CCS format: (indices start with 0 !!!)
   A   the non-zero entries
   iA   the row indices of the non-zero entries
   jA   the beginning indices of the columns in A and iA, (jA(N)=NNZ)
*/

/* CRS format: (indices start with 0 !!!)
   A   the non-zero entries
   jA   the column indices of the non-zero entries
   iA   the beginning indices of the rows in A and jA, (iA(N)=NNZ)
*/

// wurde die Matrix in der originalen Indizierung generiert, wird eine
// Umsortierung noetig. op_perm bildet die originalen Indizes auf die
// permutierten ab, po_perm ist die inverse Permutation.

// frees a 3x2 matrix
template<class T>
void freemblsH_3_2_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                       mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) freembls_recursive(blclTree, A);
  else {
    unsigned cols = blclTree->getncs();
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH0i = blclTree->getson(0, i);
        freemblsH_3_2_ND_(begp, p/2, rootH0i, A);
	if (nsons==3) {
	  blcluster* rootH2i = blclTree->getson(2,i);
	  freembls_recursive(rootH2i, A);
	}
      }
    } else {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH1i = blclTree->getson(1, i);
        freemblsH_3_2_ND_(begp+p/2, p/2, rootH1i, A);
	if (nsons==3) {
	  blcluster* rootH2i = blclTree->getson(2,i);
	  freembls_recursive(rootH2i, A);
	}
      }
    }
  }
}

// frees a 2x3 matrix
template<class T>
void freemblsH_2_3_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                       mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) freembls_recursive(blclTree, A);
  else {
    unsigned rows = blclTree->getnrs();
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi0 = blclTree->getson(i, 0);
        freemblsH_2_3_ND_(begp, p/2, rootHi0, A);
	if (nsons==3) {
	  blcluster* rootHi2 = blclTree->getson(i,2);
	  freembls_recursive(rootHi2, A);
	}
      }
    } else {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi1 = blclTree->getson(i, 1);
        freemblsH_2_3_ND_(begp+p/2, p/2, rootHi1, A);
	if (nsons==3) {
	  blcluster* rootHi2 = blclTree->getson(i,2);
	  freembls_recursive(rootHi2, A);
	}
      }
    }
  }
}

// frees a 3x3 matrix
template<class T>
void freemblsH_0_ND_(unsigned begp, unsigned p, blcluster* bl, mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p == 1 || bl->isleaf()) freembls_recursive(bl, A);
  else {
    unsigned nsons = bl->getnrs();
    if (begp <= rank && rank < begp + p/2) {
      blcluster* rootH00 = bl->getson(0, 0);
      freemblsH_0_ND_(begp, p/2, rootH00, A);
      if (nsons==3) {
	blcluster* rootH02 = bl->getson(0, 2);
	blcluster* rootH20 = bl->getson(2, 0);
	blcluster* rootH22 = bl->getson(2, 2);
	freemblsH_3_2_ND_(begp, p/2, rootH02, A);
	freemblsH_2_3_ND_(begp, p/2, rootH20, A);
      }
      freembls_recursive(bl, A);
    } else {
      blcluster* rootH11 = bl->getson(1, 1);
      freemblsH_0_ND_(begp+p/2, p/2, rootH11, A);
      if (nsons==3) {
	blcluster* rootH12 = bl->getson(1, 2);
	blcluster* rootH21 = bl->getson(2, 1);
	blcluster* rootH22 = bl->getson(2, 2);
	freemblsH_3_2_ND_(begp+p/2, p/2, rootH12, A);
	freemblsH_2_3_ND_(begp+p/2, p/2, rootH21, A);
      }
      freembls_recursive(bl, A);
    }
  }
}

// initializes a 3xcols matrix
template<class T>
void initH_3_2_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                   mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) initGeH_0_withoutAlloc(blclTree, A);
  else {
    unsigned cols = blclTree->getncs();
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH0i = blclTree->getson(0, i);
        initH_3_2_ND_(begp, p/2, rootH0i, A);
      }
      if (nsons == 3) {
	for (unsigned i=0; i<cols; ++i) {
	  blcluster* rootH2i = blclTree->getson(2,i);
	  initGeH_0_withoutAlloc(rootH2i, A);
	}
      }
    } else {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH1i = blclTree->getson(1, i);
        initH_3_2_ND_(begp+p/2, p/2, rootH1i, A);
      }
      if (nsons == 3) {
	for (unsigned i=0; i<cols; ++i) {
	  blcluster* rootH2i = blclTree->getson(2,i);
	  initGeH_0_withoutAlloc(rootH2i, A);
	}
      }
    }
  }
}

// initializes a rowsx3 matrix
template<class T>
void initH_2_3_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                   mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) initGeH_0_withoutAlloc(blclTree, A);
  else {
    unsigned rows = blclTree->getnrs();
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi0 = blclTree->getson(i, 0);
        initH_2_3_ND_(begp, p/2, rootHi0, A);
      }
      if (nsons == 3) {
	for (unsigned i = 0; i < rows; i++) {
	  blcluster* rootHi2 = blclTree->getson(i,2);
	  initGeH_0_withoutAlloc(rootHi2, A);
	}
      }
    } else {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi1 = blclTree->getson(i, 1);
        initH_2_3_ND_(begp+p/2, p/2, rootHi1, A);
      }
      if (nsons == 3) {
	for (unsigned i = 0; i < rows; i++) {
	  blcluster* rootHi2 = blclTree->getson(i,2);
	  initGeH_0_withoutAlloc(rootHi2, A);
	}
      }
    }
  }
}

// initializes a 3x3 matrix
template<class T>
void initGeH_0_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                 mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) initGeH_0_withoutAlloc(blclTree, A);
  else {
    unsigned nsons  = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      initGeH_0_ND_(begp, p/2, rootH00, A);
      if (nsons == 3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	blcluster* rootH20 = blclTree->getson(2, 0);
	blcluster* rootH22 = blclTree->getson(2, 2);
	initH_3_2_ND_(begp, p/2, rootH02, A);
	initH_2_3_ND_(begp, p/2, rootH20, A);
	initGeH_0_withoutAlloc(rootH22, A);
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      initGeH_0_ND_(begp+p/2, p/2, rootH11, A);
      if (nsons == 3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	blcluster* rootH21 = blclTree->getson(2, 1);
	blcluster* rootH22 = blclTree->getson(2, 2);
	initH_3_2_ND_(begp+p/2, p/2, rootH12, A);
	initH_2_3_ND_(begp+p/2, p/2, rootH21, A);
	initGeH_0_withoutAlloc(rootH22, A);
      }
    }
  }
}

void initGeH_0_ND(unsigned p, blcluster* blclTree, mblock<double>** &A)
{
  allocmbls(blclTree, A);
  initGeH_0_ND_(0, p, blclTree, A);
}

void initGeH_0_ND(unsigned p, blcluster* blclTree, mblock<float>** &A)
{
  allocmbls(blclTree, A);
  initGeH_0_ND_(0, p, blclTree, A);
}

void initGeH_0_ND(unsigned p, blcluster* blclTree, mblock<dcomp>** &A)
{
  allocmbls(blclTree, A);
  initGeH_0_ND_(0, p, blclTree, A);
}

void initGeH_0_ND(unsigned p, blcluster* blclTree, mblock<scomp>** &A)
{
  allocmbls(blclTree, A);
  initGeH_0_ND_(0, p, blclTree, A);
}


// initializes a 3x3 symmetric matrix
template<class T>
void initHeH_0_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                    mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) initHeH_0_withoutAlloc(blclTree, A);
  else {
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      initHeH_0_ND_(begp, p/2, rootH00, A);
      if (nsons == 3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	blcluster* rootH22 = blclTree->getson(2, 2);
	initH_3_2_ND_(begp, p/2, rootH02, A);
	initHeH_0_withoutAlloc(rootH22, A);
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      initHeH_0_ND_(begp+p/2, p/2, rootH11, A);
      if (nsons == 3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	blcluster* rootH22 = blclTree->getson(2, 2);
	initH_3_2_ND_(begp+p/2, p/2, rootH12, A);
	initHeH_0_withoutAlloc(rootH22, A);
      }
    }
  }
}

void initHeH_0_ND(unsigned p, blcluster* blclTree, mblock<double>** &A)
{
  allocmbls(blclTree, A);
  initHeH_0_ND_(0, p, blclTree, A);
}

void initHeH_0_ND(unsigned p, blcluster* blclTree, mblock<float>** &A)
{
  allocmbls(blclTree, A);
  initHeH_0_ND_(0, p, blclTree, A);
}


// initializes a lower 3x3 matrix
template<class T>
void initLtH_0_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                   mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) initLtH_0_withoutAlloc(blclTree, A);
  else {
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      initLtH_0_ND_(begp, p/2, rootH00, A);
      if (nsons == 3) {
	blcluster* rootH20 = blclTree->getson(2, 0);
	initH_2_3_ND_(begp, p/2, rootH20, A);
	if (rank == begp) {
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  initLtH_0_withoutAlloc(rootH22, A);
	}
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      initLtH_0_ND_(begp+p/2, p/2, rootH11, A);
      if (nsons == 3) {
	blcluster* rootH21 = blclTree->getson(2, 1);
	initH_2_3_ND_(begp+p/2, p/2, rootH21, A);
      }
    }
  }
}

void initLtH_0_ND(unsigned p, blcluster* blclTree, mblock<double>** &A)
{
  allocmbls(blclTree, A);
  initLtH_0_ND_(0, p, blclTree, A);
}

void initLtH_0_ND(unsigned p, blcluster* blclTree, mblock<float>** &A)
{
  allocmbls(blclTree, A);
  initLtH_0_ND_(0, p, blclTree, A);
}

// initializes an upper 3x3 matrix
template<class T>
void initUtH_0_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                   mblock<T>** &A)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) initUtH_0_withoutAlloc(blclTree, A);
  else {
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      initUtH_0_ND_(begp, p/2, rootH00, A);
      if (nsons == 3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	initH_3_2_ND_(begp, p/2, rootH02, A);
	if (rank == begp) {
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  initUtH_0_withoutAlloc(rootH22, A);
	}
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      initUtH_0_ND_(begp+p/2, p/2, rootH11, A);
      if (nsons == 3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	initH_3_2_ND_(begp+p/2, p/2, rootH12, A);
      }
    }
  }
}

void initUtH_0_ND(unsigned p, blcluster* blclTree, mblock<double>** &A)
{
  allocmbls(blclTree, A);
  initUtH_0_ND_(0, p, blclTree, A);
}

void initUtH_0_ND(unsigned p, blcluster* blclTree, mblock<float>** &A)
{
  allocmbls(blclTree, A);
  initUtH_0_ND_(0, p, blclTree, A);
}

//converts a 3x2 matrix
template<class T, class S>
void convCRS_toGeH3_2_ND_(unsigned begp, unsigned p, T* A,
                       unsigned* jA, unsigned* iA, unsigned* op_perm,
                       unsigned* po_perm, double eps, blcluster* blclTree,
                       mblock<S>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toGeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
  else {
    unsigned cols = blclTree->getncs();
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH0i = blclTree->getson(0, i);
        convCRS_toGeH3_2_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
                          rootH0i, AH);
      }
      if (rank == begp && nsons == 3) {
        for (unsigned i = 0; i < cols; i++) {
          blcluster* rootH2i = blclTree->getson(2,i);
          convCRS_toGeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps,
                                  rootH2i, AH);
        }
      }
    } else {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH1i = blclTree->getson(1, i);
        convCRS_toGeH3_2_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm, eps,
                          rootH1i, AH);
      }
    }
  }
}

//converts a 2x3 matrix
template<class T, class S>
void convCRS_toGeH2_3_ND_(unsigned begp, unsigned p, T* A,
                       unsigned* jA, unsigned* iA, unsigned* op_perm,
                       unsigned* po_perm, double eps, blcluster* blclTree,
                       mblock<S>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toGeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
  else {
    unsigned rows = blclTree->getnrs();
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi0 = blclTree->getson(i, 0);
        convCRS_toGeH2_3_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
                          rootHi0, AH);
      }
      if (rank == begp && nsons==3) {
        for (unsigned i = 0; i < rows; i++) {
          blcluster* rootHi2 = blclTree->getson(i,2);
          convCRS_toGeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps,
                                  rootHi2, AH);
        }
      }
    } else {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi1 = blclTree->getson(i, 1);
        convCRS_toGeH2_3_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm,
                          eps, rootHi1, AH);
      }
    }
  }
}

// converts a 3x3 matrix with permutation
template<class T, class S>
void convCRS_toGeH_ND_(unsigned begp, unsigned p, T* A,
                    unsigned* jA, unsigned* iA, unsigned* op_perm,
                    unsigned* po_perm, double eps, blcluster* blclTree,
                    mblock<S>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toGeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
  else {
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      convCRS_toGeH_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
                     rootH00, AH);
      if (nsons==3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	blcluster* rootH20 = blclTree->getson(2, 0);
	convCRS_toGeH3_2_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
			  rootH02, AH);
	convCRS_toGeH2_3_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
			  rootH20, AH);
	if (rank == begp)
	  convCRS_toGeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps,
				  blclTree->getson(2, 2), AH);
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      convCRS_toGeH_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm, eps,
                     rootH11, AH);
      if (nsons==3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	blcluster* rootH21 = blclTree->getson(2, 1);
	convCRS_toGeH3_2_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm, eps,
			  rootH12, AH);
	convCRS_toGeH2_3_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm,
			  eps, rootH21, AH);
      }
    }
  }
}

// converts a 3xcols matrix
template<class T, class S>
void convCRS_toGeH3_2_ND_(unsigned begp, unsigned p, T* A,
                       unsigned* jA, unsigned* iA, double eps,
                       blcluster* blclTree, mblock<S>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toGeH_withoutAlloc(A, jA, iA, eps, blclTree, AH);
  else {
    unsigned cols = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH0i = blclTree->getson(0, i);
        convCRS_toGeH3_2_ND_(begp, p/2, A, jA, iA, eps, rootH0i, AH);
      }
      unsigned nsons = blclTree->getnrs();
      if (rank == begp && nsons == 3) {
        for (unsigned i = 0; i < cols; i++) {
          blcluster* rootH2i = blclTree->getson(2,i);
          convCRS_toGeH_withoutAlloc(A, jA, iA, eps, rootH2i, AH);
        }
      }
    } else {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH1i = blclTree->getson(1, i);
        convCRS_toGeH3_2_ND_(begp+p/2, p/2, A, jA, iA, eps, rootH1i, AH);
      }
    }
  }
}

// converts a 2x3 matrix
template<class T, class S>
void convCRS_toGeH2_3_ND_(unsigned begp, unsigned p, T* A, unsigned* jA,
                       unsigned* iA, double eps, blcluster* blclTree,
                       mblock<S>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toGeH_withoutAlloc(A, jA, iA, eps, blclTree, AH);
  else {
    unsigned rows = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi0 = blclTree->getson(i, 0);
        convCRS_toGeH2_3_ND_(begp, p/2, A, jA, iA, eps, rootHi0, AH);
      }
      unsigned nsons = blclTree->getncs();
      if (rank == begp && nsons == 3) {
        for (unsigned i = 0; i < rows; i++) {
          blcluster* rootHi2 = blclTree->getson(i,2);
          convCRS_toGeH_withoutAlloc(A, jA, iA, eps, rootHi2, AH);
        }
      }
    } else {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi1 = blclTree->getson(i, 1);
        convCRS_toGeH2_3_ND_(begp+p/2, p/2, A, jA, iA, eps, rootHi1, AH);
      }
    }
  }
}

// converts a 3x3 matrix
template<class T, class S>
void convCRS_toGeH_ND_(unsigned begp, unsigned p, T* A, unsigned* jA,
                    unsigned* iA, double eps, blcluster* blclTree,
                    mblock<S>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toGeH_withoutAlloc(A, jA, iA, eps, blclTree, AH);
  else {
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      convCRS_toGeH_ND_(begp, p/2, A, jA, iA, eps, rootH00, AH);
      if (nsons==3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	blcluster* rootH20 = blclTree->getson(2, 0);
	convCRS_toGeH3_2_ND_(begp, p/2, A, jA, iA, eps, rootH02, AH);
	convCRS_toGeH2_3_ND_(begp, p/2, A, jA, iA, eps, rootH20, AH);
	if (rank == begp)
	  convCRS_toGeH_withoutAlloc(A, jA, iA, eps, blclTree->getson(2, 2), AH);
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      convCRS_toGeH_ND_(begp+p/2, p/2, A, jA, iA, eps, rootH11, AH);
      if (nsons==3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	blcluster* rootH21 = blclTree->getson(2, 1);
	convCRS_toGeH3_2_ND_(begp+p/2, p/2, A, jA, iA, eps, rootH12, AH);
	convCRS_toGeH2_3_ND_(begp+p/2, p/2, A, jA, iA, eps, rootH21, AH);
      }
    }
  }
}

template<class T> static
void addGeHId_ND_(unsigned begp, unsigned p,
                blcluster* blclTree, mblock<T>** AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) addGeHId(blclTree, AH);
  else {
    if (begp<=rank && rank<begp+p/2) {
      addGeHId_ND_(begp, p/2, blclTree->getson(0, 0), AH);
      unsigned nsons = blclTree->getncs();
      if (rank == begp && nsons == 3) addGeHId(blclTree->getson(2, 2), AH);
    } else
      addGeHId_ND_(begp+p/2, p/2, blclTree->getson(1, 1), AH);
  }
}

void addGeHId_ND(unsigned p, blcluster* blclTree, mblock<double>** AH)
{
  addGeHId_ND_(0, p, blclTree, AH);
}

void addGeHId_ND(unsigned p, blcluster* blclTree, mblock<float>** AH)
{
  addGeHId_ND_(0, p, blclTree, AH);
}

void addGeHId_ND(unsigned p, blcluster* blclTree, mblock<dcomp>** AH)
{
  addGeHId_ND_(0, p, blclTree, AH);
}

void addGeHId_ND(unsigned p, blcluster* blclTree, mblock<scomp>** AH)
{
  addGeHId_ND_(0, p, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, double* A,
                   unsigned* jA, unsigned* iA, unsigned* op_perm,
                   unsigned* po_perm, double eps, blcluster* blclTree,
                   mblock<double>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, double* A,
                   unsigned* jA, unsigned* iA, unsigned* op_perm,
                   unsigned* po_perm, double eps,blcluster* blclTree,
                   mblock<float>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, dcomp* A,
                   unsigned* jA, unsigned* iA, unsigned* op_perm,
                   unsigned* po_perm, double eps, blcluster* blclTree,
                   mblock<dcomp>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, dcomp* A,
                   unsigned* jA, unsigned* iA, unsigned* op_perm,
                   unsigned* po_perm, double eps,blcluster* blclTree,
                   mblock<scomp>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, double* A, unsigned* jA, unsigned* iA,
                   double eps, blcluster* blclTree, mblock<double>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, double* A, unsigned* jA, unsigned* iA,
                   double eps, blcluster* blclTree, mblock<float>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, dcomp* A, unsigned* jA, unsigned* iA,
                   double eps, blcluster* blclTree, mblock<dcomp>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, eps, blclTree, AH);
}

void convCRS_toGeH_ND(unsigned p, dcomp* A, unsigned* jA, unsigned* iA,
                   double eps, blcluster* blclTree, mblock<scomp>**& AH)
{
  initGeH_0_ND(p, blclTree, AH);
  convCRS_toGeH_ND_(0, p, A, jA, iA, eps, blclTree, AH);
}

// converts a 3x3 symmetric matrix
template <class T>
void convCRS_toHeH_ND_(unsigned begp, unsigned p, double* A,
                       unsigned* jA, unsigned* iA, unsigned* op_perm,
                       unsigned* po_perm, double eps, blcluster* blclTree,
                       mblock<T>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toHeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps,
                               blclTree, AH);
  else {
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      convCRS_toHeH_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
                        rootH00, AH);
      if (nsons == 3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	convCRS_toGeH3_2_ND_(begp, p/2, A, jA, iA, op_perm, po_perm, eps,
			  rootH02, AH);
	if (rank == begp)
	  convCRS_toHeH_withoutAlloc(A, jA, iA, op_perm, po_perm, eps,
				     blclTree->getson(2, 2), AH);
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      convCRS_toHeH_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm, eps,
                        rootH11, AH);
      if (nsons == 3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	convCRS_toGeH3_2_ND_(begp+p/2, p/2, A, jA, iA, op_perm, po_perm, eps,
			  rootH12, AH);
      }
    }
  }
}

void convCRS_toHeH_ND(unsigned p, double* A,
                      unsigned* jA, unsigned* iA, unsigned* op_perm,
                      unsigned* po_perm, double eps, blcluster* blclTree,
                      mblock<double>**& AH)
{
  initHeH_0_ND(p, blclTree, AH);
  convCRS_toHeH_ND_(0, p, A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
}

void convCRS_toHeH_ND(unsigned p, double* A,
                      unsigned* jA, unsigned* iA, unsigned* op_perm,
                      unsigned* po_perm, double eps, blcluster* blclTree,
                      mblock<float>**& AH)
{
  initHeH_0_ND(p, blclTree, AH);
  convCRS_toHeH_ND_(0, p, A, jA, iA, op_perm, po_perm, eps, blclTree, AH);
}

// converts a 3x3 symmetric matrix
template <class T>
void convCRS_toHeH_ND_(unsigned begp, unsigned p, double* A, unsigned* jA,
                       unsigned* iA, double eps, blcluster* blclTree,
                       mblock<T>**& AH)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf())
    convCRS_toHeH_withoutAlloc(A, jA, iA, eps, blclTree, AH);
  else {
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      convCRS_toHeH_ND_(begp, p/2, A, jA, iA, eps, rootH00, AH);
      if (nsons == 3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	convCRS_toGeH3_2_ND_(begp, p/2, A, jA, iA, eps, rootH02, AH);
	if (rank == begp)
	  convCRS_toHeH_withoutAlloc(A, jA, iA, eps, blclTree->getson(2,2), AH);
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      convCRS_toHeH_ND_(begp+p/2, p/2, A, jA, iA, eps, rootH11, AH);
      if (nsons == 3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	convCRS_toGeH3_2_ND_(begp+p/2, p/2, A, jA, iA, eps, rootH12, AH);
      }
    }
  }
}

void convCRS_toHeH_ND(unsigned p, double* A, unsigned* jA, unsigned* iA,
                      double eps, blcluster* blclTree, mblock<double>**& AH)
{
  initHeH_0_ND(p, blclTree, AH);
  convCRS_toHeH_ND_(0, p, A, jA, iA, eps, blclTree, AH);
}

void convCRS_toHeH_ND(unsigned p, double* A, unsigned* jA, unsigned* iA,
                      double eps, blcluster* blclTree, mblock<float>**& AH)
{
  initHeH_0_ND(p, blclTree, AH);
  convCRS_toHeH_ND_(0, p, A, jA, iA, eps, blclTree, AH);
}


// copys a H-Matrix
template<class T>
void copyH3_2_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                  mblock<T>** A, mblock<T>** B)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) copyH(blclTree, A, B);
  else {
    unsigned cols = blclTree->getncs();
    unsigned nsons = blclTree->getnrs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH0i = blclTree->getson(0, i);
        copyH3_2_ND_(begp, p/2, rootH0i, A, B);
        if (rank == begp && nsons == 3) {
          blcluster* rootH2i = blclTree->getson(2, i);
          copyH(rootH2i, A, B);
        }
      }
    } else {
      for (unsigned i = 0; i < cols; i++) {
        blcluster* rootH1i = blclTree->getson(1, i);
        copyH3_2_ND_(begp+p/2, p/2, rootH1i, A, B);
      }
    }
  }
}

template<class T>
void copyH2_3_ND_(unsigned begp, unsigned p, blcluster* blclTree,
                  mblock<T>** A, mblock<T>** B)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) copyH(blclTree, A, B);
  else {
    unsigned rows = blclTree->getnrs();
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi0 = blclTree->getson(i, 0);
        copyH2_3_ND_(begp, p/2, rootHi0, A, B);
        if (rank == begp && nsons == 3) {
          blcluster* rootHi2 = blclTree->getson(i, 2);
          copyH(rootHi2, A, B);
        }
      }
    } else {
      for (unsigned i = 0; i < rows; i++) {
        blcluster* rootHi1 = blclTree->getson(i, 1);
        copyH2_3_ND_(begp+p/2, p/2, rootHi1, A, B);
      }
    }
  }
}

template<class T>
void copyH_ND_(unsigned begp, unsigned p, blcluster* blclTree,
               mblock<T>** A, mblock<T>** B)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || blclTree->isleaf()) copyH(blclTree, A, B);
  else {
    unsigned nsons = blclTree->getncs();
    if (begp<=rank && rank<begp+p/2) {
      blcluster* rootH00 = blclTree->getson(0, 0);
      copyH_ND_(begp, p/2, rootH00, A, B);
      if (nsons==3) {
	blcluster* rootH02 = blclTree->getson(0, 2);
	blcluster* rootH20 = blclTree->getson(2, 0);
	copyH3_2_ND_(begp, p/2, rootH02, A, B);
	copyH2_3_ND_(begp, p/2, rootH20, A, B);
	if (rank == begp) {
	  blcluster* rootH22 = blclTree->getson(2, 2);
	  copyH(rootH22, A, B);
	}
      }
    } else {
      blcluster* rootH11 = blclTree->getson(1, 1);
      copyH_ND_(begp+p/2, p/2, rootH11, A, B);
      if (nsons==3) {
	blcluster* rootH12 = blclTree->getson(1, 2);
	blcluster* rootH21 = blclTree->getson(2, 1);
	copyH3_2_ND_(begp+p/2, p/2, rootH12, A, B);
	copyH2_3_ND_(begp+p/2, p/2, rootH21, A, B);
      }
    }
  }
}

void copyH_ND(unsigned p, blcluster* blclTree, mblock<double>** A,
              mblock<double>** B)
{
  copyH_ND_(0, p, blclTree, A, B);
}

void copyH_ND(unsigned p, blcluster* blclTree, mblock<float>** A,
              mblock<float>** B)
{
  copyH_ND_(0, p, blclTree, A, B);
}
