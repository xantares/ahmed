/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "parallel.h"

// solve L x = b for x
// forward substitution, b is destroyed
template<class T> static
void LtHVec_solve_ND_(MPI::Intracomm& comm, blcluster* blL,
		      mblock<T>** L, T*& b)
{
  unsigned p = comm.Get_size();
  if (p==1 || blL->isleaf()) LtHVec_solve(blL, L, b);
  else {
    unsigned rank = comm.Get_rank();
    unsigned nsons = blL->getnrs();

    T *b_shift1 = b + blL->getson(0,0)->getn1();
    T *b_shift2 = b + blL->getson(0,0)->getn1()+blL->getson(1,1)->getn1();

    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      LtHVec_solve_ND_(intra, blL->getson(0, 0), L, b);
      if (nsons==3) {
	mltaH1vec_ND(intra, (T)-1.0, blL->getson(2, 0), L, b, b_shift2);
	if (rank==0) {
	  unsigned length2 = blL->getson(2,2)->getn1();
	  T *y = new T[length2];
	  MPI_Recv_intracom(comm, y, length2, p/2, 20);
	  blas::add(length2, y, b_shift2);
	  LtHVec_solve(blL->getson(2, 2), L, b_shift2);
	  delete [] y;
	}
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      LtHVec_solve_ND_(intra, blL->getson(1, 1), L, b_shift1);
      if (nsons==3) {
	unsigned length2 = blL->getson(2,2)->getn1();
	blas::setzero(length2, b_shift2);
	mltaH1vec_ND(intra, (T)-1.0, blL->getson(2, 1), L, b_shift1,
		      b_shift2);
	if (rank==p/2)
	  MPI_Send_intracom(comm, b_shift2, length2, 0, 20);
      }
      intra.Free();
    }
    if (nsons==3) {
      unsigned length2 = blL->getson(2,2)->getn1();
      broadcast_intracom(comm, 0, b_shift2, length2);
    }
  }
}

// solve U x = b for x
// forward substitution, b is destroyed
template<class T> static
void UtHVec_solve_ND_(MPI::Intracomm& comm, blcluster* blU, mblock<T>** U, T* b)
{
  unsigned p = comm.Get_size();
  if (p==1 || blU->isleaf()) UtHVec_solve(blU, U, b);
  else {
    unsigned rank = comm.Get_rank();
    unsigned nsons = blU->getnrs();

    T *b_shift1 = b + blU->getson(0,0)->getn1();
    T *b_shift2 = b + blU->getson(0,0)->getn1()+blU->getson(1,1)->getn1();
    if (nsons==3) {
      unsigned length2 = blU->getson(2,2)->getn1();
      if (rank==0) UtHVec_solve(blU->getson(2,2), U, b_shift2);
      broadcast_intracom(comm, 0, b_shift2, length2);
    }
    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      if (nsons==3)
	mltaH2vec_ND(intra, (T)-1.0, blU->getson(0, 2), U, b_shift2, b);
      UtHVec_solve_ND_(intra, blU->getson(0, 0), U, b);
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      if (nsons==3)
	mltaH2vec_ND(intra, (T)-1.0, blU->getson(1, 2), U, b_shift2, b_shift1);
      UtHVec_solve_ND_(intra, blU->getson(1, 1), U, b_shift1);
      intra.Free();
    }
  }
}

// solve U^H x = b for x
// forward substitution, b contains solution on exit
template<class T> static
void UtHhVec_solve_ND_(MPI::Intracomm& comm, blcluster* blU, mblock<T>** U, T* b)
{
  unsigned p = comm.Get_size();

  if (p==1 || blU->isleaf()) UtHhVec_solve(blU, U, b);
  else {
    unsigned rank = comm.Get_rank();
    unsigned nsons = blU->getnrs();

    T *b_shift1 = b + blU->getson(0,0)->getn2();
    T *b_shift2 = b + blU->getson(0,0)->getn2()+blU->getson(1,1)->getn2();

    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      UtHhVec_solve_ND_(intra, blU->getson(0, 0), U, b);
      if (nsons==3) {
	mltaH1hvec_ND(intra, (T)-1.0, blU->getson(0, 2), U, b, b_shift2);
	if (rank==0) {
	  unsigned length2 = blU->getson(2,2)->getn2();
	  T *y1 = new T[length2];
	  MPI_Recv_intracom(comm, y1, length2, p/2, 22);
	  blas::add(length2, y1, b_shift2);
	  UtHhVec_solve(blU->getson(2, 2), U, b_shift2);
	  delete [] y1;
	}
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      UtHhVec_solve_ND_(intra, blU->getson(1, 1), U, b_shift1);
      if (nsons==3) {
	unsigned length2 = blU->getson(2,2)->getn2();
	blas::setzero(length2, b_shift2);
	mltaH1hvec_ND(intra, (T)-1.0, blU->getson(1, 2), U, b_shift1, b_shift2);
	if (rank==p/2) MPI_Send_intracom(comm, b_shift2, length2, 0, 22);
      }
      intra.Free();
    }
    if (nsons==3) {
      unsigned length2 = blU->getson(2,2)->getn2();
      broadcast_intracom(comm, 0, b_shift2, length2);
    }
  }
}

// solve U^H D x = b for x
// forward substitution, b contains solution on exit
template<class T> static
void UtHhDVec_solve_ND_(MPI::Intracomm& comm, blcluster* blU, mblock<T>** U,
			int* piv, T* b)
{
  unsigned p = comm.Get_size();
  if (p==1 || blU->isleaf()) UtHhDVec_solve(blU, U, piv, b);
  else {
    unsigned rank = comm.Get_rank();
    unsigned nsons = blU->getnrs();

    T* b_shift1 = b + blU->getson(0,0)->getn2();
    T* b_shift2 = b_shift1 + blU->getson(1,1)->getn2();

    if (rank<p/2) { // I
      MPI::Intracomm intra = comm.Split(0, rank);
      UtHhDVec_solve_ND_(intra, blU->getson(0,0), U, piv, b);
      if (nsons==3) {
	unsigned length2 = blU->getson(2,2)->getn2();
	mltaH1hDvec_ND(intra, (T)-1.0, U, blU->getson(0,2),
			blU->getson(0,0), piv, b, b_shift2);
	if (rank==0) { // 1
	  T *y1 = new T[length2];
	  MPI_Recv_intracom(comm, y1, length2, p/2, 22);
	  blas::add(length2, y1, b_shift2);
	  UtHhDVec_solve(blU->getson(2,2), U, piv, b_shift2);
	  delete [] y1;
	}
      }
      intra.Free();
    } else { // II
      MPI::Intracomm intra = comm.Split(1, rank);
      UtHhDVec_solve_ND_(intra, blU->getson(1,1), U, piv, b_shift1);
      if (nsons==3) {
	unsigned length2 = blU->getson(2,2)->getn2();
	blas::setzero(length2, b_shift2);      
	mltaH1hDvec_ND(intra, (T)-1.0, U, blU->getson(1,2), blU->getson(1,1),
			piv, b_shift1, b_shift2);
	if (rank==p/2) MPI_Send_intracom(comm, b_shift2, length2, 0, 22);
      }
      intra.Free();
    }
    if (nsons==3) {
      unsigned length2 = blU->getson(2,2)->getn2();
      broadcast_intracom(comm, 0, b_shift2, length2);
    }
  }
}

 // solve U x = b for x
 // backward substitution, b contains solution on exit
template<class T> static
void UtHVec_solve_ND_(MPI::Intracomm& comm, blcluster* blU, mblock<T>** U,
		      int* piv, T* b)
{
  unsigned p = comm.Get_size();

  if (p==1 || blU->isleaf()) UtHVec_solve(blU, U, piv, b);
  else {
    unsigned rank = comm.Get_rank();
    unsigned nsons = blU->getnrs();

    T* b_shift1 = b + blU->getson(0,0)->getn1();
    T* b_shift2 = b_shift1 + blU->getson(1,1)->getn1();

    if (nsons==3) {
      unsigned length2 = blU->getson(2,2)->getn1();
      if (rank==0) UtHVec_solve(blU->getson(2,2), U, piv, b_shift2);  // 1
      broadcast_intracom(comm, 0, b_shift2, length2);
    }

    if (rank<p/2) { // I
      MPI::Intracomm intra = comm.Split(0, rank);
      if (nsons==3)
	mltaH2vec_ND(intra, (T)-1.0, blU->getson(0,2), U, b_shift2, b);
      UtHVec_solve_ND_(intra, blU->getson(0,0), U, piv, b);
      intra.Free();
    } else { // II
      MPI::Intracomm intra = comm.Split(1, rank);
      if (nsons==3)
	mltaH2vec_ND(intra, (T)-1.0, blU->getson(1,2), U, b_shift2, b_shift1);
      UtHVec_solve_ND_(intra, blU->getson(1,1), U, piv, b_shift1);
      intra.Free();
    }
  }
}


void LtHVec_solve_ND(blcluster* blL, mblock<double>** L, double*& b)
{
  LtHVec_solve_ND_(COMM_AHMED, blL, L, b);
}

void LtHVec_solve_ND(blcluster* blL, mblock<float>** L, float*& b)
{
  LtHVec_solve_ND_(COMM_AHMED, blL, L, b);
}

void LtHVec_solve_ND(blcluster* blL, mblock<dcomp>** L, dcomp*& b)
{
  LtHVec_solve_ND_(COMM_AHMED, blL, L, b);
}

void LtHVec_solve_ND(blcluster* blL, mblock<scomp>** L, scomp*& b)
{
  LtHVec_solve_ND_(COMM_AHMED, blL, L, b);
}


void UtHVec_solve_ND(blcluster* blU, mblock<double>** U, double*& b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<float>** U, float*& b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<dcomp>** U, dcomp*& b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<scomp>** U, scomp*& b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHhVec_solve_ND(blcluster* blU, mblock<double>** U, double*& b)
{
  UtHhVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHhVec_solve_ND(blcluster* blU, mblock<float>** U, float*& b)
{
  UtHhVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHhVec_solve_ND(blcluster* blU, mblock<dcomp>** U, dcomp*& b)
{
  UtHhVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHhVec_solve_ND(blcluster* blU, mblock<scomp>** U, scomp*& b)
{
  UtHhVec_solve_ND_(COMM_AHMED, blU, U, b);
}

void UtHhDVec_solve_ND(blcluster* blU, mblock<double>** U,  int* piv, double* b)
{
  UtHhDVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHhDVec_solve_ND(blcluster* blU, mblock<float>** U, int* piv, float* b)
{
  UtHhDVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHhDVec_solve_ND(blcluster* blU, mblock<dcomp>** U, int* piv, dcomp* b)
{
  UtHhDVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHhDVec_solve_ND(blcluster* blU, mblock<scomp>** U, int* piv, scomp* b)
{
  UtHhDVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<double>** U, int* piv, double* b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<float>** U,  int* piv, float* b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<dcomp>** U, int* piv, dcomp* b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

void UtHVec_solve_ND(blcluster* blU, mblock<scomp>** U, int* piv, scomp* b)
{
  UtHVec_solve_ND_(COMM_AHMED, blU, U, piv, b);
}

