/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "parallel.h"
#include "mblock.h"
#include "H.h"

// y += alphaHx   H(2x3)
template<class T> static
void mltaH1vec_ND_(MPI::Intracomm& comm, T alpha, blcluster* rootH,
                    mblock<T>** H, T* x, T* y)
{
  unsigned p = comm.Get_size();
  unsigned nsons = rootH->getncs();
  if (p==1 || nsons<2) mltaGeHVec(alpha, rootH, H, x, y);
  else {
    unsigned rank = comm.Get_rank();
    unsigned rows = rootH->getnrs();

    T *x_shift1 = x + rootH->getson(0,0)->getn2();
    T *x_shift2 = x + rootH->getson(0,0)->getn2()+rootH->getson(0,1)->getn2();
    T **y_shift = new T*[rows];

    unsigned* length = new unsigned[rows];
    for (unsigned i = 0; i < rows; i++) {
      y_shift[i] = y + rootH->getson(i,1)->getb1() - rootH->getb1();
      length[i] = rootH->getson(i,0)->getn1();
    }

    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      for (unsigned i = 0; i < rows; i++) {
        mltaH1vec_ND_(intra, alpha, rootH->getson(i,0), H, x, y_shift[i]);
        if (rank==0) {
	  if (nsons==3)
	    mltaGeHVec(alpha, rootH->getson(i,2), H, x_shift2, y_shift[i]);
          T* temp = new T [length[i]];
          MPI_Recv_intracom(comm, temp, length[i], p/2, 31);
          blas::add(length[i], temp, y_shift[i]);
          delete [] temp;
        }
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      for (unsigned i = 0; i < rows; i++) {
	blas::setzero(length[i], y_shift[i]);
        mltaH1vec_ND_(intra, alpha, rootH->getson(i,1), H, x_shift1,
                       y_shift[i]);
        if (rank==p/2)
          MPI_Send_intracom(comm, y_shift[i], length[i], 0, 31);
      }
      intra.Free();
    }
    broadcast_intracom(comm, 0, y, rootH->getn1());
    delete [] y_shift;
    delete [] length;
  }
}

// y += Hx   H(3xcols)
template<class T> static
void mltaH2vec_ND_(MPI::Intracomm& comm, T alpha, blcluster* rootH,
                    mblock<T>** H, T* x, T* y)
{
  unsigned p = comm.Get_size();
  unsigned nsons =  rootH->getnrs();
  if (p==1 || nsons<2) mltaGeHVec(alpha, rootH, H, x, y);
  else {
    unsigned rank = comm.Get_rank();
    unsigned col = rootH->getncs();

    T *y_shift1 = y + rootH->getson(0,0)->getn1();
    T *y_shift2 = y_shift1 + rootH->getson(1,0)->getn1();
    T **x_shift = new T*[col];
    for (unsigned i = 0; i < col; i++)
      x_shift[i] = x + rootH->getson(0,i)->getb2() - rootH->getb2();

    if (rank<p/2) { // I
      MPI::Intracomm intra = comm.Split(0, rank);
      for (unsigned i = 0; i < col; i++)
        mltaH2vec_ND_(intra, alpha, rootH->getson(0,i), H, x_shift[i], y);

      if (rank==0 && nsons==3) { // 1
        for (unsigned i = 0; i < col; i++)
          mltaGeHVec(alpha, rootH->getson(2,i), H, x_shift[i], y_shift2);
      }
      intra.Free();
    } else { // II
      MPI::Intracomm intra = comm.Split(1, rank);
      for (unsigned i = 0; i < col; i++)
        mltaH2vec_ND_(intra, alpha, rootH->getson(1,i), H, x_shift[i],
		       y_shift1);
      intra.Free();
    }
    if (nsons==3) {
      unsigned length = rootH->getson(2,0)->getn1();
      broadcast_intracom(comm, 0, y_shift2, length);
    }
    delete [] x_shift;
  }
}

// y += HTx   H(3x2)
template<class T> static
void mltaH1hvec_ND_(MPI::Intracomm& comm, T alpha, blcluster* rootH,
                     mblock<T>** H, T* x, T* y)
{
  unsigned p = comm.Get_size();
  unsigned nsons = rootH->getnrs();
  if (p==1 || nsons<2) mltaGeHhVec(alpha, rootH, H, x, y);
  else {
    unsigned rank = comm.Get_rank();
    unsigned col = rootH->getncs();

    T *x_shift1 = x + rootH->getson(0,0)->getn1();
    T *x_shift2 = x + rootH->getson(0,0)->getn1()+rootH->getson(1,0)->getn1();
    T **y_shift = new T*[col];

    unsigned* length = new unsigned[col];
    for (unsigned i = 0; i < col; i++) {
      y_shift[i] = y + rootH->getson(0,i)->getb2() - rootH->getb2();
      length[i] = rootH->getson(0,i)->getn2();
    }

    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      for (unsigned i = 0; i < col; i++) {
        mltaH1hvec_ND_(intra, alpha, rootH->getson(0,i), H, x, y_shift[i]);
        if (rank==0) {
	  if (nsons==3)
	    mltaGeHhVec(alpha, rootH->getson(2,i), H, x_shift2, y_shift[i]);
          T* temp = new T [length[i]];
          MPI_Recv_intracom(comm, temp, length[i], p/2, 31);
          blas::add(length[i], temp, y_shift[i]);
          delete [] temp;
        }
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      for (unsigned i = 0; i < col; i++) {
	blas::setzero(length[i], y_shift[i]);
        mltaH1hvec_ND_(intra, alpha, rootH->getson(1,i), H, x_shift1,
                        y_shift[i]);
        if (rank==p/2)
          MPI_Send_intracom(comm, y_shift[i], length[i], 0, 31);
      }
      intra.Free();
    }
    broadcast_intracom(comm, 0, y, rootH->getn2());
    delete [] y_shift;
    delete [] length;
  }
}

// y += H^T D x   H(3xcols)
template<class T> static
void mltaH1hDvec_ND_(MPI::Intracomm& comm, T d, mblock<T>** H, blcluster* blH,
		      blcluster* blD, int* piv, T* x, T* y)
{
  unsigned p = comm.Get_size();
  unsigned nsons = blH->getnrs();
  if (p==1 || nsons<2) mltaGeHhDiHVec(d, H, blH, blD, piv, x, y);
  else {
    unsigned rank = comm.Get_rank();
    unsigned col = blH->getncs();

    T* x_shift1 = x + blH->getson(0,0)->getn1();
    T* x_shift2 = x_shift1 + blH->getson(1,0)->getn1();
    T** y_shift = new T*[col];
    unsigned* length = new unsigned[col];
    for (unsigned i=0; i<col; ++i) {
      y_shift[i] = y + blH->getson(0,i)->getb2() - blH->getb2();
      length[i] = blH->getson(0,i)->getn2();
    }

    if (rank<p/2) { // I
      MPI::Intracomm intra = comm.Split(0, rank);
      blcluster* blD00 = blD->getson(0,0);
      for (unsigned i=0; i<col; ++i) {
        mltaH1hDvec_ND_(intra, d, H, blH->getson(0,i), blD00, piv,
			 x, y_shift[i]);
        if (rank==0) { // 1
	  if (nsons==3) {
	    mltaGeHhDiHVec(d, H, blH->getson(2,i), blD->getson(2,2), piv,
			x_shift2, y_shift[i]);
	  }
	  T* temp = new T [length[i]];
          MPI_Recv_intracom(comm, temp, length[i], p/2, 31);
          blas::add(length[i], temp, y_shift[i]);
          delete [] temp;
	}
      }
      intra.Free();
    } else { // II
      MPI::Intracomm intra = comm.Split(1, rank);
      blcluster* blD11 = blD->getson(1,1);
      for (unsigned i=0; i<col; ++i) {
	blas::setzero(length[i], y_shift[i]);
        mltaH1hDvec_ND_(intra, d, H, blH->getson(1,i), blD11, piv,
			 x_shift1, y_shift[i]);
        if (rank==p/2) MPI_Send_intracom(comm, y_shift[i], length[i], 0, 31);
      }
      intra.Free();
    }
    broadcast_intracom(comm, 0, y, blH->getn2());

    delete [] y_shift;
    delete [] length;
  }
}

void mltaH1vec_ND(MPI::Intracomm& comm, double alpha, blcluster* rootH,
                   mblock<double>** H, double* x, double* y)
{
  mltaH1vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1vec_ND(MPI::Intracomm& comm, float alpha, blcluster* rootH,
                   mblock<float>** H, float* x, float* y)
{
  mltaH1vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1vec_ND(MPI::Intracomm& comm, dcomp alpha, blcluster* rootH,
                   mblock<dcomp>** H, dcomp* x, dcomp* y)
{
  mltaH1vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1vec_ND(MPI::Intracomm& comm, scomp alpha, blcluster* rootH,
                   mblock<scomp>** H, scomp* x, scomp* y)
{
  mltaH1vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1hvec_ND(MPI::Intracomm& comm, double alpha, blcluster* rootH,
                    mblock<double>** H, double* x, double* y)
{
  mltaH1hvec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1hvec_ND(MPI::Intracomm& comm, float alpha, blcluster* rootH,
                    mblock<float>** H, float* x, float* y)
{
  mltaH1hvec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1hvec_ND(MPI::Intracomm& comm, dcomp alpha, blcluster* rootH,
                    mblock<dcomp>** H, dcomp* x, dcomp* y)
{
  mltaH1hvec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1hvec_ND(MPI::Intracomm& comm, scomp alpha, blcluster* rootH,
                    mblock<scomp>** H, scomp* x, scomp* y)
{
  mltaH1hvec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH1hDvec_ND(MPI::Intracomm& comm, double d, mblock<double>** H,
		     blcluster* blH, blcluster* blD, int* piv,
		     double* x, double* y)
{
  mltaH1hDvec_ND_(comm, d, H, blH, blD, piv, x, y);
}

void mltaH1hDvec_ND(MPI::Intracomm& comm, float d, mblock<float>** H,
		     blcluster* blH, blcluster* blD, int* piv,
		     float* x, float* y)
{
  mltaH1hDvec_ND_(comm, d, H, blH, blD, piv, x, y);
}

void mltaH1hDvec_ND(MPI::Intracomm& comm, dcomp d, mblock<dcomp>** H,
		     blcluster* blH, blcluster* blD, int* piv,
		     dcomp* x, dcomp* y)
{
  mltaH1hDvec_ND_(comm, d, H, blH, blD, piv, x, y);
}

void mltaH1hDvec_ND(MPI::Intracomm& comm, scomp d, mblock<scomp>** H,
		     blcluster* blH, blcluster* blD, int* piv,
		     scomp* x, scomp* y)
{
  mltaH1hDvec_ND_(comm, d, H, blH, blD, piv, x, y);
}

void mltaH2vec_ND(MPI::Intracomm& comm, double alpha, blcluster* rootH,
                   mblock<double>** H, double* x, double* y)
{
  mltaH2vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH2vec_ND(MPI::Intracomm& comm, float alpha, blcluster* rootH,
                   mblock<float>** H, float* x, float* y)
{
  mltaH2vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH2vec_ND(MPI::Intracomm& comm, dcomp alpha, blcluster* rootH,
                   mblock<dcomp>** H, dcomp* x, dcomp* y)
{
  mltaH2vec_ND_(comm, alpha, rootH, H, x, y);
}

void mltaH2vec_ND(MPI::Intracomm& comm, scomp alpha, blcluster* rootH,
                   mblock<scomp>** H, scomp* x, scomp* y)
{
  mltaH2vec_ND_(comm, alpha, rootH, H, x, y);
}
