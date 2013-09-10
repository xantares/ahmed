/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "parallel.h"
#include "blcluster.h"
#include "blas.h"

template<class T>
void setzero_ND_(unsigned begp, unsigned p, blcluster* bl, T* par)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || bl->isleaf()) blas::setzero(bl->getn1(), par);
  else {
    T* par1 = par + bl->getson(0,0)->getn1();
    T* par2 = par1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2)
      setzero_ND_(begp, p/2, bl->getson(0,0), par);
    else
      setzero_ND_(begp+p/2, p/2, bl->getson(1,1), par1);

    blas::setzero(bl->getson(2,2)->getn1(), par2);
  }
}

template<class T>
void copy_ND_(unsigned begp, unsigned p, blcluster* bl, T* source, T* dest)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || bl->isleaf()) blas::copy(bl->getn1(), source, dest);
  else {
    T* source1 = source + bl->getson(0,0)->getn1();
    T* dest1 = dest + bl->getson(0,0)->getn1();
    T* source2 = source1 + bl->getson(1,1)->getn1();
    T* dest2 = dest1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2)
      copy_ND_(begp, p/2, bl->getson(0,0), source, dest);
    else
      copy_ND_(begp+p/2, p/2, bl->getson(1,1), source1, dest1);

    blas::copy(bl->getson(2,2)->getn1(), source2, dest2);
  }
}

template<class T>
void axpy_ND_(unsigned begp, unsigned p, blcluster* bl, T a, T* x, T* y)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || bl->isleaf()) blas::axpy(bl->getn1(), a, x, y);
  else {
    T* x1 = x + bl->getson(0,0)->getn1();
    T* y1 = y + bl->getson(0,0)->getn1();
    T* x2 = x1 + bl->getson(1,1)->getn1();
    T* y2 = y1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2)
      axpy_ND_(begp, p/2, bl->getson(0,0), a, x, y);
    else
      axpy_ND_(begp+p/2, p/2, bl->getson(1,1), a, x1, y1);

    blas::axpy(bl->getson(2,2)->getn1(), a, x2, y2);
  }
}

template<class T>
T nrm2_ND_(unsigned begp, unsigned p, blcluster* bl, T* par)
{
  unsigned rank = COMM_AHMED.Get_rank();
  T temp = 0.0;
  if (p==1 || bl->isleaf()) temp += blas::nrm2(bl->getn1(), par);
  else {
    unsigned length2 = bl->getson(2,2)->getn1();
    T* par1 = par + bl->getson(0,0)->getn1();
    T* par2 = par1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2) {
      temp += nrm2_ND_(begp, p/2, bl->getson(0,0), par);
      if (rank == begp) {
        temp += blas::nrm2(length2, par2);
        T y;
        MPI_Recv(&y, 1, begp+p/2, 41);
        temp += y;
      }
    } else {
      temp += nrm2_ND_(begp+p/2, p/2, bl->getson(1,1), par1);
      if (rank == begp+p/2) MPI_Send(&temp, 1, begp, 41);
    }
    MPI::Group group_all = COMM_AHMED.Get_group();
    const int range1[3] = {begp, begp+p-1, 1};
    MPI::Group group1 = group_all.Range_incl(1, &range1);
    MPI::Intracomm intra1 (COMM_AHMED.Create(group1));
    broadcast_intracom(intra1, 0, &temp, 1);
  }
  return temp;
}

template<class T>
T scpr_ND_(unsigned begp, unsigned p, blcluster* bl, T* v, T* w)
{
  unsigned rank = COMM_AHMED.Get_rank();
  T temp = 0.0;
  if (p==1 || bl->isleaf()) temp += blas::scpr(bl->getn1(), v, w);
  else {
    unsigned length2 = bl->getson(2,2)->getn1();
    T* v1 = v + bl->getson(0,0)->getn1();
    T* w1 = w + bl->getson(0,0)->getn1();
    T* v2 = v1 + bl->getson(1,1)->getn1();
    T* w2 = w1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2) {
      temp += scpr_ND_(begp, p/2, bl->getson(0,0), v, w);
      if (rank == begp) {
        temp += blas::scpr(length2, v2, w2);
        T y;
        MPI_Recv(&y, 1, begp+p/2, 42);
        temp += y;
      }
    } else {
      temp += scpr_ND_(begp+p/2, p/2, bl->getson(1,1), v1, w1);
      if (rank == begp+p/2) MPI_Send(&temp, 1, begp, 42);
    }
    MPI::Group group_all = COMM_AHMED.Get_group();
    const int range1[3] = {begp, begp+p-1, 1};
    MPI::Group group1 = group_all.Range_incl(1, &range1);
    MPI::Intracomm intra1 (COMM_AHMED.Create(group1));
    broadcast_intracom(intra1, 0, &temp, 1);
  }
  return temp;
}

template<class T>
void scal_ND_(unsigned begp, unsigned p, blcluster* bl, T a, T* par)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1 || bl->isleaf()) blas::scal(bl->getn1(), a, par);
  else {
    T* par1 = par + bl->getson(0,0)->getn1();
    T* par2 = par1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2)
      scal_ND_(begp, p/2, bl->getson(0,0), a, par);
    else
      scal_ND_(begp+p/2, p/2, bl->getson(1,1), a, par1);

    blas::scal(bl->getson(2,2)->getn1(), a, par2);
  }
}

namespace blas
{
void setzero_ND(unsigned nproc, blcluster* bl, double* par)
{
  setzero_ND_(0, nproc, bl, par);
}
void setzero_ND(unsigned nproc,blcluster* bl, float* par)
{
  setzero_ND_(0, nproc, bl, par);
}

void copy_ND(unsigned nproc,blcluster* bl, double* source, double* dest)
{
  copy_ND_(0, nproc, bl, source, dest);
}
void copy_ND(unsigned nproc, blcluster* bl, float* source, float* dest)
{
  copy_ND_(0, nproc, bl, source, dest);
}

void axpy_ND(unsigned nproc,blcluster* bl, double a, double* x, double* y)
{
  axpy_ND_(0, nproc, bl, a, x, y);
}
void axpy_ND(unsigned nproc,blcluster* bl, float a, float* x, float* y)
{
  axpy_ND_(0, nproc, bl, a, x, y);
}

double nrm2_ND(unsigned nproc,blcluster* bl, double* par)
{
  return nrm2_ND_(0, nproc, bl, par);
}
float nrm2_ND(unsigned nproc,blcluster* bl, float* par)
{
  return nrm2_ND_(0, nproc, bl, par);
}

double scpr_ND(unsigned nproc,blcluster* bl, double* v, double* w)
{
  return scpr_ND_(0, nproc, bl, v, w);
}
float scpr_ND(unsigned nproc,blcluster* bl, float* v, float* w)
{
  return scpr_ND_(0, nproc, bl, v, w);
}

void scal_ND(unsigned nproc,blcluster* bl, double a, double* par)
{
  scal_ND_(0, nproc, bl, a, par);
}
void scal_ND(unsigned nproc,blcluster* bl, float a, float* par)
{
  scal_ND_(0, nproc, bl, a, par);
}
}
