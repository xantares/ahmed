/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "mblock.h"
#include "blcluster.h"
#include "bllist.h"
#include <omp.h>

template<class T> static
bool mltaGeHVec_omp_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  return A[bl->getidx()]->mltaVec(d, x+bl->getb2(), y+bl->getb1());
}

template<class T> static
bool mltaGeHVec_omp(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  unsigned nblcks = bl->nleaves();
  bool changed = false;
  blcluster** BlList;
  gen_BlSequence(bl, BlList);

#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<(int) nblcks; i++)
    if (mltaGeHVec_omp_(d, BlList[i], A, x, y)) changed = true;

  return changed;
}

bool mltaGeHVec_omp(double d, blcluster* bl, mblock<double>** A, double* x,
		   double* y)
{
  return mltaGeHVec_omp_(d, bl, A, x, y);
}

bool mltaGeHVec_omp(float d, blcluster* bl, mblock<float>** A, float* x,
		   float* y)
{
  return mltaGeHVec_omp_(d, bl, A, x, y);
}

bool mltaGeHVec_omp(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
		   dcomp* y)
{
  return mltaGeHVec_omp_(d, bl, A, x, y);
}

bool mltaGeHVec_omp(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
		   scomp* y)
{
  return mltaGeHVec_omp_(d, bl, A, x, y);
}
