/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "parallel.h"
#include "H.h"

// y += alpha CRS x   CRS(2x3)
template<class T> static
void mltaCRS1vec_ND_(MPI::Intracomm& comm, const T alpha,
		      const blcluster* const rootH,
                      const CRSblock<T>* const H, T*const x, T*const y)
{
  unsigned p = comm.Get_size();
  unsigned rank = comm.Get_rank();
  unsigned nsons = rootH->getnrs();

  if (p==1 || nsons<2)
    amuxCRS(rootH->getn2(), alpha, x, y, H->iA, H->jA, H->A);
  else {
    unsigned nrs = rootH->getncs();
    if (nrs > 1) {
      std::cout << "\nFehler in mltaCRS1vec_ND_" << std::endl;
      exit(1);
    }
    T *x1 = x + rootH->getson(0,0)->getn1();
    T *x2 = x1 + rootH->getson(1,0)->getn1();
    T **y_shift = new T*[nrs];
    unsigned* length = new unsigned[nrs];
    y_shift[0] = y;
    length[0] = rootH->getson(0,0)->getn2();

    for (unsigned i = 0; i+1 < nrs; i++) {
      y_shift[i+1] = y_shift[i] + length[i];
      length[i+1] = rootH->getson(0,i+1)->getn2();
    }

    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      for (unsigned i = 0; i < nrs; i++) {
        mltaCRS1vec_ND_(intra, alpha, rootH->getson(0,i), &H->sons[0],
			 x, y_shift[i]);
        if (rank==0) {
	  if (nsons==3) {
	    amuxCRS(rootH->getson(2,i)->getn2(), alpha, x2, y_shift[i],
		    H->sons[1].iA, H->sons[1].jA, H->sons[1].A);
	  }
          T* temp = new T[length[i]];
          MPI_Recv_intracom(comm, temp, length[i], p/2, 31);
          blas::add(length[i], temp, y_shift[i]);
          delete [] temp;
        }
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      for (unsigned i = 0; i < nrs; i++) {
        blas::setzero(length[i], y_shift[i]);
        mltaCRS1vec_ND_(intra, alpha, rootH->getson(1,i), &H->sons[0],
			 x1, y_shift[i]);
        if (rank==p/2) MPI_Send_intracom(comm, y_shift[i], length[i], 0, 31);
      }
      intra.Free();
    }
    broadcast_intracom(comm, 0, y, rootH->getn2());
    delete [] y_shift;
    delete [] length;

  }
}

// y += alpha CRS x   CRS(3xcols)
template<class T> static
void mltaCRS2vec_ND_(MPI::Intracomm& comm, const T alpha,
                      const blcluster*const rootH,
                      const CRSblock<T>*const H, T*const x, T*const y)
{
  unsigned p = comm.Get_size();
  unsigned rank = comm.Get_rank();
  unsigned nsons = rootH->getnrs();
  if (p==1 || nsons<2)
    amuxCRS(rootH->getn1(), alpha, x, y, H->iA, H->jA, H->A);
  else {
    unsigned ncs = rootH->getncs();
    T *y1 = y + rootH->getson(0,0)->getn1();
    T *y2 = y1 + rootH->getson(1,0)->getn1();
    T **x_shift = new T*[ncs];

    for (unsigned i=0; i<ncs; i++)
      x_shift[i] = x + rootH->getson(0,i)->getb2() - rootH->getb2();

    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      for (unsigned i=0; i<ncs; i++) {
        mltaCRS2vec_ND_(intra, alpha, rootH->getson(0,i), &H->sons[0],
                         x_shift[i], y);
      }
      if (rank==0 && nsons==3) {
        for (unsigned i=0; i<ncs; i++) {
          amuxCRS(rootH->getson(2,i)->getn1(), alpha, x_shift[i], y2,
                  H->sons[1].iA, H->sons[1].jA, H->sons[1].A);
        }
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      for (unsigned i=0; i<ncs; i++) {
        mltaCRS2vec_ND_(intra, alpha, rootH->getson(1,i), &H->sons[0],
                         x_shift[i], y1);
      }
      intra.Free();
    }
    if (nsons==3) {
      unsigned length = rootH->getson(2,0)->getn1();
      broadcast_intracom(comm, 0, y2, length);
    }
    delete [] x_shift;
  }
}

// y += alpha CRS x   CRS having ND structure
template<class T> static
void mltaCRSvec_ND_(MPI::Intracomm& comm, T alpha,
		     const blcluster* const rootH, const CRSblock<T>* const H,
		     T* const x, T* const y)
{
  unsigned p = comm.Get_size();
  unsigned rank = comm.Get_rank();
  if (p==1 || rootH->isleaf())
    amuxCRS(rootH->getn1(), alpha, x, y, H->iA, H->jA, H->A);
  else {
    T *x1 = x + rootH->getson(0,0)->getn2();
    T *y1 = y + rootH->getson(0,0)->getn1();
    T *x2 = x1 + rootH->getson(1,1)->getn2();
    T *y2 = y1 + rootH->getson(1,1)->getn1();

    unsigned nsons = rootH->getncs();
    if (rank<p/2) {
      MPI::Intracomm intra = comm.Split(0, rank);
      mltaCRSvec_ND_(intra, alpha, rootH->getson(0,0), &H->sons[0], x, y);
      if (nsons==3) {
	unsigned length = rootH->getson(2,2)->getn1();
      	mltaCRS1vec_ND_(intra, alpha, rootH->getson(0,2), &H->sons[2], x, y2);
      	mltaCRS2vec_ND_(intra, alpha, rootH->getson(0,2), &H->sons[1], x2, y);
			
	if (rank==0) {
	  amuxCRS(length, alpha, x2, y2,
		  H->sons[3].iA, H->sons[3].jA, H->sons[3].A);

	  T* temp = new T[length];
	  MPI_Recv_intracom(comm, temp, length, p/2, 30);
	  blas::add(length, temp, y2);
	  delete [] temp;
	}
	broadcast_intracom(comm, 0, y2, length);      
      }
      intra.Free();
    } else {
      MPI::Intracomm intra = comm.Split(1, rank);
      mltaCRSvec_ND_(intra, alpha, rootH->getson(1,1), &H->sons[0], x1, y1);
      if (nsons==3) {
	unsigned length = rootH->getson(2,2)->getn1();
	mltaCRS2vec_ND_(intra, alpha, rootH->getson(1,2), &H->sons[1], x2, y1);

	blas::setzero(length, y2);
	mltaCRS1vec_ND_(intra, alpha, rootH->getson(1,2), &H->sons[2], x1, y2);

	if (rank==p/2) MPI_Send_intracom(comm, y2, length, 0, 30);
	broadcast_intracom(comm, 0, y2, length);      
      }
      intra.Free();
    }
  }
}

void mltaCRSvec_ND(float alpha, const blcluster* const rootH,
                    const CRSblock<float>* const H, float* const x,
                    float* const y)
{
  mltaCRSvec_ND_(COMM_AHMED, alpha, rootH, H, x, y);
}

void mltaCRSvec_ND(double alpha, const blcluster* const rootH,
                    const CRSblock<double>* const H, double* const x,
                    double* const y)
{
  mltaCRSvec_ND_(COMM_AHMED, alpha, rootH, H, x, y);
}

void mltaCRSvec_ND(scomp alpha, const blcluster* const rootH,
                    const CRSblock<scomp>* const H, scomp* const x,
                    scomp* const y)
{
  mltaCRSvec_ND_(COMM_AHMED, alpha, rootH, H, x, y);
}

void mltaCRSvec_ND(dcomp alpha, const blcluster* const rootH,
                    const CRSblock<dcomp>* const H, dcomp* const x,
                    dcomp* const y)
{
  mltaCRSvec_ND_(COMM_AHMED, alpha, rootH, H, x, y);
}
