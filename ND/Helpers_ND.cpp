/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <cmath>
#include "parallel.h"
#include "blcluster.h"
#include "H.h"

//#define DEBUG

MPI::Intracomm COMM_AHMED;

template<class T, class S>
void convT2S_ND_(unsigned begp, unsigned p, blcluster* bl, T* x, S* xf)
{
  unsigned rank = COMM_AHMED.Get_rank();
  if (p==1||bl->isleaf())
    for (unsigned i = 0; i < bl->getn1(); i++) xf[i] = (S) x[i];
  else {
    unsigned length2 = bl->getson(2,2)->getn1();
    T* x1 = x + bl->getson(0,0)->getn1();
    T* x2 = x1 + bl->getson(1,1)->getn1();
    S* xf1 = xf + bl->getson(0,0)->getn1();
    S* xf2 = xf1 + bl->getson(1,1)->getn1();
    if (begp<=rank && rank<begp+p/2)
      convT2S_ND_(begp, p/2, bl->getson(0,0), x, xf);
    else
      convT2S_ND_(begp+p/2, p/2, bl->getson(1,1), x1, xf1);

    for (unsigned i = 0; i < length2; i++) xf2[i] = (S) x2[i];
  }
}

void convT2S_ND(unsigned nproc, blcluster* bl, double* x, float* xf)
{
  convT2S_ND_(0, nproc, bl, x, xf);
}
void convT2S_ND(unsigned nproc,blcluster* bl, float* xf, double* x)
{
  convT2S_ND_(0, nproc, bl, xf, x);
}
void convT2S_ND(unsigned nproc,blcluster* bl, dcomp* x, scomp* xf)
{
  convT2S_ND_(0, nproc, bl, x, xf);
}
void convT2S_ND(unsigned nproc,blcluster* bl, scomp* xf, dcomp* x)
{
  convT2S_ND_(0, nproc, bl, xf, x);
}

// gathers vector par at proc 0
template<class T>
void gatherVec_ND_(unsigned begp, unsigned p, blcluster* bl, T* x)
{
  unsigned rank = COMM_AHMED.Get_rank();
  
  if(!(p==1 || bl->isleaf())){
    blcluster* son0 = bl->getson(0,0);
    blcluster* son1 = bl->getson(1,1);
    T* x1 = x + son0->getn1();
    if (begp<=rank && rank<begp+p/2) {
      gatherVec_ND_(begp,p/2,son0,x);
      if(rank==begp){
	//std::cout<<rank<<"recv"<<std::endl;
	MPI_Recv(x1,son1->getn1(),begp+p/2,26);
      }
    } else {
      gatherVec_ND_(begp+p/2,p/2,son1,x1);
      if(rank==begp+p/2){
	//std::cout<<rank<<"send"<<std::endl;
	MPI_Send(x1,son1->getn1(),begp,26);
      }
    }
  }
}

void gatherVec_ND(unsigned nproc, blcluster* bl, double* par)
{
  gatherVec_ND_(0, nproc, bl, par);
}

void gatherVec_ND(unsigned nproc, blcluster* bl, float* par)
{
  gatherVec_ND_(0, nproc, bl, par);
}

template<class T>
double error_ND_(unsigned begp, unsigned p, T* seq, blcluster* bl, T* par)
{
  unsigned rank = COMM_AHMED.Get_rank();
  double error(0.0);

  if (p==1 || bl->isleaf()) {
    unsigned length = bl->getn1();
    for (unsigned i=0; i<length; i++) error += SQR(seq[i]-par[i]);
  } else {
    blcluster* son0 = bl->getson(0,0);
    blcluster* son1 = bl->getson(1,1);
    T* seq1 = seq + son0->getn1();
    T* par1 = par + son0->getn1();
    T* seq2 = seq1 + son1->getn1();
    T* par2 = par1 + son1->getn1();

    if (begp<=rank && rank<begp+p/2) {
      error = error_ND_(begp, p/2, seq, son0, par);
      //std::cout << "Error 1 " << sqrt(error) << ", p=" << p << std::endl;

      if (rank == begp) {
        unsigned length2 = bl->getson(2,2)->getn1();
        for (unsigned i = 0; i < length2; i++) error += SQR(seq2[i]-par2[i]);
        //std::cout << "Error 3 " << sqrt(error) << ", p=" << p << std::endl;
        double temp;
        MPI_Recv(&temp, 1, begp+p/2, 22);
        error += temp;
      }

    } else {
      error = error_ND_(begp+p/2, p/2, seq1, son1, par1);
      //std::cout << "Error 2 " << sqrt(error) << ", p=" << p << std::endl;
      if (rank == begp+p/2) MPI_Send(&error, 1, begp, 22);
    }
  }
  return error;
}

double error_ND(unsigned nproc, double* seq, blcluster* bl, double* par)
{
  return sqrt(error_ND_(0, nproc, seq, bl, par));
}
double error_ND(unsigned nproc, float* seq, blcluster* bl, float* par)
{
  return sqrt(error_ND_(0, nproc, seq, bl, par));
}

// sends block bl and adds it to A from proc r1 to proc r2
template<class T>
void send_and_addH_(blcluster* bl, mblock<T>** A, unsigned r1, unsigned r2,
                    double eps, unsigned rankmax)
{
  unsigned rank = COMM_AHMED.Get_rank();
  unsigned info[8];

  if (rank==r1 || rank==r2) {

    if (bl->isleaf()) {
      if (rank==r1) {
        mblock<T>* p = A[bl->getidx()];
        info[0] = p->getn1();
        info[1] = p->getn2();
        p->get_prop(info+2);
        info[7] = p->nvals();
        COMM_AHMED.Send(info, 8, MPI::UNSIGNED, r2, 7);

        if (info[7]) MPI_Send(p->getdata(), info[7], r2, 8);

      } else {
        COMM_AHMED.Recv(info, 8, MPI::UNSIGNED, r1, 7);
        T* tmp = NULL;
        if (info[7]) {
          tmp = new T[info[7]];
          MPI_Recv(tmp, info[7], r1, 8);
          mblock<T> mbl_temp(info[0], info[1]);
          mbl_temp.cpy_mbl(info+2, tmp);
          A[bl->getidx()]->addMbl(eps, rankmax, &mbl_temp);
          delete [] tmp;
        }
      }
    } else {
      unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
      for (unsigned i=0; i<ns1; ++i)
        for (unsigned j=0; j<ns2; ++j)
          send_and_addH(bl->getson(i, j), A, r1, r2, eps, rankmax);
    }
  }
}


// sends block bl and adds it to A from proc r1 to proc r2
template<class T>
void send_and_addHSym_(blcluster* bl, mblock<T>** A, unsigned r1,
                       unsigned r2, double eps, unsigned rankmax)
{
  unsigned rank = COMM_AHMED.Get_rank();
  unsigned info[8];

  if (rank==r1 || rank==r2) {

    if (bl->isleaf()) {
      if (rank==r1) {
        mblock<T>* p = A[bl->getidx()];
        info[0] = p->getn1();
        info[1] = p->getn2();
        p->get_prop(info+2);
        info[7] = p->nvals();
        COMM_AHMED.Send(info, 8, MPI::UNSIGNED, r2, 7);

        if (info[7]) MPI_Send(p->getdata(), info[7], r2, 8);

      } else {
	COMM_AHMED.Recv(info, 8, MPI::UNSIGNED, r1, 7);
	T* tmp = NULL;
	if (info[7]) {
	  tmp = new T[info[7]];
	  MPI_Recv(tmp, info[7], r1, 8);
	  mblock<T> mbl_temp(info[0], info[1]);
	  mbl_temp.cpy_mbl(info+2, tmp);
	  A[bl->getidx()]->addMbl(eps, rankmax, &mbl_temp);
	  delete [] tmp;
	}
      }
    } else {
      unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
      for (unsigned i=0; i<ns1; ++i) {
        send_and_addHSym(bl->getson(i, i), A, r1, r2, eps, rankmax);
        for (unsigned j=i+1; j<ns2; ++j)
          send_and_addH(bl->getson(i, j), A, r1, r2, eps, rankmax);
      }
    }
  }
}

void send_and_addH(blcluster* bl, mblock<double>** A, unsigned r1,
                   unsigned r2, double eps, unsigned rankmax)
{
  send_and_addH_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addH(blcluster* bl, mblock<float>** A, unsigned r1,
                   unsigned r2, double eps, unsigned rankmax)
{
  send_and_addH_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addH(blcluster* bl, mblock<dcomp>** A, unsigned r1,
                   unsigned r2, double eps, unsigned rankmax)
{
  send_and_addH_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addH(blcluster* bl, mblock<scomp>** A, unsigned r1,
                   unsigned r2, double eps, unsigned rankmax)
{
  send_and_addH_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addHSym(blcluster* bl, mblock<double>** A, unsigned r1,
                      unsigned r2, double eps, unsigned rankmax)
{
  send_and_addHSym_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addHSym(blcluster* bl, mblock<float>** A, unsigned r1,
                      unsigned r2, double eps, unsigned rankmax)
{
  send_and_addHSym_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addHSym(blcluster* bl, mblock<dcomp>** A, unsigned r1,
                      unsigned r2, double eps, unsigned rankmax)
{
  send_and_addHSym_(bl, A, r1, r2, eps, rankmax);
}

void send_and_addHSym(blcluster* bl, mblock<scomp>** A, unsigned r1,
                      unsigned r2, double eps, unsigned rankmax)
{
  send_and_addHSym_(bl, A, r1, r2, eps, rankmax);
}
