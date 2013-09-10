/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef PARALLEL_H
#define PARALLEL_H

#include <mpi.h>
#include "bllist.h"
#include "matrix.h"
#include "H.h"

extern MPI::Intracomm COMM_AHMED;

inline void initAHMED(MPI::Intracomm comm)
{
  COMM_AHMED = comm.Clone();
  assert(COMM_AHMED != MPI::COMM_NULL);
}

extern void psoutH_MPI(std::ofstream&, blcluster**, unsigned*, unsigned,
		       mblock<double>**);
extern void transfH_MPI(mblock<double>**, mblock<float>**, unsigned*);
extern unsigned CG_MPI(const Matrix<double>&, double* const, double* const,
		       double&, unsigned&);
extern unsigned GMRes_MPI(const Matrix<double>&, double* const, double* const,
			  double&, const unsigned, unsigned&);
extern unsigned GMRes_MPI(const Matrix<dcomp>&, dcomp* const, dcomp* const,
			  double&, const unsigned, unsigned&);

extern void mltaGeHVec_MPI(double, mblock<double>**, double*, double*,
			  unsigned, unsigned*, blcluster**);
extern void mltaGeHVec_MPI(float, mblock<float>**, float*, float*,
			  unsigned, unsigned*, blcluster**);
extern void mltaGeHVec_MPI(dcomp, mblock<dcomp>**, dcomp*, dcomp*,
			  unsigned, unsigned*, blcluster**);
extern void mltaGeHVec_MPI(scomp, mblock<scomp>**, scomp*, scomp*,
			  unsigned, unsigned*, blcluster**);
extern void mltaHeHVec_MPI(double, mblock<double>**, double*, double*,
			     unsigned, unsigned*, blcluster**);

#define COUT(X) if (COMM_AHMED.Get_rank()==0) std::cout << X;

inline void MPI_Send(double* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Send(data, length, MPI::DOUBLE, i, id);
}

inline void MPI_Send(float* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Send(data, length, MPI::FLOAT, i, id);
}

inline void MPI_Send(dcomp* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Send((double*)data, 2*length, MPI::DOUBLE, i, id);
}

inline void MPI_Send(scomp* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Send((float*)data, 2*length, MPI::FLOAT, i, id);
}

inline void MPI_Send_intracom(MPI::Intracomm& intracomm, double* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Send(data, length, MPI::DOUBLE, i, id);
}

inline void MPI_Send_intracom(MPI::Intracomm& intracomm, float* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Send(data, length, MPI::FLOAT, i, id);
}

inline void MPI_Send_intracom(MPI::Intracomm& intracomm, dcomp* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Send((double*)data, 2*length, MPI::DOUBLE, i, id);
}

inline void MPI_Send_intracom(MPI::Intracomm& intracomm, scomp* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Send((float*)data, 2*length, MPI::FLOAT, i, id);
}

inline void MPI_Recv(double* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Recv(data, length, MPI::DOUBLE, i, id);
}

inline void MPI_Recv(float* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Recv(data, length, MPI::FLOAT, i, id);
}

inline void MPI_Recv(dcomp* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Recv((double*)data, 2*length, MPI::DOUBLE, i, id);
}

inline void MPI_Recv(scomp* data, unsigned length, unsigned i, unsigned id)
{
  COMM_AHMED.Recv((float*)data, 2*length, MPI::FLOAT, i, id);
}

inline void MPI_Recv_intracom(MPI::Intracomm& intracomm, double* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Recv(data, length, MPI::DOUBLE, i, id);
}

inline void MPI_Recv_intracom(MPI::Intracomm& intracomm, float* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Recv(data, length, MPI::FLOAT, i, id);
}

inline void MPI_Recv_intracom(MPI::Intracomm& intracomm, dcomp* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Recv((double*)data, 2*length, MPI::DOUBLE, i, id);
}

inline void MPI_Recv_intracom(MPI::Intracomm& intracomm, scomp* data,
			      unsigned length, unsigned i, unsigned id)
{
  intracomm.Recv((float*)data, 2*length, MPI::FLOAT, i, id);
}

inline void broadcast(unsigned root, double* x, unsigned length)
{
  COMM_AHMED.Bcast(x, length, MPI::DOUBLE, root);
}

inline void broadcast(unsigned root, float* x, unsigned length)
{
  COMM_AHMED.Bcast(x, length, MPI::FLOAT, root);
}

inline void broadcast(unsigned root, dcomp* x, unsigned length)
{
  COMM_AHMED.Bcast((double*)x, 2*length, MPI::DOUBLE, root);
}

inline void broadcast(unsigned root, scomp* x, unsigned length)
{
  COMM_AHMED.Bcast((float*)x, 2*length, MPI::FLOAT, root);
}

inline void broadcast_intracom(MPI::Intracomm& intracomm, unsigned root,
                               double* x, unsigned length)
{
  intracomm.Bcast(x, length, MPI::DOUBLE, root);
}

inline void broadcast_intracom(MPI::Intracomm& intracomm, unsigned root,
                               float* x, unsigned length)
{
  intracomm.Bcast(x, length, MPI::FLOAT, root);
}

inline void broadcast_intracom(MPI::Intracomm& intracomm, unsigned root,
                               dcomp* x, unsigned length)
{
  intracomm.Bcast((double*)x, 2*length, MPI::DOUBLE, root);
}

inline void broadcast_intracom(MPI::Intracomm& intracomm, unsigned root,
                               scomp* x, unsigned length)
{
  intracomm.Bcast((float*)x, 2*length, MPI::FLOAT, root);
}

extern unsigned long sizeH_ND(unsigned, blcluster*, mblock<double>**, char);
extern unsigned long sizeH_ND(unsigned, blcluster*, mblock<float>**, char);
extern unsigned long sizeH_ND(unsigned, blcluster*, mblock<dcomp>**, char);
extern unsigned long sizeH_ND(unsigned, blcluster*, mblock<scomp>**, char);
extern void gatherVec_ND(unsigned, blcluster*, double*);
extern void gatherVec_ND(unsigned, blcluster*, float*);
extern double error_ND(unsigned, double*, blcluster*, double*);
extern double error_ND(unsigned, float*, blcluster*, float*);
extern void initGeH_0_ND(unsigned, blcluster* , mblock<double>**&);
extern void initGeH_0_ND(unsigned, blcluster* , mblock<float>**&);
extern void initHeH_0_ND(unsigned, blcluster*, mblock<double>**&);
extern void initHeH_0_ND(unsigned, blcluster*, mblock<float>**&);
extern void initUtH_0_ND(unsigned, blcluster*, mblock<double>**&);
extern void initUtH_0_ND(unsigned, blcluster*, mblock<float>**&);
extern void initLtH_0_ND(unsigned, blcluster*, mblock<double>**&);
extern void initLtH_0_ND(unsigned, blcluster*, mblock<float>**&);
extern void copyH_ND(unsigned, blcluster*, mblock<double>**,
		     mblock<double>**);
extern void copyH_ND(unsigned, blcluster*, mblock<float>**,
		     mblock<float>**);
extern void convCRS_toGeH_ND(unsigned, double*, unsigned*,
			  unsigned*, unsigned*, unsigned*, double, blcluster*,
			  mblock<double>**&);
extern void convCRS_toGeH_ND(unsigned, double*, unsigned*,
			  unsigned*, unsigned*, unsigned*, double, blcluster*,
			  mblock<float>**&);
extern void convCRS_toGeH_ND(unsigned, dcomp*, unsigned*, unsigned*, unsigned*,
			  unsigned*, double, blcluster*, mblock<dcomp>**&);
extern void convCRS_toGeH_ND(unsigned, dcomp*, unsigned*, unsigned*, unsigned*,
			  unsigned*, double, blcluster*, mblock<scomp>**&);

extern void convCRS_toHeH_ND(unsigned, double*, unsigned*,
			     unsigned*, unsigned*, unsigned*, double,
			     blcluster*, mblock<double>**&);
extern void convCRS_toHeH_ND(unsigned, double*, unsigned*,
			     unsigned*, unsigned*, unsigned*, double,
			     blcluster*, mblock<float>**&);
extern void convCRS_toGeH_ND(unsigned, double*, unsigned*, unsigned*, double,
			  blcluster*, mblock<double>**&);
extern void convCRS_toGeH_ND(unsigned, double*, unsigned*, unsigned*, double,
			  blcluster*, mblock<float>**&);
extern void convCRS_toGeH_ND(unsigned, dcomp*, unsigned*, unsigned*, double,
			  blcluster*, mblock<dcomp>**&);
extern void convCRS_toGeH_ND(unsigned, dcomp*, unsigned*, unsigned*, double,
			  blcluster*, mblock<scomp>**&);
extern void convCRS_toHeH_ND(unsigned, double*, unsigned*, unsigned*, double,
			     blcluster*, mblock<double>**&);
extern void convCRS_toHeH_ND(unsigned, double*, unsigned*, unsigned*, double,
			     blcluster*, mblock<float>**&);

extern void freemblsH_0_ND(unsigned, blcluster*, mblock<double>**&);

void send_and_addH(blcluster*, mblock<double>**, unsigned, unsigned,
                   double, unsigned);
void send_and_addH(blcluster*, mblock<float>**, unsigned, unsigned,
                   double, unsigned);
void send_and_addH(blcluster*, mblock<dcomp>**, unsigned, unsigned,
                   double, unsigned);
void send_and_addH(blcluster*, mblock<scomp>**, unsigned, unsigned,
                   double, unsigned);
void send_and_addHSym(blcluster*, mblock<double>**, unsigned, unsigned,
                      double, unsigned);
void send_and_addHSym(blcluster*, mblock<float>**, unsigned, unsigned,
                      double, unsigned);
void send_and_addHSym(blcluster*, mblock<dcomp>**, unsigned, unsigned,
                      double, unsigned);
void send_and_addHSym(blcluster*, mblock<scomp>**, unsigned, unsigned,
                      double, unsigned);

extern void mltaH1vec_ND(MPI::Intracomm&, double, blcluster*,
			  mblock<double>**, double*, double*);
extern void mltaH1vec_ND(MPI::Intracomm&, float, blcluster*,
			  mblock<float>**, float*, float*);
extern void mltaH1vec_ND(MPI::Intracomm&, dcomp, blcluster*,
			  mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaH1vec_ND(MPI::Intracomm&, scomp, blcluster*,
			  mblock<scomp>**, scomp*, scomp*);
extern void mltaH1hvec_ND(MPI::Intracomm&, double, blcluster*,
			   mblock<double>**, double*, double*);
extern void mltaH1hvec_ND(MPI::Intracomm&, float, blcluster*,
			   mblock<float>**, float*, float*);
extern void mltaH1hvec_ND(MPI::Intracomm&, dcomp, blcluster*,
			   mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaH1hvec_ND(MPI::Intracomm&, scomp, blcluster*,
			   mblock<scomp>**, scomp*, scomp*);
extern void mltaH1hDvec_ND(MPI::Intracomm&, double, mblock<double>**,
			    blcluster*, blcluster*, int*, double*, double*);
extern void mltaH1hDvec_ND(MPI::Intracomm&, float, mblock<float>**,
			    blcluster*, blcluster*, int*, float*, float*);
extern void mltaH1hDvec_ND(MPI::Intracomm&, dcomp, mblock<dcomp>**,
			    blcluster*, blcluster*, int*, dcomp*, dcomp*);
extern void mltaH1hDvec_ND(MPI::Intracomm&, scomp, mblock<scomp>**,
			    blcluster*, blcluster*, int*, scomp*, scomp*);
extern void mltaH2vec_ND(MPI::Intracomm&, double, blcluster*,
			  mblock<double>**, double*, double*);
extern void mltaH2vec_ND(MPI::Intracomm&, float, blcluster*,
			  mblock<float>**, float*, float*);
extern void mltaH2vec_ND(MPI::Intracomm&, dcomp, blcluster*,
			  mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaH2vec_ND(MPI::Intracomm&, scomp, blcluster*,
			  mblock<scomp>**, scomp*, scomp*);

extern unsigned CG_ND(unsigned, blcluster*, const Matrix<double>&,
		      double* const, double* const, double&, unsigned&);

extern bool HLU_ND(unsigned, blcluster*, mblock<double>**, mblock<double>**,
		   mblock<double>**, double, unsigned);
extern bool HLU_ND(unsigned, blcluster*, mblock<float>**, mblock<float>**,
		   mblock<float>**, double, unsigned);
extern bool HLU_ND(unsigned, blcluster*, mblock<dcomp>**, mblock<dcomp>**,
		   mblock<dcomp>**, double, unsigned);
extern bool HLU_ND(unsigned, blcluster*, mblock<scomp>**, mblock<scomp>**,
		   mblock<scomp>**, double, unsigned);

extern bool HCholesky_ND(unsigned, blcluster*, mblock<double>**,
			 double, unsigned);
extern bool HCholesky_ND(unsigned, blcluster*, mblock<float>**,
			 double, unsigned);
extern bool HCholesky_ND(unsigned, blcluster*, mblock<dcomp>**,
			 double, unsigned);
extern bool HCholesky_ND(unsigned, blcluster*, mblock<scomp>**,
			 double, unsigned);

extern bool HUhDU_ND(unsigned, blcluster*, mblock<double>**, int*,
		     double, unsigned);
extern bool HUhDU_ND(unsigned, blcluster*, mblock<float>**, int*,
		     double, unsigned);
extern bool HUhDU_ND(unsigned, blcluster*, mblock<dcomp>**, int*,
		     double, unsigned);
extern bool HUhDU_ND(unsigned, blcluster*, mblock<scomp>**, int*,
		     double, unsigned);

extern void LtHVec_solve_ND(blcluster*, mblock<double>**, double*&);
extern void LtHVec_solve_ND(blcluster*, mblock<float>**, float*&);
extern void LtHVec_solve_ND(blcluster*, mblock<dcomp>**, dcomp*&);
extern void LtHVec_solve_ND(blcluster*, mblock<scomp>**, scomp*&);
extern void UtHVec_solve_ND(blcluster*, mblock<double>**, double*&);
extern void UtHVec_solve_ND(blcluster*, mblock<float>**, float*&);
extern void UtHVec_solve_ND(blcluster*, mblock<dcomp>**, dcomp*&);
extern void UtHVec_solve_ND(blcluster*, mblock<scomp>**, scomp*&);
extern void UtHhVec_solve_ND(blcluster*, mblock<double>**, double*&);
extern void UtHhVec_solve_ND(blcluster*, mblock<float>**, float*&);
extern void UtHhVec_solve_ND(blcluster*, mblock<dcomp>**, dcomp*&);
extern void UtHhVec_solve_ND(blcluster*, mblock<scomp>**, scomp*&);
extern void UtHVec_solve_ND(blcluster*, mblock<double>**, int*, double*);
extern void UtHVec_solve_ND(blcluster*, mblock<float>**, int*, float*);
extern void UtHVec_solve_ND(blcluster*, mblock<dcomp>**, int*, dcomp*);
extern void UtHVec_solve_ND(blcluster*, mblock<scomp>**, int*, scomp*);
extern void UtHhDVec_solve_ND(blcluster*, mblock<double>**, int*, double*);
extern void UtHhDVec_solve_ND(blcluster*, mblock<float>**, int*, float*);
extern void UtHhDVec_solve_ND(blcluster*, mblock<dcomp>**, int*, dcomp*);
extern void UtHhDVec_solve_ND(blcluster*, mblock<scomp>**, int*, scomp*);

extern void addGeHId_ND(unsigned, blcluster*, mblock<double>**);
extern void addGeHId_ND(unsigned, blcluster*, mblock<float>**);
extern void addGeHId_ND(unsigned, blcluster*, mblock<dcomp>**);
extern void addGeHId_ND(unsigned, blcluster*, mblock<scomp>**);

namespace blas
{
  extern void setzero_ND(unsigned, blcluster*, double*);
  extern void setzero_ND(unsigned, blcluster*, float*);

  extern void copy_ND(unsigned, blcluster*, double*, double*);
  extern void copy_ND(unsigned, blcluster*, float*, float*);

  extern void axpy_ND(unsigned, blcluster*, double, double*, double*);
  extern void axpy_ND(unsigned, blcluster*, float, float*, float*);

  extern double nrm2_ND(unsigned, blcluster*, double*);
  extern float nrm2_ND(unsigned, blcluster*, float*);

  extern double scpr_ND(unsigned, blcluster*, double*, double*);
  extern float scpr_ND(unsigned, blcluster*, float*, float*);

  extern void scal_ND(unsigned, blcluster*, double, double*);
  extern void scal_ND(unsigned, blcluster*, float, float*);
}

extern void subdivideCRS_ND(const unsigned, const blcluster* const,
			    CRSblock<float>*);
extern void subdivideCRS_ND(const unsigned , const blcluster* const,
			    CRSblock<double>*);
extern void subdivideCRS_ND(const unsigned, const blcluster* const,
			    CRSblock<dcomp>*);
extern void subdivideCRS_ND(const unsigned, const blcluster* const,
			    CRSblock<scomp>*);
extern void freeCRS_ND(const unsigned, const blcluster* const,
		       CRSblock<float>*);
extern void freeCRS_ND(const unsigned, const blcluster* const,
		       CRSblock<double>*);

extern void mltaCRSvec_ND(double, const blcluster* const,
			   const CRSblock<double>* const, double* const,
			   double* const);
extern void mltaCRSvec_ND(float, const blcluster* const,
			   const CRSblock<float>* const H, float* const,
			   float* const);
extern void mltaCRSvec_ND(scomp, const blcluster* const,
			   const CRSblock<scomp>* const H, scomp* const,
			   scomp* const);
extern void mltaCRSvec_ND(dcomp, const blcluster* const,
			   const CRSblock<dcomp>* const, dcomp* const,
			   dcomp* const);

extern void convT2S_ND(unsigned, blcluster*, double*, float*);
extern void convT2S_ND(unsigned, blcluster*, float*, double*);
extern void convT2S_ND(unsigned, blcluster*, dcomp*, scomp*);
extern void convT2S_ND(unsigned, blcluster*, scomp*, dcomp*);

// solve U^H U x = b for x and store x in b
template<class T>
inline void HCholesky_solve_ND(blcluster* bl, mblock<T>** U, T* b)
{
  UtHhVec_solve_ND(bl, U, b);  // Forward substitution
  UtHVec_solve_ND(bl, U, b);   // Backward substitution
}

// solve U^H D U x = b for x and store x in b
template<class T>
inline void HUhDU_solve_ND(blcluster* bl, mblock<T>** U, int* piv, T* b)
{
  UtHhDVec_solve_ND(bl, U, piv, b); // forward substitution
  UtHVec_solve_ND(bl, U, piv, b); // backward substitution
}

// solve L U x = b for x and store x in b
template<class T>
inline void HLU_solve_ND(blcluster* bl, mblock<T>** L, mblock<T>** U, T* b)
{
  LtHVec_solve_ND(bl, L, b);    // Forward substitution
  UtHVec_solve_ND(bl, U, b);    // Backward substitution
}

#endif
