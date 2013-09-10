/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef H_H
#define H_H

#include "blcluster.h"
#include "helper.h"
#include "preserveVec2.h"
#include <fstream>
#include <sstream>

template<class T> void allocmbls(unsigned n, mblock<T>** &A)
{
  A = new mblock<T>*[n];
  assert(A!=NULL);
  for (unsigned i=0; i<n; ++i) A[i] = NULL;
}

template<class T> void allocmbls(blcluster* bl, mblock<T>** &A)
{
  allocmbls(bl->nleaves(), A);
}

template<class T> void freembls(unsigned n, mblock<T>**& A)
{
  if (A!=NULL) {
    for (unsigned i=0; i<n; ++i) delete A[i];
    delete[] A;
    A = NULL;
  }
}

template<class T> void freembls(blcluster* bl, mblock<T>**& A)
{
  freembls(bl->nleaves(), A);
}

template<class T> static
void loadmbls_(const unsigned n, mblock<T>** &A, std::ifstream& is)
{
  bool flag;
  for (unsigned i=0; i<n; i++){
    is.read((char*) &flag, sizeof(bool));
    if(flag==true){
      A[i]=new mblock<T> (0,0);
      A[i]->load(is);
    }
  }
}

// loadmbls does not need allocmbls!!!
template<class T> void loadmbls(double& eta, unsigned& bmin, 
                                mblock<T>** &A, std::ifstream& is)
{
  if (!is) {
    std::cerr<<"InputStream not found"<<std::endl;
    exit(1);
  }
  is.read((char*) &eta, sizeof(double));
  is.read((char*) &bmin, sizeof(unsigned));
  unsigned n;
  is.read((char*) &n, sizeof(unsigned));
  allocmbls(n, A);
  loadmbls_(n, A, is);
}

template<class T> static
void savembls_(unsigned n, mblock<T>** &A, std::ofstream& os)
{
  assert(A!=NULL);
  bool one = true, zero = false;
  for (unsigned i=0; i<n; ++i) {
    if (A[i]!=NULL){
      os.write((char*) &one, sizeof(bool));
      A[i]->save(os);
    } else
      os.write((char*) &zero, sizeof(bool));
  }
}

template<class T> void savembls(const double eta, const unsigned bmin,
                                blcluster* bl, mblock<T>** &A, 
                                std::ofstream& os)
{
  if (!os) {
    std::cerr<<"OutputStream not found"<<std::endl;
    exit(1);
  }
  os.write((char*) &eta, sizeof(double));
  os.write((char*) &bmin, sizeof(unsigned));
  unsigned n = bl->nleaves();
  os.write((char*) &n, sizeof(unsigned));
  savembls_(n, A, os);
}

///////////////////////////////////////////////////////////////////////////
//
// double precision real
//

////initH.cpp:
extern void copyH(blcluster*, mblock<double>**, mblock<double>**);
extern void freembls_recursive(blcluster*, mblock<double>** &);
extern unsigned Hmax_rank(blcluster*, mblock<double>**, char type='H');
extern unsigned long sizeH(blcluster*, mblock<double>**, char type='H');
extern unsigned long sizeH(unsigned, mblock<double>**);
extern void initGeH_0(blcluster*, mblock<double>** &);
extern void initGeH_0_withoutAlloc(blcluster*, mblock<double>** &);
extern void initGeH_Id(blcluster*, mblock<double>** &);
extern void initHeH_0(blcluster*, mblock<double>** &);
extern void initHeH_0_withoutAlloc(blcluster*, mblock<double>** &);
extern void initHeH_Id(blcluster*, mblock<double>** &);
extern void initLtH_0(blcluster*, mblock<double>** &);
extern void initLtH_0_withoutAlloc(blcluster*, mblock<double>** &);
extern void initUtH_0(blcluster*, mblock<double>** &);
extern void initUtH_0_withoutAlloc(blcluster*, mblock<double>** &);
extern void setGeHzero(blcluster*, mblock<double>**);
extern void setHeHzero(blcluster*, mblock<double>**);
////convCS_toH.cpp:
extern void convCCS_toGeH(double*, unsigned*, unsigned*,
			  double, blcluster*, mblock<double>**&);
extern void convCCS_toGeH(double*, unsigned*, unsigned*, unsigned*,
			  unsigned*, double, blcluster*, mblock<double>**&);
extern void convCRS_toGeH(double*, unsigned*, unsigned*,
			  double, blcluster*, mblock<double>**&);
extern void convCRS_toGeH(double*, unsigned*, unsigned*, unsigned*,
			  unsigned*, blcluster*, mblock<double>**&);
extern void convCRS_toGeH(double*, unsigned*, unsigned*, unsigned*,
			  unsigned*, double, blcluster*, mblock<double>**&);
extern void convCRS_toGeH(dcomp*, unsigned*, unsigned*,
			  double, blcluster*, mblock<dcomp>**&);
extern void convCRS_toGeH(dcomp*, unsigned*, unsigned*,
			  double, blcluster*, mblock<scomp>**&);
extern void convCRS_toGeH_withoutAlloc(double*, unsigned*, unsigned*,
				       unsigned*, unsigned*, double, blcluster*,
				       mblock<double>**&);
extern void convCRS_toGeH_withoutAlloc(dcomp*, unsigned*, unsigned*,
				       unsigned*, unsigned*, double,
				       blcluster*, mblock<dcomp>**&);
extern void convCRS_toGeH_withoutAlloc(dcomp*, unsigned*, unsigned*,
				       unsigned*, unsigned*, double,
				       blcluster*, mblock<scomp>**&);
extern void convCRS_toHeH_withoutAlloc(double*, unsigned*, unsigned*,
                                       unsigned*, unsigned*, double,
                                       blcluster*, mblock<double>**&);
extern void convCRS_toGeH_withoutAlloc(double*, unsigned*, unsigned*, double,
				       blcluster*,  mblock<double>**&);
extern void convCRS_toHeH_withoutAlloc(double*, unsigned*, unsigned*, double,
                                       blcluster*, mblock<double>**&);
extern void convCRS_toGeH_withoutAlloc(dcomp*, unsigned*, unsigned*, double,
				       blcluster*, mblock<dcomp>**&);
extern void convCRS_toGeH_withoutAlloc(dcomp*, unsigned*, unsigned*, double,
				       blcluster*, mblock<scomp>**&);
extern void convCCS_toHeH(double*, unsigned*, unsigned*,
                          double, blcluster*, mblock<double>**&);
extern void convCCS_toHeH(double*, unsigned*, unsigned*, unsigned*,
                          unsigned*, double, blcluster*, mblock<double>**&);
extern void convCRS_toHeH(double*, unsigned*, unsigned*,
                          double, blcluster*, mblock<double>**&);
extern void convCRS_toHeH(double*, unsigned*, unsigned*, unsigned*,
                          unsigned*, double, blcluster*, mblock<double>**&);
////agglH.cpp:
extern void unify_sons(blcluster*, mblock<double>**, double, unsigned);
extern void fill_gaps(blcluster*, mblock<double>**, unsigned, unsigned&);
extern void agglH(blcluster*, mblock<double>**, double, unsigned,
                  contBasis<double>* haar=NULL);
extern void aggl_light_H(blcluster*, mblock<double>**, double, unsigned);
// extern void agglHSym_stab(blcluster*, mblock<double>**, double, unsigned);
////Htrunc.cpp:
extern void Htrunc_abs(blcluster*, mblock<double>**, double);
////nrmH.cpp:
extern double nrmF2GeH(blcluster*, mblock<double>**);
extern double nrmF2HeH(blcluster*, mblock<double>**);
extern double spctrlnrm_sym(blcluster*, mblock<double>**);
extern double spctrlnrm_sym(unsigned, blcluster*, mblock<double>**, blcluster*,
                            mblock<double>**);
extern double spctrlnrm(unsigned, blcluster*, mblock<double>**, blcluster*, 
			mblock<double>**);

////addHH.cpp:
extern void addGeHId(blcluster*, mblock<double>**);
extern void addGeHLrM(blcluster*, mblock<double>**, double, unsigned, unsigned,
		      double*, unsigned, double*, unsigned, 
		      contBasis<double>* haar=NULL);
extern void addGeHGeM(blcluster*, mblock<double>**, double*, unsigned,
		      double, unsigned, contBasis<double>* haar=NULL);
extern void addGeHGeH(blcluster* bl, mblock<double>** A, mblock<double>** B,
		      double eps, unsigned rankmax, 
		      contBasis<double>* haar=NULL);
extern void addHeHLrM(blcluster*, mblock<double>**, double, unsigned,
		      unsigned, double*, unsigned, double*, unsigned, 
		      contBasis<double>* haar=NULL);
extern void addHeHLrMExact(blcluster* bl, mblock<double>** A, unsigned k, 
			   double* U, unsigned ldU, double* V, unsigned ldV);
extern void addHeHLrM_stab(blcluster*, mblock<double>**, double, unsigned,
			   unsigned, double*, unsigned, double*, unsigned);
extern void addHeHGeM(blcluster*, mblock<double>**, double*, unsigned,
		      double, unsigned, contBasis<double>* haar=NULL);
////mltaGeHVec.cpp:
extern bool mltaGeHVec(double, blcluster*, mblock<double>**, double*, double*);
extern bool mltaGeHGeM(double, blcluster*, mblock<double>**, unsigned,
		       double*, unsigned, double*, unsigned);
extern bool mltaGeHhVec(double, blcluster*, mblock<double>**, double*, double*);
extern bool mltaGeHhGeM(double, blcluster*, mblock<double>**, unsigned,
			double*, unsigned, double*, unsigned);
extern bool mltaGeHhDiHVec(double, mblock<double>**, blcluster*, blcluster*, 
			   int*, double*, double*);
extern bool mltaGeHhDiHGeM(double, mblock<double>**, blcluster*, 
			   blcluster*, int*, unsigned, double*, unsigned, 
			   double*, unsigned);
extern void mltaLtHVec(double, blcluster*, mblock<double>**, double*, double*);
extern void mltaLtHGeM(double, blcluster*, mblock<double>**, unsigned,
		       double*, unsigned, double*, unsigned);
extern void mltaLtHhVec(double, blcluster*, mblock<double>**, double*, double*);
extern void mltaLtHhGeM(double, blcluster*, mblock<double>**, unsigned,
			double*, unsigned, double*, unsigned);
extern void mltaUtHVec(double, blcluster*, mblock<double>**, double*, double*);
extern void mltaUtHGeM(double, blcluster*, mblock<double>**, unsigned,
		       double*, unsigned, double*, unsigned);
extern void mltaUtHhVec(double, blcluster*, mblock<double>**, double*, double*);
extern void mltaUtHhGeM(double, blcluster*, mblock<double>**, unsigned,
                        double*, unsigned, double*, unsigned);
extern void mltaHeHVec(double, blcluster*, mblock<double>**, double*, double*);
extern void mltaHeHGeM(double, blcluster*, mblock<double>**, unsigned,
		       double*, unsigned, double*, unsigned);
////mltaGeHGeH.cpp:
extern void mltaGeHGeH(double, blcluster*, mblock<double>**, blcluster*, 
		       mblock<double>**, blcluster*, mblock<double>**, double, 
		       unsigned, contBasis<double>* haar=NULL);
extern void mltaGeHGeH_toMbl(double, blcluster*, mblock<double>**, blcluster*,
			     mblock<double>**, mblock<double>*, double, 
			     unsigned, contBasis<double>* haar=NULL);
////mltaGeHhGeH.cpp:
extern void mltaGeHhGeH(double, blcluster*, mblock<double>**, blcluster*, 
			mblock<double>**, blcluster*, mblock<double>**, double,
			unsigned, contBasis<double>* haar=NULL);
extern void mltaGeHhGeH_toMbl(double, blcluster*, mblock<double>**,
			      blcluster*, mblock<double>**, mblock<double>*, 
			      double, unsigned);
extern void mltaGeHhLrM_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double, 
                              unsigned, contBasis<double>* haar=NULL);
extern void mltaGeHhGeM_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double, 
                              unsigned, contBasis<double>* haar=NULL);
extern void mltaGeMhGeH_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double, 
                              unsigned, contBasis<double>* haar=NULL);
extern void mltaLrMhGeH_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double, 
                              unsigned, contBasis<double>* haar=NULL);
////mltaGeHhGeH_toHeH.cpp:
extern void mltaGeHhGeH_toHeH(double, blcluster*, mblock<double>**, blcluster*,
			      mblock<double>**, blcluster*, mblock<double>**, 
			      double, unsigned, contBasis<double>* haar=NULL);
////mltaGeHhDiHGeH.cpp:
extern void mltaGeHhDiHGeH(double, mblock<double>**, blcluster*, 
			   blcluster*, int*, blcluster*,
			   blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH(double, blcluster*, mblock<double>**, blcluster*,
			   mblock<double>**, int*, blcluster*, mblock<double>**,
			   blcluster*, mblock<double>**, double, unsigned);
extern void mltaGeHhDiHLrM_toMbl(double, blcluster*, mblock<double>**, 
				 blcluster*, mblock<double>**, int*, blcluster*,
				 mblock<double>**, mblock<double>*, double,
				 unsigned);
extern void mltaGeHhDiHGeM_toMbl(double, blcluster*, mblock<double>**, 
				 blcluster*, mblock<double>**, int*, blcluster*,
				 mblock<double>**, mblock<double>*, double,
				 unsigned);
extern void mltaGeMhDiHGeH_toMbl(double, blcluster*, mblock<double>**, 
				 blcluster*, mblock<double>**, int*, blcluster*,
				 mblock<double>**, mblock<double>*, double,
				 unsigned);
extern void mltaLrMhDiHGeH_toMbl(double, blcluster*, mblock<double>**, 
				 blcluster*, mblock<double>**, int*, blcluster*,
				 mblock<double>**, mblock<double>*, double,
				 unsigned);
extern void mltDiHGeM(unsigned, mblock<double>**, blcluster*, int* const,
		      double*, double*, unsigned);
////mltaGeHhDiHGeH_toHeH.cpp:
extern void mltaGeHhDiHGeH_toHeH(double, mblock<double>**, blcluster*, 
				 blcluster*, int*, blcluster*,
				 blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH_toHeH(double, blcluster*, mblock<double>**, 
				 blcluster*, mblock<double>**, int*, blcluster*,
				 mblock<double>**, blcluster*, mblock<double>**,
				 double, unsigned);
////mltaGeHGeHh.cpp:
extern void mltaGeHGeHh(double, blcluster*, mblock<double>**, blcluster*,
			mblock<double>**, blcluster*, mblock<double>**,
			double, unsigned);
extern void mltaGeHGeHh_toMbl(double, blcluster*, mblock<double>**,
			      blcluster*, mblock<double>**, mblock<double>*, 
			      double, unsigned);
extern void mltaGeHGeMh_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double,
			      unsigned);
extern void mltaGeMGeHh_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double,
			      unsigned);
extern void mltaLrMGeHh_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double, 
			      unsigned);
extern void mltaGeHLrMh_toMbl(double, blcluster*, mblock<double>**, blcluster*,
                              mblock<double>**, mblock<double>*, double, 
                              unsigned);
////mltaGeHGeHh_toHeH.cpp:
extern void mltaGeHGeHh_toHeH(double, blcluster*, mblock<double>**, blcluster*,
			      mblock<double>**, blcluster*, mblock<double>**,
			      double, unsigned);
////mltaHeHGeH.cpp:
extern void mltaHeHGeH(double, blcluster*, mblock<double>**, blcluster*,
                       mblock<double>**, blcluster*, mblock<double>**, double,
                       unsigned);
////mltaGeHHeH.cpp:
extern void mltaGeHHeH(double, blcluster*, mblock<double>**, blcluster*,
                       mblock<double>**, blcluster*, mblock<double>**,
                       double, unsigned);
////mltaGeHLtH.cpp:
extern void mltaGeHLtH(double d, blcluster* blA, mblock<double>** A, 
		       blcluster* blB, mblock<double>** B, blcluster* blC, 
		       mblock<double>** C, double eps, unsigned rankmax, 
		       contBasis<double>* haar=NULL);
////mltaUtHGeH.cpp:
extern void mltaUtHGeH(double d, blcluster* blA, mblock<double>** A, 
		       blcluster* blB, mblock<double>** B, blcluster* blC, 
		       mblock<double>** C, double eps, unsigned rankmax, 
		       contBasis<double>* haar=NULL);
////mltaUtHhGeH.cpp:
extern void mltaUtHhGeH(double, blcluster*, mblock<double>**, blcluster*,
                       mblock<double>**, blcluster*, mblock<double>**,
                       double, unsigned);
////mltaGeHUtHh.cpp:
extern void mltaGeHUtHh(double, blcluster*, mblock<double>**, blcluster*,
                       mblock<double>**, blcluster*, mblock<double>**,
                       double, unsigned);
////mltaUtHhUtH_toHeH.cpp:
extern void mltaUtHhUtH_toHeH(double, blcluster*, mblock<double>**, 
			      mblock<double>**, double, unsigned);
////mltaUtHUtHh.cpp:
extern void mltaUtHUtHh(double, blcluster*, mblock<double>**, mblock<double>**,
			double, unsigned);
extern void mltaUtHLtH(double, blcluster*, mblock<double>**, mblock<double>**,
		       mblock<double>**, double, unsigned, 
		       contBasis<double>* haar=NULL);

////invertH.cpp:
extern void invertGeH(blcluster*, mblock<double>**, mblock<double>**&,
		      double, unsigned, contBasis<double>* haar=NULL);
extern void invertHeH(blcluster*, mblock<double>**, mblock<double>**&,
		      double, unsigned);
////HLU.cpp:
extern bool HLU(blcluster*, mblock<double>**, mblock<double>**, 
		mblock<double>**, double, unsigned, 
		contBasis<double>* haar=NULL);
extern bool HCholesky(blcluster* const, mblock<double>** const, const double,
                      const unsigned, contBasis<double>* haar=NULL);
extern bool HUhDU(blcluster*, mblock<double>**, int*, double, unsigned);
extern bool genLUprecond(blcluster*, mblock<double>**, double, unsigned,
                         blcluster*&, mblock<double>**&, mblock<double>**&,
                         bool);
extern bool genCholprecond(blcluster*, mblock<double>**, double, unsigned,
                           blcluster*&, mblock<double>**&, bool);
////TU_solve.cpp:
extern void LtHGeM_solve(blcluster*, mblock<double>**, unsigned, double*,
                         unsigned);
extern void LtHhGeM_solve(blcluster*, mblock<double>**, unsigned, double*,
                          unsigned);
extern void LtHGeH_solve(blcluster*, mblock<double>**, blcluster*,
			 mblock<double>**, mblock<double>**, double, unsigned,
			 contBasis<double>* haar=NULL);
extern void UtHGeM_solve(blcluster*, mblock<double>**, unsigned, double*,
                         unsigned);
extern void UtHGeM_solve(blcluster*, mblock<double>**, int*, unsigned,
			 double*, unsigned);
extern void UtHhGeM_solve(blcluster*, mblock<double>**, unsigned, double*,
                          unsigned);
extern void UtHGeH_solve(blcluster*, mblock<double>**, blcluster*,
			 mblock<double>**, double, unsigned);
extern void UtHhGeH_solve(blcluster*, mblock<double>**, blcluster*,
			  mblock<double>**, double, unsigned, 
			  contBasis<double>* haar=NULL);
extern void GeHUtH_solve(blcluster*, mblock<double>**, blcluster*,
			 mblock<double>**, mblock<double>**, double, unsigned,
			 contBasis<double>* haar=NULL);
extern void UtHhDH_solve(mblock<double>**, blcluster*, blcluster*, int*,
			 double, unsigned);
extern void UtHhDdns_solve(blcluster*, mblock<double>**, int*, unsigned,
			   double*, unsigned);
//extern void HLtHh_solve(blcluster*, mblock<double>**, blcluster*,
// mblock<double>** B, mblock<double>**, double, unsigned);

////Hequilib.cpp:
extern void Hcolnrms(unsigned, blcluster*, mblock<double>**, double*);
extern void Hrownrms(unsigned, blcluster*, mblock<double>**, double*);
extern void Hnrms_sym(unsigned, blcluster*, mblock<double>**, double*);
extern void Hequilib_sym(unsigned, blcluster*, mblock<double>**, double*);
extern void Hequilib_cols(unsigned, blcluster*, mblock<double>**, double*);
extern void Hequilib_rows(unsigned, blcluster*, mblock<double>**, double*);
////psoutH.cpp:
extern void psoutputGeH(std::ofstream&, blcluster*, unsigned, mblock<double>**);
extern void psoutputGeH(std::ofstream&, blcluster**, unsigned, unsigned,
			mblock<double>**);
extern void psoutputHeH(std::ofstream&, blcluster**, unsigned, unsigned,
                         mblock<double>**);
extern void psoutputHeH(std::ofstream&, blcluster*, unsigned, mblock<double>**,
			unsigned level=0, cluster* cl=NULL,
			void (*f)(cluster* cl, const unsigned level, 
				  double* X, const double* const basis)=NULL);
extern void psoutputLtH(std::ofstream&, blcluster*, unsigned, mblock<double>**);
extern void psoutputUtH(std::ofstream&, blcluster*, unsigned, mblock<double>**);

////??
extern void expndHSym(blcluster*, mblock<double>**, blcluster*&,
		      mblock<double>**);
extern void convHeH_toHeM(blcluster**, unsigned, mblock<double>**, double*,
			  unsigned);
extern void initrndH(blcluster*, mblock<double>**, mblock<double>**);
extern void mlta_part_vec(double, blcluster**, unsigned, mblock<double>**, 
			  double*, double*);
extern void mlta_partT_vec(double, blcluster**, unsigned, mblock<double>**, 
			   double*, double*);
extern void mlta_part_vec(double, blcluster**, unsigned,
			  mblock<double>**, double*, double*);
extern void mlta_partT_vec(double, blcluster**, unsigned,
			   mblock<double>**, double*, double*);
extern void mlta_sympart_vec(double, blcluster**, unsigned, mblock<double>**,
			     double*, double*);


///////////////////////////////////////////////////////////////////////////
//
// single precision
//

////initH.cpp:
extern void freembls_recursive(blcluster*, mblock<float>** &);
extern void copyH(blcluster*, mblock<float>**, mblock<float>**);
extern void copyH(blcluster*, mblock<double>**, mblock<float>**);
extern void copyH(unsigned, mblock<double>**, mblock<float>**);
extern unsigned Hmax_rank(blcluster*, mblock<float>**, char type='H');
extern unsigned long sizeH(blcluster*, mblock<float>**, char type='H');
extern unsigned long sizeH(unsigned, mblock<float>**);
extern void initGeH_0(blcluster*, mblock<float>** &);
extern void initGeH_0_withoutAlloc(blcluster*, mblock<float>** &);
extern void initGeH_Id(blcluster*, mblock<float>** &);
extern void initHeH_0(blcluster*, mblock<float>** &);
extern void initHeH_0_withoutAlloc(blcluster*, mblock<float>** &);
extern void initHeH_Id(blcluster*, mblock<float>** &);
extern void initLtH_0(blcluster*, mblock<float>** &);
extern void initLtH_0_withoutAlloc(blcluster*, mblock<float>** &);
extern void initUtH_0(blcluster*, mblock<float>** &);
extern void initUtH_0_withoutAlloc(blcluster*, mblock<float>** &);
extern void setGeHzero(blcluster*, mblock<float>**);
extern void setHeHzero(blcluster*, mblock<float>**);
////convCS_toH.cpp:
extern void convCCS_toGeH(double*, unsigned*, unsigned*, double,
			  blcluster*, mblock<float>**&);
extern void convCCS_toGeH(double*, unsigned*, unsigned*, unsigned*,
			  unsigned*, double, blcluster*, mblock<float>**&);
extern void convCRS_toGeH(double*, unsigned*, unsigned*,
			  double, blcluster*, mblock<float>**&);
extern void convCRS_toGeH(double*, unsigned*, unsigned*, unsigned*,
			  unsigned*, blcluster*, mblock<float>**&);
extern void convCRS_toGeH(double*, unsigned*, unsigned*, unsigned*,
			  unsigned*, double, blcluster*, mblock<float>**&);
extern void convCCS_toHeH(double*, unsigned*, unsigned*, unsigned*,
                          unsigned*, double, blcluster*, mblock<float>**&);
extern void convCCS_toHeH(double*, unsigned*, unsigned*,
                          double, blcluster*, mblock<float>**&);
extern void convCRS_toHeH(double*, unsigned*, unsigned*,
                          double, blcluster*, mblock<float>**&);
extern void convCRS_toHeH(double*, unsigned*, unsigned*, unsigned*,
                          unsigned*, double, blcluster*, mblock<float>**&);
extern void convCRS_toGeH_withoutAlloc(double*, unsigned*, unsigned*,
				       unsigned*, unsigned*, double, blcluster*,
				       mblock<float>**&);
extern void convCRS_toHeH_withoutAlloc(double*, unsigned*, unsigned*,
                                       unsigned*, unsigned*, double,
                                       blcluster*, mblock<float>**&);
extern void convCRS_toGeH_withoutAlloc(double*, unsigned*, unsigned*,
				       double, blcluster*,  mblock<float>**&);
extern void convCRS_toHeH_withoutAlloc(double*, unsigned*, unsigned*,  double,
                                       blcluster*, mblock<float>**&);
////agglH.cpp:
extern void unify_sons(blcluster*, mblock<float>**, double, unsigned);
extern void fill_gaps(blcluster*, mblock<float>**, unsigned, unsigned&);
extern void agglH(blcluster*, mblock<float>**, double, unsigned, 
                  contBasis<float>* haar=NULL);
extern void aggl_light_H(blcluster*, mblock<float>**, double, unsigned);
// extern void agglHSym_stab(blcluster*, mblock<float>**, double, unsigned);
////nrmH.cpp:
extern double nrmF2GeH(blcluster*, mblock<float>**);
extern double nrmF2HeH(blcluster*, mblock<float>**);
extern double spctrlnrm_sym(blcluster*, mblock<float>**);
extern double spctrlnrm_sym(unsigned, blcluster*, mblock<float>**, blcluster*,
                            mblock<float>**);
extern double spctrlnrm(unsigned, blcluster*, mblock<float>**, blcluster*, 
			mblock<float>**);

////addHH.cpp:
extern void addGeHId(blcluster*, mblock<float>**);
extern void addGeHLrM(blcluster*, mblock<float>**, double, unsigned, unsigned,
		      float*, unsigned, float*, unsigned, 
		      contBasis<float>* haar=NULL);
extern void addGeHGeM(blcluster*, mblock<float>**, float*, unsigned,
		      double, unsigned, contBasis<float>* haar=NULL);
extern void addGeHGeH(blcluster* bl, mblock<float>** A, mblock<float>** B,
		      float eps, unsigned rankmax, contBasis<float>* haar=NULL);
extern void addHeHLrM(blcluster*, mblock<float>**, double, unsigned, unsigned,
		      float*, unsigned, float*, unsigned, 
		      contBasis<float>* haar=NULL);
extern void addHeHLrMExact(blcluster* bl, mblock<float>** A, unsigned k, 
			   float* U, unsigned ldU, float* V, unsigned ldV);
extern void addHeHLrM_stab(blcluster*, mblock<float>**, double, unsigned,
			   unsigned, float*, unsigned, float*, unsigned);
extern void addHeHGeM(blcluster*, mblock<float>**, float*, unsigned,
		      double, unsigned, contBasis<float>* haar=NULL);
////mltaGeHVec.cpp:
extern bool mltaGeHVec(float, blcluster*, mblock<float>**, float*, float*);
extern bool mltaGeHGeM(float, blcluster*, mblock<float>**, unsigned,
		       float*, unsigned, float*, unsigned);
extern bool mltaGeHhVec(float, blcluster*, mblock<float>**, float*, float*);
extern bool mltaGeHhGeM(float, blcluster*, mblock<float>**, unsigned,
			float*, unsigned, float*, unsigned);
extern bool mltaGeHhDiHVec(float, mblock<float>**, blcluster*, blcluster*, int*,
			   float*, float*);
extern bool mltaGeHhDiHGeM(float, mblock<float>**, blcluster*, 
			   blcluster*, int*, unsigned, float*, unsigned, 
			   float*, unsigned);
extern void mltaLtHVec(float, blcluster*, mblock<float>**, float*, float*);
extern void mltaLtHGeM(float, blcluster*, mblock<float>**, unsigned, float*,
		       unsigned, float*, unsigned);
extern void mltaLtHhVec(float, blcluster*, mblock<float>**, float*, float*);
extern void mltaLtHhGeM(float, blcluster*, mblock<float>**, unsigned, float*,
			unsigned, float*, unsigned);
extern void mltaUtHVec(float, blcluster*, mblock<float>**, float*, float*);
extern void mltaUtHGeM(float, blcluster*, mblock<float>**, unsigned, float*,
		       unsigned, float*, unsigned);
extern void mltaUtHhVec(float, blcluster*, mblock<float>**, float*, float*);
extern void mltaUtHhGeM(float, blcluster*, mblock<float>**, unsigned, float*,
			unsigned, float*, unsigned);
extern void mltaHeHVec(float, blcluster*, mblock<float>**, float*, float*);
extern void mltaHeHGeM(float, blcluster*, mblock<float>**, unsigned,
		       float*, unsigned, float*, unsigned);
////mltaGehGeH.cpp:
extern void mltaGeHGeH(float, blcluster*, mblock<float>**, blcluster*, 
		       mblock<float>**, blcluster*, mblock<float>**, double, unsigned,
		       contBasis<float>* haar=NULL);
extern void mltaGeHGeH_toMbl(float, blcluster*, mblock<float>**, blcluster*,
			     mblock<float>**, mblock<float>*, double, unsigned,
			     contBasis<float>* haar=NULL);
////mltaGeHhGeH.cpp:
extern void mltaGeHhGeH(float, blcluster*, mblock<float>**, blcluster*,
			mblock<float>**, blcluster*, mblock<float>**, double,
			unsigned, contBasis<float>* haar=NULL);
extern void mltaGeHhGeH_toMbl(float, blcluster*, mblock<float>**, blcluster*,
			      mblock<float>**, mblock<float>*, double, unsigned);
extern void mltaGeHhLrM_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned, 
                              contBasis<float>* haar=NULL);
extern void mltaGeHhGeM_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned, 
                              contBasis<float>* haar=NULL);
extern void mltaGeMhGeH_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned,
                              contBasis<float>* haar=NULL);
extern void mltaLrMhGeH_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned, 
                              contBasis<float>* haar=NULL);
////mltaGeHhGeH_toHeH.cpp:
extern void mltaGeHhGeH_toHeH(float, blcluster*, mblock<float>**, blcluster*,
			      mblock<float>**, blcluster*, mblock<float>**, 
			      double, unsigned, contBasis<float>* haar=NULL);
////mltaGeHhDiHGeH_toHeH.cpp:
extern void mltaGeHhDiHGeH_toHeH(float, mblock<float>**, blcluster*, 
				 blcluster*, int*, blcluster*,
				 blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH_toHeH(float, blcluster*, mblock<float>**, blcluster*,
				 mblock<float>**, int*, blcluster*,
				 mblock<float>**, blcluster*, mblock<float>**,
				 double, unsigned);
////mltaGeHhDiHGeH.cpp:
extern void mltaGeHhDiHGeH(float, mblock<float>**, blcluster*, 
			   blcluster*, int*, blcluster*,
			   blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH(float, blcluster*, mblock<float>**, blcluster*,
			   mblock<float>**, int*, blcluster*, mblock<float>**,
			   blcluster*, mblock<float>**, double, unsigned);
extern void mltaGeHhDiHLrM_toMbl(float, blcluster*, mblock<float>**, blcluster*,
				 mblock<float>**, int*, blcluster*,
				 mblock<float>**, mblock<float>*, double,
				 unsigned);
extern void mltaGeHhDiHGeM_toMbl(float, blcluster*, mblock<float>**, blcluster*,
				 mblock<float>**, int*, blcluster*,
				 mblock<float>**, mblock<float>*, double,
				 unsigned);
extern void mltaGeMhDiHGeH_toMbl(float, blcluster*, mblock<float>**, blcluster*,
				 mblock<float>**, int*, blcluster*,
				 mblock<float>**, mblock<float>*, double,
				 unsigned);
extern void mltaLrMhDiHGeH_toMbl(float, blcluster*, mblock<float>**, blcluster*,
				 mblock<float>**, int*, blcluster*,
				 mblock<float>**, mblock<float>*, double,
				 unsigned);
extern void mltDiHGeM(unsigned, mblock<float>**, blcluster*, int* const,
		      float*, float*, unsigned);
////mltaGeHGeHh.cpp:
extern void mltaGeHGeHh(float, blcluster*, mblock<float>**, blcluster*,
			mblock<float>**, blcluster*, mblock<float>**,
			double, unsigned);
extern void mltaGeHGeHh_toMbl(float, blcluster*, mblock<float>**,
			      blcluster*, mblock<float>**, mblock<float>*, 
			      double, unsigned);
extern void mltaGeHGeMh_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned);
extern void mltaGeMGeHh_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned);
extern void mltaGeHLrMh_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, 
			      unsigned);
extern void mltaLrMGeHh_toMbl(float, blcluster*, mblock<float>**, blcluster*,
                              mblock<float>**, mblock<float>*, double, unsigned);
////mltaGeHGeHh_toHeH.cpp:
extern void mltaGeHGeHh_toHeH(float, blcluster*, mblock<float>**, blcluster*,
			      mblock<float>**, blcluster*, mblock<float>**,
			      double, unsigned);
////mltaHeHGeH.cpp:
extern void mltaHeHGeH(float, blcluster*, mblock<float>**, blcluster*,
                       mblock<float>**, blcluster*, mblock<float>**, 
		       double, unsigned);
////mltaGeHHeH.cpp:
extern void mltaGeHHeH(float, blcluster*, mblock<float>**, blcluster*,
                       mblock<float>**, blcluster*, mblock<float>**,
                       double, unsigned);
////mltaGeHLtH.cpp:
extern void mltaGeHLtH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
		       mblock<float>** B, blcluster* blC, mblock<float>** C,
		       double eps, unsigned rankmax, contBasis<float>* haar=NULL);
////mltaUtHGeH.cpp:
extern void mltaUtHGeH(float d, blcluster* blA, mblock<float>** A, blcluster* blB,
		       mblock<float>** B, blcluster* blC, mblock<float>** C,
		       double eps, unsigned rankmax, contBasis<float>* haar=NULL);
////mltaUtHhGeH.cpp:
extern void mltaUtHhGeH(float, blcluster*, mblock<float>**, blcluster*,
                       mblock<float>**, blcluster*, mblock<float>**,
                       double, unsigned);
////mltaGeHUtHh.cpp:
extern void mltaGeHUtHh(float, blcluster*, mblock<float>**, blcluster*,
                       mblock<float>**, blcluster*, mblock<float>**,
                       double, unsigned);
////mltaUtHhUtH_toHeH.cpp:
extern void mltaUtHhUtH_toHeH(float, blcluster*, mblock<float>**, mblock<float>**,
                         double, unsigned);
////mltaUtHUtHh.cpp:
extern void mltaUtHUtHh(float, blcluster*, mblock<float>**, mblock<float>**,
                         double, unsigned);
extern void mltaUtHLtH(float, blcluster*, mblock<float>**, mblock<float>**,
			mblock<float>**, double, unsigned, 
			contBasis<float>* haar=NULL);

////invertH.cpp:
extern void invertGeH(blcluster*, mblock<float>**, mblock<float>**&,
		      double, unsigned, contBasis<float>* haar=NULL);
extern void invertHeH(blcluster*, mblock<float>**, mblock<float>**&,
		      double, unsigned);
////HLU.cpp:
extern bool HLU(blcluster*, mblock<float>**, mblock<float>**, mblock<float>**,
                double, unsigned, contBasis<float>* haar=NULL);
extern bool HCholesky(blcluster* const, mblock<float>**, double, unsigned, 
                      contBasis<float>* haar=NULL);
extern bool HUhDU(blcluster*, mblock<float>**, int*, double, unsigned);
extern bool genLUprecond(blcluster*, mblock<double>**, double, unsigned,
                         blcluster*&, mblock<float>**&, mblock<float>**&, bool);
extern bool genLUprecond(blcluster*, mblock<float>**, double, unsigned,
                         blcluster*&, mblock<float>**&, mblock<float>**&, bool);
extern bool genCholprecond(blcluster*, mblock<double>**, double, unsigned,
                           blcluster*&, mblock<float>**&, bool);
extern bool genCholprecond(blcluster*, mblock<float>**, double, unsigned,
                           blcluster*&, mblock<float>**&, bool);
////TU_solve.cpp:
extern void LtHGeM_solve(blcluster*, mblock<float>**, unsigned, float*,
			 unsigned);
extern void LtHhGeM_solve(blcluster*, mblock<float>**, unsigned, float*,
                          unsigned);
extern void LtHGeH_solve(blcluster*, mblock<float>**, blcluster*,
			 mblock<float>**, mblock<float>**, double, unsigned,
			 contBasis<float>* haar=NULL);
extern void UtHGeM_solve(blcluster*, mblock<float>**, unsigned, float*, unsigned);
extern void UtHGeM_solve(blcluster*, mblock<float>**, int*, unsigned,
			 float*, unsigned);
extern void UtHhGeM_solve(blcluster*, mblock<float>**, unsigned, float*,
                          unsigned);
extern void UtHGeH_solve(blcluster*, mblock<float>**, blcluster*,
			 mblock<float>**, double, unsigned);
extern void UtHhGeH_solve(blcluster*, mblock<float>**, blcluster*,
			  mblock<float>**, double, unsigned, 
			  contBasis<float>* haar=NULL);
extern void GeHUtH_solve(blcluster*, mblock<float>**, blcluster*,
			 mblock<float>**, mblock<float>**, double, unsigned, 
			 contBasis<float>* haar=NULL);
extern void UtHhDH_solve(mblock<float>**, blcluster*, blcluster*, int*,
			 double, unsigned);
extern void UtHhDdns_solve(blcluster*, mblock<float>**, int*, unsigned,
			   float*, unsigned);
/*
  extern void HLtHh_solve(blcluster*, mblock<float>**, blcluster*,
  mblock<float>** B, mblock<float>**, double, unsigned);
*/

////psoutH.cpp:
extern void psoutputGeH(std::ofstream&, blcluster*, unsigned, mblock<float>**);
extern void psoutputGeH(std::ofstream&, blcluster**, unsigned, unsigned,
                      mblock<float>**);
extern void psoutputHeH(std::ofstream&, blcluster**, unsigned, unsigned,
                         mblock<float>**);
extern void psoutputHeH(std::ofstream&, blcluster*, unsigned, mblock<float>**,
                         unsigned level=0, cluster* cl=NULL,
                         void (*f)(cluster* cl, const unsigned level, 
                                   double* X, const double* const basis)=NULL);
extern void psoutputLtH(std::ofstream&, blcluster*, unsigned, mblock<float>**);
extern void psoutputUtH(std::ofstream&, blcluster*, unsigned, mblock<float>**);

////??
extern void expndHSym(blcluster*, mblock<float>**, blcluster*&,
		      mblock<float>**);
extern void convHeH_toHeM(blcluster**, unsigned, mblock<float>**, float*, unsigned);
extern void initrndH(blcluster*, mblock<float>**, mblock<float>**);
extern void mlta_part_vec(float, blcluster**, unsigned, mblock<float>**, float*,
                           float*);
extern void mlta_partT_vec(float, blcluster**, unsigned, mblock<float>**, float*,
                            float*);

extern void mlta_part_vec(float, blcluster**, unsigned,
                           mblock<float>**, float*, float*);
extern void mlta_partT_vec(float, blcluster**, unsigned,
                            mblock<float>**, float*, float*);
extern void mlta_sympart_vec(float, blcluster**, unsigned, mblock<float>**,
                              float*, float*);


///////////////////////////////////////////////////////////////////////////
//
// single precision complex
//

////initH.cpp:
extern void copyH(blcluster*, mblock<scomp>**, mblock<scomp>**);
extern void copyH(blcluster*, mblock<dcomp>**, mblock<scomp>**);
extern void freembls_recursive(blcluster*, mblock<scomp>** &);
extern unsigned Hmax_rank(blcluster*, mblock<scomp>**, char type='H');
extern unsigned long sizeH(blcluster*, mblock<scomp>**, char type='H');
extern unsigned long sizeH(unsigned, mblock<scomp>**);
extern void initGeH_0(blcluster*, mblock<scomp>** &);
extern void initGeH_0_withoutAlloc(blcluster*, mblock<scomp>** &);
extern void initGeH_Id(blcluster*, mblock<scomp>** &);
extern void initHeH_0_withoutAlloc(blcluster*, mblock<scomp>** &);
extern void initHeH_0(blcluster*, mblock<scomp>** &);
extern void initHeH_Id(blcluster*, mblock<scomp>** &);
extern void initLtH_0(blcluster*, mblock<scomp>** &);
extern void initLtH_0_withoutAlloc(blcluster*, mblock<scomp>** &);
extern void initUtH_0(blcluster*, mblock<scomp>** &);
extern void initUtH_0_withoutAlloc(blcluster*, mblock<scomp>** &);
extern void setGeHzero(blcluster*, mblock<scomp>**);
extern void setHeHzero(blcluster*, mblock<scomp>**);
////agglH.cpp:
extern void unify_sons(blcluster*, mblock<scomp>**, double, unsigned);
extern void fill_gaps(blcluster*, mblock<scomp>**, unsigned, unsigned&);
extern void agglH(blcluster*, mblock<scomp>**, double, unsigned, 
                  contBasis<scomp>* haar=NULL);
extern void aggl_light_H(blcluster*, mblock<scomp>**, double, unsigned);
//extern void agglHSym_stab(blcluster*, mblock<scomp>**, double, unsigned);
////nrmH.cpp:
extern double nrmF2GeH(blcluster*, mblock<scomp>**);
extern double nrmF2HeH(blcluster*, mblock<scomp>**);
extern double spctrlnrm_sym(blcluster*, mblock<scomp>**);
extern double spctrlnrm_sym(unsigned, blcluster*, mblock<scomp>**, blcluster*,
                            mblock<scomp>**);
extern double spctrlnrm(unsigned, blcluster*, mblock<scomp>**, blcluster*, 
			mblock<scomp>**);

////addHH.cpp:
extern void addGeHId(blcluster*, mblock<scomp>**);
extern void addGeHLrM(blcluster*, mblock<scomp>**, double, unsigned, unsigned,
		      scomp*, unsigned, scomp*, unsigned, 
		      contBasis<scomp>* haar=NULL);
extern void addGeHGeM(blcluster*, mblock<scomp>**, scomp*, unsigned,
		      double, unsigned, contBasis<scomp>* haar=NULL);
extern void addGeHGeH(blcluster* bl, mblock<scomp>** A, mblock<scomp>** B,
		      scomp eps, unsigned rankmax, contBasis<scomp>* haar=NULL);
extern void addHeHLrM(blcluster*, mblock<scomp>**, double, unsigned, unsigned,
		      scomp*, unsigned, scomp*, unsigned, 
		      contBasis<scomp>* haar=NULL);
extern void addHeHLrMExact(blcluster* bl, mblock<scomp>** A, unsigned k, 
			   scomp* U, unsigned ldU, scomp* V, unsigned ldV);
extern void addHeHLrM_stab(blcluster*, mblock<scomp>**, double, unsigned,
			   unsigned, scomp*, unsigned, scomp*, unsigned);
extern void addHeHGeM(blcluster*, mblock<scomp>**, scomp*, unsigned, 
		      double, unsigned, contBasis<scomp>* haar=NULL);
////mltaGeHVec.cpp
extern bool mltaGeHVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern bool mltaGeHGeM(scomp, blcluster*, mblock<scomp>**, unsigned,
		       scomp*, unsigned, scomp*, unsigned);
extern bool mltaGeHhVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern bool mltaGeHhGeM(scomp, blcluster*, mblock<scomp>**, unsigned,
			scomp*, unsigned, scomp*, unsigned);
extern bool mltaGeHtVec(scomp, blcluster*, mblock<scomp>**, scomp*,scomp*);
extern bool mltaGeHtGeM(scomp, blcluster*, mblock<scomp>**, unsigned,
			scomp*, unsigned, scomp*, unsigned);
extern bool mltaGeHhDiHVec(scomp, mblock<scomp>**, blcluster*, blcluster*, int*,
			   scomp*, scomp*);
extern bool mltaGeHhDiHGeM(scomp, mblock<scomp>**, blcluster*, 
			   blcluster*, int*, unsigned, scomp*, unsigned, 
			   scomp*, unsigned);
extern void mltaLtHVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern void mltaLtHGeM(scomp, blcluster*, mblock<scomp>**, unsigned, scomp*,
		       unsigned, scomp*, unsigned);
extern void mltaLtHhVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern void mltaLtHhGeM(scomp, blcluster*, mblock<scomp>**, unsigned, scomp*,
			unsigned, scomp*, unsigned);
extern void mltaUtHVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern void mltaUtHGeM(scomp, blcluster*, mblock<scomp>**, unsigned, scomp*,
		       unsigned, scomp*, unsigned);
extern void mltaUtHhVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern void mltaUtHhGeM(scomp, blcluster*, mblock<scomp>**, unsigned, scomp*,
			unsigned, scomp*, unsigned);
extern void mltaHeHVec(scomp, blcluster*, mblock<scomp>**, scomp*, scomp*);
extern void mltaHeHGeM(scomp, blcluster*, mblock<scomp>**, unsigned,
		       scomp*, unsigned, scomp*, unsigned);
extern void mltaHeHtVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaSyHVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaSyHGeM(scomp, blcluster*, mblock<scomp>**, unsigned,
		       scomp*, unsigned, scomp*, unsigned);
extern void mltaSyHhVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
////mltaGeHGeH.cpp:
extern void mltaGeHGeH(scomp, blcluster*, mblock<scomp>**, blcluster*, mblock<scomp>**,
		       blcluster*, mblock<scomp>**, double, unsigned, 
		       contBasis<scomp>* haar=NULL);
extern void mltaGeHGeH_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
			     mblock<scomp>**, mblock<scomp>*, double, unsigned,
			     contBasis<scomp>* haar=NULL);
////mltaGeHhGeH.cpp:
extern void mltaGeHhGeH(scomp, blcluster*, mblock<scomp>**, blcluster*, mblock<scomp>**,
			blcluster*, mblock<scomp>**, double, unsigned, 
			contBasis<scomp>* haar=NULL);
extern void mltaGeHhGeH_toMbl(scomp, blcluster*, mblock<scomp>**,
			      blcluster*, mblock<scomp>**, mblock<scomp>*, double, 
			      unsigned);
extern void mltaGeHhLrM_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned, 
                              contBasis<scomp>* haar=NULL);
extern void mltaGeHhGeM_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned, 
                              contBasis<scomp>* haar=NULL);
extern void mltaGeMhGeH_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned, 
                              contBasis<scomp>* haar=NULL);
extern void mltaLrMhGeH_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned, 
                              contBasis<scomp>* haar=NULL);
////mltaGeHhGeH_toHeH.cpp:
extern void mltaGeHhGeH_toHeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
			      mblock<scomp>**, blcluster*, mblock<scomp>**, 
			      double, unsigned, contBasis<scomp>* haar=NULL);
////mltaGeHhDiHGeH_toHeH.cpp:
extern void mltaGeHhDiHGeH_toHeH(scomp, mblock<scomp>**, blcluster*, 
				 blcluster*, int*, blcluster*,
				 blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH_toHeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
				 mblock<scomp>**, int*, blcluster*,
				 mblock<scomp>**, blcluster*, mblock<scomp>**,
				 double, unsigned);
////mltaGeHhDiHGeH.cpp:
extern void mltaGeHhDiHGeH(scomp, mblock<scomp>**, blcluster*, 
			   blcluster*, int*, blcluster*,
			   blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
			   mblock<scomp>**, int*, blcluster*, mblock<scomp>**,
			   blcluster*, mblock<scomp>**, double, unsigned);
extern void mltaGeHhDiHLrM_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
				 mblock<scomp>**, int*, blcluster*,
				 mblock<scomp>**, mblock<scomp>*, double,
				 unsigned);
extern void mltaGeHhDiHGeM_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
				 mblock<scomp>**, int*, blcluster*,
				 mblock<scomp>**, mblock<scomp>*, double,
				 unsigned);
extern void mltaGeMhDiHGeH_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
				 mblock<scomp>**, int*, blcluster*,
				 mblock<scomp>**, mblock<scomp>*, double,
				 unsigned);
extern void mltaLrMhDiHGeH_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
				 mblock<scomp>**, int*, blcluster*,
				 mblock<scomp>**, mblock<scomp>*, double,
				 unsigned);
extern void mltDiHGeM(unsigned, mblock<scomp>**, blcluster*, int* const ,
		      scomp*, scomp*, unsigned);
////mltaGeHGeHh.cpp:
extern void mltaGeHGeHh(scomp, blcluster*, mblock<scomp>**, blcluster*,
			mblock<scomp>**, blcluster*, mblock<scomp>**,
			double, unsigned);
extern void mltaGeHGeHh_toMbl(scomp, blcluster*, mblock<scomp>**,
			      blcluster*, mblock<scomp>**, mblock<scomp>*, 
			      double, unsigned);
extern void mltaGeHGeMh_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned);
extern void mltaGeMGeHh_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned);
extern void mltaLrMGeHh_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned);
extern void mltaGeHLrMh_toMbl(scomp, blcluster*, mblock<scomp>**, blcluster*,
                              mblock<scomp>**, mblock<scomp>*, double, unsigned);
////mltaGeHGeHh_toHeH.cpp:
extern void mltaGeHGeHh_toHeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
			      mblock<scomp>**, blcluster*, mblock<scomp>**,
			      double, unsigned);
////mltaHeHGeH.cpp:
extern void mltaHeHGeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
                       mblock<scomp>**, blcluster*, mblock<scomp>**, 
		       double, unsigned);
////mltaGeHHeH.cpp:
extern void mltaGeHHeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
                       mblock<scomp>**, blcluster*, mblock<scomp>**,
                       double, unsigned);
////mltaGeHLtH.cpp:
extern void mltaGeHLtH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
		       mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
		       double eps, unsigned rankmax, contBasis<scomp>* haar=NULL);
////mltaUtHGeH.cpp:
extern void mltaUtHGeH(scomp d, blcluster* blA, mblock<scomp>** A, blcluster* blB,
		       mblock<scomp>** B, blcluster* blC, mblock<scomp>** C,
		       double eps, unsigned rankmax, contBasis<scomp>* haar=NULL);
////mltaUtHhGeH.cpp:
extern void mltaUtHhGeH(scomp, blcluster*, mblock<scomp>**, blcluster*,
                       mblock<scomp>**, blcluster*, mblock<scomp>**,
                       double, unsigned);
////mltaGeHUtHh.cpp:
extern void mltaGeHUtHh(scomp, blcluster*, mblock<scomp>**, blcluster*,
                       mblock<scomp>**, blcluster*, mblock<scomp>**,
                       double, unsigned);
////mltaUtHhUtH_toHeH.cpp:
extern void mltaUtHhUtH_toHeH(scomp, blcluster*, mblock<scomp>**, mblock<scomp>**,
                         double, unsigned);
////mltaUtHUtHh.cpp:
extern void mltaUtHUtHh(scomp, blcluster*, mblock<scomp>**, mblock<scomp>**,
                         double, unsigned);
extern void mltaUtHLtH(scomp, blcluster*, mblock<scomp>**, mblock<scomp>**,
			mblock<scomp>**, double, unsigned, 
			contBasis<scomp>* haar=NULL);

////invertH.cpp:
extern void invertGeH(blcluster*, mblock<scomp>**, mblock<scomp>**&,
		      scomp, unsigned, contBasis<scomp>* haar=NULL);
extern void invertHeH(blcluster*, mblock<scomp>**, mblock<scomp>**&,
		      scomp, unsigned);
////HLU.cpp:
extern bool HLU(blcluster*, mblock<scomp>**, mblock<scomp>**, mblock<scomp>**, double,
		unsigned, contBasis<scomp>* haar=NULL);
extern bool HCholesky(blcluster* const, mblock<scomp>** const, double, unsigned, 
                      contBasis<scomp>* haar=NULL);
extern bool HUhDU(blcluster*, mblock<scomp>**, int*, double, unsigned);
extern bool genLUprecond(blcluster*, mblock<dcomp>**, double, unsigned,
                         blcluster*&, mblock<scomp>**&, mblock<scomp>**&, bool);
extern bool genCholprecond(blcluster*, mblock<scomp>**, scomp, unsigned,
                           blcluster*&, mblock<scomp>**&, bool);
////TU_solve.cpp:
extern void LtHGeM_solve(blcluster*, mblock<scomp>**, unsigned, scomp*,
                         unsigned);
extern void LtHhGeM_solve(blcluster*, mblock<scomp>**, unsigned, scomp*,
                          unsigned);
extern void LtHGeH_solve(blcluster*, mblock<scomp>**, blcluster*,
			 mblock<scomp>**, mblock<scomp>**, double, unsigned,
			 contBasis<scomp>* haar=NULL);
extern void UtHGeM_solve(blcluster*, mblock<scomp>**, unsigned, scomp*,
                         unsigned);
extern void UtHGeM_solve(blcluster*, mblock<scomp>**, int*, unsigned,
			 scomp*, unsigned);
extern void UtHhGeM_solve(blcluster*, mblock<scomp>**, unsigned, scomp*,
                          unsigned);
extern void UtHGeH_solve(blcluster*, mblock<scomp>**, blcluster*,
			 mblock<scomp>**, double, unsigned);
extern void UtHhGeH_solve(blcluster*, mblock<scomp>**, blcluster*,
			  mblock<scomp>**, double, unsigned, 
			  contBasis<scomp>* haar=NULL);
extern void GeHUtH_solve(blcluster*, mblock<scomp>**, blcluster*,
			 mblock<scomp>**, mblock<scomp>**, double, unsigned, 
			 contBasis<scomp>* haar=NULL);
extern void UtHhDH_solve(mblock<scomp>**, blcluster*, blcluster*, int*,
			 double, unsigned);
extern void UtHhDdns_solve(blcluster*, mblock<scomp>**, int*, unsigned,
			   scomp*, unsigned);
/*
  extern void HLtHh_solve(blcluster*, mblock<scomp>**, blcluster*,
  mblock<scomp>** B, mblock<scomp>**, double, unsigned);
*/

////psoutH.cpp:
extern void psoutputGeH(std::ofstream&, blcluster*, unsigned, mblock<scomp>**);
extern void psoutputGeH(std::ofstream&, blcluster**, unsigned, unsigned,
                      mblock<scomp>**);
extern void psoutputHeH(std::ofstream&, blcluster**, unsigned, unsigned,
                         mblock<scomp>**);
extern void psoutputHeH(std::ofstream&, blcluster*, unsigned, mblock<scomp>**,
			unsigned level=0, cluster* cl=NULL,
			void (*f)(cluster* cl, const unsigned level, 
				  double* X, const double* const basis)=NULL);
extern void psoutputLtH(std::ofstream&, blcluster*, unsigned, mblock<scomp>**);
extern void psoutputUtH(std::ofstream&, blcluster*, unsigned, mblock<scomp>**);

////??
extern void expndHSym(blcluster*, mblock<scomp>**, blcluster*&,
		      mblock<scomp>**);
extern void convHeH_toHeM(blcluster**, unsigned, mblock<scomp>**, scomp*, unsigned);
extern void mlta_part_vec(scomp, blcluster**, unsigned, mblock<scomp>**, scomp*,
                           scomp*);
extern void mlta_partT_vec(scomp, blcluster**, unsigned, mblock<scomp>**, scomp*,
                            scomp*);
extern void mlta_part_vec(scomp, blcluster**, unsigned,
                           mblock<scomp>**, scomp*, scomp*);
extern void mlta_partT_vec(scomp, blcluster**, unsigned,
                            mblock<scomp>**, scomp*, scomp*);
extern void mlta_sympart_vec(scomp, blcluster**, unsigned, mblock<scomp>**,
                              scomp*, scomp*);


///////////////////////////////////////////////////////////////////////////
//
// double precision complex
//

////initH.cpp:
extern void copyH(blcluster*, mblock<dcomp>**, mblock<dcomp>**);
extern void freembls_recursive(blcluster*, mblock<dcomp>** &);
extern unsigned Hmax_rank(blcluster*, mblock<dcomp>**, char type='H');
extern unsigned long sizeH(blcluster*, mblock<dcomp>**, char type='H');
extern unsigned long sizeH(unsigned, mblock<dcomp>**);
extern void initGeH_0(blcluster*, mblock<dcomp>** &);
extern void initGeH_0_withoutAlloc(blcluster*, mblock<dcomp>** &);
extern void initGeH_Id(blcluster*, mblock<dcomp>** &);
extern void initHeH_0(blcluster*, mblock<dcomp>** &);
extern void initHeH_0_withoutAlloc(blcluster*, mblock<dcomp>** &);
extern void initHeH_Id(blcluster*, mblock<dcomp>** &);
extern void initLtH_0(blcluster*, mblock<dcomp>** &);
extern void initLtH_0_withoutAlloc(blcluster*, mblock<dcomp>** &);
extern void initUtH_0(blcluster*, mblock<dcomp>** &);
extern void initUtH_0_withoutAlloc(blcluster*, mblock<dcomp>** &);
extern void setGeHzero(blcluster*, mblock<dcomp>**);
extern void setHeHzero(blcluster*, mblock<dcomp>**);
////agglH.cpp:
extern void unify_sons(blcluster*, mblock<dcomp>**, double, unsigned);
extern void fill_gaps(blcluster*, mblock<dcomp>**, unsigned, unsigned&);
extern void agglH(blcluster*, mblock<dcomp>**, double, unsigned, 
                  contBasis<dcomp>* haar=NULL);
extern void aggl_light_H(blcluster*, mblock<dcomp>**, double, unsigned);
//extern void agglHSym_stab(blcluster*, mblock<dcomp>**, double, unsigned);
////nrmH.cpp:
extern double nrmF2GeH(blcluster*, mblock<dcomp>**);
extern double nrmF2HeH(blcluster*, mblock<dcomp>**);
extern double spctrlnrm_sym(blcluster*, mblock<dcomp>**);
extern double spctrlnrm_sym(unsigned, blcluster*, mblock<dcomp>**, blcluster*,
                            mblock<dcomp>**);
extern double spctrlnrm(unsigned, blcluster*, mblock<dcomp>**, blcluster*,
			mblock<dcomp>**);

////addHH.cpp:
extern void addGeHId(blcluster*, mblock<dcomp>**);
extern void addGeHLrM(blcluster*, mblock<dcomp>**, double, unsigned, unsigned,
		      dcomp*, unsigned, dcomp*, unsigned, 
		      contBasis<dcomp>* haar=NULL);
extern void addGeHGeM(blcluster*, mblock<dcomp>**, dcomp*, unsigned,
		      double, unsigned, contBasis<dcomp>* haar=NULL);
extern void addGeHGeH(blcluster* bl, mblock<dcomp>** A, mblock<dcomp>** B,
		      double eps, unsigned rankmax, contBasis<dcomp>* haar=NULL);
extern void addHeHLrM(blcluster*, mblock<dcomp>**, double, unsigned, unsigned,
		      dcomp*, unsigned, dcomp*, unsigned,
		      contBasis<dcomp>* haar=NULL);
extern void addHeHLrMExact(blcluster* bl, mblock<dcomp>** A, unsigned k, 
			   dcomp* U, unsigned ldU, dcomp* V, unsigned ldV);
extern void addHeHLrM_stab(blcluster*, mblock<dcomp>**, double, unsigned,
			   unsigned, dcomp*, unsigned, dcomp*, unsigned);
extern void addHeHGeM(blcluster*, mblock<dcomp>**, dcomp*, unsigned,
		      double, unsigned, contBasis<dcomp>* haar=NULL);
////mltaGeHVec.cpp
extern bool mltaGeHVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern bool mltaGeHGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned,
		       dcomp*, unsigned, dcomp*, unsigned);
extern bool mltaGeHhVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern bool mltaGeHhGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned,
			dcomp*, unsigned, dcomp*, unsigned);
extern bool mltaGeHtVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*,dcomp*);
extern bool mltaGeHtGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned,
			dcomp*, unsigned, dcomp*, unsigned);
extern bool mltaGeHhDiHVec(dcomp, mblock<dcomp>**, blcluster*, blcluster*, int*,
			   dcomp*, dcomp*);
extern bool mltaGeHhDiHVec(dcomp, mblock<dcomp>**, blcluster*, blcluster*, int*,
			   dcomp*, dcomp*);
extern bool mltaGeHhDiHGeM(dcomp, mblock<dcomp>**, blcluster*, 
			   blcluster*, int*, unsigned, dcomp*, unsigned, 
			   dcomp*, unsigned);
extern void mltaLtHVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaLtHGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned, dcomp*,
		       unsigned, dcomp*, unsigned);
extern void mltaLtHhVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaLtHhGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned, dcomp*,
			unsigned, dcomp*, unsigned);
extern void mltaUtHVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaUtHGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned, dcomp*,
		       unsigned, dcomp*, unsigned);
extern void mltaUtHhVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaUtHhGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned, dcomp*,
			unsigned, dcomp*, unsigned);
extern void mltaHeHVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaHeHGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned,
		       dcomp*, unsigned, dcomp*, unsigned);
extern void mltaHeHtVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaSyHVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
extern void mltaSyHGeM(dcomp, blcluster*, mblock<dcomp>**, unsigned,
		       dcomp*, unsigned, dcomp*, unsigned);
extern void mltaSyHhVec(dcomp, blcluster*, mblock<dcomp>**, dcomp*, dcomp*);
////mltaGeHGeH.cpp:
extern void mltaGeHGeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
		       mblock<dcomp>**, blcluster*, mblock<dcomp>**, double,
		       unsigned, contBasis<dcomp>* haar=NULL);
extern void mltaGeHGeH_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
			     mblock<dcomp>**, mblock<dcomp>*, double, unsigned,
			     contBasis<dcomp>* haar=NULL);
////mltaGeHhGeH.cpp:
extern void mltaGeHhGeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*, mblock<dcomp>**,
			blcluster*, mblock<dcomp>**, double, unsigned, 
			contBasis<dcomp>* haar=NULL);
extern void mltaGeHhGeH_toMbl(dcomp, blcluster*, mblock<dcomp>**,
			      blcluster*, mblock<dcomp>**, mblock<dcomp>*, 
			      double, unsigned);
extern void mltaGeHhLrM_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned, 
                              contBasis<dcomp>* haar=NULL);
extern void mltaGeHhGeM_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned, 
                              contBasis<dcomp>* haar=NULL);
extern void mltaGeMhGeH_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned, 
                              contBasis<dcomp>* haar=NULL);
extern void mltaLrMhGeH_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned, 
                              contBasis<dcomp>* haar=NULL);
////mltaGeHhGeH_toHeH.cpp:
extern void mltaGeHhGeH_toHeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
			      mblock<dcomp>**, blcluster*, mblock<dcomp>**, 
			      double, unsigned, contBasis<dcomp>* haar=NULL);
////mltaGeHhDiHGeH.cpp:
extern void mltaGeHhDiHGeH(dcomp, mblock<dcomp>**, blcluster*, 
			   blcluster*, int*, blcluster*,
			   blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
			   mblock<dcomp>**, int*, blcluster*, mblock<dcomp>**,
			   blcluster*, mblock<dcomp>**, double, unsigned);
extern void mltaGeHhDiHLrM_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
				 mblock<dcomp>**, int*, blcluster*,
				 mblock<dcomp>**, mblock<dcomp>*, double,
				 unsigned);
extern void mltaGeHhDiHGeM_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
				 mblock<dcomp>**, int*, blcluster*,
				 mblock<dcomp>**, mblock<dcomp>*, double,
				 unsigned);
extern void mltaGeMhDiHGeH_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
				 mblock<dcomp>**, int*, blcluster*,
				 mblock<dcomp>**, mblock<dcomp>*, double,
				 unsigned);
extern void mltaLrMhDiHGeH_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
				 mblock<dcomp>**, int*, blcluster*,
				 mblock<dcomp>**, mblock<dcomp>*, double,
				 unsigned);
extern void mltDiHGeM(unsigned, mblock<dcomp>**, blcluster*, int* const ,
		      dcomp*, dcomp*, unsigned);
////mltaGeHhDiHGeH_toHeH.cpp:
extern void mltaGeHhDiHGeH_toHeH(dcomp, mblock<dcomp>**, blcluster*, 
				 blcluster*, int*, blcluster*,
				 blcluster*, double, unsigned);
extern void mltaGeHhDiHGeH_toHeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
				 mblock<dcomp>**, int*, blcluster*,
				 mblock<dcomp>**, blcluster*, mblock<dcomp>**,
				 double, unsigned);
////mltaGeHGeHh.cpp:
extern void mltaGeHGeHh(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
			mblock<dcomp>**, blcluster*, mblock<dcomp>**,
			double, unsigned);
extern void mltaGeHGeHh_toMbl(dcomp, blcluster*, mblock<dcomp>**,
			      blcluster*, mblock<dcomp>**, mblock<dcomp>*, double, 
			      unsigned);
extern void mltaGeHGeMh_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned);
extern void mltaGeMGeHh_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned);
extern void mltaLrMGeHh_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned);
extern void mltaGeHLrMh_toMbl(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                              mblock<dcomp>**, mblock<dcomp>*, double, unsigned);
////mltaGeHGeHh_toHeH.cpp:
extern void mltaGeHGeHh_toHeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
			      mblock<dcomp>**, blcluster*, mblock<dcomp>**,
			      double, unsigned);
////mltaHeHGeH.cpp:
extern void mltaHeHGeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                       mblock<dcomp>**, blcluster*, mblock<dcomp>**, double, unsigned);
////mltaGeHHeH.cpp:
extern void mltaGeHHeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                       mblock<dcomp>**, blcluster*, mblock<dcomp>**,
                       double, unsigned);
////mltaGeHLtH.cpp:
extern void mltaGeHLtH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
		       mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
		       double eps, unsigned rankmax, contBasis<dcomp>* haar=NULL);
////mltaUtHGeH.cpp:
extern void mltaUtHGeH(dcomp d, blcluster* blA, mblock<dcomp>** A, blcluster* blB,
		       mblock<dcomp>** B, blcluster* blC, mblock<dcomp>** C,
		       double eps, unsigned rankmax, contBasis<dcomp>* haar=NULL);
////mltaGeHUtHh.cpp:
extern void mltaGeHUtHh(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                       mblock<dcomp>**, blcluster*, mblock<dcomp>**,
                       double, unsigned);
////mltaUtHhGeH.cpp:
extern void mltaUtHhGeH(dcomp, blcluster*, mblock<dcomp>**, blcluster*,
                       mblock<dcomp>**, blcluster*, mblock<dcomp>**,
                       double, unsigned);
////mltaUtHhUtH_toHeH.cpp:
extern void mltaUtHhUtH_toHeH(dcomp, blcluster*, mblock<dcomp>**, mblock<dcomp>**,
                         double, unsigned);
////mltaUtHUtHh.cpp:
extern void mltaUtHUtHh(dcomp, blcluster*, mblock<dcomp>**, mblock<dcomp>**,
                         double, unsigned);
extern void mltaUtHLtH(dcomp, blcluster*, mblock<dcomp>**, mblock<dcomp>**,
			mblock<dcomp>**, double, unsigned, 
			contBasis<dcomp>* haar=NULL);

////invertH.cpp:
extern void invertGeH(blcluster*, mblock<dcomp>**, mblock<dcomp>**&,
		      dcomp, unsigned, contBasis<dcomp>* haar=NULL);
extern void invertHeH(blcluster*, mblock<dcomp>**, mblock<dcomp>**&,
		      dcomp, unsigned);
////HLU.cpp:
extern bool HLU(blcluster*, mblock<dcomp>**, mblock<dcomp>**, mblock<dcomp>**,
                double, unsigned, contBasis<dcomp>* haar=NULL);
extern bool HCholesky(blcluster* const, mblock<dcomp>** const, double, unsigned, 
                      contBasis<dcomp>* haar=NULL);
extern bool HUhDU(blcluster*, mblock<dcomp>**, int*, double, unsigned);
extern bool genLUprecond(blcluster*, mblock<dcomp>**, double, unsigned,
                         blcluster*&, mblock<dcomp>**&, mblock<dcomp>**&, bool);
extern bool genCholprecond(blcluster*, mblock<dcomp>**, dcomp, unsigned,
                           blcluster*&, mblock<dcomp>**&, bool);
////TU_solve.cpp:
extern void LtHGeM_solve(blcluster*, mblock<dcomp>**, unsigned, dcomp*,
                         unsigned);
extern void LtHhGeM_solve(blcluster*, mblock<dcomp>**, unsigned, dcomp*,
                          unsigned);
extern void LtHGeH_solve(blcluster*, mblock<dcomp>**, blcluster*,
			 mblock<dcomp>**, mblock<dcomp>**, double, unsigned,
			 contBasis<dcomp>* haar=NULL);
extern void UtHGeM_solve(blcluster*, mblock<dcomp>**, unsigned, dcomp*,
                         unsigned);
extern void UtHGeM_solve(blcluster*, mblock<dcomp>**, int*, unsigned,
			 dcomp*, unsigned);
extern void UtHhGeM_solve(blcluster*, mblock<dcomp>**, unsigned, dcomp*,
                          unsigned);
extern void UtHGeH_solve(blcluster*, mblock<dcomp>**, blcluster*,
			 mblock<dcomp>**, double, unsigned);
extern void UtHhGeH_solve(blcluster*, mblock<dcomp>**, blcluster*,
			  mblock<dcomp>**, double, unsigned, 
			  contBasis<dcomp>* haar=NULL);
extern void GeHUtH_solve(blcluster*, mblock<dcomp>**, blcluster*,
			 mblock<dcomp>**, mblock<dcomp>**, double, unsigned,
			 contBasis<dcomp>* haar=NULL);
extern void UtHhDH_solve(mblock<dcomp>**, blcluster*, blcluster*, int*,
			 double, unsigned);

extern void UtHhDdns_solve(blcluster*, mblock<dcomp>**, int*, unsigned,
			   dcomp*, unsigned);
/*
  extern void HLtHh_solve(blcluster*, mblock<dcomp>**, blcluster*,
  mblock<dcomp>** B, mblock<dcomp>**, double, unsigned);
*/
////psoutH.cpp:
extern void psoutputGeH(std::ofstream&, blcluster*, unsigned, mblock<dcomp>**);
extern void psoutputGeH(std::ofstream&, blcluster**, unsigned, unsigned,
                      mblock<dcomp>**);
extern void psoutputHeH(std::ofstream&, blcluster**, unsigned, unsigned,
			mblock<dcomp>**);
extern void psoutputHeH(std::ofstream&, blcluster*, unsigned, mblock<dcomp>**,
			unsigned level=0, cluster* cl=NULL,
			void (*f)(cluster* cl, const unsigned level, 
				  double* X, const double* const basis)=NULL);
extern void psoutputLtH(std::ofstream&, blcluster*, unsigned, mblock<dcomp>**);
extern void psoutputUtH(std::ofstream&, blcluster*, unsigned, mblock<dcomp>**);

////??
extern void expndHSym(blcluster*, mblock<dcomp>**, blcluster*&,
		      mblock<dcomp>**);
extern void convHeH_toHeM(blcluster**, unsigned, mblock<dcomp>**, dcomp*, unsigned);
extern void mlta_part_vec(dcomp, blcluster**, unsigned, mblock<dcomp>**, dcomp*,
                           dcomp*);
extern void mlta_partT_vec(dcomp, blcluster**, unsigned, mblock<dcomp>**, dcomp*,
                            dcomp*);

extern void mlta_part_vec(dcomp, blcluster**, unsigned,
                           mblock<dcomp>**, dcomp*, dcomp*);
extern void mlta_partT_vec(dcomp, blcluster**, unsigned,
                            mblock<dcomp>**, dcomp*, dcomp*);
extern void mlta_sympart_vec(dcomp, blcluster**, unsigned, mblock<dcomp>**,
                              dcomp*, dcomp*);






// solve L x = b for x, forward substitution, b is destroyed
template<class T>
inline void LtHVec_solve(blcluster* blL, mblock<T>** L, T* b)
{
  LtHGeM_solve(blL, L, 1, b, blL->getn1());
}

// solve L^H x = b for x, forward substitution, b is destroyed
template<class T>
inline void LtHhVec_solve(blcluster* blL, mblock<T>** L, T* b)
{
  LtHhGeM_solve(blL, L, 1, b, blL->getn1());
}

// solve U x = b for x, forward substitution, b is destroyed
template<class T>
inline void UtHVec_solve(blcluster* blU, mblock<T>** U, T* b)
{
  UtHGeM_solve(blU, U, 1, b, blU->getn1());
}

// solve U^H x = b for x, forward substitution, b is destroyed
template<class T>
inline void UtHhVec_solve(blcluster* blU, mblock<T>** U, T* b)
{
  UtHhGeM_solve(blU, U, 1, b, blU->getn2());
}

// solve U^H U x = b for x and store x in b
template<class T>
inline void HCholesky_solve(blcluster* bl, mblock<T>** U, T* b)
{
  UtHhVec_solve(bl, U, b);  // Forward substitution
  UtHVec_solve(bl, U, b);   // Backward substitution
}

// solve U^H D x = b for x, forward substitution, b is destroyed
template<class T>
inline void UtHhDVec_solve(blcluster* blU, mblock<T>** U, int* piv, T* b)
{
  UtHhDdns_solve(blU, U, piv, 1, b, blU->getn2());
}

// solve U x = b for x, backward substitution, b is destroyed
template<class T>
inline void UtHVec_solve(blcluster* blU, mblock<T>** U, int* piv, T* b)
{
  UtHGeM_solve(blU, U, piv, 1, b, blU->getn1());
}

// solve U^H D U x = b for x and store x in b
template<class T>
inline void HUhDU_solve(blcluster* blU, mblock<T>** U, int* piv, T* b)
{
  UtHhDVec_solve(blU, U, piv, b); // forward substitution
  UtHVec_solve(blU, U, piv, b); // backward substitution
}

// solve L U x = b for x and store x in b
template<class T>
inline void HLU_solve(blcluster* bl, mblock<T>** L, mblock<T>** U, T* b)
{
  LtHVec_solve(bl, L, b);    // Forward substitution
  UtHVec_solve(bl, U, b);    // Backward substitution
}

// solve (L U)^H x = b for x and store x in b
template<class T>
inline void HLUh_solve(blcluster* bl, mblock<T>** L, mblock<T>** U, T* b)
{
  UtHhVec_solve(bl, U, b);    // Backward substitution
  LtHhVec_solve(bl, L, b);    // Forward substitution
}


#endif
