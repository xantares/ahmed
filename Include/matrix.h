/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <iostream>
#include <fstream>
#include <limits>

template<class T> struct Matrix {
  unsigned m, n;                                      // dimensions
  Matrix(unsigned n1, unsigned n2) : m(n1), n(n2) { }
  Matrix() : m(0), n(0) { }
  virtual void amux(T d, T *x, T *y) const = 0;       // y +=d*Ax
  
  //The default version of atmux applies for symmetric matrices only!
  virtual void atmux(T d, T *x, T *y) const
  {  amux(d, x, y);  }

  virtual void precond_apply(T*) const { }            // apply preconditioner
  virtual ~Matrix() { }
};



#include "blas.h"
/*
template<typename T> static
void unitfunc__(comp<T>& a, const double& b)
{  a.re = a.im = (T) b;  }

static void unitfunc__(float& a, const double& b)
{  a = (float) b;  }
static void unitfunc__(double& a, const double& b)
{  a = b;  }
*/

//computes spectral norm of matrix A = sqrt(sigma(A^H A))
//applies for symmetric matrices also
template<typename T> 
double spctrlnrm(Matrix<T>& A, const unsigned it_num = 10)
{
  const unsigned n = A.n, m = A.m;
  T *x = new T[n], *y = new T[m];
  assert(x!=NULL&&y!=NULL);
  blas::load(n, (T)1., x);  //start vector
  double d;

  for(unsigned i = 0; i < it_num; ++i) {
    //compute A^H A x
      blas::setzero(m, y);
      A.amux((T)1., x, y);  //y = Ax
      blas::setzero(n, x);
      A.atmux((T)1., y, x); //x = A^HAx

      d = blas::nrm2(n, x);       // normalize
      const double e = 1./d;
      blas::scal(n, (T) e, x); 
  }

  delete [] y;
  delete [] x;

  return sqrt(d); //return sqroot of max eigenvalue
}
//symmetric case: returns |sigma(A)|
template<typename T> 
double spctrlnrm_sym(Matrix<T>& A, const unsigned it_num = 10)
{
  const unsigned n(A.n);
  assert(n==A.m);
  T *x = new T[n], *y = new T[n];
  assert(x!=NULL&&y!=NULL);
  blas::load(n, (T)1., x);  //start vector
  double d;

  for(unsigned i = 0; i < it_num; ++i) {
    //compute  A x
      blas::setzero(n, y);
      A.amux((T)1., x, y);  //y = Ax

      d = blas::nrm2(n, y);       // normalize
      const double e = 1./d;
      blas::copy(n, y, x);
      blas::scal(n, (T) e, x); 
  }

  delete [] y;
  delete [] x;

  return d; //return max eigenvalue in modulus
}

//computes spectral norm of (A - B) = sqrt(sigma((A-B)^H(A-B)))
//applies for hermitian matrices also
template<typename T> 
double spctrlnrm(Matrix<T>& A, Matrix<T>& B, const unsigned it_num = 10)
{
  assert(A.n == B.n && A.m == B.m);
  const unsigned n = A.n, m = A.m;
  T *x = new T[n], *y = new T[m];
  assert(x!=NULL&&y!=NULL);
  blas::load(n, (T) 1., x);  //start vector
  double d;

  for(unsigned i = 0; i < it_num; ++i) {
    //compute (A-B)^H (A-B) x
      blas::setzero(m, y);
      A.amux((T)1., x, y);  
      B.amux((T)-1., x, y); //y = (A-B)x
      blas::setzero(n, x);
      A.atmux((T)1., y, x); 
      B.atmux((T)-1., y, x);//x = (A-B)^H (A-B)x

      d = blas::nrm2(n, x);       // normalize
      if(d==0.) break;
      const double e = 1./d;
      blas::scal(n, (T) e, x);
  }

  delete [] y;
  delete [] x;

  return sqrt(d); //return sqroot of max eigenvalue
}
//symmetric case: returns |sigma(A-B)|
template<typename T> 
double spctrlnrm_sym(Matrix<T>& A, Matrix<T>& B, const unsigned it_num = 10)
{
  const unsigned n(A.n);
  assert(n == B.n && n == A.m && n == B.m);
  T *x = new T[n], *y = new T[n];
  assert(x!=NULL&&y!=NULL);
  blas::load(n, (T) 1., x);  //start vector
  double d;

  for(unsigned i = 0; i < it_num; ++i) {
    //compute (A-B) x
      blas::setzero(n, y);
      A.amux((T)1., x, y);  
      B.amux((T)-1., x, y); //y = (A-B)x

      d = blas::nrm2(n, y);       // normalize
      const double e = 1./d;
      blas::copy(n, y, x);
      blas::scal(n, (T) e, x);
  }

  delete [] y;
  delete [] x;

  return d; //return max eigenvalue in modulus
}



template<typename T> static
void permvec(unsigned n, const unsigned* perm, const T* a, T* b)
{  while(n--) *b++ = a[*perm++];  }

template<typename T> static
void invpermvec(unsigned n, const unsigned* perm, const T* a, T* b)
{  while(n--) b[*perm++] = *a++;  }

template<typename T> static
void addinvpermvec(unsigned n, const unsigned* perm, const T* a, T* b) 
{  while(n--) b[*perm++] += *a++;  }


//for matrices with differing partitions
//applies for symmetric matrices also
template<typename T> 
double spctrlnrm(Matrix<T>& A, const unsigned* op_permrowA, const unsigned* op_permcolA,
		 Matrix<T>& B, const unsigned* op_permrowB, const unsigned* op_permcolB,
		 const unsigned it_num = 10)
{
  assert(A.n == B.n && A.m == B.m);
  const unsigned n = A.n, m = A.m;
  T *x = new T[n], *xperm = new T[n], *y = new T[m], *yperm = new T[m];
  assert(x && xperm && y && yperm);
  blas::load(n, (T)1., x);  //start vector
  double d;

  for(unsigned i = 0; i < it_num; ++i) {
    //compute (A-B)^H (A-B) x
      blas::setzero(m, yperm);
      permvec(n, op_permcolA, x, xperm);  
      A.amux((T)1., xperm, yperm);  
      invpermvec(m, op_permrowA, yperm, y);//y = Ax

      permvec(n, op_permcolB, x, xperm);
      blas::setzero(m, yperm);
      B.amux((T)-1., xperm, yperm); 
      addinvpermvec(m, op_permrowB, yperm, y);//y = (A-B)x
      
      permvec(m, op_permrowA, y, yperm);
      blas::setzero(n, xperm);
      A.atmux((T)1., yperm, xperm); 
      invpermvec(n, op_permcolA, xperm, x);//x = A^H(A-B)x

      permvec(m, op_permrowB, y, yperm);
      blas::setzero(n, xperm);
      B.atmux((T)-1., yperm, xperm);
      addinvpermvec(n, op_permcolB, xperm, x); //x = (A-B)^H (A-B)x

      d = blas::nrm2(n, x);       // normalize
      if(d==0.) break;
      const double e = 1./d;
      blas::scal(n, (T) e, x);
  }

  delete [] yperm;
  delete [] y;
  delete [] xperm;
  delete [] x;

  return sqrt(d); //return sqroot of max eigenvalue
}
//symmetric case:
template<typename T> 
double spctrlnrm_sym(Matrix<T>& A, const unsigned* op_permA, 
		     Matrix<T>& B, const unsigned* op_permB,
		     const unsigned it_num = 10)
{
  const unsigned n(A.n);
  assert(n == B.n && n == A.m && n == B.m);
  T *x = new T[n], *xperm = new T[n], *y = new T[n], *yperm = new T[n];
  assert(x && xperm && y && yperm);
  blas::load(n, (T)1., x);  //start vector
  double d;

  for(unsigned i = 0; i < it_num; ++i) {
    //compute (A-B) x
      blas::setzero(n, yperm);
      permvec(n, op_permA, x, xperm);  
      A.amux((T)1., xperm, yperm);  
      invpermvec(n, op_permA, yperm, y);//y = Ax

      permvec(n, op_permB, x, xperm);
      blas::setzero(n, yperm);
      B.amux((T)-1., xperm, yperm); 
      addinvpermvec(n, op_permB, yperm, y);//y = (A-B)x

      d = blas::nrm2(n, y);       // normalize
      const double e = 1./d;
      blas::copy(n, y, x);
      blas::scal(n, (T) e, x);
  }

  delete [] yperm;
  delete [] y;
  delete [] xperm;
  delete [] x;

  return d; //return max eigenvalue in modulus
}



//computes spectral norm of (I-AB) in order to check the contraction rate
template<typename T>
double cntrctrate(Matrix<T>& A, Matrix<T>& B, const unsigned it_num = 10)
{
  assert(A.n==B.n && A.n == A.m && A.m==B.m);
  const unsigned n = A.n;
  T* const x = new T[n];
  T* const y = new T[n];
  blas::load(n, (T) 1.0, x);   // starting vector
  double d;

  for (unsigned l=0; l < it_num; ++l) {
    // compute (I-AB)^T (I-AB) x
    blas::setzero(n, y);
    A.amux((T) 1.0, x, y);
    B.amux((T) -1.0, y, x);
    blas::setzero(n, y);
    A.atmux((T) 1.0, x, y);
    B.atmux((T) -1.0, y, x);

    d = blas::nrm2(n, x);       // normalize
    double e = 1.0/d;
    blas::scal(n, (T) e, x);
  }

  delete[] y;
  delete[] x;

  return sqrt(d);
}


#include "mblock.h"
#include "blcluster.h"
#include "H.h" //mltaGeHVec, mltaGeHhVec, mltaHeHVec, mltaSyHVec

//note that blclTree must not be deleted before call of destructor
template<class T> struct HeHMatrix : public Matrix<T> {
  mblock<T>** blcks;
  blcluster* blclTree;

  HeHMatrix(unsigned n, blcluster* tree = NULL) : Matrix<T>(n, n),
    blcks(NULL), blclTree(tree) { }

  ~HeHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(T d, T* x, T* y) const {
    mltaHeHVec(d, blclTree, blcks, x, y);
  }

  void precond_apply(T* x) const { }
};


//note that blclTree must not be deleted before call of destructor
//symmetric complex H-matrix
template<class T> struct SyHMatrix : public Matrix<T > {
  mblock<T >** blcks;
  blcluster* blclTree;

  SyHMatrix(unsigned n, blcluster* tree = NULL) : Matrix<T >(n, n),
    blcks(NULL), blclTree(tree) { }

  ~SyHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(T d, T* x, T* y) const {
    mltaSyHVec(d, blclTree, blcks, x, y);
  }

  void atmux(T d, T* x, T* y) const {
    mltaSyHhVec(d, blclTree, blcks, x, y);
  }

  void precond_apply(T* x) const { }
};


//note that blclTree must not be deleted before call of destructor
template<class T> struct GeHMatrix : public Matrix<T> {
  mblock<T>** blcks;
  blcluster* blclTree;

  GeHMatrix(unsigned n1, unsigned n2, blcluster* tree = NULL) :
      Matrix<T>(n1, n2), blcks(NULL), blclTree(tree) { }

  ~GeHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(T d, T* x, T* y) const {
    mltaGeHVec(d, blclTree, blcks, x, y);
  }
  void atmux(T d, T* x, T* y) const {
    mltaGeHhVec(conj(d), blclTree, blcks, x, y);
  }
  void precond_apply(T* x) const { }
};




#include "sparse.h"

template<class T>
struct CRSMatrix : public Matrix<T> {
  unsigned *iA, *jA;
  T* A;

  CRSMatrix(std::string &fname): Matrix<T>(), iA(NULL), jA(NULL), A(NULL) {
    std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
    if (in) {
      CS_read(in, Matrix<T>::m, iA, jA, A);
      Matrix<T>::n = Matrix<T>::m;
      in.close();
    } else {
      std::cout << "cannot open " << fname << std::endl;
      exit(1);
    }
  }

  CRSMatrix(unsigned n1): Matrix<T>(n1, n1), iA(NULL), jA(NULL), A(NULL) {}

  ~CRSMatrix() {
    delete [] A;
    delete [] jA;
    delete [] iA;
  }

  void amux(T d, T* x, T *y) const {
    amuxCRS(Matrix<T>::n, d, x, y, iA, jA, A);
  }

  void atmux(T d, T* x, T *y) const {
    atmuxCRS(Matrix<T>::n, d, x, y, iA, jA, A);
  }

  void precond_apply(T* x) const { }
};



struct CRSMatrixSym : public Matrix<double> {
  unsigned *iA, *jA;
  double* A;

  CRSMatrixSym(std::string &fname): Matrix<double>(), iA(NULL), jA(NULL), A(NULL) {
    std::ifstream in(fname.c_str(), std::ios::in | std::ios::binary);
    if (in) {
      CS_read(in, m, iA, jA, A);
      n = m;
      in.close();
    } else {
      std::cout << "cannot open " << fname << std::endl;
      exit(1);
    }
  }

  CRSMatrixSym(unsigned n1): Matrix<double>(n1, n1), iA(NULL), jA(NULL), A(NULL) {}

  ~CRSMatrixSym() {
    delete [] A;
    delete [] jA;
    delete [] iA;
  }

  void amux(double d, double* x, double *y) const {
    amuxSymmCRS(n, d, x, y, iA, jA, A);
  }
};

struct CCSMatrix : public Matrix<double> {
  double* A;
  unsigned *iA, *jA;

  CCSMatrix(unsigned n1) : Matrix<double>(n1, n1), A(NULL),
      iA(NULL), jA(NULL) { }

  ~CCSMatrix() {
    delete [] A;
    delete [] iA;
    delete [] jA;
  }

  void amux(double d, double* x, double *y) const {
    atmuxCRS(n, d, x, y, jA, iA, A);
  }

  void atmux(double d, double* x, double *y) const {
    amuxCRS(n, d, x, y, jA, iA, A);
  }

  void precond_apply(double* x) const { }
};

struct DiagPrecond : public CRSMatrix<double> {
  double *inv_diag;

  DiagPrecond (std::string &fname) : CRSMatrix<double>(fname), inv_diag(NULL) {
    inv_diag = new double[n];
    if (!generateDiagPrecond (n, A, jA, iA, inv_diag)) {
      std::cout << "Could not create diagonal preconditioner" << std::endl;
      exit(1);
    }
  }

  void precond_apply(double* x) const {
    for (unsigned k=0; k<n; ++k) x[k] = inv_diag[k]*x[k];
  }

  ~DiagPrecond() {
    delete [] inv_diag;
  }
};

#include "H.h"
struct FEMMatrix : public CRSMatrix<double> {
  mblock<float> **L, **U ;        // the preconditioner in H-Matrix format
  blcluster *bl;                  // the corresponding block cluster tree
  float* xf;
  unsigned *op_perm, *po_perm;    // permutation required for preconditioner

  FEMMatrix(unsigned n1) : CRSMatrix<double>(n1), xf(new float[n]) { }
  FEMMatrix (std::string &fname) : CRSMatrix<double>(fname), xf(new float[n]) {}

  ~FEMMatrix() {
    delete [] xf;
  }

  void precond_apply(double* x) const {
    for (unsigned i=0; i<m; ++i) xf[po_perm[i]] = (float) x[i];
    HLU_solve(bl, L, U, xf);
    for (unsigned i=0; i<m; ++i) x[i] = xf[po_perm[i]];
  }
};

struct FEMMatrixp : public CRSMatrix<double> {
  mblock<float>** U;       // the preconditioner in H-Matrix format
  mblock<float>** L;       // the preconditioner in H-Matrix format
  blcluster *blLU;         // the corresponding block cluster tree
  float* xf;

  FEMMatrixp(unsigned n1) : CRSMatrix<double>(n1) {
    xf = new float[n];
  }
  ~FEMMatrixp() {
    delete [] xf;
  }

  void precond_apply(double* x) const {
    for (unsigned i=0; i<n; ++i) xf[i] = (float) x[i];
    HLU_solve(blLU, L, U, xf);
    for (unsigned i=0; i<n; ++i) x[i] = xf[i];
  }
};

struct SymFEMMatrixCCS : public CCSMatrix {
  mblock<double>** U;       // the preconditioner in H-Matrix format
  blcluster *blU;          // the corresponding block cluster tree

  SymFEMMatrixCCS(unsigned n1) : CCSMatrix(n1) { }
  ~SymFEMMatrixCCS() {}

  void precond_apply(double* x) const {
    HCholesky_solve(blU, U, x);
  }
};

struct FEMMatrixpInvCCS : public CCSMatrix {
  mblock<double>** H;      // the H-Matrix approximation of the stiffness matrix
  mblock<double>** Inv;    // the preconditioner in H-Matrix format
  blcluster *bl;           // the corresponding block cluster tree
  double *xf;
  
 FEMMatrixpInvCCS(unsigned n1) : CCSMatrix(n1), xf(new double[n1]) 
    {
      
    }
  ~FEMMatrixpInvCCS() {
    delete [] xf;
  }

  void precond_apply(double* x) const {
    blas::copy(n, x, xf);
    blas::setzero(n, x);
    mltaHeHVec(1.0, bl, Inv, xf, x);
  }
};

template<class T> struct SymFEMMatrixCRS : public CRSMatrix<T> {
  mblock<T>** U;       // the preconditioner in H-Matrix format
  blcluster *blU;          // the corresponding block cluster tree

  SymFEMMatrixCRS(unsigned n1) : CRSMatrix<T>(n1) { }
  SymFEMMatrixCRS(std::string &fname) : CRSMatrix<T>(fname) {}

  void precond_apply(T* x) const {
    HCholesky_solve(blU, U, x);
  }
};

// S typedef for CRS Matrix
// T typedef for preconditioner
template<class S, class T> struct SymFEMMatrixCRSDiag : public CRSMatrix<S> {
  mblock<T>** U;                 // the preconditioner in H-Matrix format
  blcluster *blU;                // the corresponding block cluster tree
  T *xf;
  T* diag;                       // scaling of matrix

  SymFEMMatrixCRSDiag(unsigned n1) 
    : CRSMatrix<S>(n1), xf(new T[n1]) 
    { }
  SymFEMMatrixCRSDiag(std::string &fname) 
    : CRSMatrix<S>(fname), xf(new T[Matrix<S>::n]) 
    {}
  ~SymFEMMatrixCRSDiag() {
    delete [] xf;
  }

  void precond_apply(S* x) const {
    for (unsigned i=0; i<Matrix<S>::n; ++i) xf[i] = x[i]*diag[i];
    HCholesky_solve(blU, U, xf);
    for (unsigned i=0; i<Matrix<S>::n; ++i) x[i] = xf[i]*diag[i];
  }
};

// S typedef for CRS Matrix
// T typedef for preconditioner
template<class S, class T> struct SymFEMMatrixpCRSDiag : public CRSMatrix<S> {
  mblock<T>** U;                 // the preconditioner in H-Matrix format
  blcluster *blU;                // the corresponding block cluster tree
  T *xf;
  unsigned *op_perm, *po_perm;   // permutation required for preconditioner
  T* diag;                       // scaling of matrix

  SymFEMMatrixpCRSDiag(unsigned n1) 
    : CRSMatrix<S>(n1), xf(new T[n1]) 
    { }
  SymFEMMatrixpCRSDiag(std::string &fname) 
    : CRSMatrix<S>(fname), xf(new T[Matrix<S>::n]) 
    {}
  ~SymFEMMatrixpCRSDiag() {
    delete [] xf;
  }

  void precond_apply(S* x) const {
    for (unsigned i=0; i<Matrix<S>::n; ++i) xf[po_perm[i]] = x[i]*diag[i];
    HCholesky_solve(blU, U, xf);
    for (unsigned i=0; i<Matrix<S>::n; ++i) x[i] = xf[po_perm[i]]*diag[i];
  }
};

template<class T> struct SymFEMMatrixpCRS : public CRSMatrix<T> {
  mblock<T>** U;       // the preconditioner in H-Matrix format
  blcluster *blU;          // the corresponding block cluster tree
  T *xf;
  unsigned *op_perm, *po_perm;   // permutation required for preconditioner

  SymFEMMatrixpCRS(unsigned n1) : CRSMatrix<T>(n1), xf(new T[n1]) { }
  SymFEMMatrixpCRS(std::string &fname) : CRSMatrix<T>(fname), xf(new T[Matrix<T>::n]) {}
  ~SymFEMMatrixpCRS() {
    delete [] xf;
  }

  void precond_apply(T* x) const {
    for (unsigned i=0; i<Matrix<T>::n; ++i) xf[po_perm[i]] = x[i];
    HCholesky_solve(blU, U, xf);
    for (unsigned i=0; i<Matrix<T>::n; ++i) x[i] = xf[po_perm[i]];
  }
};

struct SymFEMMatrixp : public CRSMatrix<double> {
  mblock<float>** U;       // the preconditioner in H-Matrix format
  blcluster*blU;          // the corresponding block cluster tree
  float *xf;

  SymFEMMatrixp(unsigned n1) : CRSMatrix<double>(n1), xf(new float[n1]) { }
  ~SymFEMMatrixp() {
    delete [] xf;
  }

  void precond_apply(double* x) const {
    for (unsigned i=0; i<n; ++i) xf[i] = (float) x[i];
    HCholesky_solve(blU, U, xf);
    for (unsigned i=0; i<n; ++i) x[i] = xf[i];
  }
};

struct SymFEMMatrix : public CRSMatrix<double> {
  mblock<float>** U;             // the preconditioner in H-Matrix format
  blcluster *blU;                // the corresponding block cluster tree
  float* xf;
  unsigned *op_perm, *po_perm;   // permutation required for preconditioner

  SymFEMMatrix(unsigned n1) : CRSMatrix<double>(n1), xf(new float[n1]) { }
  SymFEMMatrix(std::string &fname) : CRSMatrix<double>(fname), xf(new float[n]) {}
  ~SymFEMMatrix() {
    delete [] xf;
  }

  void precond_apply(double* x) const {
    for (unsigned i=0; i<m; ++i) xf[po_perm[i]] = (float) x[i];
    HCholesky_solve(blU, U, xf);
    for (unsigned i=0; i<m; ++i) x[i] = xf[po_perm[i]];
  }
};

#endif
