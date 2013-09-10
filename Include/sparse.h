/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <assert.h>
#include "cmplx.h"

template<class T>
struct litem {
  unsigned i;             // the row, resp. the column
  T d;                    // the value
  litem* next;            // the next item in the list
  litem(unsigned id, T e, litem *succ=0) : i(id), d(e), next(succ) { }
};


// stores a hierarchy of CRS blocks (for ND representation of a CRS matrix)

template<class T> struct CRSblock {
  unsigned nnz, size;
  unsigned *iA, *jA;
  T *A;
  CRSblock* sons;
  CRSblock(): nnz(0), size(0), iA(NULL), jA(NULL), A(NULL), sons(NULL) { }

  ~CRSblock() {
    delete [] iA;
    delete [] jA;
    delete [] A;
    delete [] sons;
  }
};

extern void add_item(litem<double>*& list, unsigned i, double d);
extern void add_item(litem<dcomp>*& list, unsigned i, dcomp d);
extern unsigned len_list(litem<double>* p);
extern unsigned len_list(litem<dcomp>* p);
extern void pack_list(litem<double>* p, double*& pA, unsigned*& piA);
extern void pack_list(litem<dcomp>* p, dcomp*& pA, unsigned*& piA);

extern void CS_write(char*, unsigned, unsigned*, unsigned*, double*);
extern void CS_write(char*, unsigned, unsigned*, unsigned*, dcomp*);
extern void CS_read(char*, unsigned&, unsigned*&, unsigned*&, double*&);
extern void CS_read(char*, unsigned&, unsigned*&, unsigned*&, dcomp*&);

template<class T> void
CS_write(std::ostream &os, unsigned n, unsigned* iA, unsigned* jA, T* A)
{
  os.write((char*) &n, sizeof(unsigned));
  os.write((char*) iA, (n+1)*sizeof(unsigned));
  os.write((char*) jA, iA[n]*sizeof(unsigned));
  os.write((char*) A, iA[n]*sizeof(T));
}

template<class T> void
CS_read(std::istream &is, unsigned &n, unsigned* &iA, unsigned* &jA, T* &A)
{
  is.read((char*) &n, sizeof(unsigned));
  if (iA != NULL) {
    delete [] iA;
    delete [] jA;
    delete [] A;
  }
  iA = new unsigned[n+1];
  assert(iA!=NULL);
  is.read((char*) iA, (n+1)*sizeof(unsigned));

  jA = new unsigned[iA[n]];
  assert(jA!=NULL);
  is.read((char*) jA, iA[n]*sizeof(unsigned));

  A = new T[iA[n]];
  assert(A!=NULL);
  is.read((char*) A, iA[n]*sizeof(T));

#ifndef NDEBUG
  // do simple checks
  if (iA[0]!=0)
    std::cerr << std::endl
              << "CRS matrix: array iA doesn't start with 0" << std::endl;

  unsigned i = 0;
  while (i<iA[n] && jA[i]<n) ++i;
  if (i<iA[n])
    std::cerr << std::endl
              << "CRS matrix: the " << i << "th entry of jA has the value "
              << jA[i] << ", which is out of bounds." << std::endl;
#endif
}

extern void amuxCRS(unsigned, double, double*, double*,
                    unsigned*, unsigned*, double*);
extern void amuxCRS(unsigned, float, float*, float*,
                    unsigned*, unsigned*, float*);
extern void amuxCRS(unsigned, dcomp, dcomp*, dcomp*,
                    unsigned*, unsigned*, dcomp*);
extern void amuxCRS(unsigned, scomp, scomp*, scomp*,
                    unsigned*, unsigned*, scomp*);

extern void atmuxCRS(unsigned, double, double*, double*,
                     unsigned*, unsigned*, double*);
extern void atmuxCRS(unsigned, float, float*, float*,
                     unsigned*, unsigned*, float*);
extern void atmuxCRS(unsigned, dcomp, dcomp*, dcomp*,
                     unsigned*, unsigned*, dcomp*);
extern void atmuxCRS(unsigned, scomp, scomp*, scomp*,
                     unsigned*, unsigned*, scomp*);

/** function multiplies the symmetric CRS-Matrix (iA, jA, A) with vector x and
 stores the result in y */
extern void amuxSymmCRS(unsigned n, double d, double* x, double* y,
                        unsigned* iA, unsigned* jA, double* A);
/** function multiplies the symmetric CRS-Matrix (iA, jA) with vector x and
 stores the result in y. Only the structure is used. */
extern void amuxSymmCRS(unsigned, double, double*, double*,
                        unsigned*, unsigned*);

extern void genAdjM(unsigned, unsigned*&, unsigned*&);
/** converts the structure of a matrix into a symmetric structure */
extern void genFullAdjMat(unsigned, unsigned* &, unsigned* &);

extern void genAdjM(unsigned n, unsigned* &iA, unsigned* &jA, unsigned* &A);
/** converts a matrix into a symmetric matrix including the matrix entries
  (weights of edges of graph) */
extern void genFullAdjMat(unsigned n, unsigned* &iA, unsigned* &jA, unsigned* &A);
/** changes matrix entries of type double into the unsigned integer
 interval [a,b] via a linear transformation */
void transMatEntries ( unsigned a, unsigned b, unsigned nnz, double* A,
                       unsigned* idata );

extern bool isConnected(unsigned, unsigned*, unsigned*, unsigned&);

/**
 * function computes a farthest node from the start node s within a subgraph of the graph.
 * The graph is represents as an adjacency matrix, that is available in compressed row storage format.
 * @param perm the permutation of the nodes
 * @param iperm the inverse permutation of the nodes
 * @param beg the beg of the cluster
 * @param end the end of the cluster
 * @param s start node s
 * @param nrows number of nodes of the global graph = number of rows of the adjacency matrix
 * @param row_ptr the row pointer of the compressed row storage format
 * @param col_idx the column entries of the compressed row storage format
 * @param diam an approximation of the diameter
 * @return a farthest node from the start node s
 */
unsigned getRemoteNode (const unsigned* const perm, const unsigned* const iperm,
                        unsigned beg, unsigned end, unsigned s, unsigned nrows,
                        unsigned *row_ptr, unsigned *col_idx, unsigned &diam);

unsigned getRemotlyNode (unsigned s, unsigned n, unsigned *iA, unsigned *jA, unsigned &diam);

extern void CS_transp(unsigned, double*, unsigned*, unsigned*, double*,
                      unsigned*, unsigned*);
extern void CS_perm(unsigned, double*, unsigned*, unsigned*, unsigned*,
                    unsigned*, double*&, unsigned*&, unsigned*&);
extern void CS_perm(unsigned, dcomp*, unsigned*, unsigned*, unsigned*,
                    unsigned*, dcomp*&, unsigned*&, unsigned*&);

/** removes lower part from a CRS matrix */
extern void CRS2CRSSym(unsigned, double*, unsigned*, unsigned*,
                       double*&, unsigned*&, unsigned*&);

/** adds (symmetric) lower triangular part to upper triangular CRS matrix */
void CRSSym2CRS(unsigned, unsigned*&, unsigned*&, double*&);

extern "C" {
  void ilutp_(unsigned*, double*, unsigned*, unsigned*, unsigned*, double*,
  double*, unsigned*, double*, unsigned*, unsigned*, unsigned*,
  double*, unsigned*, unsigned*, int*);
  void lusol_(unsigned*, double*, double*, double*, unsigned*, unsigned*);
}

extern unsigned BiCGStabILU(unsigned, double*, unsigned*, unsigned*, double*,
                            double*, double&, unsigned&, double*, unsigned*,
                            unsigned*);
extern unsigned BiCGStabILU(unsigned, float*, unsigned*, unsigned*, float*,
                            float*, float&, unsigned&, float*, unsigned*,
                            unsigned*);

/** Function checks, if the distance between cluster t1 and t2 is smaller than max_dist.
  \param membership unsigned array of length nrows, it should be initialised as follows:
  	membership[k] = 1, if index k belongs to cluster t1;
  	membership[k] = 2, if index k belongs to cluster t2;
  	membership[k] = 0, else
  \param nbrsf unsigned array of length 2*nrows, it should be initialised as follows:
  	nbrsf[i] = t1[i], i = 1, ..., |t1|
  \param size_t1 the number of elements in cluster t1
  \param max_dist = 1/eta * min (diam t1, diam t2)
  \param nrows number of rows of adjacency matrix
  \param row_ptr the row pointer of CRS-Format
  \param col_idx the column indices of CRS-Format
  \returns true if dist (t1, t2) <= max_dist, else false
 */
bool getDist (unsigned *membership, unsigned *nbrsf, unsigned size_t1, unsigned &max_dist,
              unsigned nrows, unsigned *row_ptr, unsigned *col_idx);

/** Function cancels the membership of nodes in U(t1). */
void reverseBFS (unsigned *membership, unsigned *nbrsf, unsigned size_t1,
                 unsigned nrows, unsigned *row_ptr, unsigned *col_idx);

unsigned getDist (unsigned s, unsigned t, unsigned nrows, unsigned *row_ptr, unsigned *col_idx);

/** generate a diagonal preconditioner from a CRS matrix */
extern bool generateDiagPrecond (const unsigned n, const double* const A, 
				 const unsigned* const jA, const unsigned* const iA,
				 double* diag);

/** generate a scaled CRS matrix */
void scaleCRSSym(const unsigned n, const double* const A, const unsigned* const jA,
     const unsigned* const iA, double*& diag, double* A_New);

#endif
