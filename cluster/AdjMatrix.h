/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef ADJMATRIX_H
#define ADJMATRIX_H

#include "basmod.h"
#include "sparse.h"
//#include "blas.h"
#include <limits>

/** AdjMat implements an adjacency matrix used to represent the (possible 
   weighted) edges of a graph. Due to adjacency matrix are mostly sparse 
   the implementation makes use of compressed row storage format.
*/

class AdjMat
{
public:
  /** constructor with explicit data elements in aA
      \param s size of the quadratic matrix
      \param aiA array of length s+1 which holds pointers in jA,
        iA[k] points to the first non-zero column-entry of row k,
    	iA[k]-1 points accordingly to the last non-zero column-entry of row k,
     	the last entry of iA (iA[s]) takes the number of non zero entries(nnz)
     \param ajA array of length nnz, each entry is a colum-index
     \param aA  data-array of length nnz of type unsigned (weights)
  */
   AdjMat(unsigned s, unsigned *aiA, unsigned *ajA, unsigned *aA=NULL)
     : iA(aiA), jA(ajA), A(aA), n(s) { }

   void makeSymmetric() {
      // store upper triangular mat values
      genAdjM(n, iA, jA);
      // mirror the upper triangular part into lower
      genFullAdjMat(n, iA, jA);
   }

   int isElement(unsigned r, unsigned c) const;
   
   unsigned getElement(unsigned r, unsigned c) const;

   /** Is the matrix graph connected? */
   bool isGraphConnected(unsigned &diam) const {
     return ::isConnected(n, iA, jA, diam);
   }
   
   void getDiam(unsigned &diam) const;

   void getDiam(unsigned &diam, const unsigned* const perm,
		const unsigned* const iperm, unsigned beg, unsigned end) const;

  /** destructor */
  ~AdjMat() {
    if (A != NULL) delete [] A;
    delete [] iA;
    delete [] jA;
  }
   
  /** creates a graph and returns pointer for Metis routines */
  //  void buildGraph (unsigned &n, int *iA, int* &jA) const;

  /** getRows () returns the number of rows (value of attribute m) */
  unsigned getRows() const { return n; }

  /** getCols () returns the number of cols (value of attribute n) */
  unsigned getCols() const { return n; }

  unsigned getBeginRow(unsigned idx) const { return iA[idx]; }
   
  unsigned getCol(unsigned idx) const { return jA[idx]; }
   
  unsigned getNNZ() const { return iA[n]; }

  unsigned* getRowIdxArray() const { return iA; }
  unsigned* getColIdxArray() const { return jA; }
  unsigned* getDataArray() const { return A; }

  /** getMat returns the (irreducible) block [beg,end-1] x [beg,end-1] 
      respecting the permutation. If it is necessary to insert edges to make
      the block irreducible, the number of edges and the edges itself are returned.
      \param beg first row/column
      \param end one after last row/column  
      \param op_perm permutation -> original
      \param po_perm original -> permutation
      \param nedges number of edges
      \param edges pointer to a field of size 2*nedges
  */

  //   AdjMat* getMat(unsigned beg, unsigned end, const unsigned* const op_perm,
  //      const unsigned* const po_perm, unsigned &nedges, unsigned* &edges) const;

   /** getMat returns the (possibly reducible) block [beg,end-1] x [beg,end-1] 
      respecting the permutation. 
      \param beg first row/column, it is supposed that 0 <= beg <= n
      \param end one after last row/column , it is supposed that beg <= end <= n  
      \param op_perm permutation -> original
      \param po_perm original -> permutation
      \returns pointer to an AdjMat object
	 */
   AdjMat* getMat(unsigned beg, unsigned end, const unsigned* const op_perm,
		  const unsigned* const po_perm) const;

   /*
   AdjMat* getMat(unsigned beg, unsigned end,
		  unsigned* op_perm, unsigned* po_perm,
		  unsigned &nadj_mats, AdjMat** &adj_mats) const;
   */

   /** add the edges from list into the matrix */
   //   void addEdges(std::list<unsigned>* edges);

private:
   //   AdjMat* getSubMat(unsigned beg, unsigned end, const unsigned* const op_perm,
   //		     const unsigned* const po_perm, unsigned k,
   //		     const unsigned* const comp) const;

   /** getComponents computes the number of components, the nodes and the approx. 
      diameter of each particular component of the graph represented by this adjacency matrix.
      \param beg input: beginning index of the sub block
      \param end input: beginning index of the next sub block
      \param op_perm input: permutation of the rows and columns
      \param po_perm input: inverse permutation to op_perm
      \param ncomp output: number of components
      \param comp output: unsigned field [0, ..., end-beg] stores the belonging
         to the component 
      \param diams output: unsigned field [0, ..., ncomp] stores the approx diameters of components
    */
  void getComponents(unsigned beg, unsigned end, 
      const unsigned* const op_perm, const unsigned* const po_perm, 
      unsigned &ncomp, unsigned* comp) const; 

  /** copy constructor
      \param src another AdjMat with data for initializiation
  */
  AdjMat(const AdjMat& src) {}

  /** assignment operator
      \param rhs a instance of AdjacencyCRSMatrix
  */
  AdjMat& operator=(AdjMat& rhs) { return *this; }
  
  unsigned *iA, *jA, *A;
  unsigned n;
};

#endif
