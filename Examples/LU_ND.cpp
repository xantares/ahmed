/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <iostream>
#include "cluster_alg.h"
#include "blcluster.h"
#include "basmod.h"
#include "matrix.h"
#include "sparse.h"

#ifdef ENABLE_MPI
  #include "parallel.h"
#else
  #define COUT(X) std::cout << X;
#endif

int main(int argc, char* argv[])
{
  double eta2 = 2;
  unsigned bmin = 50, rankmax = 10000;
  double eps = 1e-14;

  // CRS matrix in full storage
  std::string file = "exmpl_spd.crs";

#ifdef ENABLE_MPI
  MPI::Init(argc, argv);
  initAHMED(MPI::COMM_WORLD);
#endif

  COUT("** Hierarchical Matrix LU Factorization." << std::endl);

  COUT("Reading '" << file << "'" << std::flush);
  CRSMatrix<double> S(file);
  COUT(" -- " << S.n << " number of unknowns." << std::endl);

#ifdef ENABLE_MPI
  unsigned nproc = COMM_AHMED.Get_size();
  unsigned rank = COMM_AHMED.Get_rank();
  COUT(std::endl << "Number of processor: "<< nproc <<std::endl);
#endif
  
  // initialize permutations
  unsigned *op_perm = new unsigned[S.n], *po_perm = new unsigned[S.n];
  for (unsigned k=0; k<S.n; ++k) po_perm[k] = op_perm[k] = k;

  // cluster variables 
  COUT("Clustering variables ..." << std::flush);
  cluster_alg* clTree = new cluster_alg(S.n, S.iA, S.jA);
  clTree->createClusterTree(bmin, op_perm, po_perm);
  COUT(" done." << std::endl);
  
  // block cluster tree
  COUT("Generating block cluster tree ..." << std::flush);
  blcluster* bl = new blcluster(0, 0, S.n, S.n);
  unsigned nblcks = 0;
  bl->subdivide(clTree, clTree, eta2, nblcks);
  COUT(" done." << std::endl);
  
  delete clTree;

  // apply permutation to coefficient matrix A
  CS_perm(S.n, S.A, S.jA, S.iA, po_perm, op_perm, S.A, S.jA, S.iA);

#ifdef ENABLE_MPI
  COMM_AHMED.Barrier();
#endif

  // double -> float fuer Cholesky-Zerlegung in float
  // CRS Matrix kann double bleiben
  mblock<double> **AP;

  // convert stiffness matrix in CRS format to H-matrix
#ifdef ENABLE_MPI
  convCRS_toGeH_ND(nproc, S.A, S.jA, S.iA, eps, bl, AP);
  COMM_AHMED.Barrier();
#else
  convCRS_toGeH(S.A, S.jA, S.iA, eps, bl, AP);
#endif

  COUT("\nFactorizing matrix (eps=" << eps << ") ... " << std::flush);
  mblock<double> **L, **U;
  initLtH_0(bl, L);
  initUtH_0(bl, U);

#ifdef ENABLE_MPI
  if (!HLU_ND(nproc, bl, AP, L, U, eps, rankmax)) {
#else
  if (!HLU(bl, AP, L, U, eps, rankmax)) {
#endif
    std::cout << "no succes." << std::endl;
    exit(1);
  }

  freembls(bl, AP);
  delete [] AP;

#ifdef ENABLE_MPI
  COMM_AHMED.Barrier();
#endif

  COUT("done." << std::endl);
  
#ifdef ENABLE_MPI
  std::cout << "Storage on proc " << rank << " is " 
            << inMB(sizeH_ND(nproc, bl, L, 'L')+sizeH_ND(nproc, bl, U, 'U'))
	    << " MB." << std::endl;
#else
  std::cout << "Storage is " << inMB(sizeH(bl, L, 'L')+sizeH(bl, U, 'U'))
	    << " MB." << std::endl;
#endif

  // --------------------------------------------------------------------------
  // check result
  double *x = new double[S.n];
  blas::load(S.n, 1.0, x);

  // forward/backward solve
#ifdef ENABLE_MPI
  HLU_solve_ND(bl, L, U, x);
#else
  HLU_solve(bl, L, U, x);
#endif
  // now x contains solution of L U x = 1

  // check whether Ax=1
  double *y = new double[S.n];
  blas::load(S.n, 1.0, y);
  
#ifdef ENABLE_MPI
  // generate hierarchy of CRS matrices (ND structure) for parallelization
  CRSblock<double> BBlock;
  BBlock.nnz =  S.iA[S.n];
  BBlock.size = S.n;
  BBlock.A = S.A;
  BBlock.jA = S.jA;
  BBlock.iA = S.iA;
  COUT("Distributing CRS matrix among the processors ...");
  subdivideCRS_ND(nproc, bl, &BBlock);
  COUT(" done." << std::endl);

  double* z = new double[S.n];
  blas::setzero(S.n,z);
  mltaCRSvec_ND(1.0, bl, &BBlock, x, z);
  gatherVec_ND(nproc, bl, z);
  blas::axpy(S.n, -1.0, z, y);
  delete [] z;

  // avoid double free
  S.A = NULL;
  S.jA = NULL;
  S.iA = NULL;
#else
  amuxCRS(S.n, -1.0, x, y, S.iA, S.jA, S.A);
#endif

  COUT("Error is " << blas::nrm2(S.n, y)/sqrt(S.n) << std::endl);

  delete [] x;
  delete [] y;

  // --------------------------------------------------------------------------

  freembls(bl, L);
  freembls(bl, U);
  delete bl;
  
#ifdef ENABLE_MPI
  MPI::Finalize();
#endif

  return 0;
}
