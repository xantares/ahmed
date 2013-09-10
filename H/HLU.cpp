/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// truncated A=LU decomposition of an H-matrix A, A is destroyed
// L and U have to be initialized with initLtH_0() and initUtH_0(), resp.
// returns true if successful, otherwise false

template<class T> static
bool HLU_(blcluster* const bl, mblock<T>** const A, mblock<T>** const L,
          mblock<T>** const U, const double eps, const unsigned rankmax,
	  contBasis<T>* haar=NULL)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    int inf = A[idx]->decomp_LU(L[idx], U[idx]);
    if (inf!=0) return false;
  } else {
    assert(bl->getnrs()==bl->getncs());
    unsigned i, j, k, ns = bl->getnrs();
    contBasis<T>* haarSon=NULL;

    for (i=0; i<ns; ++i) {
      blcluster *son = bl->getson(i, i);
      if(haar) haarSon = haar->son(i,i);

      for (k=0; k<i; ++k) {
        blcluster *son1 = bl->getson(i, k), *son2 = bl->getson(k, i);
        mltaGeHGeH((T) -1.0, son1, L, son2, U, son, A, eps, rankmax, haarSon);
      }
      if (!HLU_(son, A, L, U, eps, rankmax, haarSon)) return false;
      delete haarSon;

      for (j=i+1; j<ns; ++j) {
        blcluster* sonU = bl->getson(i, j);
	if(haar) haarSon = haar->son(i,j);
        for (k=0; k<i; ++k) {
          blcluster *son1 = bl->getson(i, k), *son2 = bl->getson(k, j);
          mltaGeHGeH((T) -1.0, son1, L, son2, U, sonU, A, eps, rankmax,
		  haarSon);
        }
        LtHGeH_solve(son, L, sonU, A, U, eps, rankmax, haarSon);
	delete haarSon;

        blcluster* sonL = bl->getson(j, i);
	if(haar) haarSon = haar->son(j,i);
        for (k=0; k<i; ++k) {
          blcluster *son1 = bl->getson(j, k), *son2 = bl->getson(k, i);
          mltaGeHGeH((T) -1.0, son1, L, son2, U, sonL, A, eps, rankmax,
		  haarSon);
        }
        GeHUtH_solve(son, U, sonL, A, L, eps, rankmax, haarSon);
	delete haarSon;
      }
    }
  }
  return true;
}

// truncated A = U^H D U decomp. of a symm. H-matrix A
// (notice: only the upper part)
// upon completion A contains the upper triangular factor
// returns true if successful, otherwise false
template<class T> static
bool HUhDU_(blcluster* const bl, mblock<T>** const A, int* const piv,
	    const double eps, const unsigned rankmax)
{
  if (bl->isleaf()) {
    if(A[bl->getidx()]->decomp_UhDU(piv+bl->getb1()))
      return false;
  } else {
    assert(bl->getnrs()==bl->getncs());
    unsigned i, j, k, ns = bl->getnrs();

    for (i=0; i<ns; ++i) {
      blcluster *son = bl->getson(i,i);
      for (k=0; k<i; ++k) {
	blcluster *son1 = bl->getson(k,i);
	blcluster *sonD = bl->getson(k,k);
	mltaGeHhDiHGeH_toHeH((T) -1.0, A, son1, sonD, piv, son1, son, eps,
			  rankmax);
      }
      if (!HUhDU_(son, A, piv, eps, rankmax)) return false;

      for (j=i+1; j<ns; ++j) {
	blcluster *sonU = bl->getson(i,j);
	for (k=0; k<i; ++k) {
	  blcluster *son1 = bl->getson(k,i), *son2 = bl->getson(k,j),
	    *sonD = bl->getson(k,k);
	  mltaGeHhDiHGeH((T) -1.0, A, son1, sonD, piv, son2, sonU, eps, rankmax);
	}
	//solve
	UtHhDH_solve(A, sonU, son, piv, eps, rankmax);
      }
    }
  }
  return true;
}

// truncated Cholesky A=U^HU decomp. of a symm. H-matrix A
// upon completion A contains the upper triangular factor
// returns true if successful, otherwise false
template<class T> static
bool HCholesky_(blcluster* const bl, mblock<T>** const A,
                const double eps, const unsigned rankmax,
                contBasis<T>* haar=NULL)
{
  if (bl->isleaf()) {
    if (A[bl->getidx()]->decomp_Cholesky())
      return false;
  } else {
    assert(bl->getnrs()==bl->getncs());
    unsigned i, j, k, ns = bl->getnrs();

    contBasis<T>* haarSon=NULL;
    for (i=0; i<ns; ++i) {
      if(haar) haarSon = haar->son(i,i);

      blcluster *son = bl->getson(i, i);
      for (k=0; k<i; ++k) {
	blcluster *son1 = bl->getson(k, i);
	mltaGeHhGeH_toHeH((T) -1.0, son1, A, son1, A, son, A, eps, rankmax,
			  haarSon);
      }
      if (!HCholesky_(son, A, eps, rankmax, haarSon)) return false;
      delete haarSon;

      for (j=i+1; j<ns; ++j) {
	blcluster* sonU = bl->getson(i, j);
	if(haar) haarSon = haar->son(i,j);

	for (k=0; k<i; ++k) {
	  blcluster *son1 = bl->getson(k, i), *son2 = bl->getson(k, j);
	  mltaGeHhGeH((T) -1.0, son1, A, son2, A, sonU, A, eps, rankmax, 
		      haarSon);
	}
	UtHhGeH_solve(son, A, sonU, A, eps, rankmax, haarSon);
	delete haarSon;
      }
    }
  }
  return true;
}

// generate LU preconditioner with precision delta
bool genLUprecond(blcluster* bl, mblock<double>** A, double delta,
                  unsigned max_rank, blcluster* &blAggl, mblock<double>** &L,
                  mblock<double>** &U, bool precoarse)
{
  mblock<double>** Aggl;
  allocmbls(bl, Aggl);
  copyH(bl, A, Aggl);

  if (precoarse) {
    blAggl = new blcluster(bl);
    agglH(blAggl, Aggl, delta, max_rank);
  } else blAggl = bl;

  allocmbls(blAggl, L);
  initLtH_0(blAggl, L);

  allocmbls(blAggl, U);
  initUtH_0(blAggl, U);

  bool inf = HLU_(blAggl, Aggl, L, U, delta, max_rank);
  freembls(blAggl, Aggl);

  if (!inf) std::cout << "genLUprecond: no success." << std::endl;

  return inf;
}

// generate LU preconditioner with precision delta
bool genLUprecond(blcluster* bl, mblock<float>** A, double delta,
                  unsigned max_rank, blcluster* &blAggl, mblock<float>** &L,
                  mblock<float>** &U, bool precoarse)
{
  mblock<float>** Aggl;
  allocmbls(bl, Aggl);
  copyH(bl, A, Aggl);

  if (precoarse) {
    blAggl = new blcluster(bl);
    agglH(blAggl, Aggl, delta, max_rank);
  } else blAggl = bl;

  allocmbls(blAggl, L);
  initLtH_0(blAggl, L);

  allocmbls(blAggl, U);
  initUtH_0(blAggl, U);

  bool inf = HLU_(blAggl, Aggl, L, U, delta, max_rank);
  freembls(blAggl, Aggl);

  if (!inf) std::cout << "genLUprecond: no success." << std::endl;

  return inf;
}

// generate LU preconditioner with precision delta
bool genLUprecond(blcluster* bl, mblock<double>** A, double delta,
                  unsigned max_rank, blcluster* &blAggl, mblock<float>** &L,
                  mblock<float>** &U, bool precoarse)
{
  mblock<float>** Aggl;
  allocmbls(bl, Aggl);
  copyH(bl, A, Aggl);

  if (precoarse) {
    blAggl = new blcluster(bl);
    agglH(blAggl, Aggl, delta, max_rank);
  } else blAggl = bl;

  allocmbls(blAggl, L);
  initLtH_0(blAggl, L);

  allocmbls(blAggl, U);
  initUtH_0(blAggl, U);

  bool inf = HLU_(blAggl, Aggl, L, U, delta, max_rank);
  freembls(blAggl, Aggl);

  if (!inf) std::cout << "genLUprecond: no success." << std::endl;

  return inf;
}


// generate LU preconditioner with precision delta
bool genLUprecond(blcluster* bl, mblock<dcomp>** A, double delta,
                  unsigned max_rank, blcluster* &blAggl, mblock<dcomp>** &L,
                  mblock<dcomp>** &U, bool precoarse)
{
  mblock<dcomp>** Aggl;
  allocmbls(bl, Aggl);
  copyH(bl, A, Aggl);

  if (precoarse) {
    blAggl = new blcluster(bl);
    agglH(blAggl, Aggl, delta, max_rank);
  } else blAggl = bl;

  allocmbls(blAggl, L);
  initLtH_0(blAggl, L);

  allocmbls(blAggl, U);
  initUtH_0(blAggl, U);

  bool inf = HLU_(blAggl, Aggl, L, U, delta, max_rank);
  freembls(blAggl, Aggl);

  if (!inf) std::cout << "genLUprecond: no success." << std::endl;

  return inf;
}

// generate LU preconditioner with precision delta
bool genLUprecond(blcluster* bl, mblock<scomp>** A, double delta,
                  unsigned max_rank, blcluster* &blAggl, mblock<scomp>** &L,
                  mblock<scomp>** &U, bool precoarse)
{
  mblock<scomp>** Aggl;
  allocmbls(bl, Aggl);
  copyH(bl, A, Aggl);

  if (precoarse) {
    blAggl = new blcluster(bl);
    agglH(blAggl, Aggl, delta, max_rank);
  } else blAggl = bl;

  allocmbls(blAggl, L);
  initLtH_0(blAggl, L);

  allocmbls(blAggl, U);
  initUtH_0(blAggl, U);

  bool inf = HLU_(blAggl, Aggl, L, U, delta, max_rank);
  freembls(blAggl, Aggl);

  if (!inf) std::cout << "genLUprecond: no success." << std::endl;

  return inf;
}



// generate LU preconditioner with precision delta
bool genLUprecond(blcluster* bl, mblock<dcomp>** A, double delta,
                  unsigned max_rank, blcluster* &blAggl, mblock<scomp>** &L,
                  mblock<scomp>** &U, bool precoarse)
{
  mblock<scomp>** Aggl;
  allocmbls(bl, Aggl);
  copyH(bl, A, Aggl);

  if (precoarse) {
    blAggl = new blcluster(bl);
    agglH(blAggl, Aggl, delta, max_rank);
  } else blAggl = bl;

  allocmbls(blAggl, L);
  initLtH_0(blAggl, L);

  allocmbls(blAggl, U);
  initUtH_0(blAggl, U);

  bool inf = HLU_(blAggl, Aggl, L, U, delta, max_rank);
  freembls(blAggl, Aggl);

  if (!inf) std::cout << "genLUprecond: no success." << std::endl;

  return inf;
}

// generate Cholesky preconditioner with precision delta, A symm. H-matrix
bool genCholprecond(blcluster* bl, mblock<double>** A, double delta,
                    unsigned max_rank, blcluster* &blU, mblock<double>** &U,
                    bool precoarse)
{
  allocmbls(bl, U);
  copyH(bl, A, U);

  if (precoarse) {
    blU = new blcluster(bl);
    agglH(blU, U, delta, max_rank);
  } else blU = bl;

  bool inf = HCholesky_(blU, U, delta, max_rank);
  if (!inf) std::cout << "genCholprecond: no success." << std::endl;

  return inf;
}

// generate Cholesky preconditioner with precision delta, A symm. H-matrix
bool genCholprecond(blcluster* bl, mblock<double>** A, double delta,
                    unsigned max_rank, blcluster* &blU, mblock<float>** &U,
                    bool precoarse)
{
  allocmbls(bl, U);
  copyH(bl, A, U);

  if (precoarse) {
    blU = new blcluster(bl);
    agglH(blU, U, delta, max_rank);
  } else blU = bl;

  bool inf = HCholesky_(blU, U, delta, max_rank);
  if (!inf) std::cout << "genCholprecond: no success." << std::endl;

  return inf;
}

// generate Cholesky preconditioner with precision delta, A symm. H-matrix
bool genCholprecond(blcluster* bl, mblock<float>** A, double delta,
                    unsigned max_rank, blcluster* &blU, mblock<float>** &U,
                    bool precoarse)
{
  allocmbls(bl, U);
  copyH(bl, A, U);

  if (precoarse) {
    blU = new blcluster(bl);
    agglH(blU, U, delta, max_rank);
  } else blU = bl;

  bool inf = HCholesky_(blU, U, delta, max_rank);
  if (!inf) std::cout << "genCholprecond: no success." << std::endl;

  return inf;
}


///////////////////////////////////////////////////////////////////////////////
// Instanzen
//


bool HLU(blcluster* const bl, mblock<double>** const A,
         mblock<double>** const L, mblock<double>** const U,
         const double eps, const unsigned rankmax,
	 contBasis<double>* haar)
{
  return HLU_(bl, A, L, U, eps, rankmax, haar);
}


bool HLU(blcluster* const bl, mblock<float>** const A,
         mblock<float>** const L, mblock<float>** const U,
         const double eps, const unsigned rankmax, 
	 contBasis<float>* haar)
{
  return HLU_(bl, A, L, U, eps, rankmax, haar);
}


bool HLU(blcluster* const bl, mblock<dcomp>** const A,
         mblock<dcomp>** const L, mblock<dcomp>** const U,
         const double eps, const unsigned rankmax,
	 contBasis<dcomp>* haar)
{
  return HLU_(bl, A, L, U, eps, rankmax, haar);
}


bool HLU(blcluster* const bl, mblock<scomp>** const A,
         mblock<scomp>** const L, mblock<scomp>** const U,
         const double eps, const unsigned rankmax,
	 contBasis<scomp>* haar)
{
  return HLU_(bl, A, L, U, eps, rankmax, haar);
}



bool HCholesky(blcluster* const bl, mblock<double>** const A,
               const double eps, const unsigned rankmax, 
               contBasis<double>* haar)
{
  return HCholesky_(bl, A, eps, rankmax, haar);
}

bool HCholesky(blcluster* const bl, mblock<float>** const A,
               const double eps, const unsigned rankmax, 
               contBasis<float>* haar)
{
  return HCholesky_(bl, A, eps, rankmax, haar);
}

bool HCholesky(blcluster* const bl, mblock<dcomp>** const A,
               const double eps, const unsigned rankmax, 
               contBasis<dcomp>* haar)
{
  return HCholesky_(bl, A, eps, rankmax, haar);
}

bool HCholesky(blcluster* const bl, mblock<scomp>** const A,
               const double eps, const unsigned rankmax, 
               contBasis<scomp>* haar)
{
  return HCholesky_(bl, A, eps, rankmax, haar);
}


bool HUhDU(blcluster* const bl, mblock<double>** const A, int* const piv,
	   const double eps, const unsigned rankmax)
{ return HUhDU_(bl, A, piv, eps, rankmax); }

bool HUhDU(blcluster* const bl, mblock<float>** const A, int* const piv,
	   const double eps, const unsigned rankmax)
{ return HUhDU_(bl, A, piv, eps, rankmax); }

bool HUhDU(blcluster* const bl, mblock<dcomp>** const A, int* const piv,
	   const double eps, const unsigned rankmax)
{ return HUhDU_(bl, A, piv, eps, rankmax); }

bool HUhDU(blcluster* const bl, mblock<scomp>** const A, int* const piv,
	   const double eps, const unsigned rankmax)
{ return HUhDU_(bl, A, piv, eps, rankmax); }
