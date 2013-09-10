/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// A += U V^H (A is an H-matrix)
template<class T> static
void addGeHLrMExact_(blcluster* bl, mblock<T>** A, unsigned k,
		     T* U, unsigned ldU, T* V, unsigned ldV)
{
  if (bl->isleaf()) {
    A[bl->getidx()]->addLrM_Exact(k, U, ldU, V, ldV);
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* Up = U + son->getb1() - bl->getb1();
        T* Vp = V + son->getb2() - bl->getb2();
        addGeHLrMExact_(son, A, k, Up, ldU, Vp, ldV);
      }
    }
  }
}

// A += U V^H (A is an H-matrix)
template<class T> static
void addGeHLrM_(blcluster* bl, mblock<T>** A, double eps, unsigned rankmax,
		unsigned k, T* U, unsigned ldU, T* V, unsigned ldV,
		contBasis<T>* haar=NULL)
{
  if (bl->isleaf()) {
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    A[bl->getidx()]->addLrM(k, U, ldU, V, ldV, eps, rankmax, haarInfo,
			    X, bl->getn2(), 
			    X+cols*bl->getn2(), bl->getn1());
    delete haarInfo;
  } else {
    contBasis<T>* haarSon = NULL;
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
	blcluster* son = bl->getson(i, j);
	T* Up = U + son->getb1() - bl->getb1();
	T* Vp = V + son->getb2() - bl->getb2();
	if (haar) haarSon = haar->son(i,j);
	addGeHLrM_(son, A, eps, rankmax, k, Up, ldU, Vp, ldV, haarSon);
	delete haarSon;
      }
    }
  }
}


// A += B  (A is an H-matrix, B is dense)
template<class T> static
void addGeHGeM_(blcluster* bl, mblock<T>** A, T* B, unsigned ldB,
		double eps, unsigned rankmax,
		contBasis<T>* haar=NULL)
{
  if (bl->isleaf()) {
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    A[bl->getidx()]->addGeM(B, ldB, eps, rankmax, haarInfo,
			    X, bl->getn2(), 
			    X+cols*bl->getn2(), bl->getn1());
    delete haarInfo;
  } else {
    contBasis<T>* haarSon = NULL;
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
	blcluster* son = bl->getson(i, j);
	if (haar) haarSon = haar->son(i,j);
	T* Bp = B + ldB*(son->getb2()-bl->getb2())
	  + son->getb1()-bl->getb1();				
	addGeHGeM_(son, A, Bp, ldB, eps, rankmax, haarSon);
	delete haarSon;
      }
    }
  }
}


// A += B (A,B are H-matrices)
template<class T> static
void addGeHGeH_(blcluster* bl, mblock<T>** A, mblock<T>** B,
		double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  if (bl->isleaf()) {
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    A[bl->getidx()]->addMbl(eps, rankmax, B[bl->getidx()], haarInfo,
			    X, bl->getn2(), 
			    X+cols*bl->getn2(), bl->getn1());
    delete haarInfo;
  } else {
    contBasis<T>* haarSon = NULL;
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i){
      for (unsigned j=0; j<ns2; ++j){
	if (bl->getson(i,j)) {
	  if (haar) haarSon = haar->son(i,j);
	  addGeHGeH_(bl->getson(i,j), A, B, eps, rankmax, haarSon);
	  delete haarSon;
	}
      }
    }
  }
}

// A += U V^H (A is a symm. H-matrix)
template<class T> static
void addHeHLrM_(blcluster* bl, mblock<T>** A, double eps, unsigned rankmax,
		unsigned k, T* U, unsigned ldU, T* V, unsigned ldV, 
		contBasis<T>* haar=NULL)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf()) {
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    A[bl->getidx()]->addLrM(k, U, ldU, V, ldV, eps, rankmax, haarInfo,
			    X, bl->getn2(), 
			    X+cols*bl->getn2(), bl->getn1());
    delete haarInfo;
  } else {
    contBasis<T>* haarSon = NULL;
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      if (haar) haarSon = haar->son(i,i);
      blcluster* son = bl->getson(i, i);
      T* Up = U + son->getb1() - bl->getb1();
      T* Vp = V + son->getb2() - bl->getb2();
      addHeHLrM_(son, A, eps, rankmax, k, Up, ldU, Vp, ldV, haarSon);
      delete haarSon;
      for (unsigned j=i+1; j<ns; ++j) {
	if (haar) haarSon = haar->son(i,j);
	son = bl->getson(i, j);
	Up = U + son->getb1() - bl->getb1();
	Vp = V + son->getb2() - bl->getb2();
	addGeHLrM_(son, A, eps, rankmax, k, Up, ldU, Vp, ldV, haarSon);
	delete haarSon;
      }
    }
  }
}

// A += U V^H (A is a symm. H-matrix), no rounding
template<class T> static
void addHeHLrMExact_(blcluster* bl, mblock<T>** A, unsigned k,
		     T* U, unsigned ldU, T* V, unsigned ldV)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf())
    A[bl->getidx()]->addLrM_Exact(k, U, ldU, V, ldV);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* Up = U + son->getb1() - bl->getb1();
      T* Vp = V + son->getb2() - bl->getb2();
      addHeHLrMExact_(son, A, k, Up, ldU, Vp, ldV);
      for (unsigned j=i+1; j<ns; ++j) {
	son = bl->getson(i, j);
	Up = U + son->getb1() - bl->getb1();
	Vp = V + son->getb2() - bl->getb2();
	addGeHLrMExact_(son, A, k, Up, ldU, Vp, ldV);
      }
    }
  }
}

// A += B (A is a symm. H-matrix, B is dense)
template<class T> static
void addHeHGeM_(blcluster* bl, mblock<T>** A, T* B, unsigned ldB,
		double eps, unsigned rankmax, contBasis<T>* haar=NULL)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf()) {
    contLowLevel<T>* haarInfo = NULL;
    unsigned cols = 0;
    T* X = NULL;
    if (haar) {
      haarInfo = new contLowLevel<T>(haar);
      cols = haar->getcols();
      X = haar->getX();
      haar->initX();
    }
    A[bl->getidx()]->addGeM(B, ldB, eps, rankmax, haarInfo, X, bl->getn2(), 
			    X+cols*bl->getn2(), bl->getn1());
    delete haarInfo;
  } else {
    contBasis<T>* haarSon = NULL;
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      if (haar)	haarSon = haar->son(i,i);
      blcluster* son = bl->getson(i, i);
      T* Bp = B + ldB*(son->getb2()-bl->getb2()) + son->getb1()-bl->getb1();
      addHeHGeM_(son, A, Bp, ldB, eps, rankmax, haarSon);
      delete haarSon;
      for (unsigned j=i+1; j<ns; ++j) {
	if (haar) haarSon = haar->son(i,j);
	son = bl->getson(i, j);
	Bp = B + ldB*(son->getb2()-bl->getb2()) + son->getb1()-bl->getb1();
	addGeHGeM_(son, A, Bp, ldB, eps, rankmax, haarSon);
	delete haarSon;
      }
    }
  }
}

template<class T> static
void addHeHLrM_diag_(blcluster* bl, mblock<T>** A, unsigned mult,
		     unsigned k, T* U, unsigned ldU)
{
  if (bl->isleaf())
    A[bl->getidx()]->addLrM_toHeM(mult, k, U, ldU, U, ldU);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* Up = U + son->getb1() - bl->getb1();
      addHeHLrM_diag_(son, A, 2*mult, k, Up, ldU);
    }
  }
}


// adds U V^H in a stable way to a symm. H-matrix
template<class T> static
void addHeHLrM_stab_(blcluster_ref* bl, mblock<T>** A, double eps,
		     unsigned rankmax, unsigned k, T* U,
		     unsigned ldU, T* V, unsigned ldV)
{
  if (bl->isleaf()) {
    T *UR, *VR;
    unsigned kR;

    A[bl->getidx()]->addLrM_rmnd(k, U, ldU, V, ldV, eps, rankmax, kR, UR, VR);

    if (kR>0) {
      addHeHLrM_diag_(bl->getccl()->getdbl(), A, 1, kR, VR, bl->getn2());
      delete [] VR;
      addHeHLrM_diag_(bl->getrcl()->getdbl(), A, 1, kR, UR, bl->getn1());
      delete [] UR;
    }
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; 
	 i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
	blcluster_ref* son = (blcluster_ref*) bl->getson(i, j);
	if (son) {
	  T* Up = U + son->getb1() - bl->getb1();
	  T* Vp = V + son->getb2() - bl->getb2();
	  addHeHLrM_stab_(son, A, eps, rankmax, k, Up, ldU, Vp, ldV);
	}
      }
    }
  }
}


template<class T> static
void addGeHId_(blcluster* bl, mblock<T>** A)
{
  assert(bl->isdbl());
  if (bl->isleaf())
    A[bl->getidx()]->add_dbl_Id();
  else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) addGeHId_(bl->getson(i, i), A);
  }
}




///////////////////////////////////////////////////////////////////////////////
// Instanzen

// double

void addGeHLrM(blcluster* bl, mblock<double>** A, double eps, unsigned rankmax,
	       unsigned k, double* U, unsigned ldU, double* V, unsigned ldV,
	       contBasis<double>* haar)
{
  addGeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addGeHGeM(blcluster* bl, mblock<double>** A, double* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<double>* haar)
{
  addGeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM(blcluster* bl, mblock<double>** A, double eps,
	       unsigned rankmax, unsigned k, double* U, unsigned ldU,
	       double* V, unsigned ldV, contBasis<double>* haar)
{
  addHeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addHeHLrMExact(blcluster* bl, mblock<double>** A,
		    unsigned k, double* U, unsigned ldU, double* V,
		    unsigned ldV)
{
  addHeHLrMExact_(bl, A, k, U, ldU, V, ldV);
}

void addHeHGeM(blcluster* bl, mblock<double>** A, double* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<double>* haar)
{
  addHeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM_stab(blcluster_ref* bl, mblock<double>** A, double eps,
		    unsigned rankmax, unsigned k, double* U, unsigned ldU,
		    double* V, unsigned ldV)
{
  addHeHLrM_stab_(bl, A, eps, rankmax, k, U, ldU, V, ldV);
}

void addGeHGeH(blcluster* bl, mblock<double>** A, mblock<double>** B,
	       double eps, unsigned rankmax, contBasis<double>* haar)
{
  addGeHGeH_(bl, A, B, eps, rankmax, haar);
}


// float

void addGeHLrM(blcluster* bl, mblock<float>** A, double eps, unsigned rankmax,
	       unsigned k, float* U, unsigned ldU, float* V, unsigned ldV,
	       contBasis<float>* haar)
{
  addGeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addGeHGeM(blcluster* bl, mblock<float>** A, float* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<float>* haar)
{
  addGeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM(blcluster* bl, mblock<float>** A, double eps,
	       unsigned rankmax, unsigned k, float* U, unsigned ldU,
	       float* V, unsigned ldV, contBasis<float>* haar)
{
  addHeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addHeHLrMExact(blcluster* bl, mblock<float>** A,
		    unsigned k, float* U, unsigned ldU, float* V,
		    unsigned ldV)
{
  addHeHLrMExact_(bl, A, k, U, ldU, V, ldV);
}

void addHeHGeM(blcluster* bl, mblock<float>** A, float* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<float>* haar)
{
  addHeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM_stab(blcluster_ref* bl, mblock<float>** A, double eps,
		    unsigned rankmax, unsigned k, float* U, unsigned ldU,
		    float* V, unsigned ldV)
{
  addHeHLrM_stab_(bl, A, eps, rankmax, k, U, ldU, V, ldV);
}

void addGeHGeH(blcluster* bl, mblock<float>** A, mblock<float>** B,
	       double eps, unsigned rankmax, contBasis<float>* haar)
{
  addGeHGeH_(bl, A, B, eps, rankmax, haar);
}

// dcomp

void addGeHLrM(blcluster* bl, mblock<dcomp>** A, double eps, unsigned rankmax,
	       unsigned k, dcomp* U, unsigned ldU, dcomp* V, unsigned ldV,
	       contBasis<dcomp>* haar)
{
  addGeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addGeHGeM(blcluster* bl, mblock<dcomp>** A, dcomp* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  addGeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM(blcluster* bl, mblock<dcomp>** A, double eps,
	       unsigned rankmax, unsigned k, dcomp* U, unsigned ldU,
	       dcomp* V, unsigned ldV, contBasis<dcomp>* haar)
{
  addHeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addHeHLrMExact(blcluster* bl, mblock<dcomp>** A,
		    unsigned k, dcomp* U, unsigned ldU, dcomp* V,
		    unsigned ldV)
{
  addHeHLrMExact_(bl, A, k, U, ldU, V, ldV);
}

void addHeHGeM(blcluster* bl, mblock<dcomp>** A, dcomp* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  addHeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM_stab(blcluster_ref* bl, mblock<dcomp>** A, double eps,
		    unsigned rankmax, unsigned k, dcomp* U, unsigned ldU,
		    dcomp* V, unsigned ldV)
{
  addHeHLrM_stab_(bl, A, eps, rankmax, k, U, ldU, V, ldV);
}

void addGeHGeH(blcluster* bl, mblock<dcomp>** A, mblock<dcomp>** B,
	       double eps, unsigned rankmax, contBasis<dcomp>* haar)
{
  addGeHGeH_(bl, A, B, eps, rankmax, haar);
}

// scomp

void addGeHLrM(blcluster* bl, mblock<scomp>** A, double eps, unsigned rankmax,
	       unsigned k, scomp* U, unsigned ldU, scomp* V, unsigned ldV,
	       contBasis<scomp>* haar)
{
  addGeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addGeHGeM(blcluster* bl, mblock<scomp>** A, scomp* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  addGeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM(blcluster* bl, mblock<scomp>** A, double eps,
	       unsigned rankmax, unsigned k, scomp* U, unsigned ldU,
	       scomp* V, unsigned ldV, contBasis<scomp>* haar)
{
  addHeHLrM_(bl, A, eps, rankmax, k, U, ldU, V, ldV, haar);
}

void addHeHLrMExact(blcluster* bl, mblock<scomp>** A,
		    unsigned k, scomp* U, unsigned ldU, scomp* V,
		    unsigned ldV)
{
  addHeHLrMExact_(bl, A, k, U, ldU, V, ldV);
}

void addHeHGeM(blcluster* bl, mblock<scomp>** A, scomp* B, unsigned ldB,
	       double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  addHeHGeM_(bl, A, B, ldB, eps, rankmax, haar);
}

void addHeHLrM_stab(blcluster_ref* bl, mblock<scomp>** A, double eps,
		    unsigned rankmax, unsigned k, scomp* U, unsigned ldU,
		    scomp* V, unsigned ldV)
{
  addHeHLrM_stab_(bl, A, eps, rankmax, k, U, ldU, V, ldV);
}

void addGeHGeH(blcluster* bl, mblock<scomp>** A, mblock<scomp>** B,
	       double eps, unsigned rankmax, contBasis<scomp>* haar)
{
  addGeHGeH_(bl, A, B, eps, rankmax, haar);
}


void addGeHId(blcluster* bl, mblock<double>** A)
{
  addGeHId_(bl, A);
}

void addGeHId(blcluster* bl, mblock<float>** A)
{
  addGeHId_(bl, A);
}

void addGeHId(blcluster* bl, mblock<dcomp>** A)
{
  addGeHId_(bl, A);
}

void addGeHId(blcluster* bl, mblock<scomp>** A)
{
  addGeHId_(bl, A);
}
