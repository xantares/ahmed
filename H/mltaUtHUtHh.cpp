/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// C += d U U^h   (U UtH-matrix, C HSym-matrix)
template<class T> static
void mltaUtHUtHh_(T d, blcluster* bl, mblock<T>** U, mblock<T>** C,
		  double eps, unsigned rankmax)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    C[idx]->mltaUtMUtMh_toHeM(d, U[idx]->getdata());
  }else {                            // C is not a leaf, hence U isn't
    unsigned ns = bl->getnrs();
    assert(ns==bl->getncs());
    
    for (unsigned i=0; i<ns; ++i) {
      blcluster* blii = bl->getson(i, i);
      mltaUtHUtHh_(d, blii, U, C, eps, rankmax);
      for (unsigned k=i+1; k<ns; ++k) {
        blcluster* blik = bl->getson(i, k);
        mltaGeHGeHh_toHeH(d, blik, U, blik, U, blii, C, eps, rankmax);
      }

      for (unsigned j=i+1; j<ns; ++j) {
        blcluster* blij = bl->getson(i, j);
        blcluster* bljj = bl->getson(j, j);
	mltaGeHUtHh(d, blij, U, bljj, U, blij, C, eps, rankmax);
        for (unsigned k=j+1; k<ns; ++k) {
          blcluster* blik = bl->getson(i, k);
          blcluster* bljk = bl->getson(j, k);
          mltaGeHGeHh(d, blik, U, bljk, U, blij, C, eps, rankmax);
        }
      }
    }
  }
}

// C += d U L   (U UtH-matrix, L LtH-matrix, C H-matrix)
template<class T> static
void mltaUtHLtH_(T d, blcluster* bl, mblock<T>** U, mblock<T>** L, 
		  mblock<T>** C, double eps, unsigned rankmax,
		  contBasis<T>* haar=NULL)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    C[idx]->mltaUtMLtM(d, U[idx]->getdata(), L[idx]->getdata());
  }else {                            // C is not a leaf, hence U isn't
    unsigned ns = bl->getnrs();
    assert(ns==bl->getncs());
    contBasis<T>* haarSon1=NULL;
    contBasis<T>* haarSon2=NULL;
    for (unsigned i=0; i<ns; ++i) {
      if(haar) haarSon1 = haar->son(i,i);
      blcluster* blii = bl->getson(i, i);
      mltaUtHLtH_(d, blii, U, L, C, eps, rankmax, haarSon1);
      for (unsigned k=i+1; k<ns; ++k) {
        blcluster* blik = bl->getson(i, k);
	blcluster* blki = bl->getson(k, i);
        mltaGeHGeH(d, blik, U, blki, L, blii, C, eps, rankmax, haarSon1);
      }
      delete haarSon1;
      
      for (unsigned j=i+1; j<ns; ++j) {
        blcluster* blij = bl->getson(i, j);
	blcluster* blji = bl->getson(j, i);
        blcluster* bljj = bl->getson(j, j);
	if(haar) haarSon1 = haar->son(i,j);
	mltaGeHLtH(d, blij, U, bljj, L, blij, C, eps, rankmax, haarSon1);
	if(haar) haarSon2 = haar->son(j,i);
	mltaUtHGeH(d, bljj, U, blji, L, blji, C, eps, rankmax, haarSon2);
        for (unsigned k=j+1; k<ns; ++k) {
          blcluster* blik = bl->getson(i, k);
	  blcluster* bljk = bl->getson(j, k);
	  blcluster* blki = bl->getson(k, i);
          blcluster* blkj = bl->getson(k, j);
          mltaGeHGeH(d, blik, U, blkj, L, blij, C, eps, rankmax, haarSon1);
	  mltaGeHGeH(d, bljk, U, blki, L, blji, C, eps, rankmax, haarSon2);
        }
	delete haarSon1;
	delete haarSon2;
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Instanzen


void mltaUtHUtHh(double d, blcluster* bl, mblock<double>** U,
                  mblock<double>** C, double eps, unsigned rankmax)
{
  mltaUtHUtHh_(d, bl, U, C, eps, rankmax);
}

void mltaUtHUtHh(float d, blcluster* bl, mblock<float>** U, mblock<float>** C,
                  double eps, unsigned rankmax)
{
  mltaUtHUtHh_(d, bl, U, C, eps, rankmax);
}

void mltaUtHUtHh(dcomp d, blcluster* bl, mblock<dcomp>** U, mblock<dcomp>** C,
                  double eps, unsigned rankmax)
{
  mltaUtHUtHh_(d, bl, U, C, eps, rankmax);
}

void mltaUtHUtHh(scomp d, blcluster* bl, mblock<scomp>** U, mblock<scomp>** C,
                  double eps, unsigned rankmax)
{
  mltaUtHUtHh_(d, bl, U, C, eps, rankmax);
}

void mltaUtHLtH(double d, blcluster* bl, mblock<double>** U, mblock<double>** L, 
		 mblock<double>** C, double eps, unsigned rankmax, 
		 contBasis<double>* haar)
{
  mltaUtHLtH_(d, bl, U, L, C, eps, rankmax, haar);
}

void mltaUtHLtH(float d, blcluster* bl, mblock<float>** U, mblock<float>** L,
		 mblock<float>** C, double eps, unsigned rankmax,
		 contBasis<float>* haar)
{
  mltaUtHLtH_(d, bl, U, L, C, eps, rankmax, haar);
}

void mltaUtHLtH(dcomp d, blcluster* bl, mblock<dcomp>** U, mblock<dcomp>** L,
		 mblock<dcomp>** C, double eps, unsigned rankmax,
		 contBasis<dcomp>* haar)
{
  mltaUtHLtH_(d, bl, U, L, C, eps, rankmax, haar);
}

void mltaUtHLtH(scomp d, blcluster* bl, mblock<scomp>** U, mblock<scomp>** L, 
		 mblock<scomp>** C, double eps, unsigned rankmax,
		 contBasis<scomp>* haar)
{
  mltaUtHLtH_(d, bl, U, L, C, eps, rankmax, haar);
}
