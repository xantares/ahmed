/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "H.h"

template<class T> static
void unify_sons_(blcluster* bl, mblock<T>** A, double eps, unsigned max_rank,
                 contBasis<T>* haar=NULL)
{
  unsigned cost_pre, ns1 = bl->getnrs(), ns2 = bl->getncs();
  
  // unify first block row
  blcluster* son = bl->getson(0, 0);
  const unsigned n0 = son->getn2();
  unsigned mi = son->getn1(), m = mi, n = n0;

  mblock<T>* R = new mblock<T>(mi, n0);
  R->copy(*A[son->getidx()]);
  cost_pre = A[son->getidx()]->size();

  contLowLevel<T>* haarInfo = NULL;
  unsigned cols = 0;
  T* X = NULL;
  if (haar!=NULL) {
    haarInfo = new contLowLevel<T>(haar);
    cols = haar->getcols();
    X = haar->getX();
    assert(haar->getccl()->size()==bl->getn2() && haar->getrcl()->size()==bl->getn1());
    haar->initX();
  }
	
  for (unsigned j=1; j<ns2; ++j) {
    mblock<T>* R1 = R;
    son = bl->getson(0, j);
    n += son->getn2();
    R = new mblock<T>(mi, n);
    R->unify_cols(eps, max_rank, *R1, *A[son->getidx()],
                  haarInfo, X, bl->getn2(), 
                  X+cols*bl->getn2(), bl->getn1());
    cost_pre += A[son->getidx()]->size();
    delete R1;
  }

  // now the other rows
  for (unsigned i=1; i<ns1; ++i) {

    son = bl->getson(i, 0);
    mi = son->getn1();
    n = n0;
    mblock<T>* Rz = new mblock<T>(mi, n0);
    Rz->copy(*A[son->getidx()]);
    cost_pre += A[son->getidx()]->size();

    for (unsigned j=1; j<ns2; ++j) {
      mblock<T>* R1 = Rz;
      son = bl->getson(i, j);
      n += son->getn2();
      Rz = new mblock<T>(mi, n);
      Rz->unify_cols(eps, max_rank, *R1, *A[son->getidx()], 
                     haarInfo, X, bl->getn2(), 
                     X+cols*bl->getn2()+son->getb1()-bl->getb1(), 
                     bl->getn1());
      cost_pre += A[son->getidx()]->size();

      delete R1;
    }

    mblock<T>* R1 = R;
    m += mi;
    R = new mblock<T>(m, n);
    R->unify_rows(eps, max_rank, *R1, *Rz,
                  haarInfo, X, bl->getn2(), 
                  X+cols*bl->getn2(), bl->getn1());
    delete R1;
    delete Rz;
  }
	
  if (cost_pre>=R->size()) {

    // let bl become a leaf; the new block gets the index of the first son
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        unsigned idx = bl->getson(i, j)->getidx();
        if (i!=0 || j!=0) {
          delete A[idx];
          A[idx] = NULL;
        }
        else {
          bl->setidx(idx);
          A[idx]->copy(*R);
        }
        delete bl->getson(i, j);
      }
    }

    bl->setsons(0, 0, NULL);
  }
  delete haarInfo;
  delete R;
}


// Agglomeration of blocks
// tries dense blocks only under certain circumstances)
template<class T> static
void aggl_light_H_(blcluster* bl, mblock<T>** A, double eps, unsigned max_rank)
{
  blcluster* son;

  if (bl->isnleaf()) {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    bool all_sons = true;              // indicates whether all sons are leaves
    bool *hsdns= new bool[ns1+ns2];    // dense blocks per block row/column

    unsigned i;
    for (i=0; i<ns1+ns2; ++i) hsdns[i] = false;

    for (unsigned k=0; k<ns1; ++k) {
      for (unsigned l=0; l<ns2; ++l) {
        son = bl->getson(k, l);
        if (son) {
          aggl_light_H_(son, A, eps, max_rank);

          if (son->isnleaf())
            all_sons = false;
          else
            if (A[son->getidx()]->isGeM()) hsdns[k] = hsdns[l+ns1] = true;
        }
        else all_sons = false;
      }
    }

    if (all_sons) {
      unsigned co1 = 0, co2 = 0;
      for (i=0; i<ns1; ++i)
        if (hsdns[i]) co1 += bl->getson(i, 0)->getn1();
      for (i=0; i<ns2; ++i)
        if (hsdns[i+ns1]) co2 += bl->getson(0, i)->getn2();

      if (co1*co2<bl->getn1()*bl->getn2()) unify_sons_(bl, A, eps, max_rank);
    }

    delete [] hsdns;
  }
}


// Agglomaration of blocks
// tries all blocks
template<class T> static
void agglH_(blcluster* bl, mblock<T>** A, double eps, unsigned max_rank,
            contBasis<T>* haar=NULL)
{
  blcluster* son;

  if (bl->isnleaf()) {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    bool all_sons = true;              // indicates whether all sons are leaves
    contBasis<T>* haarSon = NULL;
    for (unsigned k=0; k<ns1; ++k) {
      for (unsigned l=0; l<ns2; ++l) {
	son = bl->getson(k, l);
	if (son) {
	  if (haar) haarSon = haar->son(k,l);
	  agglH_(son, A, eps, max_rank, haar);
	  delete haarSon;
	  if (son->isnleaf()) all_sons = false;
	}
	else all_sons = false;
      }
    }

    if (all_sons && !bl->isdbl()) {
      if (haar)	haarSon = haar->son();
      if (all_sons && !bl->isdbl())
	unify_sons_(bl, A, eps, max_rank, haarSon);
      delete haarSon;
    }
  }
}


template<class T> static
void fill_gaps_(blcluster* bl, mblock<T>** A, unsigned nblnew, unsigned& i)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    if (idx>=nblnew) {
      while (A[i]) ++i;
      A[i] = A[idx];
      A[idx] = NULL;
      bl->setidx(i);
    }
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son) fill_gaps_(son, A, nblnew, i);
      }
    }
  }
}

// ----------------------------------------------------------------------------
// Instanzen

void unify_sons(blcluster* bl, mblock<double>** A, double eps,
                unsigned max_rank)
{
  return unify_sons_(bl, A, eps, max_rank);
}


void unify_sons(blcluster* bl, mblock<float>** A, double eps,
                unsigned max_rank)
{
  return unify_sons_(bl, A, eps, max_rank);
}


void unify_sons(blcluster* bl, mblock<dcomp>** A, double eps,
                unsigned max_rank)
{
  return unify_sons_(bl, A, eps, max_rank);
}


void unify_sons(blcluster* bl, mblock<scomp>** A, double eps,
                unsigned max_rank)
{
  return unify_sons_(bl, A, eps, max_rank);
}



void fill_gaps(blcluster* bl, mblock<double>** A, unsigned nblnew, unsigned& i)
{
  fill_gaps_(bl, A, nblnew, i);
}

void fill_gaps(blcluster* bl, mblock<dcomp>** A, unsigned nblnew, unsigned& i)
{
  fill_gaps_(bl, A, nblnew, i);
}



void agglH(blcluster* bl, mblock<double>** A, double eps, unsigned max_rank,
           contBasis<double>* haar)
{
  unsigned nblcks_old = bl->nleaves();
  contBasis<double>* haarSon=NULL;
  if(haar) haarSon = haar->son(1);
  agglH_(bl, A, eps, max_rank, haarSon);
  delete haarSon;
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}


void agglH(blcluster* bl, mblock<float>** A, double eps, unsigned max_rank,
           contBasis<float>* haar)
{
  unsigned nblcks_old = bl->nleaves();
  contBasis<float>* haarSon=NULL;
  if(haar) haarSon = haar->son(1);
  agglH_(bl, A, eps, max_rank, haarSon);
  delete haarSon;
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}

void agglH(blcluster* bl, mblock<dcomp>** A, double eps, unsigned max_rank,
           contBasis<dcomp>* haar)
{
  unsigned nblcks_old = bl->nleaves();
  contBasis<dcomp>* haarSon=NULL;
  if(haar) haarSon = haar->son(1);
  agglH_(bl, A, eps, max_rank, haarSon);
  delete haarSon;
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}

void agglH(blcluster* bl, mblock<scomp>** A, double eps, unsigned max_rank,
           contBasis<scomp>* haar)
{
  unsigned nblcks_old = bl->nleaves();
  contBasis<scomp>* haarSon=NULL;
  if(haar) haarSon = haar->son(1);
  agglH_(bl, A, eps, max_rank, haarSon);
  delete haarSon;
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}


void aggl_light_H(blcluster* bl, mblock<double>** A, double eps,
                  unsigned max_rank)
{
  unsigned nblcks_old = bl->nleaves();
  aggl_light_H_(bl, A, eps, max_rank);
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}

void aggl_light_H(blcluster* bl, mblock<float>** A, double eps,
                  unsigned max_rank)
{
  unsigned nblcks_old = bl->nleaves();
  aggl_light_H_(bl, A, eps, max_rank);
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}

void aggl_light_H(blcluster* bl, mblock<dcomp>** A, double eps,
		  unsigned max_rank)
{
  unsigned nblcks_old = bl->nleaves();
  aggl_light_H_(bl, A, eps, max_rank);
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}

void aggl_light_H(blcluster* bl, mblock<scomp>** A, double eps,
                  unsigned max_rank)
{
  unsigned nblcks_old = bl->nleaves();
  aggl_light_H_(bl, A, eps, max_rank);
  unsigned nblcks_new = bl->nleaves();

  if (nblcks_new<nblcks_old) {    // fill gaps
    unsigned i = 0;
    fill_gaps_(bl, A, nblcks_new, i);
  }
}

///////////////////////////////////////////////////////////////////////////////
/* stabilized version for symmetric matrices

   template<class T> static
   void unify_sons_stab_(blcluster* bl, mblock<T>** A, double eps,
   unsigned max_rank)
   {
   unsigned cost_pre, ns1 = bl->getnrs(), ns2 = bl->getncs();
  
   // unify first block row
   blcluster* son = bl->getson(0, 0);
   const unsigned n0 = son->getn2();
   unsigned mi = son->getn1(), m = mi, n = n0;

   mblock<T>* R = new mblock<T>(mi, n0);
   R->copy(*A[son->getidx()]);
   cost_pre = A[son->getidx()]->size();

   unsigned cols = 0;
   for (unsigned j=1; j<ns2; ++j) {
   mblock<T>* R1 = R;
   son = bl->getson(0, j);
   n += son->getn2();
   R = new mblock<T>(mi, n);
   R->unify_cols(eps, max_rank, *R1, *A[son->getidx()],
   cols, X, bl->getn2(), X+cols*bl->getn2(), bl->getn1());
   cost_pre += A[son->getidx()]->size();
   delete R1;
   }

   // now the other rows
   for (unsigned i=1; i<ns1; ++i) {

   son = bl->getson(i, 0);
   mi = son->getn1();
   n = n0;
   mblock<T>* Rz = new mblock<T>(mi, n0);
   Rz->copy(*A[son->getidx()]);
   cost_pre += A[son->getidx()]->size();

   for (unsigned j=1; j<ns2; ++j) {
   mblock<T>* R1 = Rz;
   son = bl->getson(i, j);
   n += son->getn2();
   Rz = new mblock<T>(mi, n);
   Rz->unify_cols(eps, max_rank, *R1, *A[son->getidx()], 
   cols, X, bl->getn2(), 
   X+cols*bl->getn2()+son->getb1()-bl->getb1(), 
   bl->getn1());
   cost_pre += A[son->getidx()]->size();

   delete R1;
   }

   mblock<T>* R1 = R;
   m += mi;
   R = new mblock<T>(m, n);
   R->unify_rows(eps, max_rank, *R1, *Rz,
   cols, X, bl->getn2(), X+cols*bl->getn2(), bl->getn1());
   delete R1;
   delete Rz;
   }
	
   if (cost_pre>=R->size()) {

   // let bl become a leaf; the new block gets the index of the first son
   for (unsigned i=0; i<ns1; ++i) {
   for (unsigned j=0; j<ns2; ++j) {
   unsigned idx = bl->getson(i, j)->getidx();
   if (i!=0 || j!=0) {
   delete A[idx];
   A[idx] = NULL;
   }
   else {
   bl->setidx(idx);
   A[idx]->copy(*R);
   }
   delete bl->getson(i, j);
   }
   }

   bl->setsons(0, 0, NULL);
   }
   delete R;
   }

   template<class T> static
   void agglHSym_stab_(blcluster* bl, mblock<T>** A, double eps,
   unsigned max_rank)
   {
   if (bl->isnleaf()) {
   unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
   bool all_sons = true;           // indicates whether all sons are leaves
   blcluster* son;
   for (unsigned i=0; i<ns1; ++i) {
   for (unsigned j=0; j<ns2; ++j) {
   son = bl->getson(i, j);
   if (son) {
   agglHSym_stab_(son, A, eps, max_rank);
   if (son->isnleaf()) all_sons = false;
   }
   else all_sons = false;
   }
   }

   if (all_sons)
   unify_sons_stab_(bl, A, eps, max_rank, level, X, ldX);
   }
   }


   void agglHSym_stab(blcluster* bl, mblock<double>** A, double eps,
   unsigned max_rank)
   {
   unsigned nblcks_old = bl->nleaves();
   agglHSym_stab_(bl, A, eps, max_rank);
   unsigned nblcks_new = bl->nleaves();

   if (nblcks_new<nblcks_old) {    // fill gaps
   unsigned i = 0;
   fill_gaps_(bl, A, nblcks_new, i);
   }
   }
*/

