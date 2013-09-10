/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "basmod.h"
#include "H.h"
#include "sllist.h"

/* CCS format: (indices start with 0 !!!)
   A   the non-zero entries
  iA   the row indices of the non-zero entries
  jA   the beginning indices of the columns in A and iA, (jA(N)=NNZ)
*/

/* CRS format: (indices start with 0 !!!)
   A   the non-zero entries
  jA   the column indices of the non-zero entries
  iA   the beginning indices of the rows in A and jA, (iA(N)=NNZ)
*/

static char _CNV_CRS_H[] = "Converting CRS matrix to H-matrix ... ";
static char _CNV_CCS_H[] = "Converting CCS matrix to H-matrix ... ";

template<class T1, class T2> static
void convCCS_toGeM_(T1* A, unsigned* iA, unsigned* jA, unsigned* op_perm,
                   unsigned* po_perm, blcluster* bl, T2* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  blas::setzero(n1*n2, AH);

  for (unsigned j=0; j<n2; ++j) {
    const unsigned oj = po_perm[b2+j];
    for (unsigned k=jA[oj]; k<jA[oj+1]; ++k) {
      unsigned i = op_perm[iA[k]];
      if (i>=b1 && (i-=b1)<n1) AH[i+n1*j] = (T2) A[k];
    }
  }
}

template<class T1, class T2> static
void convCCS_toGeM_(T1* A, unsigned* iA, unsigned* jA, blcluster* bl, T2* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  blas::setzero(n1*n2, AH);

  for (unsigned j=0; j<n2; ++j) {
    const unsigned oj = b2+j;
    for (unsigned k=jA[oj]; k<jA[oj+1]; ++k) {
      unsigned i = iA[k];
      if (i>=b1 && (i-=b1)<n1) AH[i+n1*j] = (T2) A[k];
    }
  }
}

template<class T1, class T2> static
void convCRS_toGeM_(T1* A, unsigned* jA, unsigned* iA, unsigned* op_perm,
                   unsigned* po_perm, blcluster* bl, T2* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  blas::setzero(n1*n2, AH);

  for (unsigned i=0; i<n1; ++i) {
    unsigned oi = po_perm[b1+i];
    for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
      unsigned j = op_perm[jA[k]];
      if (j>=b2 && (j-=b2)<n2) AH[i+n1*j] = (T2) A[k];
    }
  }
}

template<class T1, class T2> static
void convCRS_toGeM_(T1* A, unsigned* jA, unsigned* iA, blcluster* bl, T2* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  blas::setzero(n1*n2, AH);

  for (unsigned i=0; i<n1; ++i) {
    unsigned oi = b1+i;
    for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
      unsigned j = jA[k];
      if (j>=b2 && (j-=b2)<n2) AH[i+n1*j] = (T2) A[k];
    }
  }
}


template<class T1, class T2> static
void convCCS_toHeM_(T1* A, unsigned* iA, unsigned* jA, unsigned* op_perm,
                      unsigned* po_perm, blcluster* bl, T2* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  blas::setzero(n1*(n1+1)/2, AH);

  for (unsigned j=0; j<n2; ++j) {
    const unsigned oj = po_perm[b2+j];
    for (unsigned k=jA[oj]; k<jA[oj+1]; ++k) {
      unsigned i = op_perm[iA[k]];
      if (i>=b1 && (i-=b1)<n1 && j>=i) AH[j*(j+1)/2+i] = (T2) A[k];
    }
  }
}

template<class T1, class T2> static
void convCCS_toHeM_(T1* A, unsigned* iA, unsigned* jA, blcluster* bl, T2* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  blas::setzero(n1*(n1+1)/2, AH);

  for (unsigned j=0; j<n2; ++j) {
    const unsigned oj = b2+j;
    for (unsigned k=jA[oj]; k<jA[oj+1]; ++k) {
      unsigned i = iA[k];
      if (i>=b1 && (i-=b1)<n1 && j>=i) AH[j*(j+1)/2+i] = (T2) A[k];
    }
  }
}


template<class T1, class T2> static
void convCRS_toHeM_(T1* A, unsigned* jA, unsigned* iA, unsigned* op_perm,
                      unsigned* po_perm, blcluster* bl, T2* AH)
{
  assert(bl->getn1()==bl->getn2() && bl->getb1()==bl->getb2());
  unsigned n = bl->getn1(), b = bl->getb1();

  blas::setzero(n*(n+1)/2, AH);

  for (unsigned i=0; i<n; ++i) {
    unsigned oi = po_perm[b+i];
    for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
      unsigned j = op_perm[jA[k]];
      if (j>=b && (j-=b)<n && j>=i) AH[j*(j+1)/2+i] = (T2) A[k];
    }
  }
}

template<class T1, class T2> static
void convCRS_toHeM_(T1* A, unsigned* jA, unsigned* iA, blcluster* bl, T2* AH)
{
  assert(bl->getn1()==bl->getn2() && bl->getb1()==bl->getb2());
  unsigned n = bl->getn1(), b = bl->getb1();

  blas::setzero(n*(n+1)/2, AH);

  for (unsigned i=0; i<n; ++i) {
    unsigned oi = b + i;
    for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
      unsigned j = jA[k];
      if (j>=b && (j-=b)<n && j>=i) AH[j*(j+1)/2+i] = (T2) A[k];
    }
  }
}


template<class T1, class T2> static
void convCCStolwr_(T1* A, unsigned* iA, unsigned* jA, unsigned* op_perm,
                   unsigned* po_perm, double eps, blcluster* bl,
                   mblock<T2>* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  T2 **S = new T2*[n2+1];
  assert(S!=NULL);

  unsigned l = 0;
  S[0] = new T2[n1+n2];
  assert(S[0]!=NULL);
  blas::setzero(n1+n2, S[0]);

  for (unsigned j=0; j<n2; ++j) {
    bool hit = false;
    unsigned oj = po_perm[b2+j];

    for (unsigned k=jA[oj]; k<jA[oj+1]; ++k) {
      unsigned i = op_perm[iA[k]];
      if (i>=b1 && (i-=b1)<n1) {
        *(S[l]+i) = (T2) A[k];
        hit = true;
      }
    }

    if (hit) {
      *(S[l++]+j+n1) = 1.0;
      S[l] = new T2[n1+n2];
      assert(S[l]!=NULL);
      blas::setzero(n1+n2, S[l]);
    }
  }

  if (l>0) {
    T2* X = new T2[l*(n1+n2)], *Y = X + l*n1;
    assert(X!=NULL);
    for (unsigned i=0; i<l; ++i) {
      blas::copy(n1, S[i], X+i*n1);
      blas::copy(n2, S[i]+n1, Y+i*n2);
      delete[] S[i];
    }

    AH->addLrM(l, X, n1, Y, n2, eps, n1);
    delete[] X;
  }

  delete[] S[l];
  delete[] S;
}


template<class T1, class T2> static
void convCCStolwr_(T1* A, unsigned* iA, unsigned* jA, double eps,
                   blcluster* bl, mblock<T2>* AH)
{
  unsigned n1 = bl->getn1(), n2 = bl->getn2(),
                                  b1 = bl->getb1(), b2 = bl->getb2();

  T2 **S = new T2*[n2+1];
  assert(S!=NULL);

  unsigned l = 0;
  S[0] = new T2[n1+n2];
  assert(S[0]!=NULL);
  blas::setzero(n1+n2, S[0]);

  for (unsigned j=0; j<n2; ++j) {
    bool hit = false;
    unsigned oj = b2 + j;

    for (unsigned k=jA[oj]; k<jA[oj+1]; ++k) {
      unsigned i = iA[k];
      if (i>=b1 && (i-=b1)<n1) {
        *(S[l]+i) = (T2) A[k];
        hit = true;
      }
    }

    if (hit) {
      *(S[l++]+j+n1) = 1.0;
      S[l] = new T2[n1+n2];
      assert(S[l]!=NULL);
      blas::setzero(n1+n2, S[l]);
    }
  }

  if (l>0) {
    T2* X = new T2[l*(n1+n2)], *Y = X + l*n1;
    assert(X!=NULL);
    for (unsigned i=0; i<l; ++i) {
      blas::copy(n1, S[i], X+i*n1);
      blas::copy(n2, S[i]+n1, Y+i*n2);
      delete[] S[i];
    }

    AH->addLrM(l, X, n1, Y, n2, eps, n1);
    delete[] X;
  }

  delete[] S[l];
  delete[] S;
}


template<class T>
struct litem_ {
  unsigned i;             // the row, resp. the column
  T d;                    // the value
  litem_(unsigned k, T e) : i(k), d(e) { }
};


template<class T1, class T2>
void convCRStolwr_(T1* A, unsigned* jA, unsigned* iA, unsigned* op_perm,
                   unsigned* po_perm, double eps, blcluster* bl,
                   mblock<T2>* AH)
{
  const unsigned n1 = bl->getn1(), n2 = bl->getn2();
  const unsigned b1 = bl->getb1(), b2 = bl->getb2();


  if (n1<=n2) {
    
    T2 **S = new T2*[n1+1];
    assert(S!=NULL);

    unsigned l = 0;
    S[0] = new T2[n1+n2];
    assert(S[0]!=NULL);
    blas::setzero(n1+n2, S[0]);

    for (unsigned i=0; i<n1; ++i) {

      bool hit = false;
      unsigned oi = po_perm[b1+i];

      for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
	unsigned j = op_perm[jA[k]];
	if (j>=b2 && (j-=b2)<n2) {
	  *(S[l]+j+n1) = (T2) A[k];
	  hit = true;
	}
      }

      if (hit) {
	*(S[l]+i) = 1.0;
	++l;
	S[l] = new T2[n1+n2];
	assert(S[l]!=NULL);
	blas::setzero(n1+n2, S[l]);
      }
    }

    if (l>0) {
      T2* X = new T2[l*(n1+n2)], *Y = X + l*n1;
      assert(X!=NULL);
      for (unsigned i=0; i<l; ++i) {
	blas::copy(n1, S[i], X+i*n1);
	blas::copy(n2, S[i]+n1, Y+i*n2);
	delete[] S[i];
      }

      AH->addLrM(l, X, n1, Y, n2, eps, n1);
      delete[] X;
    }

    delete[] S[l];
    delete[] S;
  }else{// n2<n1 -> generate low-rank matrix column-wise

    sllist<litem_<T2> >* lists = new sllist<litem_<T2> >[n2];
    assert(lists!=NULL);
    unsigned l = 0;

    for (unsigned i=0; i<n1; ++i) {
      const unsigned oi = po_perm[b1+i];
      for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
	const int j = op_perm[jA[k]]-b2;
	if (j>=0 && j<n2) {
	  if (lists[j].empty()) ++l;
	  litem_<T2> item(i, (T2) A[k]);
	  lists[j].push_back(item);
	}
      }
    }

    if (l>0) {
      T2* X = new T2[l*(n1+n2)], *Y = X + l*n1;
      assert(X!=NULL);
      blas::setzero(l*(n1+n2), X);
      unsigned lc = 0;
      for (unsigned j=0; j<n2; ++j) {
	if (!lists[j].empty()) {
	  typename sllist<litem_<T2> >::iterator it = lists[j].begin();
	  while (!it.eol()) {
	    const litem_<T2> item = *it++;
	    X[item.i+lc*n1] = item.d;
	  }
	  Y[j+lc*n2] = 1.0;
	  ++lc;
	}
      }
      
      AH->addLrM(l, X, n1, Y, n2, eps, n1);
      delete[] X;
    }
    delete [] lists;
  }
}


template<class T1, class T2>
void convCRStolwr_(T1* A, unsigned* jA, unsigned* iA, double eps,
                   blcluster* bl, mblock<T2>* AH)
{
  const unsigned n1 = bl->getn1(), n2 = bl->getn2();
  const unsigned b1 = bl->getb1(), b2 = bl->getb2();

  if (n1<=n2) {

    T2 **S = new T2*[n1+1];
    assert(S!=NULL);

    unsigned l = 0;
    S[0] = new T2[n1+n2];
    assert(S[0]!=NULL);
    blas::setzero(n1+n2, S[0]);

    for (unsigned i=0; i<n1; ++i) {

      bool hit = false;
      const unsigned oi = b1 + i;

      for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
	unsigned j = jA[k];
	if (j>=b2 && (j-=b2)<n2) {
	  S[l][j+n1] = (T2) A[k];
	  hit = true;
	}
      }
      
      if (hit) {
	S[l][i] = 1.0;
	++l;
	S[l] = new T2[n1+n2];
	assert(S[l]!=NULL);
	blas::setzero(n1+n2, S[l]);
      }
    }

    if (l>0) {
      //T2* X = new T2[l*(n1+n2)], *Y = X + l*n1;
      if(l>n1+1){
	std::cerr<<"Error"<<std::endl;
	exit(1);
      }
	
      T2* X = new T2[l*n1];
      T2* Y = new T2[l*n2];
      assert(X!=NULL);
      for (unsigned i=0; i<l; ++i) {
	blas::copy(n1, S[i], X+i*n1);
	blas::copy(n2, S[i]+n1, Y+i*n2);
	delete[] S[i];
      }
      
      AH->cpyLrM(l, X, Y);
      //      AH->addLrM(l, X, n1, Y, n2, eps, n1);
      delete[] X;
      delete[] Y;
    }

    delete[] S[l];
    delete[] S;

  } else { // n2<n1 -> generate low-rank matrix column-wise

    sllist<litem_<T2> >* lists = new sllist<litem_<T2> >[n2];
    assert(lists!=NULL);
    unsigned l = 0;

    for (unsigned i=0; i<n1; ++i) {
      const unsigned oi = b1 + i;
      for (unsigned k=iA[oi]; k<iA[oi+1]; ++k) {
	const int j = jA[k] - b2;
	if (j>=0 && j<n2) {
	  if (lists[j].empty()) ++l;
	  litem_<T2> item(i, (T2) A[k]);
	  lists[j].push_back(item);
	}
      }
    }

    if (l>0) {
      T2* X = new T2[l*(n1+n2)], *Y = X + l*n1;
      assert(X!=NULL);
      blas::setzero(l*(n1+n2), X);
      unsigned lc = 0;
      for (unsigned j=0; j<n2; ++j) {
	if (!lists[j].empty()) {
	  typename sllist<litem_<T2> >::iterator it = lists[j].begin();
	  while (!it.eol()) {
	    const litem_<T2> item = *it++;
	    X[item.i+lc*n1] = item.d;
	  }
	  Y[j+lc*n2] = 1.0;
	  ++lc;
	}
      }
      
      AH->cpyLrM(l, X, Y);
      delete[] X;
    }
    delete [] lists;
  }
}


template<class T1, class T2> static
void convCCS_toGeH_(T1* A, unsigned* iA, unsigned* jA, unsigned* op_perm,
                 unsigned* po_perm, double eps, blcluster* bl, mblock<T2>** AH,
                 unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CCS_H, i++, nblcks, 20, true);

    if (bl->isdbl()) {
      AH[idx]->setGeM();
      convCCS_toGeM_(A, iA, jA, op_perm, po_perm, bl, AH[idx]->getdata());
    }
    else if (!bl->issep())
      convCCStolwr_(A, iA, jA, op_perm, po_perm, eps, bl, AH[idx]);

  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k)
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son)
          convCCS_toGeH_(A, iA, jA, op_perm, po_perm, eps, son, AH, i, nblcks);
      }
  }
}

template<class T1, class T2> static
void convCCS_toGeH_(T1* A, unsigned* iA, unsigned* jA,
                 double eps, blcluster* bl, mblock<T2>** AH,
                 unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CCS_H, i++, nblcks, 20, true);

    if (bl->isdbl()) {
      AH[idx]->setGeM();
      convCCS_toGeM_(A, iA, jA, bl, AH[idx]->getdata());
    }
    else if (!bl->issep())
      convCCStolwr_(A, iA, jA, eps, bl, AH[idx]);

  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k)
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son) convCCS_toGeH_(A, iA, jA, eps, son, AH, i, nblcks);
      }
  }
}

template<class T1, class T2> static
void convCRS_toGeH_withoutAlloc_(T1* A, unsigned* jA, unsigned* iA,
                              unsigned* op_perm, unsigned* po_perm, double eps,
                              blcluster* bl, mblock<T2>** AH, unsigned& i,
                              unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    assert(AH[idx]);

    if (bl->isdbl()) {
      AH[idx]->setGeM();
      convCRS_toGeM_(A, jA, iA, op_perm, po_perm, bl, AH[idx]->getdata());
    }
    else if (!bl->issep())
      convCRStolwr_(A, jA, iA, op_perm, po_perm, eps, bl, AH[idx]);

  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k)
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son)
          convCRS_toGeH_withoutAlloc_(A, jA, iA, op_perm, po_perm, eps, son, AH,
                                   i, nblcks);
      }
  }
}

template<class T1, class T2> static
void convCRS_toGeH_withoutAlloc_(T1* A, unsigned* jA, unsigned* iA, double eps,
                              blcluster* bl, mblock<T2>** AH, unsigned& i,
                              unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    assert(AH[idx]);

    if (bl->isdbl()) {
      AH[idx]->setGeM();
      convCRS_toGeM_(A, jA, iA, bl, AH[idx]->getdata());
    }
    else if (!bl->issep()) {
      convCRStolwr_(A, jA, iA, eps, bl, AH[idx]);
    }

  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k)
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son)
          convCRS_toGeH_withoutAlloc_(A, jA, iA, eps, son, AH, i, nblcks);
      }
  }
}


template<class T1, class T2> static
void convCRS_toGeH_(T1* A, unsigned* jA, unsigned* iA, unsigned* op_perm,
                 unsigned* po_perm, double eps, blcluster* bl, mblock<T2>** AH,
                 unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CRS_H, i++, nblcks, 20, true);

    if (bl->isdbl()) {
      AH[idx]->setGeM();
      convCRS_toGeM_(A, jA, iA, op_perm, po_perm, bl, AH[idx]->getdata());
    }
    else if (!bl->issep())
      convCRStolwr_(A, jA, iA, op_perm, po_perm, eps, bl, AH[idx]);

  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k)
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son)
          convCRS_toGeH_(A, jA, iA, op_perm, po_perm, eps, son, AH, i, nblcks);
      }
  }
}

template<class T1, class T2> static
void convCRS_toGeH_(T1* A, unsigned* jA, unsigned* iA,
                 double eps, blcluster* bl, mblock<T2>** AH,
                 unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CRS_H, i++, nblcks, 20, true);

    if (bl->isdbl()) {
      AH[idx]->setGeM();
      convCRS_toGeM_(A, jA, iA, bl, AH[idx]->getdata());
    }
    else if (!bl->issep())
      convCRStolwr_(A, jA, iA, eps, bl, AH[idx]);

  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k)
      for (unsigned l=0; l<ns2; ++l) {
        blcluster* son = bl->getson(k, l);
        if (son) convCRS_toGeH_(A, jA, iA, eps, son, AH, i, nblcks);
      }
  }
}


// wurde die Matrix in der originalen Indizierung generiert, wird eine
// Umsortierung noetig. op_perm bildet die originalen Indizes auf die
// permutierten ab, po_perm ist die inverse Permutation.

void convCCS_toGeH(double* A, unsigned* iA, unsigned* jA,
                unsigned* op_perm, unsigned* po_perm,
                double eps, blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCCS_toGeH_(A, iA, jA, op_perm, po_perm,
              eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCCS_toGeH(double* A, unsigned* iA, unsigned* jA,
                double eps, blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCCS_toGeH_(A, iA, jA, eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


void convCCS_toGeH(unsigned n, double* A, unsigned* iA, unsigned* jA,
                unsigned* op_perm, unsigned* po_perm,
                double eps, blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCCS_toGeH_(A, iA, jA, op_perm, po_perm,
              eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCCS_toGeH(unsigned n, double* A, unsigned* iA, unsigned* jA,
                double eps, blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCCS_toGeH_(A, iA, jA, eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


// wurde die Matrix in der originalen Indizierung generiert, wird eine
// Umsortierung noetig. op_perm bildet die originalen Indizes auf die
// permutierten ab, po_perm ist die inverse Permutation.
void convCRS_toGeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA,
                             unsigned* op_perm, unsigned* po_perm, double eps,
                             blcluster* blclTree, mblock<double>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, op_perm, po_perm,
                           eps, blclTree, AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA,
                             unsigned* op_perm, unsigned* po_perm, double eps,
                             blcluster* blclTree, mblock<float>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, op_perm, po_perm,
                           eps, blclTree, AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(dcomp* A, unsigned* jA, unsigned* iA,
                             unsigned* op_perm, unsigned* po_perm, double eps,
                             blcluster* blclTree, mblock<dcomp>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, op_perm, po_perm,
                           eps, blclTree, AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(dcomp* A, unsigned* jA, unsigned* iA,
                             unsigned* op_perm, unsigned* po_perm, double eps,
                             blcluster* blclTree, mblock<scomp>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, op_perm, po_perm,
                           eps, blclTree, AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA, double eps,
                             blcluster* blclTree, mblock<double>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, eps, blclTree,
                           AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA, double eps,
                             blcluster* blclTree, mblock<float>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, eps, blclTree,
                           AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(dcomp* A, unsigned* jA, unsigned* iA, double eps,
                             blcluster* blclTree, mblock<dcomp>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, eps, blclTree,
                           AH, i, blclTree->nleaves());
}

void convCRS_toGeH_withoutAlloc(dcomp* A, unsigned* jA, unsigned* iA, double eps,
                             blcluster* blclTree, mblock<scomp>**& AH)
{
  unsigned i = 0;
  convCRS_toGeH_withoutAlloc_(A, jA, iA, eps, blclTree,
                           AH, i, blclTree->nleaves());
}

void convCRS_toGeH(double* A, unsigned* jA, unsigned* iA,
                unsigned* op_perm, unsigned* po_perm, double eps,
                blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCRS_toGeH_(A, jA, iA, op_perm, po_perm,
              eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCRS_toGeH(double* A, unsigned* jA, unsigned* iA,
                double eps, blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCRS_toGeH_(A, jA, iA, eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCRS_toGeH(dcomp* A, unsigned* jA, unsigned* iA,
                double eps, blcluster* blclTree, mblock<dcomp>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCRS_toGeH_(A, jA, iA, eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCRS_toGeH(double* A, unsigned* jA, unsigned* iA,
                unsigned* op_perm, unsigned* po_perm, double eps,
                blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCRS_toGeH_(A, jA, iA, op_perm, po_perm,
              eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCRS_toGeH(double* A, unsigned* jA, unsigned* iA,
                double eps, blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCRS_toGeH_(A, jA, iA, eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCRS_toGeH(dcomp* A, unsigned* jA, unsigned* iA,
                double eps, blcluster* blclTree, mblock<scomp>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0;
  convCRS_toGeH_(A, jA, iA, eps, blclTree, AH, i, blclTree->nleaves());

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


template<class T1, class T2> static
void convCCS_toHeH_(T1* A, unsigned* iA, unsigned* jA, unsigned* op_perm,
                    unsigned* po_perm, double eps, blcluster* bl,
                    mblock<T2>** AH, unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CCS_H, i++, nblcks, 20, true);

    AH[idx]->setHeM();
    convCCS_toHeM_(A, iA, jA, op_perm, po_perm, bl, AH[idx]->getdata());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      convCCS_toHeH_(A, iA, jA, op_perm, po_perm,
                     eps, bl->getson(k, k), AH, i, nblcks);
      for (unsigned l=k+1; l<ns2; ++l)
        convCCS_toGeH_(A, iA, jA, op_perm, po_perm, eps,
                    bl->getson(k, l), AH, i, nblcks);
    }
  }
}

template<class T1, class T2> static
void convCCS_toHeH_(T1* A, unsigned* iA, unsigned* jA,
                    double eps, blcluster* bl, mblock<T2>** AH,
                    unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CCS_H, i++, nblcks, 20, true);

    AH[idx]->setHeM();
    convCCS_toHeM_(A, iA, jA, bl, AH[idx]->getdata());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      convCCS_toHeH_(A, iA, jA, eps, bl->getson(k, k), AH, i, nblcks);
      for (unsigned l=k+1; l<ns2; ++l)
        convCCS_toGeH_(A, iA, jA, eps, bl->getson(k, l), AH, i, nblcks);
    }
  }
}

template<class T1, class T2> static
void convCRS_toHeH_withoutAlloc_(T1* A, unsigned* jA, unsigned* iA,
                                 unsigned* op_perm, unsigned* po_perm,
                                 double eps, blcluster* bl,
                                 mblock<T2>** AH, unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    //AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CRS_H, i++, nblcks, 20, true);

    AH[idx]->setHeM();
    convCRS_toHeM_(A, jA, iA, op_perm, po_perm, bl, AH[idx]->getdata());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      convCRS_toHeH_withoutAlloc_(A, jA, iA, op_perm, po_perm, eps,
                                  bl->getson(k, k), AH, i, nblcks);
      for (unsigned l=k+1; l<ns2; ++l)
        convCRS_toGeH_withoutAlloc_(A, jA, iA, op_perm, po_perm, eps,
                                 bl->getson(k, l), AH, i, nblcks);
    }
  }
}

template<class T1, class T2> static
void convCRS_toHeH_withoutAlloc_(T1* A, unsigned* jA, unsigned* iA,
                                 double eps, blcluster* bl,
                                 mblock<T2>** AH, unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    //AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CRS_H, i++, nblcks, 20, true);

    AH[idx]->setHeM();
    convCRS_toHeM_(A, jA, iA, bl, AH[idx]->getdata());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      convCRS_toHeH_withoutAlloc_(A, jA, iA, eps, bl->getson(k, k),
                                  AH, i, nblcks);
      for (unsigned l=k+1; l<ns2; ++l)
        convCRS_toGeH_withoutAlloc_(A, jA, iA, eps, bl->getson(k, l),
                                 AH, i, nblcks);
    }
  }
}

template<class T1, class T2> static
void convCRS_toHeH_(T1* A, unsigned* jA, unsigned* iA, unsigned* op_perm,
                    unsigned* po_perm, double eps, blcluster* bl,
                    mblock<T2>** AH, unsigned& i, unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CRS_H, i++, nblcks, 20, true);

    AH[idx]->setHeM();
    convCRS_toHeM_(A, jA, iA, op_perm, po_perm, bl, AH[idx]->getdata());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      convCRS_toHeH_(A, jA, iA, op_perm, po_perm, eps,
                     bl->getson(k, k), AH, i, nblcks);
      for (unsigned l=k+1; l<ns2; ++l)
        convCRS_toGeH_(A, jA, iA, op_perm, po_perm, eps,
                    bl->getson(k, l), AH, i, nblcks);
    }
  }
}

template<class T1, class T2> static
void convCRS_toHeH_(T1* A, unsigned* jA, unsigned* iA, double eps,
                    blcluster* bl, mblock<T2>** AH, unsigned& i,
                    unsigned nblcks)
{
  if (bl->isleaf()) {
    unsigned idx = bl->getidx();
    AH[idx] = new mblock<T2>(bl->getn1(), bl->getn2());
    assert(AH[idx]);

    progressbar(std::cout, _CNV_CRS_H, i++, nblcks, 20, true);

    AH[idx]->setHeM();
    convCRS_toHeM_(A, jA, iA, bl, AH[idx]->getdata());
  } else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned k=0; k<ns1; ++k) {
      convCRS_toHeH_(A, jA, iA, eps, bl->getson(k, k), AH, i, nblcks);
      for (unsigned l=k+1; l<ns2; ++l)
        convCRS_toGeH_(A, jA, iA, eps, bl->getson(k, l), AH, i, nblcks);
    }
  }
}



// wurde die Matrix in der originalen Indizierung generiert, wird eine
// Umsortierung noetig. op_perm bildet die originalen Indizes auf die
// permutierten ab, po_perm ist die inverse Permutation.

// Ist A symmetrisch, so muss auch der untere Teil in A vorhanden sein
// eventuell st die Routine CS_CRSSymtoCRS zu verwenden

void convCCS_toHeH(double* A, unsigned* iA, unsigned* jA,
                   unsigned* op_perm, unsigned* po_perm, double eps,
                   blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCCS_toHeH_(A, iA, jA, op_perm, po_perm, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCCS_toHeH(double* A, unsigned* iA, unsigned* jA,
                   double eps, blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCCS_toHeH_(A, iA, jA, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


void convCCS_toHeH(double* A, unsigned* iA, unsigned* jA,
                   unsigned* op_perm, unsigned* po_perm, double eps,
                   blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCCS_toHeH_(A, iA, jA, op_perm, po_perm, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}

void convCCS_toHeH(double* A, unsigned* iA, unsigned* jA,
                   double eps, blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCCS_toHeH_(A, iA, jA, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CCS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


// wurde die Matrix in der originalen Indizierung generiert, wird eine
// Umsortierung noetig. op_perm bildet die originalen Indizes auf die
// permutierten ab, po_perm ist die inverse Permutation.

// Ist A symmetrisch, so muss auch der untere Teil in A vorhanden sein
// eventuell st die Routine CS_CRSSymtoCRS zu verwenden

void convCRS_toHeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA,
                                unsigned* op_perm, unsigned* po_perm, double eps,
                                blcluster* blclTree, mblock<double>**& AH)
{
  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_withoutAlloc_(A, jA, iA, op_perm, po_perm, eps, blclTree, AH, i, nblcks);
}

void convCRS_toHeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA,
                                unsigned* op_perm, unsigned* po_perm,
                                double eps, blcluster* blclTree,
                                mblock<float>**& AH)
{
  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_withoutAlloc_(A, jA, iA, op_perm, po_perm, eps, blclTree, AH, i, nblcks);
}

void convCRS_toHeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA,
                                double eps, blcluster* blclTree, mblock<double>**& AH)
{
  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_withoutAlloc_(A, jA, iA, eps, blclTree, AH, i, nblcks);
}

void convCRS_toHeH_withoutAlloc(double* A, unsigned* jA, unsigned* iA,
                                double eps, blcluster* blclTree,
                                mblock<float>**& AH)
{
  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_withoutAlloc_(A, jA, iA, eps, blclTree, AH, i, nblcks);
}

// Ist A symmetrisch, so muss auch der untere Teil in A vorhanden sein
// eventuell st die Routine CS_CRSSymtoCRS zu verwenden

void convCRS_toHeH(double* A, unsigned* jA, unsigned* iA,
                   unsigned* op_perm, unsigned* po_perm, double eps,
                   blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_(A, jA, iA, op_perm, po_perm, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


void convCRS_toHeH(double* A, unsigned* jA, unsigned* iA,
                   double eps, blcluster* blclTree, mblock<double>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_(A, jA, iA, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


void convCRS_toHeH(double* A, unsigned* jA, unsigned* iA,
                   unsigned* op_perm, unsigned* po_perm, double eps,
                   blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_(A, jA, iA, op_perm, po_perm, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}


void convCRS_toHeH(double* A, unsigned* jA, unsigned* iA,
                   double eps, blcluster* blclTree, mblock<float>**& AH)
{
  allocmbls(blclTree, AH);

  unsigned i = 0, nblcks = blclTree->nleaves();
  convCRS_toHeH_(A, jA, iA, eps, blclTree, AH, i, nblcks);

  std::cout << (char) 13 << _CNV_CRS_H << "done -- "
            << inMB(sizeH(blclTree, AH)) << " MB.                 "
            << std::endl;
}
