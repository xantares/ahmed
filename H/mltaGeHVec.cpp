/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "H.h"

// y += d A x (A is an H-matrix)
template<class T> static
bool mltaGeHVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  bool changed = false;
  if (bl->isleaf()) {
    if (A[bl->getidx()]->mltaVec(d, x, y)) changed = true;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* xp = x + son->getb2() - bl->getb2();
        T* yp = y + son->getb1() - bl->getb1();
        if (mltaGeHVec_(d, son, A, xp, yp)) changed = true;
      }
    }
  }
  return changed;
}

// Y += d A X (A is an H-matrix)
template<class T> static
bool mltaGeHGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
                T* X, unsigned ldX, T* Y, unsigned ldY)
{

  bool changed = false;
  if (bl->isleaf()) {
    if (A[bl->getidx()]->mltaGeM(d, p, X, ldX, Y, ldY)) changed = true;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* Xp = X + son->getb2() - bl->getb2();
        T* Yp = Y + son->getb1() - bl->getb1();
        if (mltaGeHGeM_(d, son, A, p, Xp, ldX, Yp, ldY)) changed = true;
      }
    }
  }
  return changed;
}

// y += d A^H x  (A is an H-matrix)
template<class T> static
bool mltaGeHhVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  bool changed = false;
  if (bl->isleaf()) {
    if (A[bl->getidx()]->mltahVec(d, x, y)) changed = true;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* xp = x + son->getb1() - bl->getb1();
        T* yp = y + son->getb2() - bl->getb2();
        if (mltaGeHhVec_(d, son, A, xp, yp)) changed = true;
      }
    }
  }
  return changed;
}

// Y += d A^H X  (A is an H-matrix)
template<class T> static
bool mltaGeHhGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
                 T* X, unsigned ldX, T* Y, unsigned ldY)
{
  bool changed = false;
  if (bl->isleaf()) {
    if (A[bl->getidx()]->mltahGeM(d, p, X, ldX, Y, ldY)) changed = true;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* Xp = X + son->getb1() - bl->getb1();
        T* Yp = Y + son->getb2() - bl->getb2();
        if (mltaGeHhGeM_(d, son, A, p, Xp, ldX, Yp, ldY)) changed = true;
      }
    }
  }
  return changed;
}


// y += d A^T x  (A is an H-matrix)
template<class T> static
bool mltaGeHtVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  bool changed = false;
  if (bl->isleaf()) {
    if (A[bl->getidx()]->mltatVec(d, x, y)) changed = true;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* xp = x + son->getb1() - bl->getb1();
        T* yp = y + son->getb2() - bl->getb2();
        if (mltaGeHtVec_(d, son, A, xp, yp)) changed = true;
      }
    }
  }
  return changed;
}

// Y += d A^T X  (A is an H-matrix)
template<class T> static
bool mltaGeHtGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
                 T* X, unsigned ldX, T* Y, unsigned ldY)
{
  bool changed = false;
  if (bl->isleaf()) {
    if (A[bl->getidx()]->mltatGeM(d, p, X, ldX, Y, ldY)) changed = true;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
        T* Xp = X + son->getb1() - bl->getb1();
        T* Yp = Y + son->getb2() - bl->getb2();
        if (mltaGeHtGeM_(d, son, A, p, Xp, ldX, Yp, ldY)) changed = true;
      }
    }
  }
  return changed;
}


// y += d A^H D x (A H-matrix)
template<class T> static
bool mltaGeHhDiHVec_(T d, mblock<T>** A, blcluster* bl, blcluster* blD, int* piv,
		  T* x, T* y)
{
  bool changed = false;
  if (bl->isleaf()) {
    if (bl->isGeM(A) || bl->rank(A)>0) {
      unsigned n = bl->getn1();
      T* tmp = new T[n];
      mltDiHGeM(1, A, blD, piv, x, tmp, n);
      if (A[bl->getidx()]->mltahVec(d, tmp, y)) changed = true;
      delete [] tmp;
    }
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
	blcluster* son = bl->getson(i, j);
	blcluster* sonD = blD->getson(i, i);
	T* xp = x + son->getb1() - bl->getb1();
	T* yp = y + son->getb2() - bl->getb2();
	if (mltaGeHhDiHVec_(d, A, son, sonD, piv, xp, yp)) changed = true;
      }
    }
  }
  return changed;
}

// Y+= d A^H D X (A H-Matrix)
template<class T> static
bool mltaGeHhDiHGeM_(T d, mblock<T>** A, blcluster* bl, blcluster* blD, int* piv,
		  unsigned p, T* X, unsigned ldX, T* Y, unsigned ldY)
{
  bool changed = false;
  if (bl->isleaf()) {
    unsigned n = bl->getn1();
    T* tmp = new T[2*n*p], *tmp1 = tmp + n*p;
    unsigned idx1=0, idx2=0;
    for (unsigned i=0; i<p; ++i) {
      blas::copy(n, X+idx1, tmp+idx2);
      idx1 += ldX; idx2 += n;
    }
    mltDiHGeM(p, A, blD, piv, tmp, tmp1, n);
    if (A[bl->getidx()]->mltahGeM(d, p, tmp1, n, Y, ldY)) changed = true;
    delete [] tmp;
  }
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i) {
      for (unsigned j=0; j<ns2; ++j) {
        blcluster* son = bl->getson(i, j);
	blcluster* sonD = blD->getson(i, i);
        T* Xp = X + son->getb1() - bl->getb1();
        T* Yp = Y + son->getb2() - bl->getb2();
        if (mltaGeHhDiHGeM_(d, A, son, sonD, piv, p, Xp, ldX, Yp, ldY))
	  changed = true;
      }
    }
  }
  return changed;
}

// y += d A x (A lower triangular H-matrix)
template<class T> static
void mltaLtHVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  if (bl->isleaf()) A[bl->getidx()]->mltaVec(d, x, y);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    blcluster* son;
    T *xp, *yp;
    for (unsigned i=0; i<ns; ++i) {
      for (unsigned j=0; j<i; ++j) {
        son = bl->getson(i, j);
        xp = x + son->getb2() - bl->getb2();
        yp = y + son->getb1() - bl->getb1();
        mltaGeHVec_(d, son, A, xp, yp);
      }
      son = bl->getson(i, i);
      xp = x + son->getb2() - bl->getb2();
      yp = y + son->getb1() - bl->getb1();
      mltaLtHVec_(d, son, A, xp, yp);
    }
  }
}


// Y += d A X (A lower triangular H-matrix)
template<class T> static
void mltaLtHGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
                  T* X, unsigned ldX, T* Y, unsigned ldY)
{
  if (bl->isleaf()) A[bl->getidx()]->mltaGeM(d, p, X, ldX, Y, ldY);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    blcluster* son;
    T *Xp, *Yp;
    for (unsigned i=0; i<ns; ++i) {
      for (unsigned j=0; j<i; ++j) {
        son = bl->getson(i, j);
        Xp = X + son->getb2() - bl->getb2();
        Yp = Y + son->getb1() - bl->getb1();
        mltaGeHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      }
      son = bl->getson(i, i);
      Xp = X + son->getb2() - bl->getb2();
      Yp = Y + son->getb1() - bl->getb1();
      mltaLtHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
    }
  }
}


// y += d A^H x (A lower triangular H-matrix)
template<class T> static
void mltaLtHhVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  if (bl->isleaf()) A[bl->getidx()]->mltahVec(d, x, y);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    blcluster* son;
    T *xp, *yp;
    for (unsigned i=0; i<ns; ++i) {
      for (unsigned j=0; j<i; ++j) {
        son = bl->getson(i, j);
        xp = x + son->getb1() - bl->getb1();
        yp = y + son->getb2() - bl->getb2();
        mltaGeHhVec_(d, son, A, xp, yp);
      }
      son = bl->getson(i, i);
      xp = x + son->getb1() - bl->getb1();
      yp = y + son->getb2() - bl->getb2();
      mltaLtHhVec_(d, son, A, xp, yp);
    }
  }
}


// Y += d A^T X (A lower triangular H-matrix)
template<class T> static
void mltaLtHhGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
		   T* X, unsigned ldX, T* Y, unsigned ldY)
{
  if (bl->isleaf()) A[bl->getidx()]->mltahGeM(d, p, X, ldX, Y, ldY);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    blcluster* son;
    T *Xp, *Yp;
    for (unsigned i=0; i<ns; ++i) {
      for (unsigned j=0; j<i; ++j) {
        son = bl->getson(i, j);
        Xp = X + son->getb1() - bl->getb1();
        Yp = Y + son->getb2() - bl->getb2();
        mltaGeHhGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      }
      blcluster* son = bl->getson(i, i);
      T* Xp = X + son->getb1() - bl->getb1();
      T* Yp = Y + son->getb2() - bl->getb2();
      mltaLtHhGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
    }
  }
}

// y += d A x (A upper triangular H-matrix)
template<class T> static
void mltaUtHVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  if (bl->isleaf()) A[bl->getidx()]->mltaVec(d, x, y);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    blcluster* son;
    T *xp, *yp;
    for (unsigned i=0; i<ns; ++i) {
      son = bl->getson(i, i);
      xp = x + son->getb2() - bl->getb2();
      yp = y + son->getb1() - bl->getb1();
      mltaUtHVec_(d, son, A, xp, yp);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        xp = x + son->getb2() - bl->getb2();
        yp = y + son->getb1() - bl->getb1();
        mltaGeHVec_(d, son, A, xp, yp);
      }
    }
  }
}


// Y += d A X (A symm. H-matrix)
template<class T> static
void mltaUtHGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
                  T* X, unsigned ldX, T* Y, unsigned ldY)
{
  if (bl->isleaf()) A[bl->getidx()]->mltaGeM(d, p, X, ldX, Y, ldY);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* Xp = X + son->getb2() - bl->getb2();
      T* Yp = Y + son->getb1() - bl->getb1();
      mltaUtHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        Xp = X + son->getb2() - bl->getb2();
        Yp = Y + son->getb1() - bl->getb1();
        mltaGeHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      }
    }
  }
}


// y += d A^T x (A upper triangular H-matrix)
template<class T> static
void mltaUtHhVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  if (bl->isleaf()) A[bl->getidx()]->mltahVec(d, x, y);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* xp = x + son->getb1() - bl->getb1();
      T* yp = y + son->getb2() - bl->getb2();
      mltaUtHhVec_(d, son, A, xp, yp);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        xp = x + son->getb1() - bl->getb1();
        yp = y + son->getb2() - bl->getb2();
        mltaGeHhVec_(d, son, A, xp, yp);
      }
    }
  }
}

// Y += d A^T X (A upper triangular H-matrix)
template<class T> static
void mltaUtHhGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
		   T* X, unsigned ldX, T* Y, unsigned ldY)
{
  if (bl->isleaf()) A[bl->getidx()]->mltahGeM(d, p, X, ldX, Y, ldY);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* Xp = X + son->getb1() - bl->getb1();
      T* Yp = Y + son->getb2() - bl->getb2();
      mltaUtHhGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        Xp = X + son->getb1() - bl->getb1();
        Yp = Y + son->getb2() - bl->getb2();
        mltaGeHhGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      }
    }
  }
}

// y += d A x (A herm. H-matrix)
template<class T> static
void mltaHeHVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf()) A[bl->getidx()]->mltaVec(d, x, y);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* xp = x + son->getb2() - bl->getb2();
      T* yp = y + son->getb1() - bl->getb1();
      mltaHeHVec_(d, son, A, xp, yp);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        xp = x + son->getb2() - bl->getb2();
        yp = y + son->getb1() - bl->getb1();
        mltaGeHVec_(d, son, A, xp, yp);
        xp = x + son->getb1() - bl->getb1();
        yp = y + son->getb2() - bl->getb2();
        mltaGeHhVec_(d, son, A, xp, yp);
      }
    }
  }
}

// Y += d A X (A symm. H-matrix)
template<class T> static
void mltaHeHGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
                   T* X, unsigned ldX, T* Y, unsigned ldY)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf()) A[bl->getidx()]->mltaGeM(d, p, X, ldX, Y, ldY);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* Xp = X + son->getb2() - bl->getb2();
      T* Yp = Y + son->getb1() - bl->getb1();
      mltaHeHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        Xp = X + son->getb2() - bl->getb2();
        Yp = Y + son->getb1() - bl->getb1();
        mltaGeHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
        Xp = X + son->getb1() - bl->getb1();
        Yp = Y + son->getb2() - bl->getb2();
        mltaGeHhGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      }
    }
  }
}

// y += d A^T x (A herm. H-matrix)
template<class T> static
void mltaHeHtVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  const unsigned n(bl->getn1());
  blas::conj(n, x);
  blas::conj(n, y);
  mltaHeHVec_(conj(d), bl, A, x, y);
  blas::conj(n, x);
  blas::conj(n, y);
}

// y += d A x (A sym. H-matrix)
template<class T> static
void mltaSyHVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf()) A[bl->getidx()]->mltaVec(d, x, y);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* xp = x + son->getb2() - bl->getb2();
      T* yp = y + son->getb1() - bl->getb1();
      mltaSyHVec_(d, son, A, xp, yp);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        xp = x + son->getb2() - bl->getb2();
        yp = y + son->getb1() - bl->getb1();
        mltaGeHVec_(d, son, A, xp, yp);
        xp = x + son->getb1() - bl->getb1();
        yp = y + son->getb2() - bl->getb2();
        mltaGeHtVec_(d, son, A, xp, yp);
      }
    }
  }
}

// Y += d A X (A symm. H-matrix)
template<class T> static
void mltaSyHGeM_(T d, blcluster* bl, mblock<T>** A, unsigned p,
		 T* X, unsigned ldX, T* Y, unsigned ldY)
{
  assert(bl->getn1()==bl->getn2());

  if (bl->isleaf()) A[bl->getidx()]->mltaGeM(d, p, X, ldX, Y, ldY);
  else {
    assert(bl->getnrs()==bl->getncs());
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      blcluster* son = bl->getson(i, i);
      T* Xp = X + son->getb2() - bl->getb2();
      T* Yp = Y + son->getb1() - bl->getb1();
      mltaSyHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      for (unsigned j=i+1; j<ns; ++j) {
        son = bl->getson(i, j);
        Xp = X + son->getb2() - bl->getb2();
        Yp = Y + son->getb1() - bl->getb1();
        mltaGeHGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
        Xp = X + son->getb1() - bl->getb1();
        Yp = Y + son->getb2() - bl->getb2();
        mltaGeHtGeM_(d, son, A, p, Xp, ldX, Yp, ldY);
      }
    }
  }
}

// y += d A^H x (A sym. H-matrix)
template<class T> static
void mltaSyHhVec_(T d, blcluster* bl, mblock<T>** A, T* x, T* y)
{
  const unsigned n(bl->getn1());
  blas::conj(n, x);
  blas::conj(n, y);
  mltaSyHVec_(conj(d), bl, A, x, y);
  blas::conj(n, x);
  blas::conj(n, y);
}




///////////////////////////////////////////////////////////////////////////////
// Instanzen (solange template 'export' noch nicht funktioniert)

//////////////
// double

bool mltaGeHVec(double d, blcluster* bl, mblock<double>** A, double* x,
               double* y)
{
  return mltaGeHVec_(d, bl, A, x, y);
}

bool mltaGeHGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
               double* X, unsigned ldX, double* Y, unsigned ldY)
{
  return mltaGeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhVec(double d, blcluster* bl, mblock<double>** A, double* x,
                double* y)
{
  return mltaGeHhVec_(d, bl, A, x, y);
}

bool mltaGeHhGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
                double* X, unsigned ldX, double* Y, unsigned ldY)
{
  return mltaGeHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

// y += d A^H D x
bool mltaGeHhDiHVec(double d, mblock<double>** A, blcluster* blA, blcluster* blD,
		 int* piv, double* x, double* y)
{
  return mltaGeHhDiHVec_(d, A, blA, blD, piv, x, y);
}

// Y += d A^H D X
bool mltaGeHhDiHGeM(double d, mblock<double>** A, blcluster* blA, 
		 blcluster* blD, int* piv, unsigned p, double* X, unsigned ldX, 
		 double* Y, unsigned ldY)
{
  return mltaGeHhDiHGeM_(d, A, blA, blD, piv, p, X, ldX, Y, ldY);
}

void mltaLtHVec(double d, blcluster* bl, mblock<double>** A, double* x,
                 double* y)
{
  mltaLtHVec_(d, bl, A, x, y);
}

void mltaLtHGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
                 double* X, unsigned ldX, double* Y, unsigned ldY)
{
  mltaLtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaLtHhVec(double d, blcluster* bl, mblock<double>** A, double* x,
                  double* y)
{
  mltaLtHhVec_(d, bl, A, x, y);
}

void mltaLtHhGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
		  double* X, unsigned ldX, double* Y, unsigned ldY)
{
  mltaLtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHVec(double d, blcluster* bl, mblock<double>** A, double* x,
                 double* y)
{
  mltaUtHVec_(d, bl, A, x, y);
}

void mltaUtHGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
                 double* X, unsigned ldX, double* Y, unsigned ldY)
{
  mltaUtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHhVec(double d, blcluster* bl, mblock<double>** A, double* x,
                  double* y)
{
  mltaUtHhVec_(d, bl, A, x, y);
}

void mltaUtHhGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
                 double* X, unsigned ldX, double* Y, unsigned ldY)
{
  mltaUtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaHeHVec(double d, blcluster* bl, mblock<double>** A, double* x,
                  double* y)
{
  mltaHeHVec_(d, bl, A, x, y);
}

void mltaHeHGeM(double d, blcluster* bl, mblock<double>** A, unsigned p,
                  double* X, unsigned ldX, double* Y, unsigned ldY)
{
  mltaHeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}



//////////////
// float

bool mltaGeHVec(float d, blcluster* bl, mblock<float>** A, float* x, float* y)
{
  return mltaGeHVec_(d, bl, A, x, y);
}

bool mltaGeHGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
               float* X, unsigned ldX, float* Y, unsigned ldY)
{
  return mltaGeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhVec(float d, blcluster* bl, mblock<float>** A, float* x, float* y)
{
  return mltaGeHhVec_(d, bl, A, x, y);
}

bool mltaGeHhGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
                float* X, unsigned ldX, float* Y, unsigned ldY)
{
  return mltaGeHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhDiHVec(float d, mblock<float>** A, blcluster* blA, blcluster* blD,
		 int* piv, float* x, float* y)
{
  return mltaGeHhDiHVec_(d, A, blA, blD, piv, x, y);
}

bool mltaGeHhDiHGeM(float d, mblock<float>** A, blcluster* blA, 
		 blcluster* blD, int* piv, unsigned p, float* X, unsigned ldX, 
		 float* Y, unsigned ldY)
{
  return mltaGeHhDiHGeM_(d, A, blA, blD, piv, p, X, ldX, Y, ldY);
}

void mltaLtHVec(float d, blcluster* bl, mblock<float>** A, float* x, float* y)
{
  mltaLtHVec_(d, bl, A, x, y);
}

void mltaLtHGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
                 float* X, unsigned ldX, float* Y, unsigned ldY)
{
  mltaLtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaLtHhVec(float d, blcluster* bl, mblock<float>** A, float* x,
                  float* y)
{
  mltaLtHhVec_(d, bl, A, x, y);
}

void mltaLtHhGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
		  float* X, unsigned ldX, float* Y, unsigned ldY)
{
  mltaLtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHVec(float d, blcluster* bl, mblock<float>** A, float* x,
                 float* y)
{
  mltaUtHVec_(d, bl, A, x, y);
}

void mltaUtHGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
                 float* X, unsigned ldX, float* Y, unsigned ldY)
{
  mltaUtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHhVec(float d, blcluster* bl, mblock<float>** A, float* x,
                  float* y)
{
  mltaUtHhVec_(d, bl, A, x, y);
}

void mltaUtHhGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
		  float* X, unsigned ldX, float* Y, unsigned ldY)
{
  mltaUtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaHeHVec(float d, blcluster* bl, mblock<float>** A, float* x,
                  float* y)
{
  mltaHeHVec_(d, bl, A, x, y);
}

void mltaHeHGeM(float d, blcluster* bl, mblock<float>** A, unsigned p,
                  float* X, unsigned ldX, float* Y, unsigned ldY)
{
  mltaHeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}



//////////////
// double complex

bool mltaGeHVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
               dcomp* y)
{
  return mltaGeHVec_(d, bl, A, x, y);
}

bool mltaGeHGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
               dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  return mltaGeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                dcomp* y)
{
  return mltaGeHhVec_(d, bl, A, x, y);
}

bool mltaGeHhGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
                dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  return mltaGeHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHtVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
               dcomp* y)
{
  return mltaGeHtVec_(d, bl, A, x, y);
}

bool mltaGeHtGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
                dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  return mltaGeHtGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhDiHVec(dcomp d, mblock<dcomp>** A, blcluster* blA, blcluster* blD,
		    int* piv, dcomp* x, dcomp* y)
{
  return mltaGeHhDiHVec_(d, A, blA, blD, piv, x, y);
}

bool mltaGeHhDiHGeM(dcomp d, mblock<dcomp>** A, blcluster* blA, 
		    blcluster* blD, int* piv, unsigned p, dcomp* X, unsigned ldX, 
		    dcomp* Y, unsigned ldY)
{
  return mltaGeHhDiHGeM_(d, A, blA, blD, piv, p, X, ldX, Y, ldY);
}

void mltaLtHVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                 dcomp* y)
{
  mltaLtHVec_(d, bl, A, x, y);
}

void mltaLtHGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
                 dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  mltaLtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaLtHhVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                  dcomp* y)
{
  mltaLtHhVec_(d, bl, A, x, y);
}

void mltaLtHhGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
		  dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  mltaLtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                 dcomp* y)
{
  mltaUtHVec_(d, bl, A, x, y);
}

void mltaUtHGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
                 dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  mltaUtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHhVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                  dcomp* y)
{
  mltaUtHhVec_(d, bl, A, x, y);
}

void mltaUtHhGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
		  dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  mltaUtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaHeHVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                  dcomp* y)
{
  mltaHeHVec_(d, bl, A, x, y);
}

void mltaHeHGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
                  dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  mltaHeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaHeHtVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                  dcomp* y)
{
  mltaHeHtVec_(d, bl, A, x, y);
}

void mltaSyHVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                  dcomp* y)
{
  mltaSyHVec_(d, bl, A, x, y);
}

void mltaSyHGeM(dcomp d, blcluster* bl, mblock<dcomp>** A, unsigned p,
		dcomp* X, unsigned ldX, dcomp* Y, unsigned ldY)
{
  mltaSyHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaSyHhVec(dcomp d, blcluster* bl, mblock<dcomp>** A, dcomp* x,
                  dcomp* y)
{
  mltaSyHhVec_(d, bl, A, x, y);
}



//////////////
// complex

bool mltaGeHVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
               scomp* y)
{
  return mltaGeHVec_(d, bl, A, x, y);
}

bool mltaGeHGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
               scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  return mltaGeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                scomp* y)
{
  return mltaGeHhVec_(d, bl, A, x, y);
}

bool mltaGeHhGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
		 scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  return mltaGeHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHtVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
               scomp* y)
{
  return mltaGeHtVec_(d, bl, A, x, y);
}

bool mltaGeHtGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
		 scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  return mltaGeHtGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

bool mltaGeHhDiHVec(scomp d, mblock<scomp>** A, blcluster* blA, blcluster* blD,
		 int* piv, scomp* x, scomp* y)
{
  return mltaGeHhDiHVec_(d, A, blA, blD, piv, x, y);
}

bool mltaGeHhDiHGeM(scomp d, mblock<scomp>** A, blcluster* blA, 
		 blcluster* blD, int* piv, unsigned p, scomp* X, unsigned ldX, 
		 scomp* Y, unsigned ldY)
{
  return mltaGeHhDiHGeM_(d, A, blA, blD, piv, p, X, ldX, Y, ldY);
}

void mltaLtHVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                 scomp* y)
{
  mltaLtHVec_(d, bl, A, x, y);
}

void mltaLtHGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
                 scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  mltaLtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaLtHhVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                  scomp* y)
{
  mltaLtHhVec_(d, bl, A, x, y);
}

void mltaLtHhGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
		  scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  mltaLtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                 scomp* y)
{
  mltaUtHVec_(d, bl, A, x, y);
}

void mltaUtHhVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                  scomp* y)
{
  mltaUtHhVec_(d, bl, A, x, y);
}

void mltaUtHGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
                 scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  mltaUtHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaUtHhGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
                 scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  mltaUtHhGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaHeHVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                  scomp* y)
{
  mltaHeHVec_(d, bl, A, x, y);
}

void mltaHeHGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
                  scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  mltaHeHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaHeHtVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                  scomp* y)
{
  mltaHeHtVec_(d, bl, A, x, y);
}

void mltaSyHVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                  scomp* y)
{
  mltaSyHVec_(d, bl, A, x, y);
}

void mltaSyHGeM(scomp d, blcluster* bl, mblock<scomp>** A, unsigned p,
                  scomp* X, unsigned ldX, scomp* Y, unsigned ldY)
{
  mltaSyHGeM_(d, bl, A, p, X, ldX, Y, ldY);
}

void mltaSyHhVec(scomp d, blcluster* bl, mblock<scomp>** A, scomp* x,
                  scomp* y)
{
  mltaSyHhVec_(d, bl, A, x, y);
}
