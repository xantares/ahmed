/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef HMATRIX_H
#define HMATRIX_H

#include "matrix.h"

template<class T> struct HeHMatrix : public Matrix<T> {
  mblock<T>** blcks;
  blcluster* blclTree;

  HeHMatrix(unsigned n, blcluster* tree=NULL) : Matrix<T>(n, n),
      blclTree(tree) {
    blcks = NULL;
  }
  ~HeHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(T d, typ* x, typ* y) const {
    mlta_HSym_vec(d, blclTree, blcks, x, y);
  }

  void precond_apply(T* x) const { }
};

//note that blclTree must not be deleted before call of destructor
//symmetric complex H-matrix
template<class T> struct SyHMatrix : public Matrix<comp<T> > {
  mblock<comp<T> >** blcks;
  blcluster* blclTree;

  SyHMatrix(unsigned n, blcluster* tree = NULL) : Matrix<comp<T> >(n, n),
    blcks(NULL), blclTree(tree) { }

  ~SyHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(comp<T> d, comp<T>* x, comp<T>* y) const {
    mltaSyHVec(d, blclTree, blcks, x, y);
  }

  void atmux(comp<T> d, comp<T>* x, comp<T>* y) const {
    mltaSyHhVec(d, blclTree, blcks, x, y);
  }

  void precond_apply(comp<T>* x) const { }
};

template<class T> struct GeHMatrix : public Matrix<T> {
  mblock<T>** blcks;
  blcluster* blclTree;

  GeHMatrix(unsigned n1, unsigned n2, blcluster* tree=NULL) : Matrix<T>(n1, n2),
      blclTree(tree) {
    blcks = NULL;
  }
  ~GeHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(typ d, typ* x, typ* y) const {
    mlta_H_vec(d, blclTree, blcks, x, y);
  }
  void atmux(typ d, typ* x, typ* y) const {
    mlta_HT_vec(d, blclTree, blcks, x, y);
  }

  void precond_apply(typ* x) const { }
};


template<class T> struct UtHMatrix : public Matrix<T> {
  mblock<T>** blcks;
  blcluster* blclTree;

  UtHMatrix(unsigned n, blcluster* tree=NULL) : Matrix<T>(n, n),
      blclTree(tree) {
    blcks = NULL;
  }
  ~UtHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(typ d, typ* x, typ* y) const {
    mlta_LtH_vec(d, blclTree, blcks, x, y);
  }

  void atmux(typ d, typ* x, typ* y) const {
    mlta_LtHT_vec(d, blclTree, blcks, x, y);
  }

  void precond_apply(typ* x) const { }
};

template<class T> struct LtHMatrix : public Matrix<T> {
  mblock<T>** blcks;
  blcluster* blclTree;

  LtHMatrix(unsigned n, blcluster* tree=NULL) : Matrix<T>(n, n),
      blclTree(tree) {
    blcks = NULL;
  }
  ~LtHMatrix() {
    freembls(blclTree, blcks);
  }

  void amux(typ d, typ* x, typ* y) const {
    mlta_LtH_vec(d, blclTree, blcks, x, y);
  }

  void atmux(typ d, typ* x, typ* y) const {
    mlta_LtHT_vec(d, blclTree, blcks, x, y);
  }

  void precond_apply(typ* x) const { }
};



#endif
