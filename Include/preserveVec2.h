/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef PRESERVEVEC2_H
#define PRESERVEVEC2_H

#include "mblock.h"
#include "blcluster.h"
#include "blas.h"

// pres vec
template<class T>
void initVec(const unsigned order, cluster* cl,
       const unsigned level, T* X, const T* const basis)
{
  blas::copy(cl->size(), basis+cl->getnbeg(), X);
}

// usual approach
template<class T>
void initVecHaar(cluster* cl, const unsigned level, T* X,
                 const T* const basis)
{
  unsigned shift(0);
  initVecHaar_(cl, cl, level, shift, X, basis);
}

template<class T>
void initVecHaar_(const cluster* const clf, cluster* cls, unsigned level, 
                  unsigned& shift, T* X, const T* const basis)
{
  const unsigned s(shift*clf->size());
  if(level==0){
    unsigned beg(cls->getnbeg()-clf->getnbeg());
    unsigned end(cls->getnend()-clf->getnbeg());
    unsigned k(0);
    for(;k<beg;k++) X[k+s] = 0.0;
    for(;k<end;k++) X[k+s] = 1.0;
    for(;k<clf->size();k++) X[k+s] = 0.0;
    shift++;
  }else{
    level--;
    for(unsigned i=0;i<cls->getns();i++)
      initVecHaar_(clf, cls->getson(i), level, shift, X, basis);
  }
}

template<class T>
void initVecHigherOrder(const unsigned order, cluster* cl, 
			const unsigned level, T* X, 
			const T* const basis)
{
  unsigned shift(0);
  if(order>1)
    initVecHigherOrder_(order, cl, cl, level, shift, X, basis);
  else
    if(level>0)
      initVecHaar_(cl, cl, level, shift, X, basis);
    else
      blas::load(cl->size(), 1.0, X);
}


template<class T>
void initVecHigherOrder_(const unsigned order, const cluster* const clf, 
			 cluster* cls, unsigned level, 
			 unsigned& shift, T* X, const T* const basis)
{
  if(level==0){
    unsigned beg(cls->getnbeg()-clf->getnbeg());
    unsigned end(cls->getnend()-clf->getnbeg());
    for(unsigned l=0;l<order;l++){
      unsigned k(0);
      for(;k<beg;k++) X[k+(shift+l)*clf->size()] = 0.0;
      for(;k<end;k++) X[k+(shift+l)*clf->size()] = pow(k-beg,l);
      for(;k<clf->size();k++) X[k+(shift+l)*clf->size()] = 0.0;
    }
    shift+=order;
  }else{
    level--;
    for(unsigned i=0;i<cls->getns();i++)
      initVecHigherOrder_(order, clf, cls->getson(i), level, shift, 
			  X, basis);
  }
}

struct contConst;

// container for information used to build the haar basis
template<class T>
class contBasis{
 protected:
  const unsigned order;
  const unsigned ldX;
  T* const X;
  T* const basis;
  void (*f)(const unsigned order, cluster* cl, const unsigned level, 
	    T* X, const T* const basis);
  unsigned level;
  unsigned c_adm;
  
  // container for constants
  contConst& constants;
  
  // blcluster is needed to determine whether admissible or not
  blcluster* bl;

  // cluster was needed if refinement up to leaves of cluster tree
  // is needed
  cluster* rcl;
  cluster* ccl;
  
 public:
 contBasis(const unsigned ord, T* const Xn, const unsigned ldXn, 
	   T* const basisn, int leveln, contConst& constantsn,
	   cluster* rcln, cluster* ccln,
	   void (*fn)(const unsigned order, cluster* cl, 
		      const unsigned level, T* X, 
		      const T* const basis))
   :order(ord), ldX(ldXn), X(Xn), basis(basisn), f(fn), constants(constantsn),
    rcl(rcln), ccl(ccln)
  {
    if(leveln<1) level=0;
    else level=leveln;
  }

 contBasis(const unsigned ord, T* const Xn, const unsigned ldXn, 
	   T* const basisn, int leveln, contConst& constantsn,
	   blcluster* bln, cluster* rcln, cluster* ccln,
	   void (*fn)(const unsigned order, cluster* cl, 
		      const unsigned level, T* X, 
		      const T* const basis))
   :order(ord), ldX(ldXn), X(Xn), basis(basisn), f(fn), constants(constantsn), 
    bl(bln), rcl(rcln), ccl(ccln)
  {
    if(leveln<1) level=0;
    else level=leveln;
  }
  
  virtual ~contBasis(){}
  T* getX() const{return X;}
  unsigned getldX() const{return ldX;}
  unsigned getorder() const{return order;}
  contConst& getconstants() const{return constants;}

  virtual contBasis<T>* son(int i=-1)=0;
  virtual contBasis<T>* son(unsigned i, unsigned j)=0;
  
 
  unsigned getNVec(unsigned l, cluster* cl){
    if(l!=0 and cl->isleaf()){
      std::cout<<"Error getNVec"<<std::endl;
      exit(1);
    }
    unsigned n(0);
    if(l==0)
      n++;
    else{
      l--;
      for(unsigned i=0;i<cl->getns();i++)
	n += getNVec(l, cl->getson(i));
    }
    return n;
  }

  unsigned getcolc(){
    return getNVec(level, ccl);
  }
  
  unsigned getcolr(){
    return getNVec(level, rcl);
  }

  unsigned getcols() {
    unsigned colc = getcolc();
    unsigned colr = getcolr();
    return order*MAX(colc,colr);
  }

  void initX(){
    unsigned col = getcols();
    f(order, ccl, level, X, basis);
    f(order, rcl, level, X+col*ccl->size(), basis);
  }

  unsigned getlevel(){return level;}
  T* getbasis(){return basis;}

  blcluster* getbl(){return bl;}

  cluster* getrcl(){return rcl;}
  cluster* getccl(){return ccl;}

  void (*getf())(const unsigned order, cluster* cl, 
		 const unsigned level, T* X, const T* const basis)
  {return f;}
};
  
  
// container for information used to build basis 
// (all blocks have same number of vectors)
template<class T>
class contVec : public contBasis<T>{  
 public:  
 contVec(T* Xn, unsigned ldXn, T* basisn,
	 int leveln, contConst& constantsn, 
	 blcluster* bln, cluster* rcln, cluster* ccln, 
	 void (*fn)(const unsigned order, cluster* cl, 
		    const unsigned level, T* X, 
		    const T* const basis), const unsigned ord)
   :contBasis<T>(ord+1, Xn, ldXn, basisn, leveln, constantsn, 
		 bln, rcln, ccln, fn){ }
  
 contVec(contBasis<T> *root, unsigned nrclr, unsigned ncclr)
   :contBasis<T>(root->getorder(),root->getX(),root->getldX(),
		 root->getbasis(),root->getlevel(),
		 root->getconstants(),
		 root->getrcl()->getson(nrclr),
		 root->getccl()->getson(ncclr),root->getf())
    {
      // can be needed for example at mltaGeHhGeH_toMbl_
      if(root->getbl()->isleaf())
	contBasis<T>::bl = root->getbl();
      else
	contBasis<T>::bl = root->getbl()->getson(nrclr,ncclr);
    }
  
 contVec(contBasis<T> *root, int i=-1)
   :contBasis<T>(root->getorder(),root->getX(),root->getldX(),
		 root->getbasis(),root->getlevel()+i,
		 root->getconstants(),
		 root->getbl(),root->getrcl(),root->getccl(),
		 root->getf()){ }
  
  contBasis<T>* son(int i=-1){
    return new contVec<T>(this, i);
  }

  virtual ~contVec(){}
 
  contBasis<T>* son(unsigned i, unsigned j){
    return new contVec<T>(this, i, j);
  }
};

// container for information used to build the haar basis
template<class T>
class contHaar : public contBasis<T>{
 public:  
 contHaar(T* Xn, unsigned ldXn, T* basisn,
	  int leveln, blcluster* bln, cluster* rcln, cluster* ccln, 
	  void (*fn)(const unsigned order, cluster* cl, 
		     const unsigned level, T* X, 
		     const T* const basis), unsigned ord)
   :contBasis<T>(ord+1,Xn,ldXn,basisn,leveln,bln,rcln,ccln,fn){ }
  
 contHaar(contBasis<T> *root, unsigned nrclr, unsigned ncclr)
   :contBasis<T>(root->getorder(),root->getX(),root->getldX(),
		 root->getbasis(),root->getlevel()-1,
		 root->getbl()->getson(nrclr,ncclr),
		 root->getrcl()->getson(nrclr),
		 root->getccl()->getson(ncclr),root->getf()){ }
  
 contHaar(contBasis<T> *root, int i=-1)
   :contBasis<T>(root->getorder(),root->getX(),root->getldX(),
		 root->getbasis(),root->getlevel()+i,
		 root->getbl(),root->getrcl(),root->getccl(),
		 root->getf()){ }
  
  virtual ~contHaar(){}
  
  contBasis<T>* son(int i=-1)
    {
      return new contHaar<T>(this, i);
    }
 
  contBasis<T>* son(unsigned i, unsigned j)
    {
      return new contHaar<T>(this, i, j);
    }
};

/*template<class T>
  void initVecHaar(cluster* cl, const unsigned level, T* X,
  const T* const basis)
  {
  initVecHaar_(cl, cl, level, 0, X, basis);
  }

  template<class T>
  void initVecHaar_(cluster* clf, cluster* cls, unsigned level, unsigned shift,
  T* X, const T* const basis)
  {
  if(level!=0 and cls->isleaf()){
  std::cout<<"Error initVecHaar_"<<std::endl;
  exit(1);
  }
  if(level==0){
  unsigned beg(cls->getnbeg()-clf->getnbeg());
  unsigned end(cls->getnend()-clf->getnbeg());
  unsigned k(0);
  for(;k<beg;k++) X[k+shift*clf->size()] = 0.0;
  for(;k<end;k++) X[k+shift*clf->size()] = 1.0;
  for(;k<clf->size();k++) X[k+shift*clf->size()] = 0.0;
  }else{
  for(unsigned i=0;i<cls->getns();i++){
  level--;
  initVecHaar_(clf, cls->getson(i), level, cls->getns()*shift+i, X, basis);
  }
  }
  }*/

/*template<class T>
  void initVecHaar_(const cluster* const clf, cluster* cls, unsigned level, 
  unsigned& shift, T* X, const T* const basis)
  {
  if(level==0){
  for(unsigned k=0;k<clf->size();k++)
  X[k+shift*clf->size()]=pow(k,shift);
  shift++;
  }else{
  level--;
  for(unsigned i=0;i<cls->getns();i++)
  initVecHaar_(clf,cls->getson(i),level,shift,X,basis);
  }
  }*/

/*template<class T>
  void initVecHaar(cluster* cl, const unsigned level, T* X,
  const T* const basis)
  {
  T* temp = X;
  blas::copy(cl->size(),basis+cl->getnbeg(),X);
  X += cl->size();
  if(level>0)
  initVecHaar_(cl, cl, 0, 0, level, X, X, basis);
  }

  template<class T>
  void initVecHaar_(cluster* clf, cluster* cls, unsigned shift, unsigned current,
  const unsigned maxlevel, T* Xorig, T* X, const T* const basis)
  {
  assert(cls->getns()!=0);
  current++;

  blas::setzero(clf->size(),X);

  //blas::load(cls->getnbeg()-clf->getnbeg(),0.0,X);

  //determine half of the cluster by sons
  double min=1.0;
  double sum=0.0;
  unsigned half=0;
  for(unsigned i=0;i<cls->getns();i++){
  sum+=cls->getson(i)->size();
  double tempd=std::abs(sum/(double)cls->size()-0.5);
  if(tempd<min){
  min=tempd;
  half=i;
  }
  }

  T* temp = X+cls->getnbeg()-clf->getnbeg();

  double alpha;
  for(unsigned i=0;i<cls->getns();i++){
  if(i>half)
  alpha=-1.0;
  else
  alpha=1.0;
  //blas::load(cls->getson(i)->size(),alpha,temp);
  blas::axpy(cls->getson(i)->size(),alpha,basis+cls->getson(i)->getnbeg(),
  temp);
  temp += cls->getson(i)->size();
  }

  //blas::load(clf->getnend()-cls->getnend(),0.0,temp);

  if(current<maxlevel){
  for(unsigned i=0;i<cls->getns();i++){
  initVecHaar_(clf, cls->getson(i), i+2*shift, current, maxlevel, Xorig,
  Xorig+(pow2(current)-1+i+2*shift)*clf->size(), basis);
  }
  }
  }

  template<class T>
  void initVecHaar2sons(cluster* cl, const unsigned level, T* X)
  {
  T* temp = X;
  blas::load(cl->size(),1.0,X);
  X += cl->size();
  if(level>0)
  initVecHaar2sons_(cl, cl, 0, 0, level, X, X);
  }

  template<class T>
  void initVecHaar2sons_(cluster* clf, cluster* cls, unsigned shift, unsigned current,
  const unsigned maxlevel, T* Xorig, T* X)
  {
  assert(cls->getns()==2);
  current++;
  cluster* son0=cls->getson(0);
  cluster* son1=cls->getson(1);

  blas::setzero(clf->size(),X);

  blas::load(son0->getnbeg()-clf->getnbeg(),0.0,X);
  T* temp = X+son0->getnbeg()-clf->getnbeg();

  blas::load(son0->size(),1.0,temp);
  temp += son0->size();

  blas::load(son1->size(),-1.0,temp);
  temp += son1->size();

  blas::load(clf->getnend()-son1->getnend(),0.0,temp);

  if(current<maxlevel){
  for(unsigned i=0;i<2;i++){
  initVecHaar_(clf, cls->getson(i), i+2*shift, current, maxlevel, Xorig,
  Xorig+(pow2(current)-1+i+2*shift)*clf->size());
  }
  }
  }

  template<class T>
  void initVecHaar_(cluster* cl, const unsigned level, T* X)
  {
  const unsigned nvec0=pow2(level);
  const unsigned nvec1=pow2(level+1);
  for(unsigned i=0;i<nvec0;i++){
  blas::load(i*getn/nvec0,0.0,X+i*getn);
  blas::load(getn/nvec1,1.0,X+i*getn/nvec0+i*getn);
  blas::load(getn/nvec1,-1.0, X+i*getn/nvec0+getn/nvec1+i*getn);
  blas::load((nvec0-1-i)*getn/nvec0,0.0,X+i*getn/nvec0+getn/nvec0+i*getn);
  }
  }*/

// enforces AX = Y, A is symm (using rounded addition)
// for pos def. can be improved (fewer summands needed)
template<class T>
void enforceVecsymAppr(blcluster* bl, mblock<T>** A, double eps,
                       unsigned rankmax, unsigned cols, T*X, T*Y)
{
  unsigned length = bl->getn1();
  // XHX
  T* XHX = new T[SQR(cols)];
  blas::setzero(SQR(cols),XHX);
  blas::symhm(length, cols, X, 1.0, XHX);
  lapack::geinv_herm(cols, XHX);
  // XXHX1
  T* XXHX1 = new T[length*cols];
  blas::setzero(length*cols, XXHX1);
  blas::gesymma(length, cols, X, XHX, 1.0, XXHX1);

  T* UVX = new T[length*cols];
  blas::setzero(length*cols, UVX);
  mltaHeHGeM(-1.0, bl, A, cols, X, length, UVX, length);

  // Y(XHX)1XH
  addHeHLrM(bl, A, eps, rankmax, cols, Y, length, XXHX1, length);
  // X(XHX)1Y
  addHeHLrM(bl, A, eps, rankmax, cols, XXHX1, length, Y, length);

  // Y_(XHX)1XH
  addHeHLrM(bl, A, eps, rankmax, cols, UVX, length, XXHX1, length);
  // X(XHX)1Y_
  addHeHLrM(bl, A, eps, rankmax, cols, XXHX1, length, UVX, length);

  // X(XHX)1XHY_(XHX)1XH
  blas::gemhm(length, cols, cols, -1.0, X, length, UVX, length, XHX, cols);
  blas::gemm(length, cols, cols, 1.0, XXHX1, length, XHX, cols, UVX, length);
  addHeHLrM(bl, A, eps, rankmax, cols, UVX, length, XXHX1, length);

  blas::gemhm(length, cols, cols, -1.0, X, length, Y, length, XHX, cols);
  blas::gemm(length, cols, cols, 1.0, XXHX1, length, XHX, cols, UVX, length);
  // X(XHX)1YHX(XHX)1XH
  addHeHLrM(bl, A, eps, rankmax, cols, UVX, length, XXHX1, length);

  delete [] UVX;
  delete [] XHX;
  delete [] XXHX1;
}

// enforces AX = Y, A is symm (using exact addition)
// for pos def. can be improved (fewer summands needed)
template<class T>
void enforceVecsym(blcluster* bl, mblock<T>** A, unsigned cols, T*X, T*Y)
{
  unsigned length = bl->getn1();
  // XHX
  T* XHX = new T[SQR(cols)];
  blas::setzero(SQR(cols),XHX);
  blas::symhm(length, cols, X, 1.0, XHX);
  lapack::geinv_herm(cols, XHX);
  // XXHX1
  T* XXHX1 = new T[length*cols];
  blas::setzero(length*cols, XXHX1);
  blas::gesymma(length, cols, X, XHX, 1.0, XXHX1);

  T* UVX = new T[length*cols];
  blas::setzero(length*cols, UVX);
  mltaHeHGeM(-1.0, bl, A, cols, X, length, UVX, length);

  // Y(XHX)1XH
  addHeHLrMExact(bl, A, cols, Y, length, XXHX1, length);
  // X(XHX)1Y
  addHeHLrMExact(bl, A, cols, XXHX1, length, Y, length);

  // Y_(XHX)1XH
  addHeHLrMExact(bl, A, cols, UVX, length, XXHX1, length);
  // X(XHX)1Y_
  addHeHLrMExact(bl, A, cols, XXHX1, length, UVX, length);

  // X(XHX)1XHY_(XHX)1XH
  blas::gemhm(length, cols, cols, -1.0, X, length, UVX, length, XHX, cols);
  blas::gemm(length, cols, cols, 1.0, XXHX1, length, XHX, cols, UVX, length);
  addHeHLrMExact(bl, A, cols, UVX, length, XXHX1, length);

  blas::gemhm(length, cols, cols, -1.0, X, length, Y, length, XHX, cols);
  blas::gemm(length, cols, cols, 1.0, XXHX1, length, XHX, cols, UVX, length);
  // X(XHX)1YHX(XHX)1XH
  addHeHLrMExact(bl, A, cols, UVX, length, XXHX1, length);

  delete [] UVX;
  delete [] XHX;
  delete [] XXHX1;
}


// enforces U^HUX = Y
template<class T>
void enforceUHUVec(blcluster* bl, mblock<T>** A, unsigned cols,
                   T* X, T* Y, unsigned rankmax)
{
  unsigned length = bl->getn1();

  T* update1 = new T[2*length*cols];
  T* update2 = new T[2*length*cols];
  // XHX
  T* XHX = new T[SQR(cols)];
  blas::setzero(SQR(cols),XHX);
  blas::symhm(length, cols, X, 1.0, XHX);
  lapack::geinv_herm(cols, XHX);
  // X(XHX)1YH and X(XHX)1Y_H
  blas::setzero(length*cols, update1);
  blas::gesymma(length, cols, X, XHX, 1.0, update1);
  blas::copy(length*cols,update1,update2+length*cols);

  T* UVX = new T[length*cols];
  blas::setzero(length*cols, UVX);
  blas::setzero(length*cols, update2);
  mltaUtHGeM(-1.0, bl, A, cols, X, length, update2, length);
  for (unsigned k=0; k<cols; k++)
    mltaUtHhVec(1.0, bl, A, update2+k*length, UVX+k*length);

  blas::copy(length*cols, Y, update2);
  blas::axpy(length*cols, 1.0, UVX, update2);
  blas::copy(length*cols, Y, update1+length*cols);
  blas::axpy(length*cols, 1.0, UVX, update1+length*cols);

  // X(XHX)1XHY_(XHX)1XH
  blas::gemhm(length, cols, cols, -1.0, X, length, UVX, length, XHX, cols);
  blas::gemma(length, cols, cols, 1.0, update1, length, XHX, cols,
              update1+length*cols, length);

  // X(XHX)1YHX(XHX)1XH
  blas::gemhm(length, cols, cols, -1.0, X, length, Y, length, XHX, cols);
  blas::gemma(length, cols, cols, 1.0, update1, length, XHX, cols,
              update1+length*cols, length);

  /*{
    double* test = new double[cols*length];
    double* temp = new double[2*cols*cols];
    double* temp2 = new double[cols*length];
    blas::setzero(cols*length,test);
    blas::setzero(cols*length,temp2);
    mltaUtHGeM(1.0, bl, A, cols, X, length, temp2, length);
    for(unsigned k=0;k<cols;k++)
    mltaUtHhVec(1.0, bl, A, temp2+k*length, test+k*length);

    //blas::gemhm(length, 2*cols, cols, 1.0, update2, length, X, length, temp,
    2*cols);
    //blas::gemma(length, 2*cols, cols, 1.0, update1, length, temp, 2*cols, test,
    length);
    double* test2 = new double[length*length];
    blas::gemmh(length,2*cols,length,1.0,update1,length,update2,length,test2,length);
    for(unsigned i=0;i<length;i++){
    for(unsigned j=i;j<length;j++){
    //test[i]+=test2[i+j*length]*X[i];
    if(i!=j){
    test[i]+=test2[i+j*length]*X[i];
    test[j]+=test2[j+i*length]*X[j];
    }else
    test[i]+=test2[i+j*length]*X[i];
    }
    }


    double errorCholwPres(0.0), nrm(0.0);
    for(unsigned i=0;i<cols*length;i++){
    errorCholwPres += SQR(Y[i]-test[i]);
    nrm += SQR(Y[i]);
    }
    std::cout<<"Error enforceUHU: "<<sqrt(errorCholwPres/nrm)<<std::endl;


    delete [] test;
    delete [] test2;
    delete [] temp;
    delete [] temp2;
    }

    //debug start
    double* u1back = new double[2*cols*length];
    blas::copy(2*cols*length,update2,u1back);
    //debug end*/

  UtHhGeM_solve(bl, A, 2*cols, update1, length);
  UtHhGeM_solve(bl, A, 2*cols, update2, length);

  /*{
    double* test1= new double[2*length*cols];
    blas::setzero(2*length*cols, test1);
    for(unsigned k=0;k<2*cols;k++)
    mltaUtHhVec(1.0, bl, A, update2+k*length, test1+k*length);
    double error(0.0), nrm(0.0);
    for(unsigned i=0;i<2*length*cols;i++){
    error+=SQR(test1[i]-u1back[i]);
    nrm+=SQR(u1back[i]);
    }
    std::cout<<"Error solve: "<<sqrt(error/nrm)<<std::endl;
    }
    delete [] u1back;

    //debug start
    double* test1 = new double[SQR(length)];
    blas::gemmh(length, 2*cols, length, 1.0, update1, length, update2, length,
    test1, length);
    for(unsigned i=0;i<length;i++) test1[i+i*length] += 1.0;
    //debug end*/

  T* Z = new T[(length-1)*2*cols];
  T* D = new T[length];
  UTUdplwr(length, 2*cols, update1, update2, D, Z);

  /*{
    double* test2 = new double[SQR(length)];
    double* U = new double[SQR(length)];
    blas::gemmh(length-1, 2*cols, length, 1.0, Z, length-1, update2, length,
    U, length);
    for(unsigned i=0;i<length;i++) U[i+i*length]=D[i];
    for(unsigned i=1;i<length;i++)
    for(unsigned j=0;j<i;j++) U[i+j*length]=0.0;
    blas::gemhm(length,length,length,1.0,U,length,U,length,test2,length);

    double error(0.0), nrm(0.0);
    for(unsigned i=0;i<SQR(length);i++){
    error+=SQR(test1[i]-test2[i]);
    nrm+=SQR(test1[i]);
    }
    std::cout<<"Error UHU: "<<sqrt(error/nrm)<<std::endl;

    delete [] U;
    delete [] test1;
    delete [] test2;
    }*/

  /*//debug start
    double* t1 = new double[length*cols];
    double* t2 = new double[2*cols*cols];
    double* t3 = new double[length*cols];
    blas::setzero(length*cols, t1);
    mltaUtHGeM(1.0, bl, A, cols, X, length, t1, length);
    mltUssGeM_(length, cols, Z, update2, t1)
    //debug end*/

  mltUssUtH(D, 2*cols, Z, update2, bl, A, 1e-14, rankmax);

  /*//debug start
    double* t4 = new double[length*cols];
    blas::setzero(length*cols, t4);
    mltaUtHGeM(1.0, bl, A, cols, X, length, t4, length);
    double error(0.0), nrm(0.0);
    for(unsigned i=0;i<SQR(length);i++){
    error+=SQR(t3[i]-t4[i]);
    nrm+=SQR(t4[i]);
    }
    std::cout<<"Error UHU: "<<sqrt(error/nrm)<<std::endl;

    delete [] t1;
    delete [] t2;
    delete [] t3;
    delete [] t4;
    //debug end*/

  delete [] Z;
  delete [] D;
  delete [] update1;
  delete [] update2;
  delete [] UVX;
  delete [] XHX;
}

#endif
