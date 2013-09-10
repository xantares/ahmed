/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <fstream>
#include "blcluster.h"

void psHeader(std::ofstream& os, unsigned N)
{
  os << "%!PS-Adobe-2.0-2.0 EPSF-2.0" << std::endl;
  const double offset = 10;
  const double scale = 550.0/N;
  const double si = N*scale+20+offset;

  os << "%%BoundingBox: 0 0 " << si << ' ' << si << std::endl;
  os << 10.0+offset << " dup translate" << std::endl;
  os << scale << " dup scale" << std::endl;

  if (scale<5.0)
    os << 0.02*scale << " setlinewidth" << std::endl;
  else
    os << "0.1 setlinewidth" << std::endl;

  os << "/Helvetica findfont " << 12.0/scale << " scalefont setfont"
  << std::endl;
  os << "0 setgray" << std::endl;
  os << "/bl{0 setgray} def" << std::endl;
  os << "/fb{0 setgray fill} def" << std::endl;
  os << "/fr{1.0 0.0 0.0 setrgbcolor fill} def" << std::endl;
  os << "/flr{1.0 0.4 0.4 setrgbcolor fill} def" << std::endl;
  os << "/fg{0.0 1.0 0.0 setrgbcolor fill} def" << std::endl;
  os << "/flg{0.4 1.0 0.4 setrgbcolor fill} def" << std::endl;
	os << "/fblue{0.0 0.0 1.0 setrgbcolor fill} def" << std::endl;
	//  os << "/fg{0.3 setgray fill} def" << std::endl;
  os << "/cs{closepath stroke} def" << std::endl;
  os << "/m{moveto} def" << std::endl;
  os << "/l{lineto} def" << std::endl;
  os << "/hf{/Helvetica findfont} def" << std::endl;
  os << "/sf{scalefont setfont} def" << std::endl;

  // draw border
  os << "bl" << std::endl;
  os << N << ' ' << 0 << " m" << std::endl;
  os << 0 << ' ' << 0 << " l" << std::endl;
  os << 0 << ' ' << N << " l" << std::endl;
  os << N << ' ' << N << " l" << std::endl;
  os << N << ' ' << 0 << " l" << std::endl;
  os << "cs" << std::endl;
}

template<class T> static
void psoutputGeH_(std::ofstream& os, blcluster* bl, unsigned N, mblock<T>** A)
{
  if (bl->isleaf())
    A[bl->getidx()]->psout(os, N, bl->getb1(), bl->getb2(), false, 0);
  else {
    unsigned ns1 = bl->getnrs(), ns2 = bl->getncs();
    for (unsigned i=0; i<ns1; ++i)
      for (unsigned j=0; j<ns2; ++j)
        psoutputGeH_(os, bl->getson(i, j), N, A);
  }
}

void psoutputGeH(std::ofstream& os, blcluster* bl, unsigned N,
               mblock<double>** A)
{
  psHeader(os, N);
  psoutputGeH_(os, bl, N, A);
  os << "showpage" << std::endl;
}


void psoutputGeH(std::ofstream& os, blcluster* bl, unsigned N,
               mblock<float>** A)
{
  psHeader(os, N);
  psoutputGeH_(os, bl, N, A);
  os << "showpage" << std::endl;
}


void psoutputGeH(std::ofstream& os, blcluster* bl, unsigned N,
               mblock<dcomp>** A)
{
  psHeader(os, N);
  psoutputGeH_(os, bl, N, A);
  os << "showpage" << std::endl;
}


void psoutputGeH(std::ofstream& os, blcluster* bl, unsigned N, mblock<scomp>** A)
{
  psHeader(os, N);
  psoutputGeH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

template<class T> static
void psoutputHeH_(std::ofstream& os, blcluster* bl, unsigned N, mblock<T>** A,
                   unsigned level=0, cluster* rcl=NULL, cluster* ccl=NULL,
                   void (*f)(cluster* cl, const unsigned level, double* X,
                             const double* const basis)=NULL)
{
  if (bl->isleaf()) {
    A[bl->getidx()]->psout(os, N, bl->getb1(), bl->getb2(), false,
                           level, ccl, f);
    A[bl->getidx()]->psout(os, N, bl->getb2(), bl->getb1(), true,
                           level, rcl, f);
  }
  else {
    unsigned nrs = bl->getnrs(), ncs = bl->getncs();
    cluster* srcl = rcl;
    cluster* sccl = ccl;
    if(level>0)level--;
    for (unsigned i=0; i<nrs; ++i)
      for (unsigned j=0; j<ncs; ++j)
        if (bl->getson(i, j)){
	  if(f!=NULL){
	    srcl = rcl->getson(i);
	    sccl = ccl->getson(j);
	  }
          psoutputHeH_(os, bl->getson(i, j), N, A, level, srcl, sccl, f);
	}
  }
}

void psoutputHeH(std::ofstream& os, blcluster* bl, unsigned N,
                  mblock<double>** A, unsigned level=0, cluster* cl=NULL,
                  void (*f)(cluster* cl, const unsigned level, double* X,
                            const double* const basis)=NULL)
{
  psHeader(os, N);
  psoutputHeH_(os, bl, N, A, level, cl, cl, f);
  os << "showpage" << std::endl;
}

void psoutputHeH(std::ofstream& os, blcluster* bl, unsigned N,
                  mblock<float>** A, unsigned level=0, cluster* cl=NULL,
                  void (*f)(cluster* cl, const unsigned level, double* X,
                            const double* const basis)=NULL)
{
  psHeader(os, N);
  psoutputHeH_(os, bl, N, A, level, cl, cl, f);
  os << "showpage" << std::endl;
}

void psoutputHeH(std::ofstream& os, blcluster* bl, unsigned N,
                  mblock<dcomp>** A, unsigned level=0, cluster* cl=NULL,
                  void (*f)(cluster* cl, const unsigned level, double* X,
                            const double* const basis)=NULL)
{
  psHeader(os, N);
  psoutputHeH_(os, bl, N, A, level, cl, cl, f);
  os << "showpage" << std::endl;
}

void psoutputHeH(std::ofstream& os, blcluster* bl, unsigned N,
                  mblock<scomp>** A, unsigned level=0, cluster* cl=NULL,
                  void (*f)(cluster* cl, const unsigned level, double* X,
                            const double* const basis)=NULL)
{
  psHeader(os, N);
  psoutputHeH_(os, bl, N, A, level, cl, cl, f);
  os << "showpage" << std::endl;
}

template<class T> static
void psoutputUtH_(std::ofstream& os, blcluster* bl, unsigned N,
                  mblock<T>** A)
{
  if (bl->isleaf())
    A[bl->getidx()]->psout(os, N, bl->getb1(), bl->getb2(), false, 0);
  else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      psoutputUtH_(os, bl->getson(i, i), N, A);
      for (unsigned j=i+1; j<ns; ++j)
        psoutputGeH_(os, bl->getson(i, j), N, A);
    }
  }
}

void psoutputUtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<double>** A)
{
  psHeader(os, N);
  psoutputUtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

void psoutputUtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<float>** A)
{
  psHeader(os, N);
  psoutputUtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

void psoutputUtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<dcomp>** A)
{
  psHeader(os, N);
  psoutputUtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

void psoutputUtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<scomp>** A)
{
  psHeader(os, N);
  psoutputUtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

template<class T> static
void psoutputLtH_(std::ofstream& os, blcluster* bl, unsigned N,
                  mblock<T>** A)
{
  if (bl->isleaf())
    A[bl->getidx()]->psout(os, N, bl->getb1(), bl->getb2(), false, 0);
  else {
    unsigned ns = bl->getnrs();
    for (unsigned i=0; i<ns; ++i) {
      for (unsigned j=0; j<i; ++j)
        psoutputGeH_(os, bl->getson(i, j), N, A);

      psoutputLtH_(os, bl->getson(i, i), N, A);
    }
  }
}

void psoutputLtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<double>** A)
{
  psHeader(os, N);
  psoutputLtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

void psoutputLtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<float>** A)
{
  psHeader(os, N);
  psoutputLtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

void psoutputLtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<dcomp>** A)
{
  psHeader(os, N);
  psoutputLtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}

void psoutputLtH(std::ofstream& os, blcluster* bl, unsigned N,
                 mblock<scomp>** A)
{
  psHeader(os, N);
  psoutputLtH_(os, bl, N, A);
  os << "showpage" << std::endl;
}
