/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <fstream>
#include <cstdlib>
#include "sparse.h"



// CRS write
template<class T> static
void CS_write_(char* fname, unsigned n, unsigned* iA, unsigned* jA, T* A)
{
  std::ofstream os(fname, std::ios::out | std::ios::binary);
  if (!os) {
    std::cerr << "Error writing '" << fname << "'." << std::endl;
    exit(1);
  }

  CS_write(os, n, iA, jA, A);
  os.close();
}

// CRS read
template<class T> static
void CS_read_(char* fname, unsigned &n, unsigned* &iA, unsigned* &jA, T* &A)
{
  std::ifstream is(fname, std::ios::in | std::ios::binary);
  if (!is) {
    std::cerr << "File '" << fname << "' not found." << std::endl;
    exit(1);
  }

  CS_read(is, n, iA, jA, A);
  is.close();
}

void CS_read(char* fname, unsigned &n, unsigned* &iA, unsigned* &jA, double* &A)
{
  CS_read_(fname, n, iA, jA, A);
}

void CS_read(char* fname, unsigned &n, unsigned* &iA, unsigned* &jA, dcomp* &A)
{
  CS_read_(fname, n, iA, jA, A);
}

void CS_write(char* fname, unsigned n, unsigned* iA, unsigned* jA,
              double* A)
{
  CS_write_(fname, n, iA, jA, A);
}

void CS_write(char* fname, unsigned n, unsigned* iA, unsigned* jA,
              dcomp* A)
{
  CS_write(fname, n, iA, jA, A);
}

