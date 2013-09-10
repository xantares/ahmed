/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef BASMOD
#define BASMOD

#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

// fuer Maschinen-Genauigkeit
//#include<limits>
// numeric_limits<double>::min() liefert kleinsten positiven Zahl doppelter
//                               Genauigkeit

extern "C" double cputime(const double);
extern "C" double realtime(const double);
extern void progressbar(std::ostream&, std::string, unsigned, unsigned,
                        unsigned, bool);

#define PRINT(X) std::cout << #X << " = " << X << std::endl;

//#define ABS(x) ((x)<0?-(x):(x))
//#define MAKEPOS(x) (x=ABS(x))
//#define SGN(x) ((x)==0?0:((x)>0?1:-1))

//#define SETBIT(x,n) ((x)|=1L<<(n))
//#define CLEARBIT(x,n) ((x)&=~(1L<<(n)))
//#define GETBIT(x,n) ((x)&1L<<(n))

//#define ODD(x) ((x)&1)
//#define EVEN(x) (!((x)&1))

inline double inMB(const double n)
{
  return n/(1024.0*1024.0);
}

template<class T> void swap(T& arg0, T& arg1)
{
  T temp(arg0);
  arg0 = arg1;
  arg1 = temp;
}

template<class T> T max(T x, T y)
{
  return (x > y) ? x : y;
}

template<class T> T min(T x, T y)
{
  return (x < y) ? x : y;
}

inline unsigned pow2(unsigned k)
{
  return (1 << k);
}


template <class T> std::string to_str(const char *style, T n)
{
  std::string tmp;

  if (style == NULL) return tmp;

  char buffer[256];

  snprintf(buffer, sizeof(buffer), style, n);
  tmp = buffer;

  return tmp;
}

#endif   // BASMOD
