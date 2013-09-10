/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef CMPLEX_H
#define CMPLEX_H

// choose standard implementation of complex class
#define STD_CMPLX 

#ifdef STD_CMPLX

#include <iostream>
#include <complex>

typedef std::complex<double> dcomp;
typedef std::complex<float> scomp;

inline double conj(const double& a) { return a; }
inline float conj(const float& a) { return a; }
template<typename T> inline void conj(std::complex<T> *a) { *a = conj(*a); }
template<typename T> inline std::complex<T> nconj(const std::complex<T>& a) {
  return std::complex<T>(-a.real(), a.imag());
}
inline double nconj(const double& a) { return -a; }
inline float nconj(const float& a) { return -a; }

template<class T> inline T abs2(const std::complex<T>& a) { return norm(a); }

template<class T> inline T Re(const std::complex<T>& a) { return a.real(); }
inline double Re(const double& a) { return a; }
inline float Re(const float& a) { return a; }

template<class T> inline T Im(const std::complex<T>& a) { return a.imag(); }

// Ausgabe
template<typename T> inline 
std::ostream& operator<<(std::ostream& os, const std::complex<T>& z)
{
  return os << '(' << z.real() << ", " << z.imag() << ')';
} 


template<typename T> inline
std::complex<T> operator*(const std::complex<T>& lhs, const unsigned& val)
{
  return lhs * static_cast<T>(val);
}
template<typename T> inline 
std::complex<T> operator*(const unsigned& val, const std::complex<T>& rhs)
{
  return static_cast<T>(val) * rhs;
}

#else  //use own implementation


#include <cmath>
#include <iostream>

template<class T> struct comp {
  T re, im;
  comp() {}
  comp(T x, T y) : re(x), im(y) {}
  comp(T x) : re(x), im(0.0) {}
  template<class S> operator comp<S>() const { return comp<S>((S)re, (S)im); }
};

typedef comp<double> dcomp;
typedef comp<float> scomp;

// Conjugate
template<class T> inline comp<T> conj(comp<T> a)
{
  a.im = -a.im;
  return a;
}

inline double conj(double a) { return a; }
inline float  conj(float a) { return a; }
template<class T> inline void conj(comp<T> *a) { a->im = -a->im; }

template<class T> inline comp<T> nconj(comp<T> a) { a.re = -a.re; return a; }
inline double nconj(double a) { return -a; }
inline float nconj(float a) { return -a; }

template<class T> inline T abs2(comp<T> a) { return a.re*a.re + a.im*a.im; }
template<class T> inline T abs(comp<T> a) { return sqrt(abs2(a)); }

template<class T> inline T Re(comp<T> a) { return a.re; }
inline double Re(double a) { return a; }
inline float  Re(float a) { return a; }

template<class T> inline T Im(comp<T> a) { return a.im; }

template<class T> inline comp<T> operator-(comp<T> a)
{
  a.re = -a.re;
  a.im = -a.im;
  return a;
}

// Additionsoperatoren
template<class T> inline comp<T> operator+=(comp<T>& a, comp<T> b)
{
  a.re += b.re;
  a.im += b.im;
  return a;
}

template<class T> inline comp<T> operator+(comp<T> a, comp<T> b)
{
  return (a += b);
}

template<class T> inline comp<T> operator+(comp<T> a, T b)
{
  a.re += b;
  return a;
}

template<class T> inline comp<T> operator+(T b, comp<T> a)
{
  a.re += b;
  return a;
}

template<class T> inline comp<T> operator+=(comp<T>& a, T b)
{
  a.re += b;
  return a;
}


// Subtraktionsoperator
template<class T> inline comp<T> operator-(comp<T> a, comp<T> b)
{
  a.re -= b.re;
  a.im -= b.im;
  return a;
}

template<class T> inline comp<T> operator-=(comp<T>& a, comp<T> b)
{
  a.re -= b.re;
  a.im -= b.im;
  return a;
}

template<class T> inline comp<T> operator-(comp<T> a, T b)
{
  a.re -= b;
  return a;
}

template<class T> inline comp<T> operator-=(comp<T>& a, T b)
{
  a.re -= b;
  return a;
}

template<class T> inline comp<T> operator-(T b, comp<T> a)
{
  comp<T> r;
  r.re = b - a.re;
  r.im = -a.im;
  return r;
}


// Multiplikationoperatoren
template<class T> inline comp<T> operator*(comp<T> a, comp<T> b)
{
  comp<T> r;
  r.re = a.re*b.re - a.im*b.im;
  r.im = a.re*b.im + a.im*b.re;
  return r;
}

template<class T> inline comp<T> operator*=(comp<T>& a, comp<T> b)
{
  T t = a.re*b.re - a.im*b.im;
  a.im *= b.re;
  a.im += a.re*b.im;
  a.re = t;
  return a;
}

template<class T> inline comp<T> operator*=(comp<T>& a, T b)
{
  a.re *= b;
  a.im *= b;
  return a;
}



template<class T> inline comp<T> operator*(comp<T> a, T b)
{
  a.re *= b;
  a.im *= b;
  return a;
}

template<class T> inline comp<T> operator*(T b, comp<T> a)
{
  a.re *= b;
  a.im *= b;
  return a;
}

template<class T> inline comp<T> operator*(int b, comp<T> a)
{
  a.re *= b;
  a.im *= b;
  return a;
}

// Divisionsoperator
template<class T> inline comp<T> operator/(comp<T> a, comp<T> b)
{
  comp<T> r;
  r.re = (a.re*b.re + a.im*b.im)/abs2(b);
  r.im = (a.im*b.re - a.re*b.im)/abs2(b);
  return r;
}

template<class T> inline comp<T> operator/(T a, comp<T> b)
{
  T d = a/abs2(b);
  b.re *=  d;
  b.im *= -d;
  return b;
}

template<class T> inline comp<T> operator/(comp<T> a, T b)
{
  T d = (T) 1.0/b;
  a.re *= d;
  a.im *= d;
  return a;
}


// Ausgabe
template<class T> inline std::ostream& operator<<(std::ostream& os, comp<T> z)
{
  return os << '(' << z.re << ", " << z.im << ')';
}

#endif //STD_CMPLX

// associate the type of the absolute value
template<typename T> struct num_traits { typedef T abs_type; };
template<> struct num_traits<double> { typedef double abs_type; };
template<> struct num_traits<float> { typedef float abs_type; };
template<> struct num_traits<dcomp> { typedef double abs_type; };
template<> struct num_traits<scomp> { typedef float abs_type; };

#endif //CMPLEX_H
