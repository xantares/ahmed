/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef BLAS_H
#define BLAS_H

#include <assert.h>
#include <stdlib.h>
#include "cmplx.h"    // Complex arthimetic

#define SIGN(a) ((a) >= 0 ? 1.0 : -1.0)
#define SQR(a) ((a)*(a))                             // Quadrat von a
inline double abs2(double a)
{
  return SQR(a);
}
inline double abs2(float a)
{
  return SQR(a);
}
#ifndef WIN32
inline double abs(double a)
{
  return fabs(a);
}
inline float abs(float a)
{
  return fabsf(a);
}
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

inline unsigned LTR(unsigned i, unsigned j, unsigned n)
{
  return ((2*n-j-1)*j)/2+i;
}
inline unsigned UTR(unsigned i, unsigned j)
{
  return (j*(j+1))/2+i;
}
inline unsigned GE(unsigned i, unsigned j, unsigned n)
{
  return i+j*n;
}

const double D_ZERO =  0.0;
const double D_ONE  =  1.0;
const double D_MONE = -1.0;
const float S_ZERO  =  0.0;
const float S_ONE   =  1.0;
const float S_MONE  = -1.0;
const scomp C_ZERO  = 0.0;
const scomp C_ONE   = 1.0;
const scomp C_MONE  = -1.0;
const dcomp Z_ZERO  = 0.0;
const dcomp Z_ONE   = 1.0;
const dcomp Z_MONE  = -1.0;

const double D_PREC = 1e-16;

const unsigned N_ONE = 1;
const int N_MONE = -1;
const char JOB_STR[] = "NTOSVULCRAF";


/*
  #ifdef WIN32
  #define dcopy_ dcopy
  #define dnrm2_ dnrm2
  #define dgemv_ dgemv
  #define dspmv_ dspmv
  #define dtpmv_ dtpmv
  #define dgemm_ dgemm
  #define daxpy_ daxpy
  #define ddot_ ddot
  #define dscal_ dscal
  #define dgesvd_ dgesvd
  #define dgeqrf_ dgeqrf
  #define dormqr_ dormqr
  #define dorgqr_ dorgqr
  #define dgetri_ dgetri
  #define dgetrf_ dgetrf
  #define dsptri_ dsptri
  #define dsptrf_ dsptrf
  #define dtptrs_ dtptrs
  #define dpptrf_ dpptrf
  #define dsytrf_ dsytrf
  #define scopy_ scopy
  #define snrm2_ snrm2
  #define sgemv_ sgemv
  #define sspmv_ sspmv
  #define stpmv_ stpmv
  #define sgemm_ sgemm
  #define saxpy_ saxpy
  #define sdot_ sdot
  #define sscal_ sscal
  #define sgesvd_ sgesvd
  #define sgeqrf_ sgeqrf
  #define sormqr_ sormqr
  #define sorgqr_ sorgqr
  #define sgetri_ sgetri
  #define sgetrf_ sgetrf
  #define ssptri_ ssptri
  #define ssptrf_ ssptrf
  #define stptrs_ stptrs
  #define spptrf_ spptrf
  #define ssytrf_ ssytrf
  #define ccopy_ ccopy
  #define scnrm2_ scnrm2
  #define cgemv_ cgemv
  #define chpmv_ chpmv
  #define ctpmv_ ctpmv
  #define cgemm_ cgemm
  #define caxpy_ caxpy
  #define cdotc_ cdotc
  #define cdotu_ cdotu
  #define cscal_ cscal
  #define cspmv_ cspmv
  #define cgesvd_ cgesvd
  #define cgeqrf_ cgeqrf
  #define cunmqr_ cunmqr
  #define cungqr_ cungqr
  #define cgetri_ cgetri
  #define cgetrf_ cgetrf
  #define chetrf_ chetrf
  #define chptri_ chptri
  #define chptrf_ chptrf
  #define csptri_ csptri
  #define csptrf_ csptrf
  #define ctptrs_ ctptrs
  #define cpptrf_ cpptr
  #define csytrf_ csytrf
  #define zcopy_ zcopy
  #define dznrm2_ dznrm2
  #define zgemv_ zgemv
  #define zhpmv_ zhpmv
  #define ztpmv_ ztpmv
  #define zgemm_ zgemm
  #define zaxpy_ zaxpy
  #define zdotc_ zdotc
  #define zdotu_ zdotu
  #define zscal_ zscal
  #define zgesvd_ zgesvd
  #define zgeqrf_ zgeqrf
  #define zunmqr_ zunmqr
  #define zungqr_ zungqr
  #define zgetri_ zgetri
  #define zgetrf_ zgetrf
  #define zhetrf_ zhetrf
  #define zhptri_ zhptri
  #define zhptrf_ zhptrf
  #define zsptri_ zsptri
  #define zsptrf_ zsptrf
  #define ztptrs_ ztptrs
  #define zpptrf_ zpptrf
  #define zspmv_ zspmv
  #define zsytrf_ zsytrf
  #endif
*/


extern "C" {
  /******************************************************************/
  //double precision real
  /******************************************************************/
  unsigned idamax_(const unsigned*, const double*, const unsigned*);
  void dcopy_(const unsigned*, const double*, const unsigned*,
              double*, const unsigned*);
  void daxpy_(const unsigned*, const double*, const double*,
              const unsigned*, double*, const unsigned*);
  void dscal_(const unsigned*, const double*, const double*,
              const unsigned*);
  double ddot_(const unsigned*, const double*, const unsigned*,
               const double*, const unsigned*);
  double dnrm2_(const unsigned*, const double* const, const unsigned*);

  void dgtsv_(const unsigned*, const unsigned*, const double*,
              const double*, const double*, const double*, const unsigned*,
              const int*);
  void dgemm_(const char*, const char*, const unsigned*, const unsigned*,
              const unsigned*, const double*, const double*, const unsigned*,
              const double*, const unsigned*, const double*, double*,
              const unsigned*);
  void dger_(const unsigned*, const unsigned*, const double*, const double*,
             const unsigned*, double*, const unsigned*, const double*,
             const unsigned*);
  void dgemv_(const char*, const unsigned*, const unsigned*, const double*,
              const double*, const unsigned*, const double*, const unsigned*,
              const double*, double*, const unsigned*);
  void dorgqr_(const unsigned*, const unsigned*, const unsigned*,
               double*, const unsigned*, double*, double*,
               const unsigned*, int*);
  void dormqr_(const char*, const char*, const unsigned*, const unsigned*,
               const unsigned*, double*, const unsigned*, double*,
               double*, const unsigned*, double*, const unsigned*,
               int*);
  void dsyev_(const char*, const char*, const unsigned*, double*,
              const unsigned*, double*, double*, const unsigned*, int*);
  void dgeqrf_(const unsigned*, const unsigned*, double*, const unsigned*,
               double*, double*, int*, int*);
  void dgeqp3_(const unsigned*, const unsigned*, const double*,
               const unsigned*, const unsigned*, const double*, const double*,
               const unsigned*, int*);
  void dgeqp3trunc_(const unsigned*, const unsigned*, const double*,
		    const unsigned*, const unsigned*, const double*, 
		    const unsigned*, const double*, const double*,
		    const double*, const unsigned*, int*);
  void dgeqpf_(const unsigned*, const unsigned*, const double*,
               const unsigned*, const unsigned*, const double*, const double*,
               const unsigned*, int*);
  void dgesvd_(const char*, const char*, const unsigned*, const unsigned*,
               double*, const unsigned*, double*, double*, const unsigned*,
               double*, const unsigned*, double*, const unsigned*, int*);
  void dgetrf_(const unsigned*, const unsigned*, double*, const unsigned*,
               unsigned*, int*);
  void dgetrs_(const char*, const unsigned*, const unsigned*, double*,
               const unsigned*, const unsigned*, double*, const unsigned*,
               int*);
  void dgetri_(const unsigned*, double*, const unsigned*, unsigned*, double*,
               const unsigned*, int*);
  void dspmv_(const char*, const unsigned*, const double*,
              const double*, const double*, const unsigned*,
              const double*, double*, const unsigned*);
  void dsptrf_(const char*, const unsigned*, const double*, int*, int*);
  void dsptri_(const char*, const unsigned*, const double*, int*, double*, int*);
  void dpotrf_(const char*, const unsigned*, const double*, const unsigned*,
               int*);
  void dpotri_(const char*, const unsigned*, const double*, const unsigned*,
               int*);
  void dpptrf_(const char*, const unsigned*, const double*, int*);
  void dpptri_(const char*, const unsigned*, const double*, int*);
  void dtptrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, double*, const unsigned*, double*,
               double*, const unsigned*, double*, const unsigned*,
               int*);
  void dtpsv_(const char*, const char*, const char*, const unsigned*,
              const double*, double*, const unsigned*);
  void dtrtrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, const double*, const unsigned*,
               double*, const unsigned*, int*);
  void dtrsm_(const char*, const char*, const char*, const char*,
              const unsigned*, const unsigned*, const double*, const double*,
              const unsigned*, double*, const unsigned*);
  void dtpmv_(const char*, const char*, const char*, const unsigned*,
              const double*, double*, const unsigned*);
  void dlacpy_(const char*, const unsigned*, const unsigned*, const double*,
               const unsigned*, double*, const unsigned*);
  void dlaset_(const char*, const unsigned*, const unsigned*, const double*,
               const double*, double*, const unsigned*);
  void dtrmm_(const char *, const char *, const char *, const char *,
              unsigned *, unsigned *, const double *, double *, unsigned *,
              double *, unsigned *);
  void dswap_(const unsigned*, double*, const unsigned*, double*,
              const unsigned*);
  void dsytrf_(const char*, const unsigned*, double*, const unsigned*,
               int*, double*, const int*, int*);
  void dlarfg_(const unsigned*, const double*, const double*, const unsigned*,
	       const double*);
  void dlarfb_(const char*, const char*, const char*, const char*,
	       const unsigned*, const unsigned*, const unsigned*,
	       const double*, const unsigned*, const double*,
	       const unsigned*, const double*, const unsigned*,
	       const double*, const unsigned*);

  /******************************************************************/
  //single precision real
  /******************************************************************/
  unsigned isamax_(const unsigned*, const float*, const unsigned*);
  void scopy_(const unsigned*, const float*, const unsigned*,
              float*, const unsigned*);
  void saxpy_(const unsigned*, const float*, const float*,
              const unsigned*, float*, const unsigned*);
  void sscal_(const unsigned*, const float*, const float*,
              const unsigned*);
  float sdot_(const unsigned*, const float*, const unsigned*,
              const float*, const unsigned*);
  float snrm2_(const unsigned*, const float* const, const unsigned*);

  void sgtsv_(const unsigned*, const unsigned*, const float*,
              const float*, const float*, const float*, const unsigned*,
              const int*);
  void sgemm_(const char*, const char*, const unsigned*, const unsigned*,
              const unsigned*, const float*, const float*, const unsigned*,
              const float*, const unsigned*, const float*, float*,
              const unsigned*);
  void sger_(const unsigned*, const unsigned*, const float*, const float*,
             const unsigned*, float*, const unsigned*, const float*,
             const unsigned*);
  void sgemv_(const char*, const unsigned*, const unsigned*, const float*,
              const float*, const unsigned*, const float*, const unsigned*,
              const float*, float*, const unsigned*);
  void sorgqr_(const unsigned*, const unsigned*, const unsigned*,
               float*, const unsigned*, float*, float*,
               const unsigned*, int*);
  void sormqr_(const char*, const char*, const unsigned*, const unsigned*,
               const unsigned*, float*, const unsigned*, float*,
               float*, const unsigned*, float*, const unsigned*,
               int*);
  void ssyev_(const char*, const char*, const unsigned*, float*,
              const unsigned*, float*, float*, const unsigned*, int*);
  void stpsv_(const char*, const char*, const char*, const unsigned*,
              const float*, float*, const unsigned*);
  void sgeqrf_(const unsigned*, const unsigned*, float*, const unsigned*,
               float*, float*, int*, int*);
  void sgesvd_(const char*, const char*, const unsigned*, const unsigned*,
               float*, const unsigned*, float*, float*, const unsigned*,
               float*, const unsigned*, float*, const unsigned*, int*);
  void sgetrf_(const unsigned*, const unsigned*, float*, const unsigned*,
               unsigned*, int*);
  void sgetri_(const unsigned*, float*, const unsigned*, unsigned*, float*,
               const unsigned*, int*);
  void sspmv_(const char*, const unsigned*, const float*,
              const float*, const float*, const unsigned*,
              const float*, float*, const unsigned*);
  void strsm_(const char*, const char*, const char*, const char*,
              const unsigned*, const unsigned*, const float*, const float*,
              const unsigned*, float*, const unsigned*);
  void ssptrf_(const char*, const unsigned*, const float*, int*, int*);
  void ssptri_(const char*, const unsigned*, const float*, int*, float*, int*);
  void spotrf_(const char*, const unsigned*, const float*, const unsigned*,
               int*);
  void spotri_(const char*, const unsigned*, const float*, const unsigned*,
               int*);
  void spptrf_(const char*, const unsigned*, const float*, int*);
  void spptri_(const char*, const unsigned*, const float*, int*);
  void stptrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, float*, float*, const unsigned*, int*);
  void strtrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, const float*, const unsigned*,
               float*, const unsigned*, int*);
  void stpmv_(const char*, const char*, const char*, const unsigned*,
              const float*, float*, const unsigned*);
  void slacpy_(const char*, const unsigned*, const unsigned*, const float*,
               const unsigned*, float*, const unsigned*);
  void sswap_(const unsigned*, float*, const unsigned*, float*,
              const unsigned*);
  void ssytrf_(const char*, const unsigned*, float*, const unsigned*,
               int*, float*, const int*, int*);
  void slarfg_(const unsigned*, const float*, const float*, const unsigned*,
	       const float*);
  void slarfb_(const char*, const char*, const char*, const char*,
	       const unsigned*, const unsigned*, const unsigned*,
	       const float*, const unsigned*, const float*,
	       const unsigned*, const float*, const unsigned*,
	       const float*, const unsigned*);

  /******************************************************************/
  //single precision complex
  /******************************************************************/

  void ccopy_(const unsigned*, const scomp*, const unsigned*,
              scomp*, const unsigned*);
  void caxpy_(const unsigned*, const scomp*, const scomp*,
              const unsigned*, scomp*, const unsigned*);
  void cscal_(const unsigned*, const scomp*, const scomp*,
              const unsigned*);
#ifdef G77CONVENTION
  void cdotc_(scomp*, const unsigned*, const scomp*, const unsigned*,
              const scomp*, const unsigned*);
  void cdotu_(scomp*, const unsigned*, const scomp*, const unsigned*,
              const scomp*, const unsigned*);
#else
  scomp cdotc_(const unsigned*, const scomp*, const unsigned*,
              const scomp*, const unsigned*);
  scomp cdotu_(const unsigned*, const scomp*, const unsigned*,
              const scomp*, const unsigned*);
#endif

  float scnrm2_(const unsigned*, const scomp* const, const unsigned*);

  void cgtsv_(const unsigned*, const unsigned*, const scomp*,
              const scomp*, const scomp*, const scomp*, const unsigned*,
              const int*);
  void cgemm_(const char*, const char*, const unsigned*, const unsigned*,
              const unsigned*, const scomp*, const scomp*, const unsigned*,
              const scomp*, const unsigned*, const scomp*, scomp*,
              const unsigned*);
  void cgerc_(const unsigned*, const unsigned*, const scomp*, const scomp*,
              const unsigned*, scomp*, const unsigned*, const scomp*,
              const unsigned*);
  void cgemv_(const char*, const unsigned*, const unsigned*, const scomp*,
              const scomp*, const unsigned*, const scomp*, const unsigned*,
              const scomp*, scomp*, const unsigned*);
  void cungqr_(const unsigned*, const unsigned*, const unsigned*,
               scomp*, const unsigned*, scomp*, scomp*,
               const unsigned*, int*);
  void cunmqr_(const char*, const char*, const unsigned*, const unsigned*,
               const unsigned*, scomp*, const unsigned*, scomp*,
               scomp*, const unsigned*, scomp*, const unsigned*,
               int*);
  void cheev_(const char*, const char*, const unsigned*, scomp*,
              const unsigned*, float*, scomp*, const unsigned*, float*, int*);
  void ctpsv_(const char*, const char*, const char*, const unsigned*,
              const scomp*, scomp*, const unsigned*);
  void ctrsm_(const char*, const char*, const char*, const char*,
              const unsigned*, const unsigned*, const scomp*, const scomp*,
              const unsigned*, scomp*, const unsigned*);
  void cgeqrf_(const unsigned*, const unsigned*, scomp*, const unsigned*,
               scomp*, scomp*, int*, int*);
  void cgesvd_(const char*, const char*, const unsigned*, const unsigned*,
               scomp*, const unsigned*, float*, scomp*, const unsigned*,
               scomp*, const unsigned*, scomp*, const unsigned*, float*, int*);
  void cgetrf_(const unsigned*, const unsigned*, scomp*, const unsigned*,
               unsigned*, int*);
  void cgetri_(const unsigned*, scomp*, const unsigned*, unsigned*, scomp*,
               const unsigned*, int*);
  void chetrf_(const char*, const unsigned*, scomp*, const unsigned*,
               int*, scomp*, const int*, int*);
  void chpmv_(const char*, const unsigned*, const scomp*,
              const scomp*, const scomp*, const unsigned*,
              const scomp*, scomp*, const unsigned*);
  void chptri_(const char*, const unsigned*, const scomp*, int*, scomp*, int*);
  void chptrf_(const char*, const unsigned*, const scomp*, int*, int*);
  void csptrf_(const char*, const unsigned*, const scomp*, int*, int*);
  void csptri_(const char*, const unsigned*, const scomp*, int*, scomp*, int*);
  void cpotrf_(const char*, const unsigned*, const scomp*, const unsigned*,
               int*);
  void cpotri_(const char*, const unsigned*, const scomp*, const unsigned*,
               int*);
  void cpptrf_(const char*, const unsigned*, const scomp*, int*);
  void cpptri_(const char*, const unsigned*, const scomp*, int*);
  void cspmv_(const char*, const unsigned*, const scomp*,
              const scomp*, const scomp*, const unsigned*,
              const scomp*, scomp*, const unsigned*);
  void cswap_(const unsigned*, scomp*, const unsigned*, scomp*,
              const unsigned*);
  void csytrf_(const char*, const unsigned*, scomp*, const unsigned*,
               int*, scomp*, const int*, int*);
  void ctptrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, scomp*, scomp*, const unsigned*, int*);
  void ctrtrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, const scomp*, const unsigned*,
               scomp*, const unsigned*, int*);
  void ctpmv_(const char*, const char*, const char*, const unsigned*,
              const scomp*, scomp*, const unsigned*);
  void clacpy_(const char*, const unsigned*, const unsigned*, const scomp*,
               const unsigned*, scomp*, const unsigned*);
 
  void clarfg_(const unsigned*, const scomp*, const scomp*, const unsigned*,
	       const scomp*);
  void clarfb_(const char*, const char*, const char*, const char*,
	       const unsigned*, const unsigned*, const unsigned*,
	       const scomp*, const unsigned*, const scomp*,
	       const unsigned*, const scomp*, const unsigned*,
	       const scomp*, const unsigned*);

  /******************************************************************/
  //double precision complex
  /******************************************************************/
  void zcopy_(const unsigned*, const dcomp*, const unsigned*,
              dcomp*, const unsigned*);
  void zaxpy_(const unsigned*, const dcomp*, const dcomp*,
              const unsigned*, dcomp*, const unsigned*);
  void zscal_(const unsigned*, const dcomp*, const dcomp*,
              const unsigned*);
#ifdef G77CONVENTION
  void zdotc_(dcomp*, const unsigned*, const dcomp*, const unsigned*,
              const dcomp*, const unsigned*);
  void zdotu_(dcomp*, const unsigned*, const dcomp*, const unsigned*,
              const dcomp*, const unsigned*);
#else
  dcomp zdotc_(const unsigned*, const dcomp*, const unsigned*,
              const dcomp*, const unsigned*);
  dcomp zdotu_(const unsigned*, const dcomp*, const unsigned*,
              const dcomp*, const unsigned*);
#endif

  double dznrm2_(const unsigned*, const dcomp* const, const unsigned*);

  void zgtsv_(const unsigned*, const unsigned*, const dcomp*,
              const dcomp*, const dcomp*, const dcomp*, const unsigned*,
              const int*);
  void zgemm_(const char*, const char*, const unsigned*, const unsigned*,
              const unsigned*, const dcomp*, const dcomp*, const unsigned*,
              const dcomp*, const unsigned*, const dcomp*, dcomp*,
              const unsigned*);
  void zgerc_(const unsigned*, const unsigned*, const dcomp*, const dcomp*,
              const unsigned*, dcomp*, const unsigned*, const dcomp*,
              const unsigned*);
  void zgemv_(const char*, const unsigned*, const unsigned*, const dcomp*,
              const dcomp*, const unsigned*, const dcomp*, const unsigned*,
              const dcomp*, dcomp*, const unsigned*);
  void zungqr_(const unsigned*, const unsigned*, const unsigned*,
               dcomp*, const unsigned*, dcomp*, dcomp*,
               const unsigned*, int*);
  void zunmqr_(const char*, const char*, const unsigned*, const unsigned*,
               const unsigned*, dcomp*, const unsigned*, dcomp*,
               dcomp*, const unsigned*, dcomp*, const unsigned*,
               int*);
  void zheev_(const char*, const char*, const unsigned*, dcomp*,
              const unsigned*, double*, dcomp*, const unsigned*, double*, int*);
  void ztpsv_(const char*, const char*, const char*, const unsigned*,
              const dcomp*, dcomp*, const unsigned*);
  void zgeqrf_(const unsigned*, const unsigned*, dcomp*, const unsigned*,
               dcomp*, dcomp*, int*, int*);
  void zgesvd_(const char*, const char*, const unsigned*, const unsigned*,
               dcomp*, const unsigned*, double*, dcomp*, const unsigned*,
               dcomp*, const unsigned*, dcomp*, const unsigned*, double*,
               int*);
  void zgetrf_(const unsigned*, const unsigned*, dcomp*, const unsigned*,
               unsigned*, int*);
  void zgetri_(const unsigned*, dcomp*, const unsigned*, unsigned*, dcomp*,
               const unsigned*, int*);
  void zhetrf_(const char*, const unsigned*, dcomp*, const unsigned*,
               int*, dcomp*, const int*, int*);
  void zhpmv_(const char*, const unsigned*, const dcomp*,
              const dcomp*, const dcomp*, const unsigned*,
              const dcomp*, dcomp*, const unsigned*);
  void zhptri_(const char*, const unsigned*, const dcomp*, int*, dcomp*, int*);
  void zhptrf_(const char*, const unsigned*, const dcomp*, int*, int*);
  void ztrsm_(const char*, const char*, const char*, const char*,
              const unsigned*, const unsigned*, const dcomp*, const dcomp*,
              const unsigned*, dcomp*, const unsigned*);
  void zsptrf_(const char*, const unsigned*, const dcomp*, int*, int*);
  void zsptri_(const char*, const unsigned*, const dcomp*, int*, dcomp*, int*);
  void zpotrf_(const char*, const unsigned*, const dcomp*, const unsigned*,
               int*);
  void zpotri_(const char*, const unsigned*, const dcomp*, const unsigned*,
               int*);
  void zpptrf_(const char*, const unsigned*, const dcomp*, int*);
  void zpptri_(const char*, const unsigned*, const dcomp*, int*);
  void zspmv_(const char*, const unsigned*, const dcomp*,
              const dcomp*, const dcomp*, const unsigned*,
              const dcomp*, dcomp*, const unsigned*);
  void zswap_(const unsigned*, dcomp*, const unsigned*, dcomp*,
              const unsigned*);
  void zsytrf_(const char*, const unsigned*, dcomp*, const unsigned*,
               int*, dcomp*, const int*, int*);
  void ztptrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, dcomp*, dcomp*, const unsigned*, int*);
  void ztrtrs_(const char*, const char*, const char*, const unsigned*,
               const unsigned*, const dcomp*, const unsigned*,
               dcomp*, const unsigned*, int*);
  void ztpmv_(const char*, const char*, const char*, const unsigned*,
              const dcomp*, dcomp*, const unsigned*);
  void zlacpy_(const char*, const unsigned*, const unsigned*, const dcomp*,
               const unsigned*, dcomp*, const unsigned*);

  void zlarfg_(const unsigned*, const dcomp*, const dcomp*, const unsigned*,
	       const dcomp*);
  void zlarfb_(const char*, const char*, const char*, const char*,
	       const unsigned*, const unsigned*, const unsigned*,
	       const dcomp*, const unsigned*, const dcomp*,
	       const unsigned*, const dcomp*, const unsigned*,
	       const dcomp*, const unsigned*);
}


namespace blas
{
inline void conj(unsigned n, double* v) {}
inline void conj(unsigned n, float* v) {}
inline void conj(unsigned n, dcomp* v)
{
  while(n--) ::conj(v++);  
}
inline void conj(unsigned n, scomp* v)
{
  while(n--) ::conj(v++); 
}
 inline void conj(unsigned n, double* v, const unsigned& inc) {}
 inline void conj(unsigned n, float* v, const unsigned& inc) {}
 inline void conj(unsigned n, dcomp* v, const unsigned& inc)
{
  while(n--) ::conj(v), v += inc;
}
 inline void conj(unsigned n, scomp* v, const unsigned& inc)
{
  while(n--) ::conj(v), v += inc;
}

inline void swap(const unsigned n, double* x, const unsigned incx,
                 double* y, const unsigned incy )
{
  dswap_(&n, x, &incx, y, &incy);
}

inline void swap(const unsigned n, float* x, const unsigned incx,
                 float* y, const unsigned incy )
{
  sswap_(&n, x, &incx, y, &incy);
}

inline void swap(const unsigned n, scomp* x, const unsigned incx,
                 scomp* y, const unsigned incy )
{
  cswap_(&n, x, &incx, y, &incy);
}

inline void swap(const unsigned n, dcomp* x, const unsigned incx,
                 dcomp* y, const unsigned incy )
{
  zswap_(&n, x, &incx, y, &incy);
}


inline void laset(const unsigned m, const unsigned n, const double a,
                  const double b, double* A, unsigned ldA)
{
  dlaset_(JOB_STR, &m, &n, &a, &b, A, &ldA);
}
inline void lasetu(const unsigned m, const unsigned n, const double a,
                   const double b, double* A, unsigned ldA)
{
  dlaset_(JOB_STR+5, &m, &n, &a, &b, A, &ldA);
}
inline void lasetl(const unsigned m, const unsigned n, const double a,
                   const double b, double* A, unsigned ldA)
{
  dlaset_(JOB_STR+6, &m, &n, &a, &b, A, &ldA);
}

inline unsigned maxi(const unsigned n, double* const v)
{
  return idamax_(&n, v, &N_ONE);
}
inline unsigned maxi(const unsigned n, float* const v)
{
  return isamax_(&n, v, &N_ONE);
}

inline void load(const unsigned n, double e, double* const v)
{
  for (unsigned i=0; i<n; ++i) v[i] = e;
}
inline void load(const unsigned n, float e, float* const v)
{
  for (unsigned i=0; i<n; ++i) v[i] = e;
}
inline void load(const unsigned n, scomp e, scomp* const v)
{
  for (unsigned i=0; i<n; ++i) v[i] = e;
}
inline void load(const unsigned n, dcomp e, dcomp* const v)
{
  for (unsigned i=0; i<n; ++i) v[i] = e;
}

inline void setzero(const unsigned n, unsigned* const v)
{
  for (unsigned i=0; i<n; ++i) v[i] = 0;
}
inline void setzero(const unsigned n, double* const v)
{
  load(n, D_ZERO, v);
}
inline void setzero(const unsigned n, float* const v)
{
  load(n, S_ZERO, v);
}
inline void setzero(const unsigned n, scomp* const v)
{
  load(n, C_ZERO, v);
}
inline void setzero(const unsigned n, dcomp* const v)
{
  load(n, Z_ZERO, v);
}

inline double nrm2(const unsigned n, const double* const v)
{
  return dnrm2_(&n, v, &N_ONE);
}
inline float nrm2(const unsigned n, const float* const v)
{
  return snrm2_(&n, v, &N_ONE);
}
inline float nrm2(const unsigned n, const scomp* const v)
{
  return scnrm2_(&n, v, &N_ONE);
}
inline double nrm2(const unsigned n, const dcomp* const v)
{
  return dznrm2_(&n, v, &N_ONE);
}

inline double diff2(const unsigned n, double* x, double* y)
{
  double s = 0.0;
  for (unsigned i=0; i<n; ++i) s += SQR(x[i]-y[i]);
  return sqrt(s);
}
inline float diff2(const unsigned n, float* x, float* y)
{
  float s = 0.0;
  for (unsigned i=0; i<n; ++i) s += SQR(x[i]-y[i]);
  return sqrtf(s);
}
inline double diff2(const unsigned n, dcomp* x, dcomp* y)
{
  double s = 0.0;
  for (unsigned i=0; i<n; ++i) s += abs2(x[i]-y[i]);
  return sqrt(s);
}
inline float diff2(const unsigned n, scomp* x, scomp* y)
{
  float s = 0.0;
  for (unsigned i=0; i<n; ++i) s += abs2(x[i]-y[i]);
  return sqrtf(s);
}

//computes ||a1 - a2||^2 in euclidean norm
template<typename T>
double dst2(unsigned n, const T *a1, const T *a2)
{
  double nrm2( 0.);
  while(n--)
    nrm2 += abs2(*a1++ - *a2++);
  
  return nrm2;
}

inline void copy(const unsigned n, const double* const orig, double* dest)
{
  dcopy_(&n, orig, &N_ONE, dest, &N_ONE);
}
inline void copy(const unsigned n, const float* const orig, float* dest)
{
  scopy_(&n, orig, &N_ONE, dest, &N_ONE);
}
inline void copy(const unsigned n, const scomp* const orig, scomp* dest)
{
  ccopy_(&n, orig, &N_ONE, dest, &N_ONE);
}
inline void copy(const unsigned n, const dcomp* const orig, dcomp* dest)
{
  zcopy_(&n, orig, &N_ONE, dest, &N_ONE);
}


inline void lacpy(const unsigned m, const unsigned n, double* A,
                  const unsigned ldA, double* B, const unsigned ldB)
{
  dlacpy_(JOB_STR, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpy(const unsigned m, const unsigned n, float* A,
                  const unsigned ldA, float* B, const unsigned ldB)
{
  slacpy_(JOB_STR, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpy(const unsigned m, const unsigned n, scomp* A,
                  const unsigned ldA, scomp* B, const unsigned ldB)
{
  clacpy_(JOB_STR, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpy(const unsigned m, const unsigned n, dcomp* A,
                  const unsigned ldA, dcomp* B, const unsigned ldB)
{
  zlacpy_(JOB_STR, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpyu(const unsigned m, const unsigned n, double* A,
                   const unsigned ldA, double* B, const unsigned ldB)
{
  dlacpy_(JOB_STR+5, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpyu(const unsigned m, const unsigned n, float* A,
                   const unsigned ldA, float* B, const unsigned ldB)
{
  slacpy_(JOB_STR+5, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpyu(const unsigned m, const unsigned n, scomp* A,
                   const unsigned ldA, scomp* B, const unsigned ldB)
{
  clacpy_(JOB_STR+5, &m, &n, A, &ldA, B, &ldB);
}
inline void lacpyu(const unsigned m, const unsigned n, dcomp* A,
                   const unsigned ldA, dcomp* B, const unsigned ldB)
{
  zlacpy_(JOB_STR+5, &m, &n, A, &ldA, B, &ldB);
}


inline void copy(const unsigned n, const double* orig, const unsigned inco,
                 double* dest, const unsigned incd)
{
  dcopy_(&n, orig, &inco, dest, &incd);
}
inline void copy(const unsigned n, const float* orig, const unsigned inco,
                 float* dest, const unsigned incd)
{
  scopy_(&n, orig, &inco, dest, &incd);
}
inline void copy(const unsigned n, const scomp* orig, const unsigned inco,
                 scomp* dest, const unsigned incd)
{
  ccopy_(&n, orig, &inco, dest, &incd);
}
inline void copy(const unsigned n, const dcomp* orig, const unsigned inco,
                 dcomp* dest, const unsigned incd)
{
  zcopy_(&n, orig, &inco, dest, &incd);
}

// Conv.Copy double2float
inline void copy(const unsigned n, double* orig, float* dest)
{
  for (unsigned i=0; i<n; i++) dest[i] = (float) orig[i];
}

// Conv.Copy float2double
inline void copy(const unsigned n, float* orig, double* dest)
{
  for (unsigned i=0; i<n; i++) dest[i] = (double) orig[i];
}

// Conv.Copy dcomp2scomp
inline void copy(const unsigned n, dcomp* orig, scomp* dest)
{
  for (unsigned i=0; i<n; i++)
    dest[i] = orig[i];
  //dest[i].re = (float) orig[i].re, dest[i].im = (float) orig[i].im;
}

// Conv.Copy scomp2dcomp
inline void copy(const unsigned n, scomp* orig, dcomp* dest)
{
  for (unsigned i=0; i<n; i++)
    dest[i] = orig[i];
  //dest[i].re = (double) orig[i].re, dest[i].im = (double) orig[i].im;
}

// Scalar product conj(x)*y
inline double scpr(const unsigned n, const double* const v1,
                   const double* const v2)
{
  return ddot_(&n, v1, &N_ONE, v2, &N_ONE);
}
inline float scpr(const unsigned n, const float* const v1,
                  const float* const v2)
{
  return sdot_(&n, v1, &N_ONE, v2, &N_ONE);
}
inline scomp scpr(unsigned n, const scomp* v1,
                  const scomp* v2)
{
#ifdef G77CONVENTION  
  scomp val;
  cdotc_(&val,&n, v1, &N_ONE, v2, &N_ONE);
  return val;
#else
  return cdotc_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
  /*
  scomp tmp(0., 0.);
  while(n--)
    tmp += conj(*v1++) * *v2++;
  return tmp;*/
}
inline dcomp scpr(unsigned n, const dcomp* v1,
                  const dcomp* v2)
{
#ifdef G77CONVENTION
  dcomp val;
  zdotc_(&val,&n, v1, &N_ONE, v2, &N_ONE);
  return val;
#else
  return zdotc_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
}

// Scalar product conj(x)*y
inline double scpr(const unsigned n, const double* const v1,
                   const unsigned ld1, const double* const v2,
		   const unsigned ld2)
{
  return ddot_(&n, v1, &ld1, v2, &ld2);
}
inline float scpr(const unsigned n, const float* const v1,
                  const unsigned ld1, const float* const v2,
		  const unsigned ld2)
{
  return sdot_(&n, v1, &ld1, v2, &ld2);
}
inline scomp scpr(unsigned n, const scomp* v1,
                  const unsigned ld1, const scomp* v2,
		  const unsigned ld2)
{
#ifdef G77CONVENTION
  scomp val;
  cdotc_(&val,&n, v1, &ld1, v2, &ld2);
  return val;
#else
  return cdotc_(&n, v1, &ld1, v2, &ld2);
#endif  
  /*
  scomp tmp(0., 0.);
  while(n--)
    tmp += conj(*v1) * *v2, v1 += ld1, v2 += ld2;
  return tmp;*/
}
inline dcomp scpr(unsigned n, const dcomp* v1,
                  const unsigned ld1, const dcomp* v2,
		  const unsigned ld2)
{
#ifdef G77CONVENTION
  dcomp val;
  zdotc_(&val,&n, v1, &ld1, v2, &ld2);
  return val;
#else
  return zdotc_(&n, v1, &ld1, v2, &ld2);
#endif
}

//computes v1^T* v2 (without conjugation) 
inline double vecpr(const unsigned n, const double* const v1,
                   const double* const v2)
{
  return ddot_(&n, v1, &N_ONE, v2, &N_ONE); //same as scpr
}
inline float vecpr(const unsigned n, const float* const v1,
                  const float* const v2)
{
  return sdot_(&n, v1, &N_ONE, v2, &N_ONE); //same as scpr
}
inline scomp vecpr(const unsigned n, const scomp* const v1,
                  const scomp* const v2)
{
#ifdef G77CONVENTION
  scomp val;
  cdotu_(&val,&n, v1, &N_ONE, v2, &N_ONE);
  return val;
#else
  return cdotu_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
}
inline dcomp vecpr(const unsigned n, const dcomp* const v1,
                  const dcomp* const v2)
{
#ifdef G77CONVENTION
  dcomp val;
  zdotu_(&val,&n, v1, &N_ONE, v2, &N_ONE);
  return val;
#else
  return zdotu_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
}
// Scalar product conj(x)*y
inline double vecpr(const unsigned n, const double* const v1,
                    const unsigned ld1, const double* const v2,
		    const unsigned ld2)
{
  return ddot_(&n, v1, &ld1, v2, &ld2);
}
inline float vecpr(const unsigned n, const float* const v1,
                   const unsigned ld1, const float* const v2,
		   const unsigned ld2)
{
  return sdot_(&n, v1, &ld1, v2, &ld2);
}
inline scomp vecpr(unsigned n, const scomp* v1,
                   const unsigned ld1, const scomp* v2,
		   const unsigned ld2)
{
#ifdef G77CONVENTION
  scomp val;
  cdotu_(&val,&n, v1, &N_ONE, v2, &N_ONE);
  return val;
#else
  return cdotu_(&n, v1, &N_ONE, v2, &N_ONE);
#endif
  /*scomp tmp(0., 0.);
  while(n--)
    tmp += *v1 * *v2, v1 += ld1, v2 += ld2;
  return tmp;*/
}
inline dcomp vecpr(unsigned n, const dcomp* v1,
                   const unsigned ld1, const dcomp* v2,
		   const unsigned ld2)
{
#ifdef G77CONVENTION
  dcomp val;
  zdotu_(&val,&n, v1, &ld1, v2, &ld2);
  return val;
#else
  return zdotu_(&n, v1, &ld1, v2, &ld2);
#endif
}


inline double sqrsum(const unsigned n, const double* const v)
{
  return ddot_(&n, v, &N_ONE, v, &N_ONE);
}

inline float sqrsum(const unsigned n, const float* const v)
{
  return sdot_(&n, v, &N_ONE, v, &N_ONE);
}

inline double sqrsum(const unsigned n, const dcomp* const v)
{
#ifdef G77CONVENTION
  dcomp z;
  zdotc_(&z,&n, v, &N_ONE, v, &N_ONE);
#else
  dcomp z = zdotc_(&n, v, &N_ONE, v, &N_ONE);
#endif
  return Re(z);
}

inline float sqrsum(const unsigned n, const scomp* const v)
{
#ifdef G77CONVENTION
  scomp z;
  cdotc_(&z,&n, v, &N_ONE, v, &N_ONE);
#else
  scomp z = cdotc_(&n, v, &N_ONE, v, &N_ONE);
#endif
  return Re(z);
}

inline void add(const unsigned n, double* const x, double* const y)
{
  daxpy_(&n, &D_ONE, x, &N_ONE, y, &N_ONE);
}
inline void add(const unsigned n, float* const x, float* const y)
{
  saxpy_(&n, &S_ONE, x, &N_ONE, y, &N_ONE);
}
inline void add(const unsigned n, scomp* const x, scomp* const y)
{
  caxpy_(&n, &C_ONE, x, &N_ONE, y, &N_ONE);
}
inline void add(const unsigned n, dcomp* const x,dcomp* const y)
{
  zaxpy_(&n, &Z_ONE, x, &N_ONE, y, &N_ONE);
}

inline void axpy(const unsigned n, const double d, const double* const x,
                 double* const y)
{
  daxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
}
inline void axpy(const unsigned n, const float d, const float* const x,
                 float* const y)
{
  saxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
}
inline void axpy(const unsigned n, const scomp d, const scomp* const x,
                 scomp* const y)
{
  caxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
}
inline void axpy(const unsigned n, const dcomp d, const dcomp* const x,
                 dcomp* const y)
{
  zaxpy_(&n, &d, x, &N_ONE, y, &N_ONE);
}

inline void scal(const unsigned n, const double d, double* const x, unsigned incx)
{
  dscal_(&n, &d, x, &incx);
}
inline void scal(const unsigned n, const float d, float* const x, unsigned incx)
{
  sscal_(&n, &d, x, &incx);
}
inline void scal(const unsigned n, const scomp d, scomp* const x, unsigned incx)
{
  cscal_(&n, &d, x, &incx);
}
inline void scal(const unsigned n, const dcomp d, dcomp* const x, unsigned incx)
{
  zscal_(&n, &d, x, &incx);
}

template<class T>
inline void scal(const unsigned n, const T d, T* const x)
{
  scal(n, d, x, N_ONE);
}

inline void scal(const unsigned n, const float d, scomp* const x)
{
  scal(n, (scomp) d, x, N_ONE);
}

inline void scal(const unsigned n, const double d, dcomp* const x)
{
  scal(n, (dcomp) d, x, N_ONE);
}



template<class T> inline void normalize(unsigned n, T* x)
{
  T s = 1.0/blas::nrm2(n, x);
  blas::scal(n, s, x);
}

template<class T> inline void mkOrth(unsigned n, const T* v1, T* v2)
{
  T s = -blas::scpr(n, v1, v2);
  blas::axpy(n, s, v1, v2);
}

template<class T>
inline void scal(const unsigned n, const T d, T* const x, T* const y)
{
  for (unsigned i=0; i<n; ++i) y[i] = d * x[i];
}

template<class T>
inline void scal(const unsigned n, T* const x, T* const d)
{
  for (unsigned i=0; i<n; ++i) x[i] *= d[i];
}





// y = d Ax
inline void gemv(const unsigned m, const unsigned n, double d,
                 const double* A, double *x, double *y)
{
  dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE);
}
inline void gemv(const unsigned m, const unsigned n, float d, const float* A,
                 float *x, float *y)
{
  sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE);
}
inline void gemv(const unsigned m, const unsigned n, scomp d, const scomp* A,
                 scomp *x, scomp *y)
{
  cgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE);
}
inline void gemv(const unsigned m, const unsigned n, dcomp d, const dcomp* A,
                 dcomp *x, dcomp *y)
{
  zgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE);
}

inline void gemv(const unsigned m, const unsigned n, double d,
                 const double* A, const unsigned ldA, double *x, double *y)
{
  dgemv_(JOB_STR, &m, &n, &d, A, &ldA, x, &N_ONE, &D_ZERO, y, &N_ONE);
}
inline void gemv(const unsigned m, const unsigned n, float d, const float* A,
                  const unsigned ldA, float *x, float *y)
{
  sgemv_(JOB_STR, &m, &n, &d, A, &ldA, x, &N_ONE, &S_ZERO, y, &N_ONE);
}
inline void gemv(const unsigned m, const unsigned n, scomp d, const scomp* A,
                  const unsigned ldA, scomp *x, scomp *y)
{
  cgemv_(JOB_STR, &m, &n, &d, A, &ldA, x, &N_ONE, &C_ZERO, y, &N_ONE);
}
inline void gemv(const unsigned m, const unsigned n, dcomp d, const dcomp* A,
                  const unsigned ldA, dcomp *x, dcomp *y)
{
  zgemv_(JOB_STR, &m, &n, &d, A, &ldA, x, &N_ONE, &Z_ZERO, y, &N_ONE);
}


// y += d Ax
inline void gemva(const unsigned m, const unsigned n, double d, const double* A,
                  const double *x, double *y)
{
  dgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);
}
inline void gemva(const unsigned m, const unsigned n, float d, const float* A,
                  const float *x, float *y)
{
  sgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);
}
inline void gemva(const unsigned m, const unsigned n, scomp d, const scomp* A,
                  const scomp *x, scomp *y)
{
  cgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);
}
inline void gemva(const unsigned m, const unsigned n, dcomp d, const dcomp* A,
                  const dcomp *x, dcomp *y)
{
  zgemv_(JOB_STR, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);
}

// y = d A^H x
inline void gemhv(const unsigned m, const unsigned n, double d, const double* A,
                  const double *x, double *y)
{
  dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ZERO, y, &N_ONE);
}
inline void gemhv(const unsigned m, const unsigned n, float d, const float* A,
                  const float *x, float *y)
{
  sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ZERO, y, &N_ONE);
}
inline void gemhv(const unsigned m, const unsigned n, scomp d, const scomp* A,
                  const scomp *x, scomp *y)
{
  cgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &C_ZERO, y, &N_ONE);
}
inline void gemhv(const unsigned m, const unsigned n, dcomp d, const dcomp* A,
                  const dcomp *x, dcomp *y)
{
  zgemv_(JOB_STR+7, &m, &n, &d, A, &m, x, &N_ONE, &Z_ZERO, y, &N_ONE);
}

// y += d A^H x
inline void gemhva(const unsigned m, const unsigned n, double d,
                   const double* A, unsigned ldA, const double *x, unsigned incx,
                   double *y, unsigned incy)
{
  dgemv_(JOB_STR+1, &m, &n, &d, A, &ldA, x, &incx, &D_ONE, y, &incy);
}
inline void gemhva(const unsigned m, const unsigned n, float d,
                   const float* A, unsigned ldA, const float *x, unsigned incx,
                   float *y, unsigned incy)
{
  sgemv_(JOB_STR+1, &m, &n, &d, A, &ldA, x, &incx, &S_ONE, y, &incy);
}
inline void gemhva(const unsigned m, const unsigned n, scomp d,
                   const scomp* A, unsigned ldA, const scomp *x, unsigned incx,
                   scomp *y, unsigned incy)
{
  cgemv_(JOB_STR+7, &m, &n, &d, A, &ldA, x, &incx, &C_ONE, y, &incy);
}
inline void gemhva(const unsigned m, const unsigned n, dcomp d,
                   const dcomp* A, unsigned ldA, const dcomp *x, unsigned incx,
                   dcomp *y, unsigned incy)
{
  zgemv_(JOB_STR+7, &m, &n, &d, A, &ldA, x, &incx, &Z_ONE, y, &incy);
}

inline void gemhva(const unsigned m, const unsigned n, double d, const double* A,
                   const double *x, double *y)
{
  gemhva(m, n, d, A, m, x, N_ONE, y, N_ONE);
}
inline void gemhva(const unsigned m, const unsigned n, float d, const float* A,
                   const float *x, float *y)
{
  gemhva(m, n, d, A, m, x, N_ONE, y, N_ONE);
}
inline void gemhva(const unsigned m, const unsigned n, scomp d, const scomp* A,
                   const scomp *x, scomp *y)
{
  gemhva(m, n, d, A, m, x, N_ONE, y, N_ONE);
}
inline void gemhva(const unsigned m, const unsigned n, dcomp d, const dcomp* A,
                   const dcomp *x, dcomp *y)
{
  gemhva(m, n, d, A, m, x, N_ONE, y, N_ONE);
}


// y += d A^T x
inline void gemtva(const unsigned m, const unsigned n, const double d, 
		   const double* A, const double *x, double *y)
{
  dgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &D_ONE, y, &N_ONE);
}
inline void gemtva(const unsigned m, const unsigned n, const float d, 
		   const float* A, const float *x, float *y)
{
  sgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &S_ONE, y, &N_ONE);
}
inline void gemtva(const unsigned m, const unsigned n, const scomp d, 
		   const scomp* A, const scomp *x, scomp *y)
{
  cgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &C_ONE, y, &N_ONE);
}
inline void gemtva(const unsigned m, const unsigned n, const dcomp d, 
		   const dcomp* A, const dcomp *x, dcomp *y)
{
  zgemv_(JOB_STR+1, &m, &n, &d, A, &m, x, &N_ONE, &Z_ONE, y, &N_ONE);
}



// y += d A x (A symm. dense in packed format)
inline void symva(const unsigned n, double d, const double* A,
                   const double *x, double *y)
{
  dspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &D_ONE, y, &N_ONE);
}
inline void symva(const unsigned n, float d, const float* A,
                   const float *x, float *y)
{
  sspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &S_ONE, y, &N_ONE);
}
inline void symva(const unsigned n, scomp d, const scomp* A,
                   const scomp *x, scomp *y)
{
  cspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &C_ONE, y, &N_ONE);
}
inline void symva(const unsigned n, dcomp d, const dcomp* A,
                   const dcomp *x, dcomp *y)
{
  zspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &Z_ONE, y, &N_ONE);
}


// y += d A x (A herm. dense in packed format)
inline void hemva(const unsigned n, double d, const double* A,
                   const double *x, double *y)
{
  dspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &D_ONE, y, &N_ONE);
}
inline void hemva(const unsigned n, float d, const float* A,
                   const float *x, float *y)
{
  sspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &S_ONE, y, &N_ONE);
}
inline void hemva(const unsigned n, scomp d, const scomp* A,
                   const scomp *x, scomp *y)
{
  chpmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &C_ONE, y, &N_ONE);
}
inline void hemva(const unsigned n, dcomp d, const dcomp* A,
                   const dcomp *x, dcomp *y)
{
  zhpmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &Z_ONE, y, &N_ONE);
}




// y += d A^T x (A herm. dense in packed format)
inline void hemtva(const unsigned n, double d, const double* A,
		    const double *x, double *y)
{
  dspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &D_ONE, y, &N_ONE);
}
inline void hemtva(const unsigned n, float d, const float* A,
		    const float *x, float *y)
{
  
  sspmv_(JOB_STR+5, &n, &d, A, x, &N_ONE, &S_ONE, y, &N_ONE);
}
inline void hemtva(const unsigned n, scomp d, const scomp* A,
		    const scomp *x, scomp *y)
{
  scomp tmp;
  for(unsigned i=0; i<n; ++i) {
    y[i] += d * A[i+i*(i+1)/2] * x[i];
    for(unsigned j=i+1; j<n; ++j) {
      tmp = A[i+j*(j+1)/2];
      y[i] += d * conj(tmp) * x[j];
      y[j] += d* conj(tmp) * x[i];
    }
  }
}
inline void hemtva(const unsigned n, dcomp d, const dcomp* A,
		    const dcomp *x, dcomp *y)
{
  dcomp tmp;
  for(unsigned i=0; i<n; ++i) {
    y[i] += d * A[i+i*(i+1)/2] * x[i];
    for(unsigned j=i+1; j<n; ++j) {
      tmp = A[i+j*(j+1)/2];
      y[i] += d * conj(tmp) * x[j];
      y[j] += d* conj(tmp) * x[i];
    }
  }
}


// sovles Ax=B, A is a triangluar Matrix
inline void gtsv(const unsigned* n, const double* DiagLower,
                 const double* Diag, const double* DiagUpper,
                 const double* B, const int* INFO)
{
  dgtsv_(n, &N_ONE, DiagLower, Diag, DiagUpper, B, n, INFO);
}
inline void gtsv(const unsigned* n, const float* DiagLower,
                 const float* Diag, const float* DiagUpper,
                 const float* B, const int* INFO)
{
  sgtsv_(n, &N_ONE, DiagLower, Diag, DiagUpper, B, n, INFO);
}
inline void gtsv(const unsigned* n, const scomp* DiagLower,
                 const scomp* Diag, const scomp* DiagUpper,
                 const scomp* B, const int* INFO)
{
  cgtsv_(n, &N_ONE, DiagLower, Diag, DiagUpper, B, n, INFO);
}
inline void gtsv(const unsigned* n, const dcomp* DiagLower,
                 const dcomp* Diag, const dcomp* DiagUpper,
                 const dcomp* B, const int* INFO)
{
  zgtsv_(n, &N_ONE, DiagLower, Diag, DiagUpper, B, n, INFO);
}



// C = d A B, A is m x p, B is p x n
inline void gemm(const unsigned m, const unsigned p, const unsigned n,
                 const double d, const double* const A, const unsigned ldA,
                 const double* const B, const unsigned ldB,
                 double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ZERO, C, &ldC);
}
inline void gemm(const unsigned m, const unsigned p, const unsigned n,
                 const float d, const float* const A, const unsigned ldA,
                 const float* const B, const unsigned ldB,
                 float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ZERO, C, &ldC);
}
inline void gemm(const unsigned m, const unsigned p, const unsigned n,
                 const scomp d, const scomp* const A, const unsigned ldA,
                 const scomp* const B, const unsigned ldB,
                 scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ZERO, C, &ldC);
}
inline void gemm(const unsigned m, const unsigned p, const unsigned n,
                 const dcomp d, const dcomp* const A, const unsigned ldA,
                 const dcomp* const B, const unsigned ldB,
                 dcomp* C, unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ZERO, C, &ldC);
}


// C += d A B, A is m x p, B is p x n
inline void gemma(const unsigned m, const unsigned p, const unsigned n,
                  const double d, const double* const A, const unsigned ldA,
                  const double* const B, const unsigned ldB,
                  double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ONE, C, &ldC);
}
inline void gemma(const unsigned m, const unsigned p, const unsigned n,
                  const float d, const float* const A, const unsigned ldA,
                  const float* const B, const unsigned ldB,
                  float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ONE, C, &ldC);
}
inline void gemma(const unsigned m, const unsigned p, const unsigned n,
                  const scomp d, const scomp* const A, const unsigned ldA,
                  const scomp* const B, const unsigned ldB,
                  scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ONE, C, &ldC);
}
inline void gemma(const unsigned m, const unsigned p, const unsigned n,
                  const dcomp d, const dcomp* const A, const unsigned ldA,
                  const dcomp* const B, const unsigned ldB,
                  dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ONE, C, &ldC);
}

// C = d A^H B, A is m x p, B is m x n
inline void gemhm(const unsigned m, const unsigned p, const unsigned n,
                  const double d, const double* A, const unsigned ldA,
                  const double *B, const unsigned ldB,
                  double* C, const unsigned ldC)
{
  dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &D_ZERO, C, &ldC);
}
inline void gemhm(const unsigned m, const unsigned p, const unsigned n,
                  const float d, const float* const A, const unsigned ldA,
                  const float* const B, const unsigned ldB,
                  float* C, const unsigned ldC)
{
  sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &S_ZERO, C, &ldC);
}
inline void gemhm(const unsigned m, const unsigned p, const unsigned n,
                  const scomp d, const scomp* const A, const unsigned ldA,
                  const scomp* const B, const unsigned ldB,
                  scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &C_ZERO, C, &ldC);
}
inline void gemhm(const unsigned m, const unsigned p, const unsigned n,
                  const dcomp d, const dcomp* const A, const unsigned ldA,
                  const dcomp* const B, const unsigned ldB,
                  dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &Z_ZERO, C, &ldC);
}

// C += d A^H B, A is m x p, B is m x n
inline void gemhma(unsigned m, unsigned p, unsigned n, double d,
                   const double* const A, const unsigned ldA, const double* const B,
                   const unsigned ldB, double* C, unsigned ldC)
{
  dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &D_ONE, C, &ldC);
}
inline void gemhma(unsigned m, unsigned p, unsigned n, float d,
                   const float* const A, const unsigned ldA, const float* const B,
                   const unsigned ldB, float* C, unsigned ldC)
{
  sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &S_ONE, C, &ldC);
}
inline void gemhma(unsigned m, unsigned p, unsigned n, scomp d,
                   const scomp* const A, const unsigned ldA,
                   const scomp* const B, const unsigned ldB,
                   scomp* C, unsigned ldC)
{
  cgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &C_ONE, C, &ldC);
}
inline void gemhma(unsigned m, unsigned p, unsigned n, dcomp d,
                   const dcomp* const A, const unsigned ldA,
                   const dcomp* const B, const unsigned ldB,
                   dcomp* C, unsigned ldC)
{
  zgemm_(JOB_STR+7, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &Z_ONE, C, &ldC);
}


// C += d A^T B, A is m x p, B is m x n
inline void gemtma(unsigned m, unsigned p, unsigned n, double d,
                   const double* const A, const unsigned ldA, const double* const B,
                   const unsigned ldB, double* C, unsigned ldC)
{
  dgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &D_ONE, C, &ldC);
}
inline void gemtma(unsigned m, unsigned p, unsigned n, float d,
                   const float* const A, const unsigned ldA, const float* const B,
                   const unsigned ldB, float* C, unsigned ldC)
{
  sgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &S_ONE, C, &ldC);
}
inline void gemtma(unsigned m, unsigned p, unsigned n, scomp d,
                   const scomp* const A, const unsigned ldA,
                   const scomp* const B, const unsigned ldB,
                   scomp* C, unsigned ldC)
{
  cgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &C_ONE, C, &ldC);
}
inline void gemtma(unsigned m, unsigned p, unsigned n, dcomp d,
                   const dcomp* const A, const unsigned ldA,
                   const dcomp* const B, const unsigned ldB,
                   dcomp* C, unsigned ldC)
{
  zgemm_(JOB_STR+1, JOB_STR, &p, &n, &m, &d, A, &ldA, B, &ldB,
         &Z_ONE, C, &ldC);
}


// C = d A B^H, A is m x p, B is n x p
inline void gemmh(const unsigned m, const unsigned p, const unsigned n,
                  const double d, const double* const A, const unsigned ldA,
                  const double* const B, const unsigned ldB,
                  double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ZERO, C, &ldC);
}
inline void gemmh(const unsigned m, const unsigned p, const unsigned n,
                  const float d, const float* const A, const unsigned ldA,
                  const float* const B, const unsigned ldB,
                  float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ZERO, C, &ldC);
}
inline void gemmh(const unsigned m, const unsigned p, const unsigned n,
                  const scomp d, const scomp* const A, const unsigned ldA,
                  const scomp *B, const unsigned ldB,
                  scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ZERO, C, &ldC);
}
inline void gemmh(const unsigned m, const unsigned p, const unsigned n,
                  const dcomp d, const dcomp* A, const unsigned ldA,
                  const dcomp *B, const unsigned ldB,
                  dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ZERO, C, &ldC);
}

// C = d A B^T, A is m x p, B is n x p
inline void gemmt(const unsigned m, const unsigned p, const unsigned n,
                  const double d, const double* const A, const unsigned ldA,
                  const double* const B, const unsigned ldB,
                  double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ZERO, C, &ldC);
}
inline void gemmt(const unsigned m, const unsigned p, const unsigned n,
                  const float d, const float* const A, const unsigned ldA,
                  const float* const B, const unsigned ldB,
                  float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ZERO, C, &ldC);
}
inline void gemmt(const unsigned m, const unsigned p, const unsigned n,
                  const scomp d, const scomp* const A, const unsigned ldA,
                  const scomp *B, const unsigned ldB,
                  scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ZERO, C, &ldC);
}
inline void gemmt(const unsigned m, const unsigned p, const unsigned n,
                  const dcomp d, const dcomp* A, const unsigned ldA,
                  const dcomp *B, const unsigned ldB,
                  dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ZERO, C, &ldC);
}


// C += d A B^H, A is m x p, B is n x p
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const double d, const double* const A, const unsigned ldA,
                   const double* const B, const unsigned ldB,
                   double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const float d, const float* const A, const unsigned ldA,
                   const float* const B, const unsigned ldB,
                   float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const scomp d, const scomp* const A, const unsigned ldA,
                   const scomp* const B, const unsigned ldB,
                   scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const dcomp d, const dcomp* const A, const unsigned ldA,
                   const dcomp* const B, const unsigned ldB,
                   dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const double* const A, const unsigned ldA,
                   const double* const B, const unsigned ldB,
                   double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &D_ONE, A, &ldA, B, &ldB,
         &D_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const float* const A, const unsigned ldA,
                   const float* const B, const unsigned ldB,
                   float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &S_ONE, A, &ldA, B, &ldB,
         &S_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const scomp* const A, const unsigned ldA,
                   const scomp* const B, const unsigned ldB,
                   scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &C_ONE, A, &ldA, B, &ldB,
         &C_ONE, C, &ldC);
}
inline void gemmha(const unsigned m, const unsigned p, const unsigned n,
                   const dcomp* const A, const unsigned ldA,
                   const dcomp* const B, const unsigned ldB,
                   dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR+7, &m, &n, &p, &Z_ONE, A, &ldA, B, &ldB,
         &Z_ONE, C, &ldC);
}


// C += d A B^T, A is m x p, B is n x p
inline void gemmta(const unsigned m, const unsigned p, const unsigned n,
                   const double d, const double* const A, const unsigned ldA,
                   const double* const B, const unsigned ldB,
                   double* C, const unsigned ldC)
{
  dgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ONE, C, &ldC);
}
inline void gemmta(const unsigned m, const unsigned p, const unsigned n,
                   const float d, const float* const A, const unsigned ldA,
                   const float* const B, const unsigned ldB,
                   float* C, const unsigned ldC)
{
  sgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ONE, C, &ldC);
}
inline void gemmta(const unsigned m, const unsigned p, const unsigned n,
                   const scomp d, const scomp* const A, const unsigned ldA,
                   const scomp* const B, const unsigned ldB,
                   scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ONE, C, &ldC);
}
inline void gemmta(const unsigned m, const unsigned p, const unsigned n,
                   const dcomp d, const dcomp* const A, const unsigned ldA,
                   const dcomp* const B, const unsigned ldB,
                   dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ONE, C, &ldC);
}


// C = d A^H B^H, A is p x m, B is n x p
inline void gemhmh(const unsigned m, const unsigned p, const unsigned n,
                   const double d, const double* const A, const unsigned ldA,
                   const double* const B, const unsigned ldB,
                   double* C, const unsigned ldC)
{
  dgemm_(JOB_STR+1, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &D_ZERO, C, &ldC);
}
inline void gemhmh(const unsigned m, const unsigned p, const unsigned n,
                   const float d, const float* const A, const unsigned ldA,
                   const float* const B, const unsigned ldB,
                   float* C, const unsigned ldC)
{
  sgemm_(JOB_STR+1, JOB_STR+1, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &S_ZERO, C, &ldC);
}
inline void gemhmh(const unsigned m, const unsigned p, const unsigned n,
                   const scomp d, const scomp* const A, const unsigned ldA,
                   const scomp* const B, const unsigned ldB,
                   scomp* C, const unsigned ldC)
{
  cgemm_(JOB_STR+7, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &C_ZERO, C, &ldC);
}
inline void gemhmh(const unsigned m, const unsigned p, const unsigned n,
                   const dcomp d, const dcomp* const A, const unsigned ldA,
                   const dcomp* const B, const unsigned ldB,
                   dcomp* C, const unsigned ldC)
{
  zgemm_(JOB_STR+7, JOB_STR+7, &m, &n, &p, &d, A, &ldA, B, &ldB,
         &Z_ZERO, C, &ldC);
}

//C += d*AB, A is mxm (packed upper half is stored), B is mxn and regular matrix
inline void sygemma(const unsigned m, const unsigned n,
                    const double* const A, const double* const B,
                    const double d, double* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=i; k<m; k++) {
        if (i==k) {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
        } else {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
          C[j*m+k] += d*A[i+k*(k+1)/2]*B[i+j*m];
        }
      }
    }
  }
}
inline void sygemma(const unsigned m, const unsigned n,
                    const float* const A, const float* const B,
                    const float d, float* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=i; k<m; k++) {
        if (i==k) {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
        } else {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
          C[j*m+k] += d*A[i+k*(k+1)/2]*B[i+j*m];
        }
      }
    }
  }
}
inline void sygemma(const unsigned m, const unsigned n,
                    const dcomp* const A, const dcomp* const B,
                    const dcomp d, dcomp* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=i; k<m; k++) {
        if (i==k) {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
        } else {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
          C[j*m+k] += d*A[i+k*(k+1)/2]*B[i+j*m];
        }
      }
    }
  }
}
inline void sygemma(const unsigned m, const unsigned n,
                    const scomp* const A, const scomp* const B,
                    const scomp d, scomp* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=i; k<m; k++) {
        if (i==k) {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
        } else {
          C[j*m+i] += d*A[i+k*(k+1)/2]*B[k+j*m];
          C[j*m+k] += d*A[i+k*(k+1)/2]*B[i+j*m];
        }
      }
    }
  }
}

//C += d*AB, A is mxn and regular matrix, B is nxn (packed upper half is stored)
inline void gesymma(const unsigned m, const unsigned n,
                    const double* const A, const double* const B,
                    const double d, double* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=j; k<n; k++) {
        if (j==k)
          C[j*m+i] += d*A[i+k*m]*B[k+j*(j+1)/2];
        else {
          C[j*m+i] += d*A[i+k*m]*B[j+k*(k+1)/2];
          C[k*m+i] += d*A[i+j*m]*B[j+k*(k+1)/2];
        }
      }
    }
  }
}
inline void gesymma(const unsigned m, const unsigned n,
                    const float* const A, const float* const B,
                    const float d, float* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=j; k<n; k++) {
        if (j==k)
          C[j*m+i] += d*A[i+k*m]*B[k+j*(j+1)/2];
        else {
          C[j*m+i] += d*A[i+k*m]*B[j+k*(k+1)/2];
          C[k*m+i] += d*A[i+j*m]*B[j+k*(k+1)/2];
        }
      }
    }
  }
}
inline void gesymma(const unsigned m, const unsigned n,
                    const scomp* const A, const scomp* const B,
                    const scomp d, scomp* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=j; k<n; k++) {
        if (j==k)
          C[j*m+i] += d*A[i+k*m]*B[k+j*(j+1)/2];
        else {
          C[j*m+i] += d*A[i+k*m]*B[j+k*(k+1)/2];
          C[k*m+i] += d*A[i+j*m]*B[j+k*(k+1)/2];
        }
      }
    }
  }
}
inline void gesymma(const unsigned m, const unsigned n,
                    const dcomp* const A, const dcomp* const B,
                    const dcomp d, dcomp* const C)
{
  for (unsigned i=0; i<m; i++) {
    for (unsigned j=0; j<n; j++) {
      for (unsigned k=j; k<n; k++) {
        if (j==k)
          C[j*m+i] += d*A[i+k*m]*B[k+j*(j+1)/2];
        else {
          C[j*m+i] += d*A[i+k*m]*B[j+k*(k+1)/2];
          C[k*m+i] += d*A[i+j*m]*B[j+k*(k+1)/2];
        }
      }
    }
  }
}


// C += d A^H A, C is a symm. matrix (packed upper half is stored), A is mxn
inline void symhm(const unsigned m, const unsigned n, const double* const A,
                  const double d, double* C)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<=j; ++i) {
      double sum = 0.0;
      for (unsigned k=0; k<m; ++k) sum += A[k+i*m] * A[k+j*m];
      C[i+j*(j+1)/2] += d * sum;
    }
  }
}
inline void symhm(const unsigned m, const unsigned n, const float* const A,
                  const float d, float* C)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<=j; ++i) {
      float sum = 0.0;
      for (unsigned k=0; k<m; ++k) sum += A[k+i*m] * A[k+j*m];
      C[i+j*(j+1)/2] += d * sum;
    }
  }
}
inline void symhm(unsigned m, unsigned n, scomp* A, scomp d, scomp* C)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<=j; ++i) {
      scomp sum = C_ZERO;
      for (unsigned k=0; k<m; ++k) sum += conj(A[k+i*m]) * A[k+j*m];
      C[i+j*(j+1)/2] += d * sum;
    }
  }
}
inline void symhm(unsigned m, unsigned n, dcomp* A, dcomp d, dcomp* C)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<=j; ++i) {
      dcomp sum = Z_ZERO;
      for (unsigned k=0; k<m; ++k) sum += conj(A[k+i*m]) * A[k+j*m];
      C[i+j*(j+1)/2] += d * sum;
    }
  }
}

// C += d A A^H, C is a symm. matrix (packed upper half is stored), A is mxn
inline void symmh(unsigned m, unsigned n, double* A, double d, double* C)
{
  for (unsigned k=0; k<n; ++k) {
    for (unsigned j=0; j<=n; ++j) {
      double e = d * A[j+k*m];
      for (unsigned i=0; i<j; ++i) C[i+j*(j+1)/2] += e * A[i+k*m];
    }
  }
}
inline void symmh(unsigned m, unsigned n, float* A, float d, float* C)
{
  for (unsigned k=0; k<n; ++k) {
    for (unsigned j=0; j<n; ++j) {
      float e = d * A[j+k*m];
      for (unsigned i=0; i<=j; ++i) C[i+j*(j+1)/2] += e * A[i+k*m];
    }
  }
}
inline void symmh(unsigned m, unsigned n, scomp* A, scomp d, scomp* C)
{
  for (unsigned k=0; k<n; ++k) {
    for (unsigned j=0; j<n; ++j) {
      scomp e = d * conj(A[j+k*m]);
      for (unsigned i=0; i<=j; ++i) C[i+j*(j+1)/2] += e * A[i+k*m];
    }
  }
}
inline void symmh(unsigned m, unsigned n, dcomp* A, dcomp d, dcomp* C)
{
  for (unsigned k=0; k<n; ++k) {
    for (unsigned j=0; j<n; ++j) {
      dcomp e = d * conj(A[j+k*m]);
      for (unsigned i=0; i<=j; ++i) C[i+j*(j+1)/2] += e * A[i+k*m];
    }
  }
}
// Singular Value Decomposition
inline int gesvdS(unsigned m, unsigned n, double* A, double* S,
                  double* U, unsigned ldU, double* VT, unsigned ldVT,
                  unsigned nwk, double* wk)
{
  int INF;
  dgesvd_(JOB_STR+3, JOB_STR+3, &m, &n, A, &m, S, U, &ldU, VT, &ldVT,
          wk, &nwk, &INF);
  return INF;
}
inline int gesvd(unsigned m, unsigned n, double* A, double* S,
                 double* VT, unsigned ldVT, unsigned nwk, double* wk)
{
  int INF;
  dgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
          wk, &nwk, &INF);
  return INF;
}
inline int gesvd(unsigned m, unsigned n, float* A, float* S,
                 float* VT, unsigned ldVT, unsigned nwk, float* wk)
{
  int INF;
  sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
          wk, &nwk, &INF);
  return INF;
}
inline int gesvd(unsigned m, unsigned n, float* A, double* S,
                 float* VT, unsigned ldVT, unsigned nwk, float* wk)
{
  int INF;
  // workaround (needs to be improved)
  float* Sf = new float[MIN(m,n)];
  sgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, Sf, A, &m, VT, &ldVT,
          wk, &nwk, &INF);
  for (unsigned i=0; i<MIN(m,n); i++)
    S[i] = (double)Sf[i];
  delete [] Sf;
  return INF;
}
inline int gesvd(unsigned m, unsigned n, scomp* A, float* S,
                 scomp* VT, unsigned ldVT, unsigned nwk, scomp* wk)
{
  int INF;
  float* rwk = new float[5*MIN(m,n)];
  cgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
          wk, &nwk, rwk, &INF);
  delete [] rwk;
  return INF;
}
inline int gesvd(unsigned m, unsigned n, scomp* A, double* S,
                 scomp* VT, unsigned ldVT, unsigned nwk, scomp* wk)
{
  int INF;
  float* rwk = new float[nwk];
  // workaround (needs to be improved)
  float* Sf = new float[MIN(m,n)];
  cgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, Sf, A, &m, VT, &ldVT,
          wk, &nwk, rwk, &INF);
  for (unsigned i=0; i<MIN(m,n); i++)
    S[i] = Sf[i];
  delete [] Sf;
  delete [] rwk;
  return INF;
}
inline int gesvd(unsigned m, unsigned n, dcomp* A, double* S,
                 dcomp* VT, unsigned ldVT, unsigned nwk, dcomp* wk)
{
  int INF;
  double* rwk = new double[5*MIN(m,n)];
  zgesvd_(JOB_STR+2, JOB_STR+3, &m, &n, A, &m, S, A, &m, VT, &ldVT,
          wk, &nwk, rwk, &INF);
  delete [] rwk;
  return INF;
}
inline int gesvd(unsigned m, unsigned n, double* A, double* S,
                 double* U, unsigned ldU, double* VT, unsigned ldVT,
                 unsigned nwk, double* wk)
{
  int INF;
  dgesvd_(JOB_STR+9, JOB_STR+9, &m, &n, A, &m, S, U, &ldU, VT, &ldVT,
          wk, &nwk, &INF);
  return INF;
}
inline int gesvd(unsigned m, unsigned n, float* A, double* S,
                 float* U, unsigned ldU, float* VT, unsigned ldVT,
                 unsigned nwk, float* wk)
{
  int INF;
  // workaround (needs to be improved)
  float* Sf = new float[MIN(m,n)];
  sgesvd_(JOB_STR+9, JOB_STR+9, &m, &n, A, &m, Sf, U, &ldU, VT, &ldVT,
          wk, &nwk, &INF);
  for (unsigned i=0; i<MIN(m,n); i++)
    S[i] = Sf[i];
  delete [] Sf;
  return INF;
}
inline int gesvd(unsigned m, unsigned n, scomp* A, double* S,
                 scomp* U, unsigned ldU, scomp* VT, unsigned ldVT,
                 unsigned nwk, scomp* wk)
{
  int INF;
  float* rwk = new float[nwk];
  // workaround (needs to be improved)
  float* Sf = new float[MIN(m,n)];
  cgesvd_(JOB_STR+9, JOB_STR+9, &m, &n, A, &m, Sf, U, &ldU, VT, &ldVT,
          wk, &nwk, rwk, &INF);
  for (unsigned i=0; i<MIN(m,n); i++)
    S[i] = Sf[i];
  delete [] Sf;
  delete [] rwk;
  return INF;
}
inline int gesvd(unsigned m, unsigned n, dcomp* A, double* S,
                 dcomp* U, unsigned ldU, dcomp* VT, unsigned ldVT,
                 unsigned nwk, dcomp* wk)
{
  int INF;
  double* rwk = new double[5*MIN(n,m)];
  zgesvd_(JOB_STR+9, JOB_STR+9, &m, &n, A, &m, S, U, &ldU, VT, &ldVT,
          wk, &nwk, rwk, &INF);
  delete [] rwk;
  return INF;
}

// compute singular values
inline int svals(unsigned m, unsigned n, double* A, double* S,
                 unsigned nwk, double* wk)
{
  int INF;
  dgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n, wk, &nwk, &INF);
  return INF;
}
inline int svals(unsigned m, unsigned n, float* A, float* S,
                 unsigned nwk, float* wk)
{
  int INF;
  sgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n, wk, &nwk, &INF);
  return INF;
}
inline int svals(unsigned m, unsigned n, scomp* A, float* S,
                 unsigned nwk, scomp* wk)
{
  int INF;
  float *rwk = new float[5*MIN(m,n)];
  cgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n,
          wk, &nwk, rwk, &INF);
  delete [] rwk;
  return INF;
}
inline int svals(unsigned m, unsigned n, dcomp* A, double* S,
                 unsigned nwk, dcomp* wk)
{
  int INF;
  double *rwk = new double[5*MIN(m,n)];
  zgesvd_(JOB_STR, JOB_STR, &m, &n, A, &m, S, A, &m, A, &n,
          wk, &nwk, rwk, &INF);
  delete [] rwk;
  return INF;
}


// triangular factorisation
inline int getrf(const unsigned n, double* A, unsigned* ipiv)
{
  int INF;
  dgetrf_(&n, &n, A, &n, ipiv, &INF);
  return INF;
}
inline int getrf(const unsigned n, float* A, unsigned* ipiv)
{
  int INF;
  sgetrf_(&n, &n, A, &n, ipiv, &INF);
  return INF;
}
inline int getrf(const unsigned n, scomp* A, unsigned* ipiv)
{
  int INF;
  cgetrf_(&n, &n, A, &n, ipiv, &INF);
  return INF;
}
inline int getrf(const unsigned n, dcomp* A, unsigned* ipiv)
{
  int INF;
  zgetrf_(&n, &n, A, &n, ipiv, &INF);
  return INF;
}
// triangular factorisation A is m by n
  inline int getrf(const unsigned m, const unsigned n, double* A, unsigned* ipiv)
  {
    int INF;
    dgetrf_(&m, &n, A, &m, ipiv, &INF); 
    return INF;
  }
  inline int getrf(const unsigned m, const unsigned n, float* A, unsigned* ipiv)
  {
    int INF;
    sgetrf_(&m, &n, A, &m, ipiv, &INF);
    return INF;
  }
  inline int getrf(const unsigned m, const unsigned n, scomp* A, unsigned* ipiv)
  {
    int INF;
    cgetrf_(&m, &n, A, &m, ipiv, &INF);
    return INF;
  }
  inline int getrf(const unsigned m, const unsigned n, dcomp* A, unsigned* ipiv)
  {
    int INF;
    zgetrf_(&m, &n, A, &m, ipiv, &INF);
    return INF;
  }

// upper triangular packed MV
inline void utrpv(const unsigned n, double* A, double* x)
{
  dtpmv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, x, &N_ONE);
}
inline void utrpv(const unsigned n, float* A, float* x)
{
  stpmv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, x, &N_ONE);
}
inline void utrpv(const unsigned n, scomp* A, scomp* x)
{
  ctpmv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, x, &N_ONE);
}
inline void utrpv(const unsigned n, dcomp* A, dcomp* x)
{
  ztpmv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, x, &N_ONE);
}

// lower triangular packed transpose MV
inline void ltrphv(const unsigned n, double* A, double* x)
{
  dtpmv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, x, &N_ONE);
}
inline void ltrphv(const unsigned n, float* A, float* x)
{
  stpmv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, x, &N_ONE);
}
inline void ltrphv(const unsigned n, scomp* A, scomp* x)
{
  ctpmv_(JOB_STR+6, JOB_STR+7, JOB_STR+5, &n, A, x, &N_ONE);
}
inline void ltrphv(const unsigned n, dcomp* A, dcomp* x)
{
  ztpmv_(JOB_STR+6, JOB_STR+7, JOB_STR+5, &n, A, x, &N_ONE);
}

// QR factorisation
inline int geqrf(const unsigned m, const unsigned n, double* A,
                 double* tau, int nwk, double* wk)
{
  int INF;
  dgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, float* A,
                 float* tau, int nwk, float* wk)
{
  int INF;
  sgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, scomp* A,
                 scomp* tau, int nwk, scomp* wk)
{
  int INF;
  cgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, dcomp* A,
                 dcomp* tau, int nwk, dcomp* wk)
{
  int INF;
  zgeqrf_(&m, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, double* A,
                 const unsigned ldA, double* tau, int nwk, double* wk)
{
  int INF;
  dgeqrf_(&m, &n, A, &ldA, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, float* A,
                 const unsigned ldA, float* tau, int nwk, float* wk)
{
  int INF;
  sgeqrf_(&m, &n, A, &ldA, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, scomp* A,
                 const unsigned ldA, scomp* tau, int nwk, scomp* wk)
{
  int INF;
  cgeqrf_(&m, &n, A, &ldA, tau, wk, &nwk, &INF);
  return INF;
}
inline int geqrf(const unsigned m, const unsigned n, dcomp* A,
                 const unsigned ldA, dcomp* tau, int nwk, dcomp* wk)
{
  int INF;
  zgeqrf_(&m, &n, A, &ldA, tau, wk, &nwk, &INF);
  return INF;
}

// QR factorisation (single vector)
inline void geqrfs(const unsigned n, const double* x, const double* tau)
{ 
  dlarfg_(&n, x, x+1, &N_ONE, tau);
}

inline void geqrfs(const unsigned n, const float* x, const float* tau)
{ 
  slarfg_(&n, x, x+1, &N_ONE, tau);
}

inline void geqrfs(const unsigned n, const scomp* x, const scomp* tau)
{ 
  clarfg_(&n, x, x+1, &N_ONE, tau);
}

inline void geqrfs(const unsigned n, const dcomp* x, const dcomp* tau)
{ 
  zlarfg_(&n, x, x+1, &N_ONE, tau);
}


//truncated QR-factorisation(RRQR) - truncated version of dgeqp3(...)
//with absolute/relative accuracy atrunc and rtrunc
//column permutation stored in perm
//set perm[i] = 0 at the beginnnig for all i, if columns are not permutated
inline int geqp3trunc(const unsigned m, const unsigned n, double* A,
		      unsigned* perm, double* tau, unsigned& ntrunc,
		      double atrunc, double rtrunc, unsigned nwk, 
		      double* wk)
{
  int INF;
  dgeqp3trunc_(&m, &n, A, &m, perm, tau, &ntrunc, &atrunc, &rtrunc,
	       wk, &nwk, &INF);
  return INF;
}


// Multiply a general Matrix with the Q-Matrix (QR factorization), Q C
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 double* A, double* tau, double* C,
                 unsigned nwk, double* wk)
{
  int INF;
  dormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 float* A, float* tau, float* C,
                 unsigned nwk, float* wk)
{
  int INF;
  sormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 scomp* A, scomp* tau, scomp* C,
                 unsigned nwk, scomp* wk)
{
  int INF;
  cunmqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 dcomp* A, dcomp* tau, dcomp* C,
                 unsigned nwk, dcomp* wk)
{
  int INF;
  zunmqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 double* A, const unsigned ldA, double* tau, double* C,
                 const unsigned ldC, unsigned nwk, double* wk)
{
  int INF;
  dormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 float* A, const unsigned ldA, float* tau, float* C,
                 const unsigned ldC, unsigned nwk, float* wk)
{
  int INF;
  sormqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 scomp* A, const unsigned ldA, scomp* tau, scomp* C,
                 const unsigned ldC, unsigned nwk, scomp* wk)
{
  int INF;
  cunmqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqr(const unsigned m, const unsigned n, const unsigned p,
                 dcomp* A, const unsigned ldA, dcomp* tau, dcomp* C,
                 const unsigned ldC, unsigned nwk, dcomp* wk)
{
  int INF;
  zunmqr_(JOB_STR+6, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}

// Multiply a general Matrix with the Q-Matrix (QR factorization), Q C
// just a single householder reflector
inline void ormqrs(const unsigned m, const unsigned n, 
		  double* A, double* tau, double* C, 
		  const unsigned ldC, const unsigned nwk, double* wk)
{
  dlarfb_(JOB_STR+6, JOB_STR, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline void ormqrs(const unsigned m, const unsigned n, 
		   float* A, float* tau, float* C,
		   const unsigned ldC, const unsigned nwk, float* wk)
{
  slarfb_(JOB_STR+6, JOB_STR, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline void ormqrs(const unsigned m, const unsigned n, 
		  scomp* A, scomp* tau, scomp* C, 
		  const unsigned ldC, const unsigned nwk, scomp* wk)
{
  clarfb_(JOB_STR+6, JOB_STR, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline void ormqrs(const unsigned m, const unsigned n, 
		   dcomp* A, dcomp* tau, dcomp* C,
		   const unsigned ldC, const unsigned nwk, dcomp* wk)
{
  zlarfb_(JOB_STR+6, JOB_STR, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}


// Q^H C
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  double* A, const unsigned ldA, double* tau, double* C,
                  const unsigned ldC, unsigned nwk, double* wk)
{
  int INF;
  dormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  float* A, const unsigned ldA, float* tau, float* C,
                  const unsigned ldC, unsigned nwk, float* wk)
{
  int INF;
  sormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  scomp* A, const unsigned ldA, scomp* tau, scomp* C,
                  const unsigned ldC, unsigned nwk, scomp* wk)
{
  int INF;
  cunmqr_(JOB_STR+6, JOB_STR+7, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  dcomp* A, const unsigned ldA, dcomp* tau, dcomp* C,
                  const unsigned ldC, unsigned nwk, dcomp* wk)
{
  int INF;
  zunmqr_(JOB_STR+6, JOB_STR+7, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}

inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  double* A, const unsigned ldA, double* tau, double* C,
                  unsigned nwk, double* wk)
{
  int INF;
  dormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &m, wk, &nwk,
          &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  float* A, const unsigned ldA, float* tau, float* C,
                  unsigned nwk, float* wk)
{
  int INF;
  sormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &ldA, tau, C, &m, wk, &nwk,
          &INF);
  return INF;
}

inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  double* A, double* tau, double* C,
                  unsigned nwk, double* wk)
{
  int INF;
  dormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  float* A, float* tau, float* C,
                  unsigned nwk, float* wk)
{
  int INF;
  sormqr_(JOB_STR+6, JOB_STR+1, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  scomp* A, scomp* tau, scomp* C,
                  unsigned nwk, scomp* wk)
{
  int INF;
  cunmqr_(JOB_STR+6, JOB_STR+7, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}
inline int ormqrh(const unsigned m, const unsigned n, const unsigned p,
                  dcomp* A, dcomp* tau, dcomp* C,
                  unsigned nwk, dcomp* wk)
{
  int INF;
  zunmqr_(JOB_STR+6, JOB_STR+7, &m, &n, &p, A, &m, tau, C, &m, wk, &nwk, &INF);
  return INF;
}

// Multiply a general Matrix with the Q-Matrix (QR factorization), Q^H C
// just a single householder reflector
inline void ormqrsh(const unsigned m, const unsigned n, 
		    double* A, double* tau, double* C, 
		    const unsigned ldC, const unsigned nwk, double* wk)
{
  dlarfb_(JOB_STR+6, JOB_STR+1, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline void ormqrsh(const unsigned m, const unsigned n, 
		    float* A, float* tau, float* C,
		    const unsigned ldC, const unsigned nwk, float* wk)
{
  slarfb_(JOB_STR+6, JOB_STR+1, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline void ormqrsh(const unsigned m, const unsigned n, 
		    scomp* A, scomp* tau, scomp* C, 
		    const unsigned ldC, const unsigned nwk, scomp* wk)
{
  clarfb_(JOB_STR+6, JOB_STR+1, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline void ormqrsh(const unsigned m, const unsigned n, 
		    dcomp* A, dcomp* tau, dcomp* C,
		    const unsigned ldC, const unsigned nwk, dcomp* wk)
{
  zlarfb_(JOB_STR+6, JOB_STR+1, JOB_STR+10, JOB_STR+7, &m, &n, 
	  &N_ONE, A, &m, tau, &N_ONE, C, &ldC, wk, &nwk);
}

inline int morqr(const unsigned m, const unsigned n, const unsigned p,
                 double* A, const unsigned ldA, double* tau, double* C,
                 const unsigned ldC, unsigned nwk, double* wk)
{
  int INF;
  dormqr_(JOB_STR+8, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int morqr(const unsigned m, const unsigned n, const unsigned p,
                 float* A, const unsigned ldA, float* tau, float* C,
                 const unsigned ldC, unsigned nwk, float* wk)
{
  int INF;
  sormqr_(JOB_STR+8, JOB_STR, &m, &n, &p, A, &ldA, tau, C, &ldC, wk, &nwk,
          &INF);
  return INF;
}
inline int morqr(const unsigned m, const unsigned n, const unsigned p,
                 scomp* A, const unsigned ldA, scomp* tau, scomp* C,
                 const unsigned ldC, unsigned nwk, scomp* wk)
{
  exit(1);
  return 0;
}
inline int morqr(const unsigned m, const unsigned n, const unsigned p,
                 dcomp* A, const unsigned ldA, dcomp* tau, dcomp* C,
                 const unsigned ldC, unsigned nwk, dcomp* wk)
{
  exit(1);
  return 0;
}

inline void ger(unsigned M, unsigned N, double d, double* X, unsigned INCX,
                double* y, unsigned INCY, double* A, unsigned LDA)
{
  dger_(&M, &N, &d, X, &INCX, y, &INCY, A, &LDA);
}

inline void ger(unsigned M, unsigned N, float d, float* X, unsigned INCX,
                float* y, unsigned INCY, float* A, unsigned LDA)
{
  sger_(&M, &N, &d, X, &INCX, y, &INCY, A, &LDA);
}

inline void ger(unsigned M, unsigned N, scomp d, scomp* X, unsigned INCX,
                scomp* y, unsigned INCY, scomp* A, unsigned LDA)
{
  cgerc_(&M, &N, &d, X, &INCX, y, &INCY, A, &LDA);
}

inline void ger(unsigned M, unsigned N, dcomp d, dcomp* X, unsigned INCX,
                dcomp* y, unsigned INCY, dcomp* A, unsigned LDA)
{
  zgerc_(&M, &N, &d, X, &INCX, y, &INCY, A, &LDA);
}

// return Q-Matrix (QR factorization) in A
inline int orgqr(const unsigned m, const unsigned n, double* A, double* tau,
                 unsigned nwk, double* wk)
{
  int INF;
  dorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int orgqr(const unsigned m, const unsigned n, float* A, float* tau,
                 unsigned nwk, float* wk)
{
  int INF;
  sorgqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int orgqr(const unsigned m, const unsigned n, scomp* A, scomp* tau,
                 unsigned nwk, scomp* wk)
{
  int INF;
  cungqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}
inline int orgqr(const unsigned m, const unsigned n, dcomp* A, dcomp* tau,
                 unsigned nwk, dcomp* wk)
{
  int INF;
  zungqr_(&m, &n, &n, A, &m, tau, wk, &nwk, &INF);
  return INF;
}

// B=A^H, A in R^(m,n)
inline void transpose(const unsigned m, const unsigned n, const double* A, 
		      double* B)
{
  for (unsigned i=0; i<m; ++i) dcopy_(&n, A+i, &m, B+i*n, &N_ONE);
}
inline void transpose(const unsigned m, const unsigned n, const float* A, 
		      float* B)
{
  for (unsigned i=0; i<m; ++i) scopy_(&n, A+i, &m, B+i*n, &N_ONE);
}
inline void transpose(const unsigned m, const unsigned n, const scomp* A, 
		      scomp* B)
{
  for (unsigned i=0; i<m; i++)
    for (unsigned j=0; j<n; j++) B[j+i*n] = conj(A[i+j*m]);
}
inline void transpose(const unsigned m, const unsigned n, const dcomp* A, 
		      dcomp* B)
{
  for (unsigned i=0; i<m; i++)
    for (unsigned j=0; j<n; j++) B[j+i*n] = conj(A[i+j*m]);
}

//same as above with leading dimension
inline void transpose(const unsigned m, const unsigned n, const double* A, 
		      const unsigned ldA, double* B, const unsigned ldB)
{
  for (unsigned i=0; i<m; ++i) dcopy_(&n, A+i, &ldA, B+i*ldB, &N_ONE);
}
inline void transpose(const unsigned m, const unsigned n, const float* A, 
		      const unsigned ldA, float* B, const unsigned ldB)
{
  for (unsigned i=0; i<m; ++i) scopy_(&n, A+i, &ldA, B+i*ldB, &N_ONE);
}
inline void transpose(const unsigned m, const unsigned n, const scomp* A, 
		      const unsigned ldA, scomp* B, const unsigned ldB)
{
  for (unsigned i=0; i<m; i++)
    for (unsigned j=0; j<n; j++) B[j+i*ldB] = conj(A[i+j*ldA]);
}
inline void transpose(const unsigned m, const unsigned n, const dcomp* A, 
		      const unsigned ldA, dcomp* B, const unsigned ldB)
{
  for (unsigned i=0; i<m; i++)
    for (unsigned j=0; j<n; j++) B[j+i*ldB] = conj(A[i+j*ldA]);
}

// product of upper triangular and regular Matrix  M = R A, R mxp, A nxp
inline void utrgemmh(const unsigned m, const unsigned p, const unsigned n,
                     const double* const R, const unsigned ldR,
                     const double* const A, const unsigned ldA,
                     double* const M, const unsigned ldM)
{
  for (unsigned i=0; i<m; ++i) {
    for (unsigned j=0; j<n; ++j) {
      double d = D_ZERO;
      for (unsigned l=i; l<p; ++l)
        d += R[i+l*ldR]*A[j+l*ldA];
      M[i+j*ldM] = d;
    }
  }
}
inline void utrgemmh(const unsigned m, const unsigned p, const unsigned n,
                     const float* const R, const unsigned ldR,
                     const float* const A, const unsigned ldA,
                     float* const M, const unsigned ldM)
{
  for (unsigned i=0; i<m; ++i) {
    for (unsigned j=0; j<n; ++j) {
      float d = S_ZERO;
      for (unsigned l=i; l<p; ++l)
        d += R[i+l*ldR]*A[j+l*ldA];
      M[i+j*ldM] = d;
    }
  }
}
inline void utrgemmh(const unsigned m, const unsigned p, const unsigned n,
                     const scomp* const R, const unsigned ldR,
                     const scomp* const A, const unsigned ldA,
                     scomp* const M, const unsigned ldM)
{
  for (unsigned i=0; i<m; ++i) {
    for (unsigned j=0; j<n; ++j) {
      scomp d = C_ZERO;
      for (unsigned l=i; l<p; ++l)
        d += R[i+l*ldR]*A[j+l*ldA];
      M[i+j*ldM] = d;
    }
  }
}
inline void utrgemmh(const unsigned m, const unsigned p, const unsigned n,
                     const dcomp* const R, const unsigned ldR,
                     const dcomp* const A, const unsigned ldA,
                     dcomp* const M, const unsigned ldM)
{
  for (unsigned i=0; i<m; ++i) {
    for (unsigned j=0; j<n; ++j) {
      dcomp d = Z_ZERO;
      for (unsigned l=i; l<p; ++l)
        d += R[i+l*ldR]*A[j+l*ldA];
      M[i+j*ldM] = d;
    }
  }
}

// product of two upper triangular matrices M = R1 R2^H, R1 mxp, R2 nxp
inline void utrmmh(const unsigned m, const unsigned p, const unsigned n,
                   const double* const R1, const unsigned ldR1,
                   const double* const R2, const unsigned ldR2,
                   double* const M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      double d = D_ZERO;
      unsigned ij = MAX(i,j);
      for (unsigned l=ij; l<p; ++l)
        d += R1[i+l*ldR1]*R2[j+l*ldR2];
      M[i+j*m] = d;
    }
  }
}
inline void utrmmh(const unsigned m, const unsigned p, const unsigned n,
                   const float* const R1, const unsigned ldR1,
                   const float* const R2, const unsigned ldR2,
                   float* const M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      float d = S_ZERO;
      unsigned ij = MAX(i,j);
      for (unsigned l=ij; l<p; ++l)
        d += R1[i+l*ldR1]*R2[j+l*ldR2];
      M[i+j*m] = d;
    }
  }
}
inline void utrmmh(const unsigned m, const unsigned p, const unsigned n,
                   const scomp* const R1, const unsigned ldR1,
                   const scomp* const R2, const unsigned ldR2,
                   scomp* const M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      scomp d = C_ZERO;
      unsigned ij = MAX(i,j);
      for (unsigned l=ij; l<p; ++l)
        d += R1[i+l*ldR1]*conj(R2[j+l*ldR2]);
      M[i+j*m] = d;
    }
  }
}
inline void utrmmh(const unsigned m, const unsigned p, const unsigned n,
                   const dcomp* const R1, const unsigned ldR1,
                   const dcomp* const R2, const unsigned ldR2,
                   dcomp* const M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      dcomp d = Z_ZERO;
      unsigned ij = MAX(i,j);
      for (unsigned l=ij; l<p; ++l)
        d += R1[i+l*ldR1]*conj(R2[j+l*ldR2]);
      M[i+j*m] = d;
    }
  }
}
// product of two upper triangular matrices M += R1 R2^T, R1 mxp, R2 nxp
inline void utrmmha(const unsigned m, const unsigned p, const unsigned n,
                    const double* const R1, const unsigned ldR1,
                    const double* const R2, const unsigned ldR2,
                    double* const M)
{
  for (unsigned l=0; l<p; l++) {
    for (unsigned j=p; j<n; j++) {
      for (unsigned i=p; i<m; i++) {
        M[i+j*m] += R1[i+l*ldR1]*R2[j+l*ldR2];
      }
    }
  }
}
inline void utrmmha(const unsigned m, const unsigned p, const unsigned n,
                    const float* const R1, const unsigned ldR1,
                    const float* const R2, const unsigned ldR2,
                    float* const M)
{
  for (unsigned l=0; l<p; l++) {
    for (unsigned j=p; j<n; j++) {
      for (unsigned i=p; i<m; i++) {
        M[i+j*m] += R1[i+l*ldR1]*R2[j+l*ldR2];
      }
    }
  }
}
inline void utrmmha(const unsigned m, const unsigned p, const unsigned n,
                    const scomp* const R1, const unsigned ldR1,
                    const scomp* const R2, const unsigned ldR2,
                    scomp* const M)
{
  for (unsigned l=0; l<p; l++) {
    for (unsigned j=p; j<n; j++) {
      for (unsigned i=p; i<m; i++) {
        M[i+j*m] += R1[i+l*ldR1]*conj(R2[j+l*ldR2]);
      }
    }
  }
}
inline void utrmmha(const unsigned m, const unsigned p, const unsigned n,
                    const dcomp* const R1, const unsigned ldR1,
                    const dcomp* const R2, const unsigned ldR2,
                    dcomp* const M)
{
  for (unsigned l=0; l<p; l++) {
    for (unsigned j=p; j<n; j++) {
      for (unsigned i=p; i<m; i++) {
        M[i+j*m] += R1[i+l*ldR1]*conj(R2[j+l*ldR2]);
      }
    }
  }
}

// product of an upper triangular matrix U and a matrix A, A:=U A
inline void utrgemm(unsigned m, unsigned n, double* U, unsigned ldU,
                    double* A, unsigned ldA)
{
  dtrmm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR, &m, &n, &D_ONE, U, &ldU,
         A, &ldA);
}

// A:=A U^H
inline void geutrmm(unsigned m, unsigned n, double* A, unsigned ldA,
                    double* U, unsigned ldU)
{
  dtrmm_(JOB_STR+8, JOB_STR+5, JOB_STR+1, JOB_STR, &m, &n, &D_ONE, U, &ldU,
         A, &ldA);
}

// C += d A U^T, U upper triangular matrix in packed storage
inline void geputrmmh(double d, unsigned m, unsigned n, double* A,
                      double* U, double* C)
{
  for (unsigned j=0; j<n; ++j)
    for (unsigned l=j; l<n; ++l) {
      double e = d*U[j+l*(l+1)/2];
      for (unsigned i=0; i<m; ++i) C[i+j*m] += e*A[i+l*m];
    }
}

inline void geputrmmh(float d, unsigned m, unsigned n, float* A,
                      float* U, float* C)
{
  for (unsigned j=0; j<n; ++j)
    for (unsigned l=j; l<n; ++l) {
      float e = d*U[j+l*(l+1)/2];
      for (unsigned i=0; i<m; ++i) C[i+j*m] += e*A[i+l*m];
    }
}

inline void geputrmmh(dcomp d, unsigned m, unsigned n, dcomp* A,
                      dcomp* U, dcomp* C)
{
  for (unsigned j=0; j<n; ++j)
    for (unsigned l=j; l<n; ++l) {
      dcomp e = d*conj(U[j+l*(l+1)/2]);
      for (unsigned i=0; i<m; ++i) C[i+j*m] += e*A[i+l*m];
    }
}

inline void geputrmmh(scomp d, unsigned m, unsigned n, scomp* A,
                      scomp* U, scomp* C)
{
  for (unsigned j=0; j<n; ++j)
    for (unsigned l=j; l<n; ++l) {
      scomp e = d*conj(U[j+l*(l+1)/2]);
      for (unsigned i=0; i<m; ++i) C[i+j*m] += e*A[i+l*m];
    }
}

// C += A U^T, U upper triangular matrix
inline void geutrTmm(unsigned m, unsigned n, double* U, unsigned ldU,
                     double* A, unsigned ldA)
{
  dtrmm_(JOB_STR+8, JOB_STR+5, JOB_STR+1, JOB_STR, &m, &n, &D_ONE, U, &ldU,
         A, &ldA);
}

// product of two upper triangular matrices M = R1 R2^T, R1 mxp, R2 nxp
inline void utrmmh(unsigned m, unsigned p, unsigned n, double* R1,
                   unsigned ldR1, double* R2, unsigned ldR2, double* M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      double d = D_ZERO;
      for (unsigned l=MAX(i,j); l<p; ++l) d += R1[i+l*ldR1] * R2[j+l*ldR2];
      M[i+j*m] = d;
    }
  }
}
inline void utrmmh(unsigned m, unsigned p, unsigned n, float* R1,
                   unsigned ldR1, float* R2, unsigned ldR2, float* M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      float d = S_ZERO;
      for (unsigned l=MAX(i,j); l<p; ++l) d += R1[i+l*ldR1] * R2[j+l*ldR2];
      M[i+j*m] = d;
    }
  }
}
inline void utrmmh(unsigned m, unsigned p, unsigned n, scomp* R1,
                   unsigned ldR1, scomp* R2, unsigned ldR2, scomp* M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      scomp d = C_ZERO;
      for (unsigned l=MAX(i,j); l<p; ++l)
        d += R1[i+l*ldR1] * conj(R2[j+l*ldR2]);
      M[i+j*m] = d;
    }
  }
}
inline void utrmmh(unsigned m, unsigned p, unsigned n, dcomp* R1,
                   unsigned ldR1, dcomp* R2, unsigned ldR2, dcomp* M)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<m; ++i) {
      dcomp d = Z_ZERO;
      for (unsigned l=MAX(i,j); l<p; ++l)
        d += R1[i+l*ldR1] * conj(R2[j+l*ldR2]);
      M[i+j*m] = d;
    }
  }
}

// C += d U A, where U is upper triangular in packed storage
template<class T>
inline void putrgemm(T d, unsigned n, T* U, unsigned p, T* A, T* C)
{
  for (unsigned j=0; j<p; ++j) {
    for (unsigned l=0; l<n; ++l) {
      T e = d * A[l+j*n];
      unsigned idl = l*(l+1)/2;
      for (unsigned i=0; i<=l; ++i) C[i+j*n] += e * U[i+idl];
    }
  }
}

// A += d U L, where U is upper triangular in packed storage
// and L is lower triangular in packed storage
inline void putrltrmm(double d, unsigned n, double* U, double* L, double* A)
{
  unsigned* ip = new unsigned[n];
  for (unsigned i=0; i<n; i++) ip[i] = (unsigned) L[((2*n-i+1)*i)/2];

  for (unsigned j=0; j<n; j++) {
    const unsigned ipj = ip[j];
    const unsigned idUj = (ipj*(ipj+1))/2;
    unsigned idL = ((2*n-j+1)*j)/2;
    for (unsigned i=0; i<=ipj; i++) A[i]+=d*U[idUj+i];
    for (unsigned k=j+1; k<n; k++) {
      const unsigned ipk = ip[k];
      const double t = d*L[++idL];
      const unsigned idUk = (ipk*(ipk+1))/2;
      for (unsigned i=0; i<=ipk; i++) A[i]+=U[idUk+i]*t;
    }
    A+=n;
  }
  delete [] ip;
}

inline void putrltrmm(float d, unsigned n, float* U, float* L, float* A)
{
  unsigned* ip = new unsigned[n];
  for (unsigned i=0; i<n; i++) ip[i] = (unsigned) L[((2*n-i+1)*i)/2];

  for (unsigned j=0; j<n; j++) {
    const unsigned ipj = ip[j];
    const unsigned idUj = (ipj*(ipj+1))/2;
    unsigned idL = ((2*n-j+1)*j)/2;
    for (unsigned i=0; i<=ipj; i++)
      A[i]+=d*U[idUj+i];
    for (unsigned k=j+1; k<n; k++) {
      const unsigned ipk = ip[k];
      const float t = d*L[++idL];
      const unsigned idUk = (ipk*(ipk+1))/2;
      for (unsigned i=0; i<=ipk; i++)
        A[i]+=U[idUk+i]*t;
    }
    A+=n;
  }
  delete [] ip;
}

inline void putrltrmm(dcomp d, unsigned n, dcomp* U, dcomp* L, dcomp* A)
{
  unsigned* ip = new unsigned[n];
  for (unsigned i=0; i<n; i++) ip[i]= (unsigned) abs(L[((2*n-i+1)*i)/2]);

  for (unsigned j=0; j<n; j++) {
    const unsigned ipj = ip[j];
    const unsigned idUj = (ipj*(ipj+1))/2;
    unsigned idL = ((2*n-j+1)*j)/2;
    for (unsigned i=0; i<=ipj; i++) A[i]+=d*U[idUj+i];
    for (unsigned k=j+1; k<n; k++) {
      const unsigned ipk = ip[k];
      const dcomp t = d*L[++idL];
      const unsigned idUk = (ipk*(ipk+1))/2;
      for (unsigned i=0; i<=ipk; i++) A[i]+=U[idUk+i]*t;
    }
    A+=n;
  }
  delete [] ip;
}

inline void putrltrmm(scomp d, unsigned n, scomp* U, scomp* L, scomp* A)
{
  unsigned* ip = new unsigned[n];
  for (unsigned i=0; i<n; i++) ip[i]= (unsigned) abs(L[((2*n-i+1)*i)/2]);

  for (unsigned j=0; j<n; j++) {
    const unsigned ipj = ip[j];
    const unsigned idUj = (ipj*(ipj+1))/2;
    unsigned idL = ((2*n-j+1)*j)/2;
    for (unsigned i=0; i<=ipj; i++) A[i]+=d*U[idUj+i];
    for (unsigned k=j+1; k<n; k++) {
      const unsigned ipk = ip[k];
      const scomp t = d*L[++idL];
      const unsigned idUk = (ipk*(ipk+1))/2;
      for (unsigned i=0; i<=ipk; i++) A[i]+=U[idUk+i]*t;
    }
    A+=n;
  }
  delete [] ip;
}

// A += d U^H U, where U is upper triangular in packed storage
// only the upper triangular part of A is computed
inline void putrmhm(const double d, const unsigned n, 
		    const double* const U, double* A)
{
  for (unsigned i=0; i<n; ++i) {
    const unsigned idi = i*(i+1)/2;
    for (unsigned j=i; j<n; ++j) {
      const unsigned idj = j*(j+1)/2;
      double t(0.0);
      const unsigned mij = MIN(i,j);
      for (unsigned l=0; l<=mij; ++l) t += U[l+idi] * U[l+idj];
      
      A[i+idj] += d * t;
    }
  }
}

inline void putrmhm(const float d, const unsigned n, 
		    const float* const U, float* A)
{
  for (unsigned i=0; i<n; ++i){
    const unsigned idi = i*(i+1)/2;
    for (unsigned j=i; j<n; ++j) {
      const unsigned idj = j*(j+1)/2;
      float t(0.0);
      for (unsigned l=0; l<=i; ++l) {
	t += U[l+idi] * U[l+idj];
      }
      A[i+idj] += d * t;
    }
  }
}

inline void putrmhm(const dcomp d, const unsigned n, 
		    const dcomp* const U, dcomp* A)
{
  for (unsigned i=0; i<n; ++i){
    const unsigned idi = i*(i+1)/2;
    for (unsigned j=i; j<n; ++j) {
      const unsigned idj = j*(j+1)/2;
      dcomp t(0.0);
      for (unsigned l=0; l<=i; ++l) {
	t += conj(U[l+idi]) * U[l+idj];
      }
      A[i+idj] += d * t;
    }
  }
}

inline void putrmhm(const scomp d, const unsigned n, 
		    const scomp* const U, scomp* A)
{
  for (unsigned i=0; i<n; ++i){
    const unsigned idi = i*(i+1)/2;
    for (unsigned j=i; j<n; ++j) {
      const unsigned idj = j*(j+1)/2;
      scomp t(0.0);
      for (unsigned l=0; l<=i; ++l) {
	t += conj(U[l+idi]) * U[l+idj];
      }
      A[i+idj] += d * t;
    }
  }
}

// A += d U U^H, where U is upper triangular in packed storage
// only the upper triangular part of A is computed
inline void putrmmh(double d, unsigned n, double* U, double* A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned l=j; l<n; ++l) {
      unsigned idl = l*(l+1)/2;
      double e = d * U[j+idl];
      for (unsigned i=0; i<=j; ++i) A[i+j*(j+1)/2] += e * U[i+idl];
    }
  }
}

inline void putrmmh(float d, unsigned n, float* U, float* A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned l=j; l<n; ++l) {
      unsigned idl = l*(l+1)/2;
      float e = d * U[j+idl];
      for (unsigned i=0; i<=j; ++i) A[i+j*(j+1)/2] += e * U[i+idl];
    }
  }
}

inline void putrmmh(dcomp d, unsigned n, dcomp* U, dcomp* A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned l=j; l<n; ++l) {
      unsigned idl = l*(l+1)/2;
      dcomp e = d * conj(U[j+idl]);
      for (unsigned i=0; i<=j; ++i) A[i+j*(j+1)/2] += e * U[i+idl];
    }
  }
}

inline void putrmmh(scomp d, unsigned n, scomp* U, scomp* A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned l=j; l<n; ++l) {
      unsigned idl = l*(l+1)/2;
      scomp e = d * conj(U[j+idl]);
      for (unsigned i=0; i<=j; ++i) A[i+j*(j+1)/2] += e * U[i+idl];
    }
  }
}

inline void fill0_ltr(unsigned n, double* A)
{
  for (unsigned j=0; j<n; ++j) {
    *A++ = (double) j;
    for (unsigned i=j+1; i<n; ++i) *A++ = D_ZERO;
  }
}
inline void fill0_ltr(unsigned n, float* A)
{
  for (unsigned j=0; j<n; ++j) {
    *A++ = (float) j;
    for (unsigned i=j+1; i<n; ++i) *A++ = S_ZERO;
  }
}
inline void fill0_ltr(unsigned n, scomp* A)
{
  for (unsigned j=0; j<n; ++j) {
    *A++ = scomp( (float)j, S_ZERO);
    for (unsigned i=j+1; i<n; ++i) *A++ = C_ZERO;
  }
}
inline void fill0_ltr(unsigned n, dcomp* A)
{
  for (unsigned j=0; j<n; ++j) {
    *A++ = dcomp( (double) j, D_ZERO);
    for (unsigned i=j+1; i<n; ++i) *A++ = Z_ZERO;
  }
}

// fill nxn matrix A with identity
inline void fillId(unsigned n, double *A)
{
  for (unsigned j=0; j<n; ++j) {
    unsigned i = 0;
    for (; i<j; ++i) *A++ = D_ZERO;
    *A++ = D_ONE;
    ++i;
    for (; i<n; ++i) *A++ = D_ZERO;
  }
}
inline void fillId(unsigned n, float *A)
{
  for (unsigned j=0; j<n; ++j) {
    unsigned i = 0;
    for (; i<j; ++i) *A++ = S_ZERO;
    *A++ = S_ONE;
    ++i;
    for (; i<n; ++i) *A++ = S_ZERO;
  }
}
inline void fillId(unsigned n, scomp *A)
{
  for (unsigned j=0; j<n; ++j) {
    unsigned i = 0;
    for (; i<j; ++i) *A++ = C_ZERO;
    *A++ = C_ONE;
    ++i;
    for (; i<n; ++i) *A++ = C_ZERO;
  }
}
inline void fillId(unsigned n, dcomp *A)
{
  for (unsigned j=0; j<n; ++j) {
    unsigned i = 0;
    for (; i<j; ++i) *A++ = Z_ZERO;
    *A++ = Z_ONE;
    ++i;
    for (; i<n; ++i) *A++ = Z_ZERO;
  }
}

// fill nxn upper triang. matrix A with identity (packed storage)
inline void fillId_utr(unsigned n, double *A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<j; ++i) *A++ = D_ZERO;
    *A++ = D_ONE;
  }
}
inline void fillId_utr(unsigned n, float *A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<j; ++i) *A++ = S_ZERO;
    *A++ = S_ONE;
  }
}
inline void fillId_utr(unsigned n, scomp *A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<j; ++i) *A++ = C_ZERO;
    *A++ = C_ONE;
  }
}
inline void fillId_utr(unsigned n, dcomp *A)
{
  for (unsigned j=0; j<n; ++j) {
    for (unsigned i=0; i<j; ++i) *A++ = Z_ZERO;
    *A++ = Z_ONE;
  }
}

// fill nxn normalized lower triang. matrix A with identity (packed storage)
inline void fillId_ltr(unsigned n, double *A)
{
  for (unsigned i=0; i<n; ++i) {
    for (unsigned j=0; j<i; ++j) *A++ = D_ZERO;
    // for pivoting, a ltr is assumed to have ones on the diagonal
    *A++ = (double) i;
  }
}
inline void fillId_ltr(unsigned n, float *A)
{
  for (unsigned i=0; i<n; ++i) {
    for (unsigned j=0; j<i; ++j) *A++ = S_ZERO;
    // for pivoting, a ltr is assumed to have ones on the diagonal
    *A++ = (float) i;
  }
}
inline void fillId_ltr(unsigned n, scomp *A)
{
  for (unsigned i=0; i<n; ++i) {
    for (unsigned j=0; j<i; ++j) *A++ = C_ZERO;
    // for pivoting, a ltr is assumed to have ones on the diagonal
    *A++ = scomp( (float) i, S_ZERO);
  }
}
inline void fillId_ltr(unsigned n, dcomp *A)
{
  for (unsigned i=0; i<n; ++i) {
    for (unsigned j=0; j<i; ++j) *A++ = Z_ZERO;
    // for pivoting, a ltr is assumed to have ones on the diagonal
    *A++ = dcomp( (double)i, D_ZERO);
  }
}


}


namespace lapack
{

// general inversion
inline void geinv(const unsigned n, double* const A)
{
  unsigned* const ipiv = new unsigned[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  double* const work = new double[lwork];
  assert(work!=NULL);

  int INFO = blas::getrf(n, A, ipiv);
  assert(INFO==0);

  dgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv(const unsigned n, float* const A)
{
  unsigned* const ipiv = new unsigned[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  float* const work = new float[lwork];
  assert(work!=NULL);

  int INFO = blas::getrf(n, A, ipiv);
  assert(INFO==0);

  sgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv(const unsigned n, scomp* const A)
{
  unsigned* const ipiv = new unsigned[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  scomp* const work = new scomp[lwork];
  assert(work!=NULL);

  int INFO = blas::getrf(n, A, ipiv);
  assert(INFO==0);

  cgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv(const unsigned n, dcomp* const A)
{
  unsigned* const ipiv = new unsigned[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  dcomp* const work = new dcomp[lwork];
  assert(work!=NULL);

  int INFO = blas::getrf(n, A, ipiv);
  assert(INFO==0);

  zgetri_(&n, A, &n, ipiv, work, &lwork, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}



inline void geinv_sym(const unsigned n, double* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  double* const work = new double[lwork];
  assert(work!=NULL);

  int INFO;
  dsptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  dsptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv_sym(const unsigned n, float* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  float* const work = new float[lwork];
  assert(work!=NULL);

  int INFO;
  ssptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  ssptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv_sym(const unsigned n, scomp* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  scomp* const work = new scomp[lwork];
  assert(work!=NULL);

  int INFO;
  csptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  csptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv_sym(const unsigned n, dcomp* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  dcomp* const work = new dcomp[lwork];
  assert(work!=NULL);

  int INFO;
  zsptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  zsptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}



inline void geinv_herm(const unsigned n, double* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  double* const work = new double[lwork];
  assert(work!=NULL);

  int INFO;
  dsptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  dsptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv_herm(const unsigned n, float* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  float* const work = new float[lwork];
  assert(work!=NULL);

  int INFO;
  ssptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  ssptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv_herm(const unsigned n, scomp* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  scomp* const work = new scomp[lwork];
  assert(work!=NULL);

  int INFO;
  chptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  chptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}
inline void geinv_herm(const unsigned n, dcomp* const A)
{
  int* const ipiv = new int[n];
  assert(ipiv!=NULL);

  const unsigned lwork = 4*n;
  dcomp* const work = new dcomp[lwork];
  assert(work!=NULL);

  int INFO;
  zhptrf_(JOB_STR+5, &n, A, ipiv, &INFO);
  assert(INFO==0);

  zhptri_(JOB_STR+5, &n, A, ipiv, work, &INFO);
  assert(INFO==0);

  delete [] work;
  delete [] ipiv;
}

//computation of eigenvalues and (orthonormal) eigenvectors of symmetric or
//hermitian  upper triangular matrix A, A is n by n
//the real eigenvalues are stored in diag
inline int utreigv(const unsigned n, float* A, const unsigned ldA, float *diag,
                   const unsigned nwk, float* wk)
{
  int inf;
  ssyev_(JOB_STR+4, JOB_STR+5, &n, A, &ldA, diag, wk, &nwk, &inf);
  return inf;
}
inline int utreigv(const unsigned n, double* A, const unsigned ldA, double *diag,
                   const unsigned nwk, double* wk)
{
  int inf;
  dsyev_(JOB_STR+4, JOB_STR+5, &n, A, &ldA, diag, wk, &nwk, &inf);
  return inf;
}
inline int utreigv(const unsigned n, scomp* A, const unsigned ldA, scomp *diag,
                   const unsigned nwk, scomp* wk)
{
  int inf;
  float *float_diag = new float[n], *rwk = new float[MAX(1,3*n-2)];
  assert(float_diag!=NULL && rwk!=NULL);

  cheev_(JOB_STR+4, JOB_STR+5, &n, A, &ldA, float_diag, wk, &nwk, rwk, &inf);

  for (unsigned i = 0; i < n; ++i) diag[i] = float_diag[i];
  delete [] float_diag;
  delete [] rwk;

  return inf;
}
inline int utreigv(const unsigned n, dcomp* A, const unsigned ldA, dcomp *diag,
                   const unsigned nwk, dcomp* wk)
{
  int inf;
  double *double_diag = new double[n], *rwk = new double[MAX(1,3*n-2)];
  assert(double_diag!=NULL && rwk!=NULL);

  zheev_(JOB_STR+4, JOB_STR+5, &n, A, &ldA, double_diag, wk, &nwk, rwk, &inf);

  for (unsigned i = 0; i < n; ++i) diag[i] = double_diag[i];
  delete [] double_diag;
  delete [] rwk;

  return inf;
}


// packed triangular factorisation of positive definite matrix
inline int pptrf(const unsigned n, double* A)
{
  int inf;
  dpptrf_(JOB_STR+5, &n, A, &inf);
  return inf;
}
inline int pptrf(const unsigned n, float* A)
{
  int inf;
  spptrf_(JOB_STR+5, &n, A, &inf);
  return inf;
}
inline int pptrf(const unsigned n, scomp* A)
{
  int inf;
  cpptrf_(JOB_STR+5, &n, A, &inf);
  return inf;
}
inline int pptrf(const unsigned n, dcomp* A)
{
  int inf;
  zpptrf_(JOB_STR+5, &n, A, &inf);
  return inf;
}

// packed triangular factorization of symmetric(!) matrix
inline int sptrf(const unsigned n, double* A, int* ipiv)
{
  int inf;
  dsptrf_(JOB_STR+5, &n, A, ipiv, &inf);
  return inf;
}
inline int sptrf(const unsigned n, float* A, int* ipiv)
{
  int inf;
  ssptrf_(JOB_STR+5, &n, A, ipiv, &inf);
  return inf;
}
inline int sptrf(const unsigned n, scomp* A, int* ipiv)
{
  int inf;
  csptrf_(JOB_STR+5, &n, A, ipiv, &inf);
  return inf;
}
inline int sptrf(const unsigned n, dcomp* A, int* ipiv)
{
  int inf;
  zsptrf_(JOB_STR+5, &n, A, ipiv, &inf);
  return inf;
}



// triangular factorization of hermitian matrix
inline int hetrf(const unsigned n, double* A, int* ipiv, double* work,
                 const int lwork)
{
  int inf;
  dsytrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
inline int hetrf(const unsigned n, float* A, int* ipiv, float* work,
                 const int lwork)
{
  int inf;
  ssytrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
inline int hetrf(const unsigned n, scomp* A, int* ipiv, scomp* work,
                 const int lwork)
{
  int inf;
  chetrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
inline int hetrf(const unsigned n, dcomp* A, int* ipiv, dcomp* work,
                 const int lwork)
{
  int inf;
  zhetrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}



/*
// triangular factorization of symmetric(!) matrix
inline int sytrf(const unsigned n, double* A, int* ipiv, double* work,
                 const int lwork)
{
  int inf;
  dsytrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
inline int sytrf(const unsigned n, float* A, int* ipiv, float* work,
                 const int lwork)
{
  int inf;
  ssytrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
inline int sytrf(const unsigned n, scomp* A, int* ipiv, scomp* work,
                 const int lwork)
{
  int inf;
  csytrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
inline int sytrf(const unsigned n, dcomp* A, int* ipiv, dcomp* work,
                 const int lwork)
{
  int inf;
  zsytrf_(JOB_STR+6, &n, A, &n, ipiv, work, &lwork, &inf);
  return inf;
}
*/

// lower triangular solve
inline void ltrs(const unsigned n, double* A,
                 const unsigned p, double* B, const unsigned ldB)
{
  //  dtptrs_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    dtpsv_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
inline void ltrs(const unsigned n, float* A,
                 const unsigned p, float* B, const unsigned ldB)
{
  // stptrs_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    stpsv_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
inline void ltrs(const unsigned n, scomp* A,
                 const unsigned p, scomp* B, const unsigned ldB)
{
  //ctptrs_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ctpsv_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
inline void ltrs(const unsigned n, dcomp* A,
                 const unsigned p, dcomp* B, const unsigned ldB)
{
  //ztptrs_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ztpsv_(JOB_STR+6, JOB_STR, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
// lower triangular transpose solve
inline void ltrhs(const unsigned n, double* A,
                  const unsigned p, double* B, const unsigned ldB)
{
  //  dtptrs_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    dtpsv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
inline void ltrhs(const unsigned n, float* A,
                  const unsigned p, float* B, const unsigned ldB)
{
  //  stptrs_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    stpsv_(JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
inline void ltrhs(const unsigned n, scomp* A,
                  const unsigned p, scomp* B, const unsigned ldB)
{
  //ctptrs_(JOB_STR+6, JOB_STR+7, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ctpsv_(JOB_STR+6, JOB_STR+7, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}
inline void ltrhs(const unsigned n, dcomp* A,
                  const unsigned p, dcomp* B, const unsigned ldB)
{
  //ztptrs_(JOB_STR+6, JOB_STR+7, JOB_STR+5, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ztpsv_(JOB_STR+6, JOB_STR+7, JOB_STR+5, &n, A, B+i*ldB, &N_ONE);
}


// unit upper triangular solve (with L and R stored in one matrix)
// XR=B, R is pxp, B is nxp
inline void utrcs(const unsigned p, const double* LR, const unsigned ldLR,
                  const unsigned n, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}

inline void utrcs(const unsigned p, const float* LR, const unsigned ldLR,
                  const unsigned n, float* X, const unsigned ldX)
{
  strsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
         LR, &ldLR, X, &ldX);
}

inline void utrcs(const unsigned p, const scomp* LR, const unsigned ldLR,
                  const unsigned n, scomp* X, const unsigned ldX)
{
  ctrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &C_ONE,
         LR, &ldLR, X, &ldX);
}

inline void utrcs(const unsigned p, const dcomp* LR, const unsigned ldLR,
                  const unsigned n, dcomp* X, const unsigned ldX)
{
  ztrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &Z_ONE,
         LR, &ldLR, X, &ldX);
}

// unit upper triangular solve (with L and R stored in one matrix)
// RX=B, R is nxn, B is nxp
inline void utlcs(const unsigned n, const float* LR, const unsigned ldLR,
                  const unsigned p, float* X, const unsigned ldX)
{
  strsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
         LR, &ldLR, X, &ldX);
}
inline void utlcs(const unsigned n, const double* LR, const unsigned ldLR,
                  const unsigned p, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}
inline void utlcs(const unsigned n, const scomp* LR, const unsigned ldLR,
                  const unsigned p, scomp* X, const unsigned ldX)
{
  ctrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &C_ONE,
         LR, &ldLR, X, &ldX);
}
inline void utlcs(const unsigned n, const dcomp* LR, const unsigned ldLR,
                  const unsigned p, dcomp* X, const unsigned ldX)
{
  ztrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR+5, &n, &p, &Z_ONE,
         LR, &ldLR, X, &ldX);
}

// unit lower triangular solve (with L and R stored in one matrix)
// XL=B, L is pxp, B is nxp
inline void ltrcs(const unsigned p, const float* LR, const unsigned ldLR,
                  const unsigned n, float* X, const unsigned ldX)
{
  strsm_(JOB_STR+8, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
         LR, &ldLR, X, &ldX);
}
inline void ltrcs(const unsigned p, const double* LR, const unsigned ldLR,
                  const unsigned n, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+8, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}
inline void ltrcs(const unsigned p, const scomp* LR, const unsigned ldLR,
                  const unsigned n, scomp* X, const unsigned ldX)
{
  ctrsm_(JOB_STR+8, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &C_ONE,
         LR, &ldLR, X, &ldX);
}
inline void ltrcs(const unsigned p, const dcomp* LR, const unsigned ldLR,
                  const unsigned n, dcomp* X, const unsigned ldX)
{
  ztrsm_(JOB_STR+8, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &Z_ONE,
         LR, &ldLR, X, &ldX);
}


// unit lower triangular transposed solve (with L and R stored in one matrix)
// XL^T=B, L is pxp, B is nxp
inline void ltrhcs(const unsigned p, const double* LR, const unsigned ldLR,
                   const unsigned n, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+8, JOB_STR+6, JOB_STR+1, JOB_STR+5, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}

// unit lower triangular solve (with L and R stored in one matrix)
// LX=B, L is nxn, B is nxp
inline void ltlcs(const unsigned n, const float* LR, const unsigned ldLR,
                  const unsigned p, float* X, const unsigned ldX)
{
  strsm_(JOB_STR+6, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &S_ONE,
         LR, &ldLR, X, &ldX);
}
inline void ltlcs(const unsigned n, const double* LR, const unsigned ldLR,
                  const unsigned p, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+6, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}
inline void ltlcs(const unsigned n, const scomp* LR, const unsigned ldLR,
                  const unsigned p, scomp* X, const unsigned ldX)
{
  ctrsm_(JOB_STR+6, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &C_ONE,
         LR, &ldLR, X, &ldX);
}
inline void ltlcs(const unsigned n, const dcomp* LR, const unsigned ldLR,
                  const unsigned p, dcomp* X, const unsigned ldX)
{
  ztrsm_(JOB_STR+6, JOB_STR+6, JOB_STR, JOB_STR+5, &n, &p, &Z_ONE,
         LR, &ldLR, X, &ldX);
}

// upper triangular solve
inline void utrs(const unsigned n, const double* A,
                 const unsigned p, double* B, const unsigned ldB)
{
  //  dtptrs_(JOB_STR+5, JOB_STR, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    dtpsv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}
inline void utrs(const unsigned n, const float* A,
                 const unsigned p, float* B, const unsigned ldB)
{
  //stptrs_(JOB_STR+5, JOB_STR, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    stpsv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}
inline void utrs(const unsigned n, const scomp* A,
                 const unsigned p, scomp* B, const unsigned ldB)
{
  //ctptrs_(JOB_STR+5, JOB_STR, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ctpsv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}
inline void utrs(const unsigned n, const dcomp* A,
                 const unsigned p, dcomp* B, const unsigned ldB)
{
  //ztptrs_(JOB_STR+5, JOB_STR, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ztpsv_(JOB_STR+5, JOB_STR, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}



// upper triangluar transpose solve
inline void utrhs(const unsigned n, double* A,
                  const unsigned p, double* B, const unsigned ldB)
{
  //dtptrs_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    dtpsv_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}
inline void utrhs(const unsigned n, float* A,
                  const unsigned p, float* B, const unsigned ldB)
{
  //stptrs_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    stpsv_(JOB_STR+5, JOB_STR+1, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}
inline void utrhs(const unsigned n, scomp* A,
                  const unsigned p, scomp* B, const unsigned ldB)
{
  //  ctptrs_(JOB_STR+5, JOB_STR+7, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ctpsv_(JOB_STR+5, JOB_STR+7, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}
inline void utrhs(const unsigned n, dcomp* A,
                  const unsigned p, dcomp* B, const unsigned ldB)
{
  //ztptrs_(JOB_STR+5, JOB_STR+7, JOB_STR, &n, &p, A, B, &ldB, &inf);
  for (unsigned i=0; i<p; ++i)
    ztpsv_(JOB_STR+5, JOB_STR+7, JOB_STR, &n, A, B+i*ldB, &N_ONE);
}



// upper triangular solve (with L and R stored in one matrix)
// XR=B, R is pxp, B is nxp
inline void utrs(const unsigned p, const double* LR, const unsigned ldLR,
                  const unsigned n, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}

inline void utrs(const unsigned p, const float* LR, const unsigned ldLR,
                  const unsigned n, float* X, const unsigned ldX)
{
  strsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &S_ONE,
         LR, &ldLR, X, &ldX);
}

inline void utrs(const unsigned p, const scomp* LR, const unsigned ldLR,
                  const unsigned n, scomp* X, const unsigned ldX)
{
  ctrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &C_ONE,
         LR, &ldLR, X, &ldX);
}

inline void utrs(const unsigned p, const dcomp* LR, const unsigned ldLR,
                  const unsigned n, dcomp* X, const unsigned ldX)
{
  ztrsm_(JOB_STR+8, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &Z_ONE,
         LR, &ldLR, X, &ldX);
}


// upper triangular solve (with L and R stored in one matrix)
// RX=B, R is nxn, B is nxp
inline void utls(const unsigned n, const float* LR, const unsigned ldLR,
                  const unsigned p, float* X, const unsigned ldX)
{
  strsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &S_ONE,
         LR, &ldLR, X, &ldX);
}
inline void utls(const unsigned n, const double* LR, const unsigned ldLR,
                  const unsigned p, double* X, const unsigned ldX)
{
  dtrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &D_ONE,
         LR, &ldLR, X, &ldX);
}
inline void utls(const unsigned n, const scomp* LR, const unsigned ldLR,
                  const unsigned p, scomp* X, const unsigned ldX)
{
  ctrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &C_ONE,
         LR, &ldLR, X, &ldX);
}
inline void utls(const unsigned n, const dcomp* LR, const unsigned ldLR,
                  const unsigned p, dcomp* X, const unsigned ldX)
{
  ztrsm_(JOB_STR+6, JOB_STR+5, JOB_STR, JOB_STR, &n, &p, &Z_ONE,
         LR, &ldLR, X, &ldX);
}
}




#endif
