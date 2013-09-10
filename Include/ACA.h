/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef ACA_H
#define ACA_H

#include <cmath>
#include <assert.h>
#include "blas.h"
#include "basmod.h"
#include "mblock.h"
#include "bemcluster.h"
#include "cmplx.h"


// on entry: eps contains wanted accuracy, k max rank
//           U, V of size m*(k+1), n*(k+1)
// on exit:  eps contains achieved accuracy, k the rank

// Adaptive Cross Approximation
inline unsigned nexti(unsigned n, int* Z, unsigned j)
{
  if (Z[j]!=-1) return j;
  else {
    if (j+1<n) return nexti(n, Z, j+1);
    else return nexti(n, Z, 0);
  }
}

// Z[i]==-1 means that either the i-th row has been used as a pivot
//                         or the size of its entries is below 1e-14
// in either case the i-th row of the remainder can be treated as zero
// Z[i]>=0 is the number of successful rank-1 approximations
//         that have been applied to the i-th row


// returns false if no large enough pivot can be found
template<class T, class MATGEN_T> static
bool ACA_row_step(MATGEN_T& MatGen, unsigned b1, unsigned n1, unsigned b2, unsigned n2,
		  unsigned& klast, unsigned& i0, unsigned k, unsigned& no, int* Z, int* S,
		  T* U, T* V, typename num_traits<T>::abs_type& nrmlsk2, typename num_traits<T>::abs_type scale)
{
  typedef typename num_traits<T>::abs_type abs_T;
  unsigned l, j0;
  abs_T nrmu2, nrmv2, absmax;
  T *const pu = U + k*n1, *const pv = V + k*n2;

  // find next non-zero row, start with i0
  do {
    MatGen.cmpbl(b1+i0, 1, b2, n2, pv);
    blas::conj(n2, pv);

    // subtract previously generated outer products from this row
    for (l=0; l<k; ++l) {
      const T d = nconj(U[l*n1+i0]);
      blas::axpy(n2, d, V + l*n2, pv);
    }

    // find maximum
    absmax = 0.0;
    for (l=0; l<n2; ++l) {
      if (S[l]>=0) {
        const abs_T eabs = abs(pv[l]);
        if (eabs>absmax) {
          absmax = eabs;
          j0 = l;
        }
      }
    }

    // pivot too small ? get new row
    if (absmax < 1e-14 * scale) {
      //      std::cout << "pivot too small " << no+1 << ' ' << n1 << std::endl;
      Z[i0] = -1;
      if (++no<n1) i0 = nexti(n1, Z, i0);
    }
  } while (absmax < 1e-14 * scale && no<n1);

  if (absmax < 1e-14 * scale) return false;    // remainder is zero, return

  // else: we have found a reasonable pivot (i0,j0)
  //  std::cout << "i0=" << i0 << ", j0=" << j0 << std::endl;
  abs_T sqrtpiv = sqrt(absmax);
  nrmv2 = 0.0;
  T sca = sqrtpiv/pv[j0];        // scale \tilde{v}_{k+1} -> v_{k+1}
  for (l=0; l<n2; ++l) {
    if (S[l]>=0) {
      const T e = (pv[l] *= sca);
      nrmv2 += abs2(e);
      if (abs(e) > 1e-8*scale*sqrtpiv) ++S[l];  // col l will be succ. approx
    }
  }

  MatGen.cmpbl(b1, n1, b2+j0, 1, pu);

  for (l=0; l<k; ++l) {
    const T d = nconj(V[l*n2 + j0]);
    blas::axpy(n1, d, U + l*n1, pu);
  }

  klast = Z[i0];
  Z[i0] = S[j0] = -1;
  ++no;

  abs_T dsca = 1.0/sqrtpiv;

  absmax = nrmu2 = 0.0;                 // compute norm of u_{k+1}, new i0
  for (l=0; l<n1; ++l) {
    const abs_T eabs = abs(pu[l] *= dsca);
    if (Z[l]>=0 && eabs > 1e-8*scale*sqrtpiv) ++Z[l];
    nrmu2 += eabs*eabs;

    // find new pivot i0
    if (Z[l]>=0 && eabs>absmax) {
      i0 = l;
      absmax = eabs;
    }
  }

  nrmlsk2 = nrmu2 * nrmv2;
  return true;
}


// returns true if successful
template<class T, class MATGEN_T> static
bool ACA_col_step(MATGEN_T& MatGen, unsigned b1, unsigned n1, unsigned b2, unsigned n2,
		  unsigned& j0, unsigned k, unsigned& no, int* Z, int* S, T* U, T* V,
		  typename num_traits<T>::abs_type& nrmlsk2, typename num_traits<T>::abs_type scale)
{
  typedef typename num_traits<T>::abs_type abs_T;
  abs_T nrmu2, nrmv2, absmax;
  T *const pu = U + k*n1, *const pv = V + k*n2;
  unsigned l, i0;

  // find next non-zero row, start with j0
  do {
    MatGen.cmpbl(b1, n1, b2+j0, 1, pu);

    for (l=0; l<k; ++l) {
      const T d = nconj(V[l*n2+j0]);
      blas::axpy(n1, d, U+l*n1, pu);
    }

    absmax = 0.0;
    for (l=0; l<n1; ++l) {
      if (Z[l]>=0) {
        const abs_T eabs = abs(pu[l]);
        if (eabs>absmax) {
          absmax = eabs;
          i0 = l;
        }
      }
    }

    if (absmax < 1e-14*scale) {
      S[j0] = -1;
      if (++no<n2) j0 = nexti(n2, S, j0);
    }
  } while (absmax < 1e-14*scale && no<n2);

  if (absmax < 1e-14*scale) return false;

  abs_T sqrtpiv = sqrt(absmax);
  nrmu2 = 0.0;
  T sca = sqrtpiv/pu[i0];           // scale \tilde{u}_{k+1} -> u_{k+1}
  for (l=0; l<n1; ++l) {
    if (Z[l]>=0) {
      const T e = (pu[l] *= sca);
      nrmu2 += abs2(e);
      if (abs(e) > 1e-8 * scale*sqrtpiv) ++Z[l];
    }
  }

  MatGen.cmpbl(b1+i0, 1, b2, n2, pv);
  blas::conj(n2, pv);

  for (l=0; l<k; ++l) {
    const T d = nconj(U[l*n1+i0]);
    blas::axpy(n2, d, V+l*n2, pv);
  }

  Z[i0] = S[j0] = -1;
  ++no;

  abs_T dsca = 1.0/sqrtpiv;

  absmax = nrmv2 = 0.0;                 // compute norm of v_{k+1}, new j0
  for (l=0; l<n2; ++l) {

    const abs_T eabs = abs(pv[l] *= dsca);
    if (S[l]>=0 && eabs > 1e-8*scale*sqrtpiv) ++S[l];
    nrmv2 += eabs*eabs;
    
    // find new pivot imax
    if (eabs>absmax) {
      j0 = l;
      absmax = eabs;
    }
  }


  nrmlsk2 = nrmu2 * nrmv2;
  return true;
}



template<class T, class MATGEN_T>
bool ACA(MATGEN_T& MatGen, unsigned b1, unsigned n1, unsigned b2, unsigned n2,
         double eps, unsigned kmax, unsigned i0, unsigned& k, T* &U, T* &V)
{
  typedef typename num_traits<T>::abs_type abs_T;
  unsigned l, no = 0, j0, inext, klast;
  // klast is the number of successful steps applied to the last pivotal row
  int kmin = -2, keps = -1;
  abs_T nrmlsk2, nrms2 = 0.0;
  T sum;
  bool ok, stpcrit;

  abs_T scale = MatGen.scale(b1, n1, b2, n2); // set initial scale

  U = new T[(kmax+1)*n1];
  V = new T[(kmax+1)*n2];
  assert(U!=NULL && V!=NULL);

  int *Z = new int[n1], *S = new int[n2];
  assert(Z!=NULL && S!=NULL);

  for (l=0; l<n1; ++l) Z[l] = 0;
  for (l=0; l<n2; ++l) S[l] = 0;

  k = 0;

  do {

    // compute a cross
    inext = i0;
    ok = ACA_row_step(MatGen, b1, n1, b2, n2, klast, inext, k, no, Z, S,
                      U, V, nrmlsk2, scale);

    // std::cout << "no=" << no << ", n1=" << n1 << ", k=" << k << std::endl;
    if (!ok) { // in last step no non-zero row could be found
      nrmlsk2 = 0.0;
      goto SUCC;
    }

    // else: step was successful

    // ------------------------------------------------------------------------
    // check stopping criterion

    sum = (T) 0.0;                            // update nrms2
    for (l=0; l<k; ++l)
      sum += blas::scpr(n1, U+l*n1, U+k*n1) * blas::scpr(n2, V+l*n2, V+k*n2);

    nrms2 += 2.0 * Re(sum) + nrmlsk2;

    stpcrit = (nrmlsk2<eps*eps*nrms2);
    if (stpcrit && keps==-1) keps = klast+1;  // keps is the rank for acc. eps

    //    std::cout << "keps=" << keps << std::endl;

    ++k;

    // adjust scale (estimated entry size of the next remainder)
    scale = sqrt(nrmlsk2/(n1*n2));

    // -----------------------------------------------------------------------
    // find row with smallest number of succ. subtractions

    if (keps>0) {  // use minimal succ row for next pivot

      //      for (unsigned j=0; j<n1; ++j) std::cout << Z[j] << ' ';
      //      std::cout << std::endl;

      for (i0=0; i0<n1 && Z[i0]<0; ++i0);
      if (i0==n1) {
        nrmlsk2 = 0.0;
        goto SUCC;  // last step was the last possible step
      }

      // else: Z[i0]>=0
      kmin = Z[i0];
      for (unsigned i=i0+1; i<n1; ++i)
        if (Z[i]>=0 && Z[i]<kmin) kmin = Z[i0=i];

      if ((unsigned)kmin==k) i0 = inext;
    } else i0 = inext; // use max of last column as new pivot

    //    std::cout << "kmin=" << kmin << ", keps=" << keps << std::endl;
  } while (no<n1 && k<kmax && kmin<keps);

  //  for (unsigned j=0; j<n2; ++j) std::cout << S[j] << ' ';
  //  std::cout << std::endl;

  // -----------------------------------------------------------------------
  // continue with column oriented ACA

  if ((no==n1 || k==kmax) && kmin<keps) goto NSUCC;

  // guarantee kmin>=keps for the columns, too
  do {
    for (j0=0; j0<n2 && S[j0]<0; ++j0);
    if (j0==n2) {
      nrmlsk2 = 0.0;
      goto SUCC;              // last step was the last possible step
    }

    // else: S[j0]>=0
    kmin = S[j0];
    for (unsigned j=j0+1; j<n2; ++j)
      if (S[j]>=0 && S[j]<kmin) kmin = S[j0=j];

    if (kmin>=keps) goto SUCC;   // if all were succ., stop iteration

    // else: use j0 as new pivot
    ok = ACA_col_step(MatGen, b1, n1, b2, n2, j0, k, no, Z, S, U, V, nrmlsk2, scale);

    if (!ok) {
      nrmlsk2 = 0.0;
      goto SUCC;
    }
    ++k;    // next rank-1 approximation

  } while (no<n2 && k<kmax);

NSUCC:
  //std::cout << "NSUCC " << n1 << ", " << n2 << ", k=" << k << std::endl;
  delete [] S;
  delete [] Z;
  return false;

SUCC:
  double epsa = sqrt(nrmlsk2/nrms2);
  //  std::cout << "SUCC: " << epsa << ' ' << eps << std::endl;
  if (epsa>=eps) goto NSUCC;

  delete [] S;
  delete [] Z;
  return true;
}

// ----------------------------------------------------------------------------
// new version of ACA
// 
// klast not required!

template<class T, class MATGEN_T>
bool ACAn(MATGEN_T& MatGen, unsigned b1, unsigned n1, unsigned b2, unsigned n2,
	  double eps, unsigned kmax, unsigned i0, unsigned& k, T* &U, T* &V)
{
  typedef typename num_traits<T>::abs_type abs_T;
  unsigned l, no = 0, j0, inext;
  // klast is the number of successful steps applied to the last pivotal row
  unsigned kmin, keps, k_old;
  abs_T nrmlsk2, nrms2 = 0.0;
  T sum;
  bool ok;

  abs_T scale = MatGen.scale(b1, n1, b2, n2); // set initial scale

  U = new T[(kmax+1)*n1];
  V = new T[(kmax+1)*n2];
  assert(U!=NULL && V!=NULL);

  // global row/column counters
  int *Zgl = new int[n1], *Sgl = new int[n2];
  assert(Zgl!=NULL && Sgl!=NULL);
  for (l=0; l<n1; ++l) Zgl[l] = 0;
  for (l=0; l<n2; ++l) Sgl[l] = 0;

  // counters for each step
  int *Z = new int[n1], *S = new int[n2];
  assert(Z!=NULL && S!=NULL);

  keps = k = 0;

#ifdef CHECK_ACA_ERROR
  std::cout << std::endl << "new block " << n1 << " x " << n2 << std::endl;
#endif

  do {
    for (l=0; l<n1; ++l) Z[l] = 0;
    for (l=0; l<n2; ++l) S[l] = 0;

    // start with a loop of row-steps until criterion is satisfied

    k_old = k;
    
    do {

      // compute a cross
      inext = i0;

      unsigned klast;
      ok = ACA_row_step(MatGen, b1, n1, b2, n2, klast, inext, k, no, Z, S,
			U, V, nrmlsk2, scale);
      
      if (!ok) goto NSUCC;        // in last step no non-zero row could be found
      
      // else: step was successful, check stopping criterion      
      sum = (T) 0.0;                            // update nrms2
      for (l=0; l<k; ++l)
	sum += blas::scpr(n1, U+l*n1, U+k*n1) * blas::scpr(n2, V+l*n2, V+k*n2);

      nrms2 += 2.0 * Re(sum) + nrmlsk2;

      ++k;

      // adjust scale (estimated entry size of the next remainder)
      scale = sqrt(nrmlsk2/(n1*n2));
      
    } while (no<n1 && k<kmax+1 && (nrmlsk2>=eps*eps*nrms2));

#ifdef CHECK_ACA_ERROR
    std::cout << "row loop: no=" << no << " k=" << k
	      << " kmax=" << kmax << std::endl;
#endif

    // a row that has been successfully approximated the maximum number of times
    // is marked as done (-1), all other approx. are added to the global counter
    for (unsigned i=0; i<n1; ++i) {
      if (Z[i]+k_old >= k || Z[i] == -1) Zgl[i] = -1;
      else if (Zgl[i]!=-1) Zgl[i] += Z[i];
    }

    for (unsigned j=0; j<n2; ++j) {
      if (S[j]+k_old >= k || S[j] == -1) Sgl[j] = -1;
      else if (Sgl[j]!=-1) Sgl[j] += S[j];
    }

    // set the number of required successful apprx. steps
    if (k-k_old > keps) keps = k-k_old;
#ifdef CHECK_ACA_ERROR
    std::cout << "keps=" << keps << std::endl;
#endif

    // -----------------------------------------------------------------------
    // find col with smallest number of succ. subtractions

    for (j0=0; j0<n2 && (Sgl[j0]<0 || Sgl[j0]>=(int)keps); ++j0);
    if (j0==n2) {                  // all columns successfully approximated
#ifdef CHECK_ACA_ERROR
      std::cout << "cols ok" << std::endl;
#endif
      goto PREP_ROW;               // check rows
    } else {                       // error estimator is nor reliable
      if (no==n1 || k==kmax+1) goto NSUCC; // no further step can be performed?
    }
    
    // else: perpare a column-oriented loop

    // find minimal k and corresponding j0
    kmin = Sgl[j0];
    for (unsigned j=j0+1; j<n2; ++j)
      if (Sgl[j]>=0 && Sgl[j]<(int)kmin) kmin = Sgl[j0=j];

#ifdef CHECK_ACA_ERROR
    std::cout << "treating cols: kmin=" << kmin << std::endl;
#endif

    for (l=0; l<n1; ++l) Z[l] = 0;
    for (l=0; l<n2; ++l) S[l] = 0;

    // continue with a loop of column-steps until criterion is satisfied

    k_old = k;
  
    do {

      // compute a cross
      unsigned jnext = j0;
      k_old = k;
      ok = ACA_col_step(MatGen, b1, n1, b2, n2, jnext, k, no, Z, S,
			U, V, nrmlsk2, scale);
      
      if (!ok) goto NSUCC;     // in last step no non-zero row could be found
      
      // else: step was successful, check stopping criterion
      
      sum = (T) 0.0;                            // update nrms2
      for (l=0; l<k; ++l)
	sum += blas::scpr(n1, U+l*n1, U+k*n1) * blas::scpr(n2, V+l*n2, V+k*n2);

      nrms2 += 2.0 * Re(sum) + nrmlsk2;

      ++k;
      
      // adjust scale (estimated entry size of the next remainder)
      scale = sqrt(nrmlsk2/(n1*n2));
      
    } while (no<n2 && k<kmax+1 && (nrmlsk2>=eps*eps*nrms2));

#ifdef CHECK_ACA_ERROR
    std::cout << "column loop: no=" << no << " k=" << k
	      << " kmax=" << kmax << std::endl;
#endif

    // update Zgl and Sgl

    for (unsigned i=0; i<n1; ++i) {
      if (Z[i]+k_old >= k || Z[i] == -1) Zgl[i] = -1;
      else if (Zgl[i]!=-1) Zgl[i] += Z[i];
    }

    for (unsigned j=0; j<n2; ++j) {
      if (S[j]+k_old >= k || S[j] == -1) Sgl[j] = -1;
      else if (Sgl[j]!=-1) Sgl[j] += S[j];
    }

    // set the number of required successful apprx. steps
    if (k-k_old > keps) keps = k-k_old;
#ifdef CHECK_ACA_ERROR
    std::cout << "keps=" << keps << std::endl;
#endif

    // -----------------------------------------------------------------------
    // find row with smallest number of succ. subtractions

  PREP_ROW:
    for (i0=0; i0<n1 && (Zgl[i0]<0 || Zgl[i0]>=(int)keps); ++i0);
    if (i0==n1) {         // all rows are approximated
#ifdef CHECK_ACA_ERROR
      std::cout << "rows ok" << std::endl;
#endif
      goto SUCC;
    } else {                 // error estimator not reliable
      if (no==n2 || k==kmax+1) goto NSUCC;
    }

    // else: prepare another row-oriented loop
    kmin = Zgl[i0];
    for (unsigned i=i0+1; i<n1; ++i)
      if (Zgl[i]>=0 && Zgl[i]<(int)kmin) kmin = Zgl[i0=i];

#ifdef CHECK_ACA_ERROR
    std::cout << "treating rows: kmin=" << kmin << std::endl;
#endif

  } while (no<n1 && kmin<keps);

  if (kmin>=keps) goto SUCC;

NSUCC:
  delete [] S;
  delete [] Z;
  delete [] Sgl;
  delete [] Zgl;
  return false;

SUCC:
  double epsa = sqrt(nrmlsk2/nrms2);
#ifdef CHECK_ACA_ERROR
  std::cout << "SUCC: " << epsa << ' ' << eps << std::endl;
#endif
  if (epsa>=eps) goto NSUCC;

  delete [] S;
  delete [] Z;
  delete [] Sgl;
  delete [] Zgl;
  return true;
}

// ----------------------------------------------------------------------------
// simplified ACA procedure

// returns false if no large enough pivot can be found
template<class T, class MATGEN_T> static
bool ACAr_step(MATGEN_T& MatGen, unsigned b1, unsigned n1, unsigned b2, unsigned n2,
	       unsigned& i0, unsigned k, unsigned& no, int* Z, T* U, T* V,
	       typename num_traits<T>::abs_type& nrmlsk2, typename num_traits<T>::abs_type scale)
{
  typedef typename num_traits<T>::abs_type abs_T;
  unsigned l, j0;
  abs_T nrmu2, nrmv2, absmax;
  T *const pu = U + k*n1, *const pv = V + k*n2;

  // find next non-zero row, start with i0
  do {
    MatGen.cmpbl(b1+i0, 1, b2, n2, pv);
    blas::conj(n2, pv);

    // subtract previously generated outer products from this row
    for (l=0; l<k; ++l) {
      const T d = nconj(U[l*n1+i0]);
      blas::axpy(n2, d, V + l*n2, pv);
    }

    // find maximum
    absmax = 0.0;
    for (l=0; l<n2; ++l) {
      const abs_T eabs = abs(pv[l]);
      if (eabs>absmax) {
        absmax = eabs;
        j0 = l;
      }
    }

    // pivot too small ? get new row
    if (absmax < 1e-14 * scale) {
      //      std::cout << "pivot too small " << no+1 << ' ' << n1 << std::endl;
      Z[i0] = -1;
      if (++no<n1) i0 = nexti(n1, Z, i0);
    }
  } while (absmax < 1e-14 * scale && no<n1);

  if (absmax < 1e-14 * scale) return false;    // remainder is zero, return

  // else: we have found a reasonable pivot (i0,j0)
  //  std::cout << "i0=" << i0 << ", j0=" << j0 << std::endl;
  abs_T sqrtpiv = sqrt(absmax);
  nrmv2 = 0.0;
  T sca = sqrtpiv/pv[j0];        // scale \tilde{v}_{k+1} -> v_{k+1}
  for (l=0; l<n2; ++l) {
    const T e = (pv[l] *= sca);
    nrmv2 += abs2(e);
  }

  MatGen.cmpbl(b1, n1, b2+j0, 1, pu);

  for (l=0; l<k; ++l) {
    const T d = nconj(V[l*n2 + j0]);
    blas::axpy(n1, d, U + l*n1, pu);
  }

  Z[i0] = -1;
  ++no;

  abs_T dsca = 1.0/sqrtpiv;

  absmax = nrmu2 = 0.0;                 // compute norm of u_{k+1}, new i0
  for (l=0; l<n1; ++l) {
    const abs_T eabs = abs(pu[l] *= dsca);
    nrmu2 += eabs*eabs;

    // find new pivot i0
    if (Z[l]>=0 && eabs>absmax) {
      i0 = l;
      absmax = eabs;
    }
  }

  nrmlsk2 = nrmu2 * nrmv2;
  return true;
}

template<class T, class MATGEN_T>
bool ACAr(MATGEN_T& MatGen, unsigned b1, unsigned n1, unsigned b2, unsigned n2,
          double eps, unsigned kmax, unsigned i0, unsigned& k, T* &U, T* &V)
{
  typedef typename num_traits<T>::abs_type abs_T;
  unsigned l, no = 0;
  abs_T nrmlsk2, nrms2 = 0.0;
  T sum;
  bool ok;

  abs_T scale = MatGen.scale(b1, n1, b2, n2); // set initial scale

  U = new T[(kmax+1)*n1];
  V = new T[(kmax+1)*n2];
  assert(U!=NULL && V!=NULL);

  int *Z = new int[n1];
  assert(Z!=NULL);

  for (l=0; l<n1; ++l) Z[l] = 0;

  k = 0;

  do {

    // compute a cross
    ok = ACAr_step(MatGen, b1, n1, b2, n2, i0, k, no, Z,
                   U, V, nrmlsk2, scale);

    //    std::cout << "Norm last cross " << sqrt(nrmlsk2) << ", ok=" << (int) ok << std::endl;
    //    std::cout << "no=" << no << ", n1=" << n1 << ", k=" << k << std::endl;
    if (!ok) { // in last step no non-zero row could be found
      nrmlsk2 = 0.0;
      goto SUCC;
    }

    // else: step was successful

    sum = (T) 0.0;                            // update nrms2
    for (l=0; l<k; ++l)
      sum += blas::scpr(n1, U+l*n1, U+k*n1) * blas::scpr(n2, V+l*n2, V+k*n2);

    nrms2 += 2.0 * Re(sum) + nrmlsk2;

    ++k;

    // adjust scale (estimated entry size of the next remainder)
    scale = sqrt(nrmlsk2/(n1*n2));

  } while (no<n2 && k<kmax && nrmlsk2>=eps*eps*nrms2);

SUCC:
  double epsa = sqrt(nrmlsk2/nrms2);
  //  std::cout << "SUCC: " << epsa << ' ' << eps << std::endl;
  if (epsa>=eps) goto NSUCC;

  delete [] Z;
  return true;

NSUCC:
  //std::cout << "NSUCC " << n1 << ", " << n2 << ", k=" << k << std::endl;
  delete [] Z;
  return false;
}

#endif

