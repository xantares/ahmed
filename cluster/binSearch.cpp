/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "binSearch.h"

bool binSearch(unsigned b, unsigned e, unsigned* arr, unsigned k, unsigned &p)
{
  if (k < arr[b] || k > arr[e-1]) return false;
  unsigned m = (b+e) / 2;
  if (arr[m] == k) {
  	p = m;
  	return true;
  }
  if (b == e-1) return false;
  if (arr[m] < k) return binSearch (m, e, arr, k, p);
  else return binSearch (b, m, arr, k, p);
}