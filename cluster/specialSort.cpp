/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "basmod.h"
#include "specialSort.h"

unsigned part (unsigned* arr1, unsigned beg, unsigned end, unsigned *arr2)
{
  unsigned i = beg+1;
  unsigned j = end-1;
  unsigned m = arr1[beg];
  unsigned n = arr2[beg];

  for (;;) {
    while ((i<end) && ((arr1[i] < m) || (arr1[i] == m && arr2[i] < n))) i++;
    while ((j>beg) && ((arr1[j] > m) || (arr1[j] == m && arr2[i] > n))) j--;

    if (i >= j) break;
    swap(arr1[i], arr1[j]);
    swap(arr2[i], arr2[j]);
    // if the set is not unique
    if (arr1[i] == arr1[j] && arr2[i] == arr2[j]) i++;
  }

  swap(arr1[beg], arr1[j]);
  swap(arr2[beg], arr2[j]);
  return j;
}

void specialSortRek (unsigned* arr1, unsigned beg, unsigned end, unsigned* arr2)
{
  if (beg < end) {
    unsigned p = part (arr1, beg, end, arr2);
    specialSortRek(arr1, beg, p, arr2);
    specialSortRek(arr1, p+1, end, arr2);
  }
}

void specialSort (unsigned* arr1, unsigned beg, unsigned end, unsigned* arr2)
{
  specialSortRek(arr1, beg, end, arr2);
}

