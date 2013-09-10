/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef QUICKSORT_H
#define QUICKSORT_H
#include "basmod.h"

template <class T> static
unsigned partition_(T* array, unsigned beg, unsigned end)
{
  unsigned i = beg+1;
  unsigned j = end-1;
  T m = array[beg];

  for (;;) {
    while ((i<end) && (array[i] < m)) i++;
    while ((j>beg) && !(array[j] < m)) j--;

    if (i >= j) break;
    swap(array[i], array[j]);
  }

  swap(array[beg], array[j]);
  return j;
}



template <class T> static
unsigned partition_(T* array, unsigned beg, unsigned end,
                    bool (*less)(T, T))
{
  unsigned i = beg+1;
  unsigned j = end-1;
  T m = array[beg];

  for (;;) {
    while ((i<end) && less(array[i], m)) i++;
    while ((j>beg) && !less(array[j], m)) j--;

    if (i >= j) break;
    swap(array[i], array[j]);
  }
  swap(array[beg], array[j]);
  return j;
}


template <class T> void quickSort(T* array, unsigned beg, unsigned end)
{
  if (beg < end) {
    unsigned p = partition_(array, beg, end);
    quickSort(array, beg, p);
    quickSort(array, p+1, end);
  }
}

template <class T> void quickSort(T* array, unsigned beg, unsigned end,
                                  bool (*less)(T, T))
{
  if (beg < end) {
    unsigned p = partition_(array, beg, end, less);
    quickSort(array, beg, p, less);
    quickSort(array, p+1, end, less);
  }
}


/** version of partition_ that additionaly updates the permutation vector */
template <class T>
unsigned partition_(T* array, unsigned beg, unsigned end, unsigned *perm)
{
  unsigned i = beg+1;
  unsigned j = end-1;
  T m = array[beg];

  for (;;) {
    while ((i<end) && (array[i] < m)) i++;
    while ((j>beg) && !(array[j] < m)) j--;

    if (i >= j) break;
    swap(array[i], array[j]);
    swap(perm[i], perm[j]);
  }

  swap(array[beg], array[j]);
  swap(perm[beg], perm[j]);
  return j;
}

/** version of quickSort that stores the permutation */
template <class T>
void quickSort(T* array, unsigned beg, unsigned end, unsigned* perm)
{
  if (beg < end) {
    unsigned p = partition_(array, beg, end, perm);
    quickSort(array, beg, p, perm);
    quickSort(array, p+1, end, perm);
  }
}

#endif

