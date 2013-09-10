/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "sparse.h"

template<class T>
void add_item_(litem<T>*& list, unsigned i, T d)
{
  if (list==0 || i<list->i)
    list = new litem<T>(i, d, list);              // prepend
  else {
    litem<T> *p = list, *last = 0;
    while (p!=0 && i>p->i) {
      last = p;
      p = p->next;
    }

    if (p==0) last->next = new litem<T>(i, d);    // end of list, append
    else if (i==p->i) p->d += d;               // items coincide
    else last->next = new litem<T>(i, d, p);      // put item into list
  }
}

void add_item(litem<double>*& list, unsigned i, double d)
{
  add_item_(list, i, d);
}

void add_item(litem<dcomp>*& list, unsigned i, dcomp d)
{
  add_item_(list, i, d);
}


template<class T>
unsigned len_list_(litem<T>* p)
{
  unsigned len = 0;
  while (p!=0) {
    ++len;
    p = p->next;
  }
  return len;
}

unsigned len_list(litem<double>* p)
{
  return len_list_(p);
}

unsigned len_list(litem<dcomp>* p)
{
  return len_list_(p);
}


template<class T>
void pack_list_(litem<T>* p, T*& pA, unsigned*& piA)
{
  litem<T>* next;

  while (p!=0) {
    *pA++ = p->d;
    *piA++ = p->i;
    next = p->next;
    delete p;
    p = next;
  }
}

void pack_list(litem<double>* p, double*& pA, unsigned*& piA)
{
  pack_list_(p, pA, piA);
}

void pack_list(litem<dcomp>* p, dcomp*& pA, unsigned*& piA)
{
  pack_list_(p, pA, piA);
}
