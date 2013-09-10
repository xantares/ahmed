/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef DLLIST_H
#define DLLIST_H

#include<stdlib.h>
#include<assert.h>
#ifndef NDEBUG
#include<iostream>
#endif

template<class T> class dllist
{
public:
  struct item {
    item *prev, *next;        // pointers to previous and next item
    T data;                   // data stored in the item

    item(item* p, item* n, T x) : prev(p), next(n), data(x) { }
  };

  // class for an iterator for the previous list
  class iterator
  {
    friend class dllist<T>;

    item* pitem;              // currently selected item
#ifndef NDEBUG
    void eolerr() {
      std::cerr << "dllist: end of list" << std::endl;
      exit(1);
    }
#endif

  public:
    // constructor and destructor
    iterator() : pitem(NULL) {}
    iterator(item *p) : pitem(p) {}
    iterator(const iterator &it) : pitem(it.pitem) {}

    // get status (end-of-list), i.e., iterator is in front or behind the list
    bool eol() const {
      return (pitem == NULL);
    }

    // get status "is on list"
    bool isol() const {
      return (pitem != NULL);
    }

    // increment
    iterator &next() {
      return ++*this;
    }

    // decrement
    iterator &prev() {
      return --*this;
    }

    // successor
    iterator &succ() {
      iterator tmp(*this);
      return ++tmp;
    }

    // predecessor
    iterator &pred() {
      iterator tmp(*this);
      return --tmp;
    }

    // access current element
    T &data() {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      return pitem->data;
    }

    const T &data() const {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      return pitem->data;
    }

    // copy
    iterator &operator=(const iterator &it) {
      pitem = it.pitem;
      return *this;
    }

    // access items
    T &operator()() {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      return pitem->data;
    }
    const T &operator()() const {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      return pitem->data;
    }

    // iterate (prefix/postfix)
    iterator &operator++() {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      pitem = pitem->next;
      return *this;
    }
    iterator operator++(int) {
      iterator tmp(*this);
      ++*this;
      return tmp;
    }

    iterator &operator--() {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      pitem = pitem->prev;
      return *this;
    }
    iterator operator--(int) {
      iterator tmp(*this);
      --*this;
      return tmp;
    }

    // dereference
    T &operator*() {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      return pitem->data;
    }

    T *operator->() {
#ifndef NDEBUG
      if (pitem == NULL) eolerr();
#endif
      return &pitem->data;
    }

    // compare
    bool operator==(const iterator &it) const {
      return (pitem == it.pitem);
    }
    bool operator!=(const iterator &it) const {
      return (pitem != it.pitem);
    }
  };


protected:
  // pointers to beginning and end of list
  item *first_, *last_;

  // size of the list aka number of stored items
  unsigned size_;

public:
  dllist() : first_(NULL), last_(NULL), size_(0) { }
  dllist(const dllist<T>&);

  ~dllist() {
    clear();
  }

  // access member variables

  iterator begin() {
    return first_;
  }
  const iterator begin() const {
    return first_;
  }

  iterator end() {
    return last_;
  }
  const iterator end() const {
    return last_;
  }

  // same as above but access real data
  T& front() {
    assert(first_ != NULL);
    return first_->data;
  }
  const T &front() const {
    assert(first_ != NULL);
    return first_->data;
  }

  T& back() {
    assert(last_ != NULL);
    return last_->data;
  }
  const T& back() const {
    assert(last_ != NULL);
    return last_->data;
  }

  unsigned size() const {
    return size_;
  }
  bool empty() const {
    return (size_==0);
  }

  // add a new object to the beginning/end of list
  void push_front(const T&);
  void push_back(const T&);

  // insert x before it
  iterator insert(iterator it, const T& x) {
    item *p = it.pitem;
    assert(p != NULL && size_ > 0);
    item *tmp = new item(p->prev, p, x);
    if (p->prev == NULL) first_ = tmp;
    else p->prev->next = tmp;
    it.pitem = p->prev = tmp;
    ++size_;
    return it;
  }

  // remove one element from a list
  iterator erase(iterator it) {
    item *p = it.pitem;
    assert(p != NULL && size_ > 0);
    item *np = p->next;

    if (p == first_) pop_front();
    else if (p == last_) pop_back();
    else {
      p->prev->next = np;
      p->next->prev = p->prev;
      --size_;
      delete p;
    }
    return np;
  }

  // remove all elements from list
  void clear();

  // erases the first (or last) element from a list
  void pop_front() {
    assert(size_ > 0);
    item* tmp = first_;
    first_ = first_->next;
    if (first_ != NULL) first_->prev = NULL;
    else last_ = NULL;
    delete tmp;
    --size_;
  }
  void pop_back() {
    assert(size_ > 0);
    item* tmp = last_;
    last_ = last_->prev;
    if (last_ != NULL) last_->next = NULL;
    else first_ = NULL;
    delete tmp;
    --size_;
  }

  // check whether element is in list
  bool is_in(const T& x) const {
    item *tmp = first_;

    while (tmp != NULL) {
      if (tmp->data == x) return true;
      tmp = tmp->next;
    }

    return false;
  }
};


template<class T> dllist<T>::dllist(const dllist<T> &list)
{
  item *tmp = list.first_, *oitem;

  if (tmp != NULL) {
    oitem = first_ = new item(NULL, NULL, tmp->data);
    tmp = tmp->next;

    while (tmp != NULL) {
      oitem->next = new item(oitem, NULL, tmp->data);
      oitem = oitem->next;
      tmp = tmp->next;
    }

    last_ = oitem;
  } else first_ = last_ = NULL;

  size_ = list.size_;
}


template<class T> void dllist<T>::push_front(const T& x)
{
  item *tmp = new item(NULL, first_, x);

  if (first_ != NULL) first_->prev = tmp;
  first_ = tmp;

  if (last_ == NULL) last_ = tmp;
  ++size_;
}


template<class T> void dllist<T>::push_back(const T& x)
{
  item *tmp = new item(last_, NULL, x);

  if (last_ != NULL) last_->next = tmp;
  last_ = tmp;

  if (first_ == NULL) first_ = tmp;
  ++size_;
}


template<class T> void dllist<T>::clear()
{
  item *next = first_, *tmp;
  while ((tmp=next)!= NULL) {
    next = tmp->next;
    delete tmp;
  }

  first_ = last_ = NULL;
  size_ = 0;
}


#endif
