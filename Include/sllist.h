/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef SLL_H
#define SLL_H

#include<stdlib.h>
#include<assert.h>
#include<iostream>

// Template class for a single concatenated list

template<class T> class sllist
{

public:
  // class for a single item in the list
  struct item {
    T data;
    item* next;
    item(item* n, const T& e) : data(e), next(n) {}
  };

  // class for an iterator for the previous list
  class iterator
  {
    friend class sllist<T>;

    item *pitem;                                // currently selected item
    void eolerr() {
      std::cerr << "dllist: end of list" << std::endl;
      exit(1);
    }

  public:
    // constructor and destructor
    iterator() : pitem(NULL) {}
    iterator(item *item) : pitem(item) {}
    iterator(const iterator &iter) : pitem(iter.pitem) {}

    // get status (end-of-list)
    bool eol() const {
      return (pitem == NULL);
    }

    // increment
    iterator &next() {
      return ++*this;
    }

    // successor
    iterator &succ() {
      iterator tmp(*this);
      return ++tmp;
    }

    // access current element
    T &data() {
      assert(pitem != NULL);
      return pitem->_data;
    }
    const T &data() const {
      assert(pitem != NULL);
      return pitem->_data;
    }

    // copy
    iterator &operator=(const iterator& it) {
      pitem = it.pitem;
      return *this;
    }

    // access items
    T &operator()() {
      if (pitem == NULL) eolerr();
      return pitem->data;
    }
    const T &operator()() const {
      if (pitem == NULL) eolerr();
      return pitem->data;
    }

    // iterate (prefix/postfix)
    iterator &operator++() {
      if (pitem) pitem = pitem->next;
      return *this;
    }
    iterator operator++(int) {
      iterator tmp(*this);
      ++*this;
      return tmp;
    }

    // dereference
    T &operator*() {
      if (pitem == NULL) eolerr();
      return pitem->data;
    }
    T *operator->() {
      if (pitem == NULL) eolerr();
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
  // first and last element of list
  item *first_;
  item *last_;

  // size of the list aka number of stored items
  unsigned size_;

public:
  // constructors and destructors
  sllist() : first_(NULL), last_(NULL), size_(0) {}
  sllist(const sllist<T>&);

  ~sllist() {
    clear();
  }

  // access member variables
  iterator begin() {
    return first_;
  }
  const iterator begin() const {
    return first_;
  }

  // same as above but access real data
  T &front() {
    assert(first_ != NULL);
    return first_->data;
  }
  const T &front() const {
    assert(first_ != NULL);
    return first_->data;
  }

  unsigned size() const {
    return size_;
  }
  bool empty() const {
    return (size_==0);
  }

  // add a new object to the beginning/end of list
  void push_front(const T&);
  void push_front_unique(const T& x);
  void push_back(const T&);

  // insert x after it
  iterator insert(iterator it, const T& x) {
    item *p = it.pitem;
    assert(p != NULL && size_ > 0);
    item *tmp = new item(p->next, x);
    it.pitem = p->next = tmp;
    if (it == last_) last_ = tmp;
    ++size_;
    return it;
  }

  // erases the first element from a list
  void pop_front();

  // remove next element from list
  void remove(iterator &it) {
    item *tmp = it.pitem, *nitem = tmp->next;
    // assert(nitem != NULL);
    tmp->next = nitem->next;
    if (nitem == last_) last_ = tmp;
    delete nitem;
    --size_;
  }

  // remove all elements from list
  void clear();

  // look if element is in list
  bool is_in(const T& x) const {
    item *tmp = first_;

    while (tmp) {
      if (tmp->data == x) return true;
      tmp = tmp->next;
    }

    return false;
  }

};


template <class T> sllist<T>::sllist(const sllist<T> &list)
{
  if (list.first_ != NULL) {
    item *tmp = list.first_, // tmp ... run through src list
                *nitem; // item in new list

    nitem = first_ = new item(NULL, tmp->data);
    tmp = tmp->next;

    while (tmp != NULL) {
      last_ = nitem;
      nitem->next = new item(NULL, tmp->data);
      nitem = nitem->next;
      tmp = tmp->next;
    }

  } else first_ = last_ = NULL;

  size_ = list.size_;
}

/*
template <class T> sllist<T>::sllist(const sllist<T> &list)
{
  item *tmp = list.first_, *oitem, *nitem;

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
*/

template<class T> void sllist<T>::push_front(const T& x)
{
  item *tmp = new item(first_, x);
  first_ = tmp;

  if (last_ == NULL) last_ = tmp;
  ++size_;
}

template<class T> void sllist<T>::push_front_unique(const T& x)
{
  if (!first_ || x < first_->data) push_front(x);
  else {
    item *p = first_;
    while (p->next && x>p->next->data) p = p->next;
    if (x>p->data && (p->next==NULL || x!=p->next->data)) {
      item *tmp = new item(p->next, x);
      p->next = tmp;
	  
      if (last_ == p) last_ = tmp;
      ++size_;
    }
  }
}

template<class T> void sllist<T>::push_back(const T& x)
{
  if (size_ != 0) {
    last_->next = new item (NULL, x);
    last_ = last_->next;
  } else {
    first_ = last_ = new item (NULL, x);
  }

  ++size_;
}

template<class T> void sllist<T>::pop_front()
{
  assert(size_ > 0);
  item* tmp = first_->next;
  delete first_;
  first_ = tmp;
  --size_;
  if (size_ == 0) last_ = NULL;
}

template<class T> void sllist<T>::clear()
{
  item *next = first_, *tmp;
  while ((tmp=next)!= NULL) {
    next = tmp->next;
    delete tmp;
  }

  first_ = NULL;
  last_ = NULL;
  size_ = 0;
}


#endif   // SLL_H
