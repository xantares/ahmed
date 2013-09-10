/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef BINSEARCH_H_
#define BINSEARCH_H_

/** binary search in the domain [b,e] of the sorted array.
	\param b begin index
	\param e end index
	\param arr sorted array
	\param k key which should be searched for
	\param p position of element with key k
	\returns true iff element with key k is in the array arr, else false
 */
bool binSearch(unsigned b, unsigned e, unsigned* arr, unsigned k, unsigned &p);

#endif /*BINSEARCH_H_*/
