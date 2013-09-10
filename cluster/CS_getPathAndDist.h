/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef GET_PATH_AND_DIST_H
#define GET_PATH_AND_DIST_H

unsigned getPathAndDist (unsigned s, unsigned t, int* tree,
   unsigned nrows, unsigned *row_ptr, unsigned *col_idx);

#endif
