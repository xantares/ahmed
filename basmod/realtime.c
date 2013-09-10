/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <stdio.h>
#include <stdlib.h>

#ifndef WIN32
  #include <sys/time.h>
  #include <sys/resource.h>
#else
  #include <time.h>
#endif


double realtime(double t0)
{

#ifndef WIN32
  struct timeval curtim;

  if (gettimeofday(&curtim, NULL) == -1) {
    perror ("gettimeofday");
    exit (1);
  }
  
  return ((double)curtim.tv_sec) + (((double)curtim.tv_usec) * 1.0e-6) - t0;

#else  // for WIN32

  time_t curtime;
  time(&curtime);
  return difftime(curtime, 0) - t0;

#endif

}
