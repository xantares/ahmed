/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#ifndef WIN32

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

double cputime(double x)
{
  struct rusage rsrc;
  double usr, sys;

  if (getrusage(RUSAGE_SELF, &rsrc) == -1) {
    perror("times");
    exit(1);
  }

  usr = rsrc.ru_utime.tv_sec + 1.0e-6 * rsrc.ru_utime.tv_usec;
  sys = rsrc.ru_stime.tv_sec + 1.0e-6 * rsrc.ru_stime.tv_usec;

  return usr + sys - x;
}

#else

double cputime(double x)
{
  return 0.0;
}

#endif
