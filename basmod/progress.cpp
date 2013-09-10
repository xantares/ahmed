/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <iostream>
#include <sstream>
#include <string>

void progressbar(std::ostream& os, std::string s, unsigned i, unsigned n,
                 unsigned nast, bool perc)
{
  const unsigned k = (i+1) * nast;
  std::stringstream sstr(std::stringstream::in | std::stringstream::out);
  std::string mystr;
  if (i==0 || (100*(i+1))/n > (100*i)/n) {
    unsigned j, l = k/n;
    sstr << (char) 13 << s << '[';
    for (j=0; j<l; ++j) sstr << '#';
    for (; j<nast; ++j) sstr << ' ';
    sstr << ']';
    if (perc) {
      float per = (float)((i+1) * 100.0)/n;
      unsigned nz;
      if (per<10.0) nz = 1;
      else if (per<100.0) nz = 2;
      else nz = 3;
      for (j=nz; j<4; ++j) sstr << ' ';
      sstr << (unsigned) per << '%';
    }
    mystr = sstr.str();
    os << mystr << std::flush;
  }
}

