/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "blcluster.h"
#include "mblock.h"
#include "parallel.h"

// psout all blocks associated with this node


extern void psHeader(std::ofstream& os, unsigned N);

template<class T> static
void psoutH_MPI_(std::ofstream& os, blcluster** blList,
		 unsigned* seq_part, unsigned N, mblock<T>** A)
{
  unsigned rank = COMM_AHMED.Get_rank();

  for (unsigned i=seq_part[rank]; i<seq_part[rank+1]; ++i) {
    blcluster* bl = (blcluster*)blList[i];
    const unsigned idx = i - seq_part[rank];
    A[idx]->psout(os, N, bl->getb1(), bl->getb2(), false, 0);
  }
}

void psoutH_MPI(std::ofstream& os, blcluster** blList,
		unsigned* seq_part, unsigned N, mblock<double>** A)
{
  psHeader(os, N);
  psoutH_MPI_(os, blList, seq_part, N, A);
  os << "showpage" << std::endl;
}

