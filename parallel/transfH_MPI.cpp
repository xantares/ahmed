/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <mpi.h>
#include "mblock.h"
#include "H.h"

extern MPI::Intracomm COMM_AHMED;

void transfH_MPI(mblock<double>** Ap, mblock<float>** A, unsigned* seq_part)
{
  unsigned rank = COMM_AHMED.Get_rank(),
                  nproc = COMM_AHMED.Get_size();

  unsigned info[8];

  if (rank==0) {
    copyH(seq_part[1]-seq_part[0], Ap, A);

    for (unsigned j=1; j<nproc; ++j) {
      for (unsigned i=seq_part[j]; i<seq_part[j+1]; ++i) {

        COMM_AHMED.Recv(info, 8, MPI::UNSIGNED, j, 7);
        A[i] = new mblock<float>(info[0], info[1]);
        float* tmp = NULL;
        if (info[7]) {
          tmp = new float[info[7]];
          COMM_AHMED.Recv(tmp, info[7], MPI::FLOAT, j, 8);
        }
        A[i]->cpy_mbl(info+2, tmp);
        delete [] tmp;
      }
    }
  } else {
    for (unsigned i=seq_part[rank]; i<seq_part[rank+1]; ++i) {
      mblock<double>* p = Ap[i-seq_part[rank]];
      info[0] = p->getn1();
      info[1] = p->getn2();
      p->get_prop(info+2);
      info[7] = p->nvals();
      COMM_AHMED.Send(info, 8, MPI::UNSIGNED, 0, 7);

      if (info[7]) {
        float* tmp = new float[info[7]];
        for (unsigned i=0; i<info[7]; ++i) tmp[i] = (float) p->getdata()[i];
        COMM_AHMED.Send(tmp, info[7], MPI::FLOAT, 0, 8);
        delete [] tmp;
      }
    }
  }
}
