/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "parallel.h"
#include "blcluster.h"
#include "H.h"

#ifndef NDEBUG
#define DEBUG
#endif

// H +=HH
template<class T>
void mltaGeHGeH_ND_(unsigned begp, unsigned p, T alpha, blcluster* rootA,
                 mblock<T>** A, blcluster* rootB, mblock<T>** B,
                 blcluster* &rootTmp, mblock<T>** &Tmp,
                 double eps, unsigned rankmax)
{
#ifdef DEBUG
  unsigned rank = COMM_AHMED.Get_rank();
  std::cout<< "\n" << rank << "Starting Multa: begp(" << begp << "), p(" << p
           << ")" << std::flush;
#endif
  unsigned nsons = rootA->getncs();
  if (p==1 || nsons<2) {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) start Multa" << std::flush;
#endif
    mltaGeHGeH(alpha, rootA, A, rootB, B, rootTmp, Tmp, eps, rankmax);
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) beendet Multa"
              << std::flush;
#endif
  } else {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t!(p==1 || rootU->isleaf()) Multa"
              << std::flush;
#endif
    unsigned rank = COMM_AHMED.Get_rank();
    unsigned m = rootA->getnrs();
    unsigned n = rootB->getncs();

    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < m; i++) {
	blcluster* rootAi0 = rootA->getson(i, 0);
        for (unsigned j = 0; j < n; ++j) {
          blcluster* rootB0j = rootB->getson(0, j);
          blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i, j);
          mltaGeHGeH_ND_(begp, p/2, alpha, rootAi0, A, rootB0j, B, rootTij,
                      Tmp, eps, rankmax);
        }
      }
      if (rank==begp) {
	if (nsons==3) {
	  for (unsigned i = 0; i < m; i++) {
	    blcluster* rootAi2 = rootA->getson(i, 2);
	    for (unsigned j = 0; j < n; ++j) {
	      blcluster* rootB2j = rootB->getson(2, j);
	      blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i, j);
	      mltaGeHGeH(alpha, rootAi2, A, rootB2j, B, rootTij, Tmp, eps,
		      rankmax);
	    }
	  }
	}
#ifdef DEBUG
        std::cout<< "\n" << rank << "Waiting for Data... Multa"
                 << std::flush;
#endif
        send_and_addH(rootTmp, Tmp, begp+p/2, begp, eps, rankmax);
      }
    } else {
      for (unsigned i=0; i<m; i++) {
	blcluster* rootAi1 = rootA->getson(i, 1);
        for (unsigned j=0; j<n; ++j) {
          blcluster* rootB1j = rootB->getson(1, j);
          blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i, j);
          setGeHzero(rootTij, Tmp);
          mltaGeHGeH_ND_(begp+p/2, p/2, alpha, rootAi1, A, rootB1j, B,
                      rootTij, Tmp, eps, rankmax);
        }
      }

      if (rank==begp+p/2) {
#ifdef DEBUG
        std::cout<< "\n" << rank << "Sending Data... Multa" << std::flush;
#endif
        send_and_addH(rootTmp, Tmp, begp+p/2, begp, eps, rankmax);
      }
    }
#ifdef DEBUG
    std::cout << "\n" << rank << "\t!(p==1 || rootU->isleaf()) beendet Multa"
              << std::flush;
#endif
  }
}


template<class T>
void GeHUtH_solve_ND_(unsigned begp, unsigned p, blcluster* rootU, mblock<T>** U,
                    blcluster* rootA, mblock<T>** A, mblock<T>** L,
                    double eps, unsigned rankmax)
{
#ifdef DEBUG
  unsigned rank = COMM_AHMED.Get_rank();
  std::cout<< "\n" << rank << "Starting Backwds: begp(" << begp << "), p(" << p
           << ")" << std::flush;
#endif
  if (p==1 || rootU->isleaf()) {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) Backwds"
              << std::flush;
#endif
    GeHUtH_solve(rootU, U, rootA, A, L, eps, rankmax);
  } else {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t!(p==1 || rootU->isleaf()) Backwds"
              << std::flush;
#endif
    unsigned m = rootA->getnrs();
    unsigned rank = COMM_AHMED.Get_rank();
    unsigned nsons = rootU->getnrs();

    if (begp<=rank && rank<begp+p/2) {
#ifdef DEBUG
      std::cout << "\n" << rank << "\tErste Hälfte Backwds" << std::flush;
#endif
      blcluster* sonU00 = rootU->getson(0, 0);

      for (unsigned i=0; i<m; i++) {
        blcluster* sonAi0 = rootA->getson(i, 0);
        GeHUtH_solve_ND_(begp, p/2, sonU00, U, sonAi0, A, L, eps, rankmax);
	if (nsons==3) {
	  blcluster* sonAi2 = rootA->getson(i, 2);
	  blcluster* sonU02 = rootU->getson(0, 2);
	  mltaGeHGeH_ND_(begp, p/2, (T)-1.0, sonAi0, L, sonU02, U, sonAi2, A,
		      eps, rankmax);
	  if (rank==begp) {
#ifdef DEBUG
	    std::cout<< "\n" << rank << "Waiting for Data... Backwds"
		     << std::flush;
#endif
	    blcluster* sonU22 = rootU->getson(2, 2);
	    send_and_addH(sonAi2, A, begp+p/2, begp, eps, rankmax);
	    GeHUtH_solve(sonU22, U, sonAi2, A, L, eps, rankmax);
	  }
	}
      }
#ifdef DEBUG
      std::cout << "\n" << rank << "Erste Hälfte beendet Backwds"
                << std::flush;
#endif
    } else {
#ifdef DEBUG
      std::cout << "\n" << rank << "\tZweite Hälfte Backwds" << std::flush;
#endif
      blcluster* sonU11 = rootU->getson(1, 1);

      for (unsigned i=0; i<m; i++) {
        blcluster* sonAi1 = rootA->getson(i, 1);
        GeHUtH_solve_ND_(begp+p/2, p/2, sonU11, U, sonAi1, A, L, eps, rankmax);
        if (nsons==3) {
	  blcluster* sonAi2 = rootA->getson(i, 2);
	  blcluster* sonU12 = rootU->getson(1, 2);
	  setGeHzero(sonAi2, A);
	  mltaGeHGeH_ND_(begp+p/2, p/2, (T)-1.0, sonAi1, L, sonU12, U, sonAi2, A,
		      eps, rankmax);
#ifdef DEBUG
	  std::cout << "\n" << rank << "Sending Data... Backwds" << std::flush;
#endif
	  send_and_addH(sonAi2, A, begp+p/2, begp, eps, rankmax);
	}
      }
#ifdef DEBUG
      std::cout << "\n" << rank << "Zweite Hälfte beendet Backwds"
                << std::flush;
#endif
    }
  }
}

template<class T>
void LtHGeH_solve_ND_(unsigned begp, unsigned p, blcluster* rootL, mblock<T>** L,
                    blcluster* rootA, mblock<T>** A, mblock<T>** U,
                    double eps, unsigned rankmax)
{
#ifdef DEBUG
  unsigned rank = COMM_AHMED.Get_rank();
  std::cout<< "\n" << rank << "Starting Forwds: begp(" << begp << "), p(" << p
           << ")" << std::flush;
#endif
  if (p==1 || rootL->isleaf()) {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) Forwds"
              << std::flush;
#endif
    LtHGeH_solve(rootL, L, rootA, A, U, eps, rankmax);
  } else {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t!(p==1 || rootU->isleaf()) Forwds"
              << std::flush;
#endif
    unsigned n = rootA->getncs();
    unsigned rank = COMM_AHMED.Get_rank();
    unsigned nsons = rootL->getncs();

    if (begp<=rank && rank<begp+p/2) {
#ifdef DEBUG
      std::cout << "\n" << rank << "\tErste Hälfte Forwds" << std::flush;
#endif
      blcluster* sonL00 = rootL->getson(0, 0);

      for (unsigned j=0; j<n; j++) {
        blcluster *sonA0j = rootA->getson(0, j);
        LtHGeH_solve_ND_(begp, p/2, sonL00, L, sonA0j, A, U, eps, rankmax);
	if (nsons==3) {
	  blcluster* sonL20 = rootL->getson(2, 0);
	  blcluster* sonA2j = rootA->getson(2, j);
	  mltaGeHGeH_ND_(begp, p/2, (T)-1.0, sonL20, L, sonA0j, U, sonA2j, A,
		      eps, rankmax);
	  if (rank==begp) {
#ifdef DEBUG
	    std::cout<< "\n" << rank << "Waiting for Data... Forwds"
		     << std::flush;
#endif
	    blcluster* sonL22 = rootL->getson(2, 2);
	    send_and_addH(sonA2j, A, begp+p/2, begp, eps, rankmax);
	    LtHGeH_solve(sonL22, L, sonA2j, A, U, eps, rankmax);
	  }
	}
      }
#ifdef DEBUG
      std::cout << "\n" << rank << "\tErste Hälfte beendet Forwds" << std::flush;
#endif
    } else {
#ifdef DEBUG
      std::cout << "\n" << rank << "\tZweite Hälfte Forwds" << std::flush;
#endif
      blcluster* sonL11 = rootL->getson(1, 1);

      for (unsigned j=0; j<n; j++) {
        blcluster *sonA1j = rootA->getson(1, j);
        LtHGeH_solve_ND_(begp+p/2, p/2, sonL11, L, sonA1j, A, U, eps, rankmax);
        if (nsons==3) {
	  blcluster* sonL21 = rootL->getson(2, 1);
	  blcluster* sonA2j = rootA->getson(2, j);
	  setGeHzero(sonA2j, A);
	  mltaGeHGeH_ND_(begp+p/2, p/2, (T)-1.0, sonL21, L, sonA1j, U, sonA2j, A,
		      eps, rankmax);
#ifdef DEBUG
	  std::cout << "\n" << rank << "Sending Data... Forwds" << std::flush;
#endif
	  send_and_addH(sonA2j, A, begp+p/2, begp, eps, rankmax);
	}
      }
#ifdef DEBUG
      std::cout << "\n" << rank << "\tZweite Hälfte beendet Forwds" << std::flush;
#endif

    }
  }
}

template<class T>
bool HLU_ND_(unsigned begp, unsigned p, blcluster* root, mblock<T>** A,
             mblock<T>** L, mblock<T>** U, double eps, unsigned rankmax)
{
  unsigned rank = COMM_AHMED.Get_rank();

#ifdef DEBUG
  std::cout<< "\n" << rank << "Decompose: begp(" << begp << "), p(" << p
           << "), rank(" << rank << ")" << std::flush;
#endif
  if (p==1 || root->isleaf()) {
    int inf = HLU(root, A, L, U, eps, rankmax);
    if (!inf) {
      std::cout << "\n" << rank << " HLU_ND Fehler" << std::flush;
      return false;
    }
  } else {
    unsigned nsons = root->getncs();
    if (begp<=rank && rank<begp+p/2) {
#ifdef DEBUG
      std::cout<< "\n" << rank << "\tErste Hälfte" << std::flush;
#endif
      blcluster *son00 = root->getson(0, 0);
#ifdef DEBUG
      std::cout<< "\n" << rank << "\t\tDecompose 1: begp(" << begp
               << "), p(" << p/2 << "), rank(" << rank << ")" << std::flush;
#endif
      HLU_ND_(begp, p/2, son00, A, L, U, eps, rankmax);
#ifdef DEBUG
      std::cout<< "\n" << rank << "*********************************"
               << "\n" << rank << "Decompose 1 erledigt" <<  std::endl;
#endif
      if (nsons==3) {
	blcluster* son02 = root->getson(0, 2), *son20 = root->getson(2, 0);
	blcluster* son22 = root->getson(2, 2);
	GeHUtH_solve_ND_(begp, p/2, son00, U, son20, A, L, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Backwds 1 erledigt" << std::flush;
#endif
	LtHGeH_solve_ND_(begp, p/2, son00, L, son02, A, U, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Forwds 1 erledigt" << std::flush;
#endif

	mltaGeHGeH_ND_(begp, p/2, (T)-1.0, son20, L, son02, U, son22, A, eps, rankmax);
	if (rank==begp) {
	  send_and_addH(son22, A, begp+p/2, begp, eps, rankmax);
	  if (!HLU(son22, A, L, U, eps, rankmax)) {
	    std::cout << "\n" << rank << " HLU_ND Fehler" << std::flush;
	    exit(1);
	  }
	}
      }
    } else {
#ifdef DEBUG
      std::cout << "\n" << rank << "\tZweite Hälfte" << std::flush;
#endif
      blcluster *son11 = root->getson(1, 1);
#ifdef DEBUG
      std::cout<< "\n" << rank << "\t\tDecompose 2: begp(" << begp+p/2
               << "), p(" << p/2 << "), rank(" << rank << ")" << std::flush;
#endif
      HLU_ND_(begp+p/2, p/2, son11, A, L, U, eps, rankmax);
#ifdef DEBUG
      std::cout<< "\n" << rank << "**********************************"
               << "\n" << rank << "Decompose 2 erledigt"
               << std::endl;
#endif
      if (nsons==3) {
	blcluster* son12 = root->getson(1, 2), *son21 = root->getson(2, 1);
	blcluster* son22 = root->getson(2, 2);
	GeHUtH_solve_ND_(begp+p/2, p/2, son11, U, son21, A, L, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Backwds 2 erledigt"
		 << std::flush;
#endif
	LtHGeH_solve_ND_(begp+p/2, p/2, son11, L, son12, A, U, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Forwds 2 erledigt"
		 << std::flush;
#endif
	setGeHzero(son22, A);

	mltaGeHGeH_ND_(begp+p/2, p/2, (T)-1.0, son21, L, son12, U, son22, A,
		    eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Multa 2 erledigt" << std::flush;
	std::cout<< "\n" << rank << "Sending Data..." << std::flush;
#endif
	send_and_addH(son22, A, begp+p/2, begp, eps, rankmax);
      }
#ifdef DEBUG
      std::cout << "\n" << rank << "\tZweite Hälfte beendet" << std::flush;
#endif
    }
  }
  return true;
}

bool HLU_ND(unsigned p, blcluster* root,
            mblock<double>** A, mblock<double>** L, mblock<double>** U,
            double eps, unsigned rankmax)
{
  return HLU_ND_(0, p, root, A, L, U, eps, rankmax);
}

bool HLU_ND(unsigned p,blcluster* root,
            mblock<float>** A, mblock<float>** L, mblock<float>** U,
            double eps, unsigned rankmax)
{
  return HLU_ND_(0, p, root, A, L, U, eps, rankmax);
}

bool HLU_ND(unsigned p, blcluster* root,
            mblock<dcomp>** A, mblock<dcomp>** L, mblock<dcomp>** U,
            double eps, unsigned rankmax)
{
  return HLU_ND_(0, p, root, A, L, U, eps, rankmax);
}

bool HLU_ND(unsigned p, blcluster* root,
            mblock<scomp>** A, mblock<scomp>** L, mblock<scomp>** U,
            double eps, unsigned rankmax)
{
  return HLU_ND_(0, p, root, A, L, U, eps, rankmax);
}

