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
#include "preserveVec2.h"

#ifndef NDEBUG
#define DEBUG
#endif

// H +=HTH
template<class T> static
void mltaGeHhGeH_ND_(unsigned begp, unsigned p, T alpha, blcluster* rootA,
                  mblock<T>** A, blcluster* rootB, mblock<T>** B,
                  blcluster* &rootTmp, mblock<T>** &Tmp,
                  double eps, unsigned rankmax)
{
#ifdef DEBUG
  unsigned rank = COMM_AHMED.Get_rank();
  std::cout << "\n" << rank << "Starting Multa: begp(" << begp << "), p(" << p
            << ")" << std::flush;
#endif
  unsigned nsons = rootA->getnrs();
  if (p==1 || nsons<2) {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) start Multa"
              << std::flush;
#endif
    mltaGeHhGeH(alpha, rootA, A, rootB, B, rootTmp, Tmp, eps, rankmax);
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
    unsigned m = rootA->getncs();
    unsigned n = rootB->getncs();

    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < m; i++) {
	blcluster* rootA0i = rootA->getson(0, i);
        for (unsigned j = 0; j < n; ++j) {
          blcluster* rootB0j = rootB->getson(0, j);
          blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i, j);
          mltaGeHhGeH_ND_(begp, p/2, alpha, rootA0i, A, rootB0j, B, rootTij,
                       Tmp, eps, rankmax);
        }
      }
      if (rank==begp) {
	if (nsons==3) {
	  for (unsigned i = 0; i < m; i++) {
	    blcluster* rootA2i = rootA->getson(2, i);
	    for (unsigned j = 0; j < n; ++j) {
	      blcluster* rootB2j = rootB->getson(2, j);
	      blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i, j);
	      mltaGeHhGeH(alpha, rootA2i, A, rootB2j, B, rootTij, Tmp, eps,
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
	blcluster* rootA1i = rootA->getson(1, i);
        for (unsigned j=0; j<n; ++j) {
          blcluster* rootB1j = rootB->getson(1, j);
          blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,j);
          setGeHzero(rootTij, Tmp);
          mltaGeHhGeH_ND_(begp+p/2, p/2, alpha, rootA1i, A, rootB1j, B,
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
    std::cout<< "\n" << rank << "\t!(p==1 || rootU->isleaf()) beendet Multa"
             << std::flush;
#endif
  }
}

// HSym += HTH
template<class T> static
void mltaGeHhGeH_toHeH_ND_(unsigned begp, unsigned p, T alpha, blcluster* rootA,
                         mblock<T>** A, blcluster* rootB, mblock<T>** B,
                         blcluster* &rootTmp, mblock<T>** &Tmp,
                         double eps, unsigned rankmax)
{
#ifdef DEBUG
  unsigned rank = COMM_AHMED.Get_rank();
  std::cout<< "\n" << rank << "Starting Multa: begp(" << begp << "), p(" << p
           << ")" << std::flush;
#endif
  unsigned nsons = rootA->getnrs();
  if (p==1 || nsons<2) {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) start Multa" << std::flush;
#endif
    mltaGeHhGeH_toHeH(alpha, rootA, A, rootB, B, rootTmp, Tmp, eps, rankmax);
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) beendet Multa" << std::flush;
#endif
  } else {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t!(p==1 || rootU->isleaf()) Multa" << std::flush;
#endif
    unsigned rank = COMM_AHMED.Get_rank();
    unsigned m = rootA->getncs();
    unsigned n = rootB->getncs();

    if (begp<=rank && rank<begp+p/2) {
      for (unsigned i = 0; i < m; i++) {
        blcluster* rootA0i = rootA->getson(0, i);
        blcluster* rootB0i = rootB->getson(0, i);
        blcluster* rootTii = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,i);
        mltaGeHhGeH_toHeH_ND_(begp, p/2, alpha, rootA0i, A, rootB0i, B,
                            rootTii, Tmp, eps, rankmax);
        for (unsigned j = i+1; j < n; ++j) {
          blcluster* rootB0j = rootB->getson(0, j);
          blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,j);
          mltaGeHhGeH_ND_(begp, p/2, alpha, rootA0i, A, rootB0j, B, rootTij,
                       Tmp, eps, rankmax);
        }
      }
      if (rank==begp) {
	if (nsons==3) {
	  for (unsigned i = 0; i < m; i++) {
	    blcluster* rootA2i = rootA->getson(2, i);
	    blcluster* rootB2i = rootB->getson(2, i);
	    blcluster* rootTii = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,i);
	    mltaGeHhGeH_toHeH(alpha, rootA2i, A, rootB2i, B, rootTii, Tmp, eps,
			      rankmax);
	    for (unsigned j = i+1; j < n; ++j) {
	      blcluster* rootB2j = rootB->getson(2, j);
	      blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,j);
	      mltaGeHhGeH(alpha, rootA2i, A, rootB2j, B, rootTij, Tmp, eps,
		       rankmax);
	    }
	  }
	}
#ifdef DEBUG
        std::cout<< "\n" << rank << "Waiting for Data... Multa"
                 << std::flush;
#endif
        send_and_addHSym(rootTmp, Tmp, begp+p/2, begp, eps, rankmax);
      }
    } else {
      for (unsigned i=0; i<m; i++) {
        blcluster* rootA1i = rootA->getson(1, i);
        blcluster* rootB1i = rootB->getson(1, i);
        blcluster* rootTii = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,i);
        setHeHzero(rootTii, Tmp);
        mltaGeHhGeH_toHeH_ND_(begp+p/2, p/2, alpha, rootA1i, A, rootB1i, B,
                            rootTii, Tmp, eps, rankmax);
        for (unsigned j=i+1; j<n; ++j) {
          blcluster* rootB1j = rootB->getson(1, j);
          blcluster* rootTij = rootTmp->isleaf() ? rootTmp : rootTmp->getson(i,j);
          setGeHzero(rootTij, Tmp);
          mltaGeHhGeH_ND_(begp+p/2, p/2, alpha, rootA1i, A, rootB1j, B,
                       rootTij, Tmp, eps, rankmax);
        }
      }

      if (rank==begp+p/2) {
#ifdef DEBUG
        std::cout<< "\n" << rank << "Sending Data... Multa" << std::flush;
#endif
        send_and_addHSym(rootTmp, Tmp, begp+p/2, begp, eps, rankmax);
      }
    }
#ifdef DEBUG
    std::cout<< "\n" << rank << "\t!(p==1 || rootU->isleaf()) beendet Multa"
             << std::flush;
#endif
  }
}

template<class T> static
void LtHGeH_solve_ND_(unsigned begp, unsigned p, blcluster* rootU,
                    blcluster* rootA, mblock<T>** A, double eps,
                    unsigned rankmax)
{
#ifdef DEBUG
  unsigned rank = COMM_AHMED.Get_rank();
  std::cout<< "\n" << rank << "Starting Forwds: begp(" << begp << "), p(" << p
           << ")" << std::flush;
#endif
  if (p==1 || rootU->isleaf()) {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t(p==1 || rootU->isleaf()) Forwds"
              << std::flush;
#endif
    UtHhGeH_solve(rootU, A, rootA, A, eps, rankmax);
  } else {
#ifdef DEBUG
    std::cout << "\n" << rank << "\t!(p==1 || rootU->isleaf()) Forwds" << std::flush;
#endif
    unsigned n = rootA->getncs();
    unsigned rank = COMM_AHMED.Get_rank();
    unsigned nsons = rootU->getncs();

    if (begp<=rank && rank<begp+p/2) {
#ifdef DEBUG
      std::cout << "\n" << rank << "\tErste Hälfte Forwds" << std::flush;
#endif
      blcluster* sonU00 = rootU->getson(0, 0);
      for (unsigned j=0; j<n; j++) {
        blcluster *sonA0j = rootA->getson(0, j);
        LtHGeH_solve_ND_(begp, p/2, sonU00, sonA0j, A, eps, rankmax);
	if (nsons==3) {
	  blcluster* sonA2j = rootA->getson(2, j), *sonU02 = rootU->getson(0, 2);
	  mltaGeHhGeH_ND_(begp, p/2, (T)-1.0, sonU02, A, sonA0j, A, sonA2j, A,
		       eps, rankmax);
	  if (rank==begp) {
	    blcluster* sonU22 = rootU->getson(2, 2);
#ifdef DEBUG
	    std::cout<< "\n" << rank << "Waiting for Data... Forwds"
		     << std::flush;
#endif
	    send_and_addH(sonA2j, A, begp+p/2, begp, eps, rankmax);
	    UtHhGeH_solve(sonU22, A, sonA2j, A, eps, rankmax);
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
      blcluster* sonU11 = rootU->getson(1, 1);

      for (unsigned j=0; j<n; j++) {
        blcluster *sonA1j = rootA->getson(1, j);
	LtHGeH_solve_ND_(begp+p/2, p/2, sonU11, sonA1j, A, eps, rankmax);
	if (nsons==3) {
	  blcluster*sonU12 = rootU->getson(1, 2), *sonA2j = rootA->getson(2, j);
	  setGeHzero(sonA2j, A);
	  mltaGeHhGeH_ND_(begp+p/2, p/2, (T)-1.0, sonU12, A, sonA1j, A, sonA2j, A,
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
bool HCholesky_ND_(unsigned begp, unsigned p, blcluster* root,
                   mblock<T>** A, double eps, unsigned rankmax)
{
  unsigned rank = COMM_AHMED.Get_rank();

#ifdef DEBUG
  std::cout<< "\n" << rank << "Decompose: begp(" << begp << "), p(" << p
           << "), rank(" << rank << ")" << std::flush;
#endif

  if (p==1 || root->isleaf()) {
#ifdef DEBUG
    std::cout<< "\n" << rank << "\t(p == 1 || root is leaf) decompose"
             << "\n" << rank << "\tAufruf HCholesky"
             << std::flush;
#endif
    int inf = HCholesky(root, A, eps, rankmax);
    if (!inf) {
      std::cout << "\n" << rank << "HCholesky_ND Fehler" << std::flush;
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
      HCholesky_ND_(begp, p/2, son00, A, eps, rankmax);
#ifdef DEBUG
      std::cout<< "\n" << rank << "*********************************"
               << "\n" << rank << "Decompose 1 erledigt"
               <<  std::endl;
#endif
      if (nsons == 3) {
	blcluster*son22 = root->getson(2, 2), *son02 = root->getson(0, 2);
	LtHGeH_solve_ND_(begp, p/2, son00, son02, A, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Forwds 1 erledigt"
		 << std::flush;
#endif

	mltaGeHhGeH_toHeH_ND_(begp, p/2, (T)-1.0, son02, A, son02, A, son22, A,
			    eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "mltaSym beendet" << std::flush;
#endif
	if (rank==begp) {
#ifdef DEBUG
	  std::cout<< "\n" << rank << "Empfangen..." << std::flush;
#endif
	  send_and_addHSym(son22, A, begp+p/2, begp, eps, rankmax);
	  if (!HCholesky(son22, A, eps, rankmax)) {
	    std::cout << "\n" << rank << " HCholesky_ND Fehler" << std::flush;
	    exit(1);
	  }
	}
#ifdef DEBUG
	std::cout<< "\n" << rank << "Empfangen beendet..." << std::flush;
#endif
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
      HCholesky_ND_(begp+p/2, p/2, son11, A, eps, rankmax);
#ifdef DEBUG
      std::cout<< "\n" << rank << "**********************************"
               << "\n" << rank << "Decompose 2 erledigt" << std::endl;
#endif
      if (nsons == 3) {
	blcluster *son12 = root->getson(1, 2), *son22 = root->getson(2, 2);
	LtHGeH_solve_ND_(begp+p/2, p/2, son11, son12, A, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Forwds 2 erledigt" << std::flush;
#endif
	setHeHzero(son22, A);
	mltaGeHhGeH_toHeH_ND_(begp+p/2, p/2, (T)-1.0, son12, A, son12, A, son22,
			    A, eps, rankmax);
#ifdef DEBUG
	std::cout<< "\n" << rank << "Multa 2 erledigt" << std::flush;
	std::cout<< "\n" << rank << "Sending Data..." << std::flush;
#endif
	send_and_addHSym(son22, A, begp+p/2, begp, eps, rankmax);
      }
#ifdef DEBUG
      std::cout << "\n" << rank << "\tZweite Hälfte beendet" << std::flush;
#endif
    }
  }
  return true;
}

bool HCholesky_ND(unsigned p, blcluster* root, mblock<double>** A,
                  double eps, unsigned rankmax)
{
  return HCholesky_ND_(0, p, root, A, eps, rankmax);
}

bool HCholesky_ND(unsigned p, blcluster* root, mblock<float>** A,
                  double eps, unsigned rankmax)
{
  return HCholesky_ND_(0, p, root, A, eps, rankmax);
}

bool HCholesky_ND(unsigned p, blcluster* root, mblock<dcomp>** A,
                  double eps, unsigned rankmax)
{
  return HCholesky_ND_(0, p, root, A, eps, rankmax);
}

bool HCholesky_ND(unsigned p, blcluster* root, mblock<scomp>** A,
                  double eps, unsigned rankmax)
{
  return HCholesky_ND_(0, p, root, A, eps, rankmax);
}
