/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include "matrix.h"

extern unsigned GMRes(const Matrix<double>&, double* const, double* const,
                      double&, const unsigned, unsigned&);
extern unsigned FGMRes(const Matrix<double>&, double* const, double* const,
                       double&, const unsigned, unsigned&);
extern unsigned GMRes(const Matrix<dcomp>&, dcomp* const, dcomp* const,
                      double&, const unsigned, unsigned&);
extern unsigned FGMRes(const Matrix<dcomp>&, dcomp* const, dcomp* const,
                       double&, const unsigned, unsigned&);


extern unsigned BiCGStab(const Matrix<double>&, double* const,
                         double* const, double&, unsigned&);
extern unsigned BiCGStab(const Matrix<dcomp>&, dcomp* const,
                         dcomp* const, double&, unsigned&);

extern unsigned CG(const Matrix<double>&, double* const, double* const,
                   double&, unsigned&);
extern unsigned CG(const Matrix<dcomp>&, dcomp* const, dcomp* const,
                   double&, unsigned&);

extern unsigned MinRes(const Matrix<double>&, double* const, double* const,
                       double&, unsigned&);

