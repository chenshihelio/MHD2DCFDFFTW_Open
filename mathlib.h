//Defines the basic mathematical functions
#ifndef _MATHLIB_
#define _MATHLIB_

#include <stdlib.h>
#include "parameter.h"

void tridiagonalSolver(int n, real* A, real *B, real*C, real *R, real *X);

void derivCFD1st(size_t N, real DX, real *A, real *H0, real *H1, real *B, real *X,
    int boundaryType, int iProc, int nProc);

void derivCFD2nd(size_t N, real DX, real *A, real *H0, real *H1, real *B, real *X,
    int boundaryType, int iProc, int nProc);

void filterCFD(size_t N, real *A, real *H0, real *H1, real *B, real *X,
    int boundaryType, int iProc, int nProc);

#endif