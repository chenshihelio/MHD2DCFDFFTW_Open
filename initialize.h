//Declare the fields and variables like the grid points

#ifndef _INITIALIZE_
#define _INITIALIZE_

#include <stdio.h>
#include <mpi.h>
#include <fftw3.h>
#include "parameter.h"

extern int NVAR;
extern int boundaryTypeX;
extern real cfl,Tmax,dtOut;
extern size_t nx, ny, iOffset, iSize, 
    jOffset, jSize, arrSize, nky;
extern real Lx, dx,Ly,dy; 
extern real adiabaticIdx,resistivity,viscosity;
extern int Hall, divBcleaning, outDivB;
extern real Ch_divBcleaning, Cp_divBcleaning;// divB cleaning parameters
extern real constDiffDivB;
extern real di;
extern real *uu;  //main array
extern real *dudx, *dudy; //1st order derivative
extern real *d2udx2, *d2udy2; // 2nd order derivs of velocity
extern real *d2Bdx2, *d2Bdy2; // magnetic field 2nd-deriv
extern real *divB, *ddivBdx, *ddivBdy;
extern real *Efield, *dEdx, *dEdy;
extern real *dudt, *dudtRK; //RHS
extern real *xgrid, *ygrid;
extern real *tukeyWindowX; // used for smoothing the Hall term in the open boundary case


//CFD schemes, 1st order derivatives
extern real alphaInner1,aInner1,bInner1,cInner1; 
extern real alphaBound1,aBound1,bBound1,cBound1,dBound1;
extern real alphaSub1,aSub1,bSub1,cSub1,dSub1,eSub1;
//CFD schemes, 2nd order derivatives
extern real alphaInner2,aInner2,bInner2,cInner2;
extern real alphaBound2,aBound2,bBound2,cBound2,dBound2,eBound2;
extern real alphaSub2,aSub2,bSub2,cSub2,dSub2,eSub2;
//CFD schemes, filter
extern real alphaInnerF,aInnerF,bInnerF,cInnerF,dInnerF;
extern real alphaSubF,aSubF,bSubF,cSubF,dSubF,eSubF;
extern real alphaBoundF0,aBoundF0,bBoundF0,cBoundF0,dBoundF0,eBoundF0;
extern real alphaBoundF1,aBoundF1,bBoundF1,cBoundF1,dBoundF1,eBoundF1;
extern real alphaBoundF2,aBoundF2,bBoundF2,cBoundF2;


extern int filteringOptionY;
extern real afy;
extern real *filtY;

extern real *diagLow, *diagCent, *diagUp,*cfdRhs,*cfdSol,
    *cfdArr, *cfdHalo0, *cfdHalo1;

extern int npe, myRank, iRank;
extern real *uuHaloX0, *uuHaloX1;
extern real *divBHaloX0, *divBHaloX1;
extern real *EfieldHaloX0, *EfieldHaloX1;

extern MPI_Datatype typeHaloSendX0, typeHaloSendX1;
extern MPI_Datatype typeHalodivBSendX0, typeHalodivBSendX1;
extern MPI_Datatype typeHaloEfieldSendX0, typeHaloEfieldSendX1;
extern MPI_Datatype typeMainArr, typeDivBArr; //for MPI I/O

/* variables for FFTW */
extern fftw_plan fftwPlanR2C, fftwPlanC2R;
extern real *fftwIn;
extern fftw_complex *fftwOut;


void initialize();
void finalize();

void initTriMatrixCFD1st(int n, int boundaryType, int rank, int nProc);
void initTriMatrixCFD2nd(int n, int boundaryType, int rank, int nProc);
void initTriMatrixCFDFilter(int n, int boundaryType, int rank, int nProc);

#endif