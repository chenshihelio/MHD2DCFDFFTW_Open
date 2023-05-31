// Basic math functions like the solver for linear equation set, etc.
#include "mathlib.h"
#include "initialize.h"
#include "parameter.h"
#include <string.h>
#include <stdio.h>

void tridiagonalSolver(int n, real* A, real *B, real*C, real *R, real *X)
{
    // A = real[0,1,...,n-2] :lower diagonal
    // B = real[0,1,...,n-1] :diagonal
    // C = real[0,1,...,n-2] :upperdiagonal
    // R = real[0,1,...,n-1] :rhs
    // X = real[0,1,...,n-1] :solution

    real m;
    for(int i=1;i<n;i++)
    {
        m = A[i-1]/B[i-1];
        B[i] = B[i] - m * C[i-1];
        R[i] = R[i] - m * R[i-1];
    }

    X[n-1] = R[n-1] / B[n-1];

    for (int i=n-2;i>=0;i--)
    {
        X[i] = (R[i] - C[i] * X[i+1]) / B[i];
    }
}


void derivCFD1st(size_t N, real DX, real *A, real *H0, real *H1, real *B, real *X,
    int boundaryType, int iProc, int nProc)
{
    // general function to calculate derivative

    // N: length of the array
    // DX: grid size
    // A: 1D array to be derived
    // H0: left Halo points
    // H1: right Halo points
    // B: array to store the right-hand-side for derivative
    // X: solution to return

    if ( boundaryType==0 || ((boundaryType==1) && 
        (iProc!=0) && (iProc!=(nProc-1) )) )
    {
        // i = 0
        B[0] = 1/DX * ( aSub1 * H0[0] + bSub1 * H0[1] + 
                        cSub1 * A[0] + dSub1 * A[1] + 
                        eSub1 * A[2]);

        // i = 1 
        B[1] = 1/DX * ( aInner1/2 * (A[2]-A[0]) + 
                        bInner1/4 * (A[3]-H0[1]) + 
                        cInner1/6 * (A[4]-H0[0])) ;

        B[2] = 1/DX * ( aInner1/2 * (A[3]-A[1]) + 
                        bInner1/4 * (A[4]-A[0]) + 
                        cInner1/6 * (A[5]-H0[1]));
                    
        for(size_t i=3;i<N-3;i++)
        {
            B[i] = 1/DX * (aInner1/2 * (A[i+1]-A[i-1]) + 
                    bInner1/4 * (A[i+2]-A[i-2]) + 
                    cInner1/6 * (A[i+3]-A[i-3]));
        }

        B[N-3] = 1/DX * (aInner1/2 * (A[N-2]-A[N-4]) + 
                        bInner1/4 * (A[N-1]-A[N-5]) + 
                        cInner1/6 * (H1[0]-A[N-6]));

        // i = N-2
        B[N-2] = 1/DX * (aInner1/2 * (A[N-1]-A[N-3]) + 
                        bInner1/4 * (H1[0]-A[N-4]) + 
                        cInner1/6 * (H1[1]-A[N-5]));

        // i = N-1
        B[N-1] = -1/DX * (aSub1 * H1[1] + bSub1 * H1[0] + 
                        cSub1 * A[N-1] + dSub1 * A[N-2] + 
                        eSub1 * A[N-3]);

        initTriMatrixCFD1st(N,boundaryType,iProc,nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
    else if ((boundaryType==1) && (iProc==0))
    {
        // i = 0
        B[0] = 1/DX * (aBound1 * A[0] + bBound1 * A[1] + 
                    cBound1 * A[2] + dBound1 * A[3]) ;

        // i = 1 
        B[1] = 1/DX * ( A1P1/2 * (A[2]-A[0])) ;

        B[2] = 1/DX * ( A1P2/2 * (A[3]-A[1]) + 
                        B1P2/4 * (A[4]-A[0]));
                
        for(size_t i=3;i<N-3;i++)
        {
            B[i] = 1/DX * (aInner1/2 * (A[i+1]-A[i-1]) + 
                bInner1/4 * (A[i+2]-A[i-2]) + 
                cInner1/6 * (A[i+3]-A[i-3]));
        }

        B[N-3] = 1/DX * (aInner1/2 * (A[N-2]-A[N-4]) + 
                        bInner1/4 * (A[N-1]-A[N-5]) + 
                        cInner1/6 * (H1[0]-A[N-6]));

        // i = N-2
        B[N-2] = 1/DX * (aInner1/2 * (A[N-1]-A[N-3]) + 
                        bInner1/4 * (H1[0]-A[N-4]) + 
                        cInner1/6 * (H1[1]-A[N-5]));

        // i = N-1
        B[N-1] = -1/DX * (aSub1 * H1[1] + bSub1 * H1[0] + 
                        cSub1 * A[N-1] + dSub1 * A[N-2] + 
                        eSub1 * A[N-3]);

        initTriMatrixCFD1st(N,boundaryType, iProc, nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
    else // boundaryType==1 and iProc==nProc-1
    {
        // i = 0
        B[0] = 1/DX * ( aSub1 * H0[0] + bSub1 * H0[1] + 
                        cSub1 * A[0] + dSub1 * A[1] + 
                        eSub1 * A[2]);

        // i = 1 
        B[1] = 1/DX * ( aInner1/2 * (A[2]-A[0]) + 
                        bInner1/4 * (A[3]-H0[1]) + 
                        cInner1/6 * (A[4]-H0[0])) ;

        B[2] = 1/DX * ( aInner1/2 * (A[3]-A[1]) + 
                        bInner1/4 * (A[4]-A[0]) + 
                        cInner1/6 * (A[5]-H0[1]));
                
        for(size_t i=3;i<N-3;i++)
        {
            B[i] = 1/DX * (aInner1/2 * (A[i+1]-A[i-1]) + 
                bInner1/4 * (A[i+2]-A[i-2]) + 
                cInner1/6 * (A[i+3]-A[i-3]));
        }


        B[N-3] = 1/DX * ( A1P2/2 * (A[N-2]-A[N-4]) + 
                        B1P2/4 * (A[N-1]-A[N-5]));

        // i = N-2
        B[N-2] = 1/DX * ( A1P1/2 * (A[N-1]-A[N-3]));

        // i = N-1
        B[N-1] = -1/DX * (aBound1 * A[N-1] + bBound1 * A[N-2] + 
                    cBound1 * A[N-3] + dBound1 * A[N-4]) ;

        initTriMatrixCFD1st(N,boundaryType, iProc, nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
}


void derivCFD2nd(size_t N, real DX, real *A, real *H0, real *H1, real *B, real *X,
    int boundaryType, int iProc, int nProc)
{
    // general function to calculate derivative

    // N: length of the array
    // DX: grid size
    // A: 1D array to be derived
    // H0: left Halo points
    // H1: right Halo points
    // B: array to store the right-hand-side for derivative
    // X: solution to return

    real DX2 = DX*DX;

    if ( boundaryType==0 || ((boundaryType==1) && 
        (iProc!=0) && (iProc!=(nProc-1) )) )
    {
        // i = 0
        B[0] = 1/DX2 * ( aSub2 * H0[0] + bSub2 * H0[1] + 
                        cSub2 * A[0] + dSub2 * A[1] + 
                        eSub2 * A[2]);

        // i = 1 
        B[1] = 1/DX2 * ( aInner2/1 * (A[2]-2*A[1]+A[0]) + 
                        bInner2/4 * (A[3]-2*A[1]+H0[1]) + 
                        cInner2/9 * (A[4]-2*A[1]+H0[0])) ;

        B[2] = 1/DX2 * ( aInner2 * (A[3]-2*A[2]+A[1]) + 
                        bInner2/4 * (A[4]-2*A[2]+A[0]) + 
                        cInner2/9 * (A[5]-2*A[2]+H0[1]));
                    
        for(size_t i=3;i<N-3;i++)
        {
            B[i] = 1/DX2 * (aInner2 * (A[i+1]-2*A[i]+A[i-1]) + 
                    bInner2/4 * (A[i+2]-2*A[i]+A[i-2]) + 
                    cInner2/9 * (A[i+3]-2*A[i]+A[i-3]));
        }

        B[N-3] = 1/DX2 * (aInner2 * (A[N-2]-2*A[N-3]+A[N-4]) + 
                        bInner2/4 * (A[N-1]-2*A[N-3]+A[N-5]) + 
                        cInner2/9 * (H1[0]-2*A[N-3]+A[N-6]));

        // i = N-2
        B[N-2] = 1/DX2 * (aInner2 * (A[N-1]-2*A[N-2]+A[N-3]) + 
                        bInner2/4 * (H1[0]-2*A[N-2]+A[N-4]) + 
                        cInner2/9 * (H1[1]-2*A[N-2]+A[N-5]));

        // i = N-1
        B[N-1] = 1/DX2 * (aSub2 * H1[1] + bSub2 * H1[0] + 
                        cSub2 * A[N-1] + dSub2 * A[N-2] + 
                        eSub2 * A[N-3]);

        initTriMatrixCFD2nd(N,boundaryType,iProc,nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
    else if ((boundaryType==1) && (iProc==0))
    {
        // i = 0
        B[0] = 1/DX2 * (aBound2 * A[0] + bBound2 * A[1] + 
                    cBound2 * A[2] + dBound2 * A[3] + 
                    eBound2 * A[4]) ;

        // i = 1 
        B[1] = 1/DX2 * ( A2P1 * (A[2]-2*A[1]+A[0])) ;

        B[2] = 1/DX2 * ( A2P2 * (A[3]-2*A[2]+A[1]) + 
                        B2P2/4 * (A[4]-2*A[2]+A[0]));
                
        for(size_t i=3;i<N-3;i++)
        {
            B[i] = 1/DX2 * (aInner2 * (A[i+1]-2*A[i]+A[i-1]) + 
                bInner2/4 * (A[i+2]-2*A[i]+A[i-2]) + 
                cInner2/9 * (A[i+3]-2*A[i]+A[i-3]));
        }

        B[N-3] = 1/DX2 * (aInner2 * (A[N-2]-2*A[N-3]+A[N-4]) + 
                        bInner2/4 * (A[N-1]-2*A[N-3]+A[N-5]) + 
                        cInner2/9 * (H1[0]-2*A[N-3]+A[N-6]));

        // i = N-2
        B[N-2] = 1/DX2 * (aInner2 * (A[N-1]-2*A[N-2]+A[N-3]) + 
                        bInner2/4 * (H1[0]-2*A[N-2]+A[N-4]) + 
                        cInner2/9 * (H1[1]-2*A[N-2]+A[N-5]));

        // i = N-1
        B[N-1] = 1/DX2 * (aSub2 * H1[1] + bSub2 * H1[0] + 
                        cSub2 * A[N-1] + dSub2 * A[N-2] + 
                        eSub2 * A[N-3]);

        initTriMatrixCFD2nd(N,boundaryType, iProc, nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
    else // boundaryType==1 and iProc==nProc-1
    {
        // i = 0
        B[0] = 1/DX2 * ( aSub2 * H0[0] + bSub2 * H0[1] + 
                        cSub2 * A[0] + dSub2 * A[1] + 
                        eSub2 * A[2]);

        // i = 1 
        B[1] = 1/DX2 * ( aInner2 * (A[2]-2*A[1]+A[0]) + 
                        bInner2/4 * (A[3]-2*A[1]+H0[1]) + 
                        cInner2/9 * (A[4]-2*A[1]+H0[0])) ;

        B[2] = 1/DX2 * ( aInner2 * (A[3]-2*A[2]+A[1]) + 
                        bInner2/4 * (A[4]-2*A[2]+A[0]) + 
                        cInner2/9 * (A[5]-2*A[2]+H0[1]));
                
        for(size_t i=3;i<N-3;i++)
        {
            B[i] = 1/DX2 * (aInner2 * (A[i+1]-2*A[i]+A[i-1]) + 
                bInner2/4 * (A[i+2]-2*A[i]+A[i-2]) + 
                cInner2/9 * (A[i+3]-2*A[i]+A[i-3]));
        }


        B[N-3] = 1/DX2 * ( A2P2 * (A[N-2]-2*A[N-3]+A[N-4]) + 
                        B2P2/4 * (A[N-1]-2*A[N-3]+A[N-5]));

        // i = N-2
        B[N-2] = 1/DX2 * ( A2P1 * (A[N-1]-2*A[N-2]+A[N-3]));

        // i = N-1
        B[N-1] = 1/DX2 * (aBound2 * A[N-1] + bBound2 * A[N-2] + 
                    cBound2 * A[N-3] + dBound2 * A[N-4] + 
                    eBound2 * A[N-5]);

        initTriMatrixCFD2nd(N,boundaryType, iProc, nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
}

void filterCFD(size_t N, real *A, real *H0, real *H1, real *B, real *X,
    int boundaryType, int iProc, int nProc)
{
    // general function to calculate filtering

    // N: length of the array
    // DX: grid size
    // A: 1D array to be filtered
    // H0: left Halo points
    // H1: right Halo points
    // B: array to store the right-hand-side for derivative
    // X: solution to return

    if ( boundaryType==0 || ((boundaryType==1) && 
        (iProc!=0) && (iProc!=(nProc-1) )) )
    {
        // i = 0
        B[0] =  aSubF * H0[0] + bSubF * H0[1] + 
                cSubF * A[0] + dSubF * A[1] + 
                eSubF * A[2];

        // i = 1 
        B[1] =  aInnerF * A[1] + 
                bInnerF/2 * (A[2]+A[0]) + 
                cInnerF/2 * (A[3]+H0[1]) +
                dInnerF/2 * (A[4]+H0[0]);

        B[2] =  aInnerF * A[2] + 
                bInnerF/2 * (A[3]+A[1]) + 
                cInnerF/2 * (A[4]+A[0]) +
                dInnerF/2 * (A[5]+H0[1]);
                    
        for(size_t i=3;i<N-3;i++)
        {
            B[i] =  aInnerF * A[i] + 
                    bInnerF/2 * (A[i+1]+A[i-1]) + 
                    cInnerF/2 * (A[i+2]+A[i-2]) + 
                    dInnerF/2 * (A[i+3]+A[i-3]);
        }

        B[N-3] =  aInnerF * A[N-3] + 
                bInnerF/2 * (A[N-2]+A[N-4]) + 
                cInnerF/2 * (A[N-1]+A[N-5]) + 
                dInnerF/2 * (H1[0]+A[N-6]);

        // i = N-2
        B[N-2] =  aInnerF * A[N-2] + 
                bInnerF/2 * (A[N-1]+A[N-3]) + 
                cInnerF/2 * (H1[0]+A[N-4]) + 
                dInnerF/2 * (H1[1]+A[N-5]);

        // i = N-1
        B[N-1] = aSubF * H1[1] + bSubF * H1[0] + 
                cSubF * A[N-1] + dSubF * A[N-2] + 
                eSubF * A[N-3];

        initTriMatrixCFDFilter(N,boundaryType,iProc,nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
    else if ((boundaryType==1) && (iProc==0))
    {
        // i = 0
        B[0] =  aBoundF0 * A[0] + bBoundF0 * A[1] + 
                cBoundF0 * A[2] + dBoundF0 * A[3] +
                eBoundF0 * A[4];

        // i = 1 
        B[1] =  aBoundF1 * A[0] + bBoundF1 * A[1] + 
                cBoundF1 * A[2] + dBoundF1 * A[3] +
                eBoundF1 * A[4];

        B[2] =  aBoundF2 * A[2] + bBoundF2/2*(A[3]+A[1]) +
                cBoundF2/2 * (A[4]+A[0]);
                
        for(size_t i=3;i<N-3;i++)
        {
            B[i] =  aInnerF * A[i] + 
                    bInnerF/2 * (A[i+1]+A[i-1]) + 
                    cInnerF/2 * (A[i+2]+A[i-2]) + 
                    dInnerF/2 * (A[i+3]+A[i-3]);
        }

        B[N-3] =  aInnerF * A[N-3] + 
                bInnerF/2 * (A[N-2]+A[N-4]) + 
                cInnerF/2 * (A[N-1]+A[N-5]) + 
                dInnerF/2 * (H1[0]+A[N-6]);

        // i = N-2
        B[N-2] =  aInnerF * A[N-2] + 
                bInnerF/2 * (A[N-1]+A[N-3]) + 
                cInnerF/2 * (H1[0]+A[N-4]) + 
                dInnerF/2 * (H1[1]+A[N-5]);

        // i = N-1
        B[N-1] = aSubF * H1[1] + bSubF * H1[0] + 
                cSubF * A[N-1] + dSubF * A[N-2] + 
                eSubF * A[N-3];

        initTriMatrixCFDFilter(N,boundaryType, iProc, nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
    else // boundaryType==1 and iProc==nProc-1
    {
        // i = 0
        B[0] =  aSubF * H0[0] + bSubF * H0[1] + 
                cSubF * A[0] + dSubF * A[1] + 
                eSubF * A[2];

        // i = 1 
        B[1] =  aInnerF * A[1] + 
                bInnerF/2 * (A[2]+A[0]) + 
                cInnerF/2 * (A[3]+H0[1]) +
                dInnerF/2 * (A[4]+H0[0]);

        B[2] =  aInnerF * A[2] + 
                bInnerF/2 * (A[3]+A[1]) + 
                cInnerF/2 * (A[4]+A[0]) +
                dInnerF/2 * (A[5]+H0[1]);
                
        for(size_t i=3;i<N-3;i++)
        {
            B[i] =  aInnerF * A[i] + 
                    bInnerF/2 * (A[i+1]+A[i-1]) + 
                    cInnerF/2 * (A[i+2]+A[i-2]) + 
                    dInnerF/2 * (A[i+3]+A[i-3]);
        }


        B[N-3] =  aBoundF2 * A[N-3] + 
                bBoundF2/2 * (A[N-4]+A[N-2]) + 
                cBoundF2/2 * (A[N-5]+A[N-1]);

        // i = N-2
        B[N-2] =  aBoundF1 * A[N-1] + bBoundF1 * A[N-2] + 
                cBoundF1 * A[N-3] + dBoundF1 * A[N-4] +
                eBoundF1 * A[N-5];

        // i = N-1
        B[N-1] = aBoundF0 * A[N-1] + bBoundF0 * A[N-2] + 
                cBoundF0 * A[N-3] + dBoundF0 * A[N-4] +
                eBoundF0 * A[N-5];

        initTriMatrixCFDFilter(N,boundaryType, iProc, nProc);
        tridiagonalSolver(N, diagLow, diagCent, diagUp, B, X);
    }
}