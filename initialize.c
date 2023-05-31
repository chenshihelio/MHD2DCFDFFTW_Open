#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include <mpi.h>
#include "macros.h"
#include "parameter.h"

//macros
#define MAXBUFFSIZE 200

//declare all the variables 

/* variables that are read from input file */
int boundaryTypeX;
int divBcleaning = 0, outDivB = 0;
real constDiffDivB = 2.0;
real cfl,Tmax,dtOut;
size_t nx, ny, iOffset, iSize, jOffset, jSize, arrSize, nky;
real Lx,dx,Ly,dy;
real adiabaticIdx,resistivity,viscosity;
int Hall = 0;
real di;

/* variables related to numerical schemes */
//CFD schemes, 1st order derivatives
real alphaInner1 = ALPHAINNER1,aInner1,bInner1,cInner1; 
real alphaBound1 = ALPHABOUND1,aBound1,bBound1,cBound1,dBound1;
real alphaSub1 = ALPHASUB1,aSub1,bSub1,cSub1,dSub1,eSub1;
//CFD schemes, 2nd order derivatives
real alphaInner2 = ALPHAINNER2,aInner2,bInner2,cInner2;
real alphaBound2 = ALPHABOUND2,aBound2,bBound2,cBound2,dBound2,eBound2;
real alphaSub2 = ALPHASUB1,aSub2,bSub2,cSub2,dSub2,eSub2;
//CFD schemes, filter
real alphaInnerF = ALPHAINNERF,aInnerF,bInnerF,cInnerF,dInnerF;
real alphaSubF = ALPHASUBF,aSubF,bSubF,cSubF,dSubF,eSubF;
real alphaBoundF0 = ALPHABOUNDF0,aBoundF0,bBoundF0,cBoundF0,dBoundF0,eBoundF0;
real alphaBoundF1 = ALPHABOUNDF1,aBoundF1,bBoundF1,cBoundF1,dBoundF1,eBoundF1;
real alphaBoundF2 = ALPHABOUNDF2,aBoundF2,bBoundF2,cBoundF2;

// 0: cut-off at nz/3
// 1: CFD-like filtering use (Lele, 1992) Equation (C.2.4) with d = beta = 0 
int filteringOptionY = 1;
real afy = 0.45;

/* main arrays */
int NVAR;
real *uu, *dudx, *dudy,*d2udx2, *d2udy2, 
    *d2Bdx2, *d2Bdy2, *dudt, *dudtRK;
real *divB, *ddivBdx, *ddivBdy;
real *Efield, *dEdx, *dEdy;
real *xgrid, *ygrid;
real *filtY; // for filtering along Y
real *tukeyWindowX; // used for smoothing the Hall term in the open boundary case

/* divB cleaning parameters*/
real Ch_divBcleaning, Cp_divBcleaning;

/* halo points */
real *uuHaloX0, *uuHaloX1; // 0: left, 1: right
real *divBHaloX0, *divBHaloX1;
real *EfieldHaloX0, *EfieldHaloX1;

/* variables for compact finite difference schemes */
real *diagLow, *diagCent, *diagUp,
    *cfdRhs, *cfdSol, *cfdArr, *cfdHalo0, *cfdHalo1,*filterDeltaX;


/* variables related to MPI */
int npe, myRank, iRank;
MPI_Datatype typeHaloSendX0, typeHaloSendX1;
MPI_Datatype typeHalodivBSendX0, typeHalodivBSendX1;
MPI_Datatype typeHaloEfieldSendX0, typeHaloEfieldSendX1;
MPI_Datatype typeMainArr, typeDivBArr;

/* variables for FFTW */
fftw_plan fftwPlanR2C, fftwPlanC2R;
real *fftwIn;
fftw_complex *fftwOut;

// declare functions in this file
void readInput(const char *fileName);
void decomposeGrid1D(size_t n, int nProc, int iProc, size_t *iOffset, size_t *iSize);
void initializeCFDSchemes();

//initialize all the arrays and plans
void initialize()
{
    int dims[2], periods[2], mpiCoords[2];
    int bigSize[3], smallSize[3], startIdx[3];

    MPI_Init(NULL, NULL);
    
    MPI_Comm_size(MPI_COMM_WORLD, &npe);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    iRank = myRank;

    if (myRank==0)
    {
        printf("2D MHD CFD X FFTW code......\n");
        readInput("./input");
    }


    MPI_Bcast(&boundaryTypeX, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cfl,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&Tmax,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&dtOut,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&nx,1,mpiSIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&ny,1,mpiSIZE_T,0,MPI_COMM_WORLD);
    MPI_Bcast(&Lx,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&Ly,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&adiabaticIdx,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&Hall,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&di,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&resistivity,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&viscosity,1,mpi_real,0,MPI_COMM_WORLD);
    MPI_Bcast(&divBcleaning,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&outDivB,1,MPI_INT,0,MPI_COMM_WORLD);


    nky = (int)(floor(ny/2)+1);

    if(boundaryTypeX==0)
    {
        dx = Lx/nx;
    }
    else
    {
        dx = Lx/(nx-1);
    }


    dy = Ly/ny;

    decomposeGrid1D(nx, npe, iRank, &iOffset, &iSize);

    jOffset = 0;
    jSize = ny;

    xgrid = fftw_malloc(iSize * sizeof(real));
    ygrid = fftw_malloc(jSize * sizeof(real));

    for(int i=0;i<iSize;i++)
    {
        xgrid[i] = (real)(i+iOffset) * dx;
    }
    for(int j=0;j<jSize;j++)
    {
        ygrid[j] = (real)(j+jOffset) * dy;
    }

    constDiffDivB = constDiffDivB * MAX(dx,dy) * MAX(dx,dy);

    // initialize the CFD schemes
    initializeCFDSchemes();


    // initialize the filter along Z
    if(filteringOptionY==1)
    {
        filtY = fftw_malloc(nky * sizeof(real));
        memset(filtY, 0, nky * sizeof(real));

        real aj = (5.0 + 6.0 * afy ) / 8.0;
        real bj = (1.0 + 2.0 * afy ) / 2.0;
        real cj =-(1.0 - 2*afy) / 8.0;

        real wav_len;
        for(size_t k=0;k<nky;k++)
        { 
            wav_len = 2. * PI * (real)k/ny;
            filtY[k] = (aj + bj * cos(wav_len) + cj * cos(2*wav_len)) /
                (1 + 2 * afy * cos(wav_len));
        }
    }


    // initialize the tukey window
    tukeyWindowX = fftw_malloc(iSize * sizeof(real));
    memset(tukeyWindowX, 0, iSize * sizeof(real));

    real normX;
    real tukeyWindowWidthX = 0.2;
    for (size_t i=0;i<iSize;i++)
    {
        normX = xgrid[i] / Lx;

        if ((normX >= tukeyWindowWidthX) && 
            (normX <= (1.0-tukeyWindowWidthX)))
        {
            tukeyWindowX[i] = 1.0;
        }
        else if (normX < tukeyWindowWidthX)
        {
            tukeyWindowX[i] = (1 - cos(PI * normX/tukeyWindowWidthX))*0.5;
        }
        else 
        {
            tukeyWindowX[i] = (1 - cos(PI * (1-normX)/tukeyWindowWidthX))*0.5;
        }
    }

    // determine NVAR: if divBcleaning == 2 (hyperbolic cleaning (Dedner+2002)),
    // we need an additional variable
    if(divBcleaning==2)
    {
        NVAR = 9;
    }
    else 
    {
        NVAR = 8;
    }

    arrSize = NVAR * iSize * jSize;

    //allocate the main array and the derivative array
    // Recommend to memset: otherwise strange things may happen
    uu = fftw_malloc( arrSize * sizeof(real)); 
    memset(uu, 0, arrSize * sizeof(real)); 

    dudx = fftw_malloc( arrSize * sizeof(real));
    memset(dudx, 0, arrSize * sizeof(real));

    dudy = fftw_malloc( arrSize * sizeof(real));
    memset(dudy, 0, arrSize * sizeof(real));

    if(viscosity>0)
    {
        d2udx2 = fftw_malloc( 3 * iSize * jSize * sizeof(real));
        memset(d2udx2, 0, 3 * iSize * jSize * sizeof(real));

        d2udy2 = fftw_malloc( 3 * iSize * jSize * sizeof(real));
        memset(d2udy2, 0, 3 * iSize * jSize * sizeof(real));
    }

    if(resistivity>0 || Hall==1)
    {
        d2Bdx2 = fftw_malloc( 3 * iSize * jSize * sizeof(real));
        memset(d2Bdx2, 0, 3 * iSize * jSize * sizeof(real));

        d2Bdy2 = fftw_malloc( 3 * iSize * jSize * sizeof(real));
        memset(d2Bdy2, 0, 3 * iSize * jSize * sizeof(real));
    }

    divB = fftw_malloc( iSize * jSize * sizeof(real));
    memset(divB, 0, iSize * jSize * sizeof(real));

    if(divBcleaning == 1)
    {
        ddivBdx = fftw_malloc( iSize * jSize * sizeof(real));
        memset(ddivBdx, 0, iSize * jSize * sizeof(real));

        ddivBdy = fftw_malloc( iSize * jSize * sizeof(real));
        memset(ddivBdy, 0, iSize * jSize * sizeof(real));

        divBHaloX0 = fftw_malloc( NHALO * jSize * sizeof(real));
        memset(divBHaloX0, 0, NHALO * jSize * sizeof(real));

        divBHaloX1 = fftw_malloc( NHALO * jSize * sizeof(real));
        memset(divBHaloX1, 0, NHALO * jSize * sizeof(real));
    }


    // allocate Electric field
    Efield = fftw_malloc( 3 * iSize * jSize * sizeof(real));
    memset(Efield, 0, 3 * iSize * jSize * sizeof(real));

    dEdx = fftw_malloc( 3 * iSize * jSize * sizeof(real));
    memset(dEdx, 0, 3 * iSize * jSize * sizeof(real));

    dEdy = fftw_malloc( 3 * iSize * jSize * sizeof(real));
    memset(dEdy, 0, 3 * iSize * jSize * sizeof(real));

    EfieldHaloX0 = fftw_malloc( 3 * NHALO * jSize * sizeof(real));
    memset(EfieldHaloX0, 0, NHALO * jSize * sizeof(real));

    EfieldHaloX1 = fftw_malloc( 3 * NHALO * jSize * sizeof(real));
    memset(EfieldHaloX1, 0, NHALO * jSize * sizeof(real));



    dudt = fftw_malloc( arrSize * sizeof(real));
    memset(dudt, 0, arrSize * sizeof(real));

    dudtRK = fftw_malloc( arrSize * sizeof(real));
    memset(dudtRK, 0, arrSize * sizeof(real));

    uuHaloX0 = fftw_malloc( NVAR * NHALO * jSize * sizeof(real));
    memset(uuHaloX0, 0, NVAR * NHALO * jSize * sizeof(real));

    uuHaloX1 = fftw_malloc( NVAR * NHALO * jSize * sizeof(real));
    memset(uuHaloX1, 0, NVAR * NHALO * jSize * sizeof(real));

    
    // arrays needed for compact difference scheme
    diagLow = fftw_malloc( (nx-1) * sizeof(real) );
    memset(diagLow, 0, (nx-1) * sizeof(real));

    diagCent = fftw_malloc( (nx) * sizeof(real) );
    memset(diagCent, 0, (nx) * sizeof(real));

    diagUp = fftw_malloc( (nx-1) * sizeof(real) );
    memset(diagUp, 0, (nx-1) * sizeof(real));

    cfdRhs = fftw_malloc( (nx) * sizeof(real) );
    memset(cfdRhs, 0, (nx) * sizeof(real));

    cfdSol = fftw_malloc( (nx) * sizeof(real) );
    memset(cfdSol, 0, (nx) * sizeof(real));

    cfdArr = malloc( (nx) * sizeof(real) );
    memset(cfdArr, 0, (nx) * sizeof(real));

    cfdHalo0 = malloc( NHALO * sizeof(real) );
    memset(cfdHalo0, 0, NHALO * sizeof(real));

    cfdHalo1 = malloc( NHALO * sizeof(real) );
    memset(cfdHalo1, 0, NHALO * sizeof(real));

    //fftw plans and arrays
    fftwIn = fftw_malloc( ny * sizeof(real) );
    fftwOut = (fftw_complex*) fftw_malloc(  nky * sizeof(fftw_complex));
    fftwPlanR2C = fftw_plan_dft_r2c_1d(ny, fftwIn, fftwOut, FFTW_MEASURE);
    fftwPlanC2R = fftw_plan_dft_c2r_1d(ny, fftwOut, fftwIn, FFTW_MEASURE);


    // create datatypes for MPI
    // Type of the main array: for writing files
    bigSize[0] = NVAR; bigSize[1] = nx; bigSize[2] = ny;
    smallSize[0] = NVAR; smallSize[1] = iSize; smallSize[2] = jSize;
    startIdx[0] = 0; startIdx[1] = iOffset; startIdx[2] = jOffset;

    MPI_Type_create_subarray(3, bigSize, smallSize, startIdx, 
        MPI_ORDER_C, mpi_real, &typeMainArr);
    MPI_Type_commit(&typeMainArr);

    // Type of the divB array: for writing files
    bigSize[0] = 1; bigSize[1] = nx; bigSize[2] = ny;
    smallSize[0] = 1; smallSize[1] = iSize; smallSize[2] = jSize;
    startIdx[0] = 0; startIdx[1] = iOffset; startIdx[2] = jOffset;

    MPI_Type_create_subarray(3, bigSize, smallSize, startIdx, 
        MPI_ORDER_C, mpi_real, &typeDivBArr);
    MPI_Type_commit(&typeDivBArr);


    // Type of Halo points of main array along X
    bigSize[0] = NVAR; bigSize[1] = iSize; bigSize[2] = jSize;
    smallSize[0] = NVAR; smallSize[1] = NHALO; smallSize[2] = jSize;

    startIdx[0] = 0; startIdx[1] = 0; startIdx[2] = 0;
    MPI_Type_create_subarray(3, bigSize, smallSize, startIdx,
        MPI_ORDER_C, mpi_real, &typeHaloSendX0);
    MPI_Type_commit(&typeHaloSendX0);

    startIdx[0] = 0; startIdx[1] = iSize-NHALO; startIdx[2] = 0;
    MPI_Type_create_subarray(3, bigSize, smallSize, startIdx,
        MPI_ORDER_C, mpi_real, &typeHaloSendX1);
    MPI_Type_commit(&typeHaloSendX1);

    if(divBcleaning == 1)
    {
        // Type of Halo points of ddivBdx along X
        bigSize[0] = 1; bigSize[1] = iSize; bigSize[2] = jSize;
        smallSize[0] = 1; smallSize[1] = NHALO; smallSize[2] = jSize;

        startIdx[0] = 0; startIdx[1] = 0; startIdx[2] = 0;
        MPI_Type_create_subarray(3, bigSize, smallSize, startIdx,
            MPI_ORDER_C, mpi_real, &typeHalodivBSendX0);
        MPI_Type_commit(&typeHalodivBSendX0);

        startIdx[0] = 0; startIdx[1] = iSize-NHALO; startIdx[2] = 0;
        MPI_Type_create_subarray(3, bigSize, smallSize, startIdx,
            MPI_ORDER_C, mpi_real, &typeHalodivBSendX1);
        MPI_Type_commit(&typeHalodivBSendX1);
    }


    // MPI_type for Efield halo points
    bigSize[0] = 3; bigSize[1] = iSize; bigSize[2] = jSize;
    smallSize[0] = 3; smallSize[1] = NHALO; smallSize[2] = jSize;

    startIdx[0] = 0; startIdx[1] = 0; startIdx[2] = 0;
    MPI_Type_create_subarray(3, bigSize, smallSize, startIdx,
        MPI_ORDER_C, mpi_real, &typeHaloEfieldSendX0);
    MPI_Type_commit(&typeHaloEfieldSendX0);

    startIdx[0] = 0; startIdx[1] = iSize-NHALO; startIdx[2] = 0;
    MPI_Type_create_subarray(3, bigSize, smallSize, startIdx,
        MPI_ORDER_C, mpi_real, &typeHaloEfieldSendX1);
    MPI_Type_commit(&typeHaloEfieldSendX1);
    
}


// destroy and free all the arrays
void finalize()
{
    free(uu); free(dudx); free(dudy); 
    fftw_free(d2udx2); fftw_free(d2udy2);
    fftw_free(d2Bdx2); fftw_free(d2Bdy2);

    free(xgrid); free(ygrid);
    free(dudt); free(dudtRK);

    free(diagLow); free(diagCent);
    free(diagUp); free(cfdRhs);
    free(cfdSol);

    free(uuHaloX0); free(uuHaloX1);

    fftw_destroy_plan(fftwPlanR2C);
    fftw_destroy_plan(fftwPlanC2R);
    fftw_free(fftwIn); fftw_free(fftwOut);

    MPI_Finalize();
}


// function that reads the value from one single line of the file
// char *readLine(FILE *fp)
// {
//     char fileLine[MAXBUFFSIZE];
//     char *token, *buf, *temp;

//     temp = fgets(fileLine, MAXBUFFSIZE, fp);

//     buf = fileLine;
//     //printf("buf=%s\n",buf);
//     token = strsep(&buf, ";");

//     // do not understand how the original string is changed...
//     // printf("buf=%s\n",buf);
//     // printf("fileLine=%s\n",fileLine); 

//     return token;
// }

//the above version does not work on COMET(strsep)
void readLine(FILE *fp, char* returnLine) 
{
    char fileLine[MAXBUFFSIZE];
    char  *temp;

    temp = fgets(fileLine, MAXBUFFSIZE, fp);

    for(int i=0;i<16;i++)
    {
        returnLine[i] = fileLine[i];
    }
}

// Function that reads the input file
void readInput(const char *fileName)
{
    FILE *fpRead=fopen(fileName, "r");

    char *valueRead;

    valueRead = malloc(16*sizeof(char));

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    boundaryTypeX = atoi(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    cfl = atof(valueRead);

    // readLine(fpRead,valueRead);
    // dealisingOption = (int)atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    Tmax = atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    dtOut = atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    nx = atoi(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    ny = atoi(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    Lx = atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    Ly = atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    adiabaticIdx = atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    Hall = atoi(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    di = atof(valueRead);

    // valueRead = readLine(fpRead);
    readLine(fpRead,valueRead);
    resistivity = atof(valueRead);

    readLine(fpRead,valueRead);
    viscosity = atof(valueRead);

    readLine(fpRead,valueRead);
    divBcleaning = atoi(valueRead);

    readLine(fpRead,valueRead);
    outDivB = atoi(valueRead);

    int closeStatus = fclose(fpRead);
}


// // old version: not good: if n mod nProc /=0, the last processor will usually have
// // much more points than others
// void decomposeGrid1D(size_t n, int nProc, int iProc, size_t *iOffset, size_t *iSize)
// {
//     size_t normalSize = n/nProc;

//     if (iProc < nProc-1)
//     {
//         *iSize = normalSize;
//         *iOffset = iProc * normalSize;
//     }
//     else
//     {
//         *iSize = n - normalSize * (nProc-1);
//         *iOffset = (nProc-1) * normalSize;
//     }
// }


void decomposeGrid1D(size_t n, int nProc, int iProc, size_t *iOffset, size_t *iSize)
{
    size_t normalSize = n/nProc;
    size_t remaining = n%nProc;
    size_t arrayN[nProc], arrayOffset[nProc];

    for(int i=0;i<nProc;i++)
    {
        arrayN[i] = normalSize;
    }
    for(int i=0;i<remaining;i++)
    {
        arrayN[i] = arrayN[i] + 1;
    }
    arrayOffset[0] = 0;
    for(int i=1;i<nProc;i++)
    {
        arrayOffset[i] = arrayOffset[i-1] + arrayN[i-1];
    }

    *iSize = arrayN[iProc];
    *iOffset = arrayOffset[iProc];
}


void initializeCFDSchemes()
{
    aInner1 = (alphaInner1 + 9.)/6.;
    bInner1 = (32*alphaInner1 - 9.)/15.;
    cInner1 = (-3.*alphaInner1 + 1)/10.;

    aInner2 = (6. - 9. * alphaInner2)/4.;
    bInner2 = (-3. + 24.*alphaInner2)/5.;
    cInner2 = (2. - 11. * alphaInner2)/20.;

    aBound1 = -(11. + 2. * alphaBound1)/6.;
    bBound1 = (6. - alphaBound1)/2.;
    cBound1 = (2. * alphaBound1-3.)/2.;
    dBound1 = (2. - alphaBound1)/6.;

    aBound2 = (11.*alphaBound2 + 35.)/12.;
    bBound2 = -(5.*alphaBound2 + 26.)/3.;
    cBound2 = (alphaBound2 + 19.)/2.;
    dBound2 = (alphaBound2 - 14.)/3.;
    eBound2 = (11. - alphaBound2)/12.;

    aSub1 = (1. - alphaSub1)/12.;
    bSub1 = (3. * alphaSub1 - 4.)/6.;
    cSub1 = -3./2.*alphaSub1;
    dSub1 = (4. + 5. * alphaSub1) / 6.;
    eSub1 = (3. * alphaSub1 - 1.) / 12.;
    
    aSub2 = -(1. + alphaSub2) / 12.;
    bSub2 = (alphaSub2 + 4.) / 3.;
    cSub2 = (alphaSub2 - 5.) / 2.;
    dSub2 = (4. - 5. * alphaSub2) / 3.;
    eSub2 = (11. * alphaSub2 -  1.) / 12.;

    aInnerF = (11. + 10.*alphaInnerF) / 16.;
    bInnerF = (15. + 34.*alphaInnerF) / 32.;
    cInnerF = (-3. + 6. *alphaInnerF) / 16.;
    dInnerF = (1. - 2. * alphaInnerF) / 32.;

    aSubF = (alphaSubF - 1.)/16.;
    bSubF = (1 - alphaSubF)/4.;
    cSubF = (3. * alphaSubF + 5.)/8.;
    dSubF = (3. * alphaSubF + 1.)/4.;
    eSubF = (alphaSubF - 1.)/16.;

    aBoundF0 = (15.+alphaBoundF0)/16.;
    bBoundF0 = (1.+3.*alphaBoundF0)/4.;
    cBoundF0 = 3.*(alphaBoundF0-1)/8.;
    dBoundF0 = (1.-alphaBoundF0)/4.;
    eBoundF0 = (alphaBoundF0-1.)/16.;

    aBoundF1 = (1.+14.*alphaBoundF1)/16.;
    bBoundF1 = (3.+2.*alphaBoundF1)/4.;
    cBoundF1 = (3.+2.*alphaBoundF1)/8.;
    dBoundF1 = (-1.+2.*alphaBoundF1)/4.;
    eBoundF1 = (1.-2.*alphaBoundF1)/16.;

    aBoundF2 = (5.+6.*alphaBoundF2)/8.;
    bBoundF2 = (1.+2.*alphaBoundF2)/2.;
    cBoundF2 = -(1.-2.*alphaBoundF2)/8.;
}

void initTriMatrixCFD1st(int n, int boundaryType, int rank, int nProc) 
//Be cautious about the size n: should not be larger than nMax
{
    for(int i=0;i<n;i++)
    {
        diagCent[i] = 1.0;
    }

    if (boundaryType == 0)
    {
        for(int i=0;i<n-2;i++)
        {
            diagLow[i] = alphaInner1;
        }
        diagLow[n-2] = alphaSub1;

        diagUp[0] = alphaSub1;
        for(int i=1;i<n-1;i++)
        {
            diagUp[i] = alphaInner1;
        }
    }
    else if (boundaryType == 1)
    {
        if ((rank != 0) && (rank != nProc-1))
        {
            for(int i=0;i<n-2;i++)
            {
                diagLow[i] = alphaInner1;
            }
            diagLow[n-2] = alphaSub1;

            diagUp[0] = alphaSub1;
            for(int i=1;i<n-1;i++)
            {
                diagUp[i] = alphaInner1;
            }
        }
        else if (rank==0)
        {
            diagLow[0] = ALPHA1P1;
            diagLow[1] = ALPHA1P2;
            for(int i=2;i<n-2;i++)
            {
                diagLow[i] = alphaInner1;
            }
            diagLow[n-2] = alphaSub1;

            diagUp[0] = alphaBound1;
            diagUp[1] = ALPHA1P1;
            diagUp[2] = ALPHA1P2;
            for(int i=3;i<n-1;i++)
            {
                diagUp[i] = alphaInner1;
            }
        }
        else
        {
            for(int i=0;i<n-4;i++)
            {
                diagLow[i] = alphaInner1;
            }
            diagLow[n-4] = ALPHA1P2;
            diagLow[n-3] = ALPHA1P1;
            diagLow[n-2] = alphaBound1;

            diagUp[0] = alphaSub1;
            for(int i=1;i<n-3;i++)
            {
                diagUp[i] = alphaInner1;
            }
            diagUp[n-3] = ALPHA1P2;
            diagUp[n-2] = ALPHA1P1;
        }
    }
}

void initTriMatrixCFD2nd(int n, int boundaryType, int rank, int nProc) 
//Be cautious about the size n: should not be larger than nMax
{
    for(int i=0;i<n;i++)
    {
        diagCent[i] = 1.0;
    }

    if (boundaryType == 0)
    {
        for(int i=0;i<n-2;i++)
        {
            diagLow[i] = alphaInner2;
        }
        diagLow[n-2] = alphaSub2;

        diagUp[0] = alphaSub2;
        for(int i=1;i<n-1;i++)
        {
            diagUp[i] = alphaInner2;
        }
    }
    else if (boundaryType == 1)
    {
        if ((rank != 0) && (rank != nProc-1))
        {
            for(int i=0;i<n-2;i++)
            {
                diagLow[i] = alphaInner2;
            }
            diagLow[n-2] = alphaSub2;

            diagUp[0] = alphaSub2;
            for(int i=1;i<n-1;i++)
            {
                diagUp[i] = alphaInner2;
            }
        }
        else if (rank==0)
        {
            diagLow[0] = ALPHA2P1;
            diagLow[1] = ALPHA2P2;
            for(int i=2;i<n-2;i++)
            {
                diagLow[i] = alphaInner2;
            }
            diagLow[n-2] = alphaSub2;

            diagUp[0] = alphaBound2;
            diagUp[1] = ALPHA2P1;
            diagUp[2] = ALPHA2P2;
            for(int i=3;i<n-1;i++)
            {
                diagUp[i] = alphaInner2;
            }
        }
        else
        {
            for(int i=0;i<n-4;i++)
            {
                diagLow[i] = alphaInner2;
            }
            diagLow[n-4] = ALPHA2P2;
            diagLow[n-3] = ALPHA2P1;
            diagLow[n-2] = alphaBound2;

            diagUp[0] = alphaSub2;
            for(int i=1;i<n-3;i++)
            {
                diagUp[i] = alphaInner2;
            }
            diagUp[n-3] = ALPHA2P2;
            diagUp[n-2] = ALPHA2P1;
        }
    }
}


void initTriMatrixCFDFilter(int n, int boundaryType, int rank, int nProc)
{
    for(int i=0;i<n;i++)
    {
        diagCent[i] = 1.0;
    }

    if (boundaryType == 0)
    {
        for(int i=0;i<n-2;i++)
        {
            diagLow[i] = alphaInnerF;
        }
        diagLow[n-2] = alphaSubF;

        diagUp[0] = alphaSubF;
        for(int i=1;i<n-1;i++)
        {
            diagUp[i] = alphaInnerF;
        }
    }
    else if (boundaryType == 1)
    {
        if ((rank != 0) && (rank != nProc-1))
        {
            for(int i=0;i<n-2;i++)
            {
                diagLow[i] = alphaInnerF;
            }
            diagLow[n-2] = alphaSubF;

            diagUp[0] = alphaSubF;
            for(int i=1;i<n-1;i++)
            {
                diagUp[i] = alphaInnerF;
            }
        }
        else if (rank==0)
        {
            diagLow[0] = alphaBoundF1;
            diagLow[1] = alphaBoundF2;
            for(int i=2;i<n-2;i++)
            {
                diagLow[i] = alphaInnerF;
            }
            diagLow[n-2] = alphaSubF;

            diagUp[0] = alphaBoundF0;
            diagUp[1] = alphaBoundF1;
            diagUp[2] = alphaBoundF2;
            for(int i=3;i<n-1;i++)
            {
                diagUp[i] = alphaInnerF;
            }
        }
        else
        {
            for(int i=0;i<n-4;i++)
            {
                diagLow[i] = alphaInnerF;
            }
            diagLow[n-4] = alphaBoundF2;
            diagLow[n-3] = alphaBoundF1;
            diagLow[n-2] = alphaBoundF0;

            diagUp[0] = alphaSubF;
            for(int i=1;i<n-3;i++)
            {
                diagUp[i] = alphaInnerF;
            }
            diagUp[n-3] = alphaBoundF2;
            diagUp[n-2] = alphaBoundF1;
        }
    }
}