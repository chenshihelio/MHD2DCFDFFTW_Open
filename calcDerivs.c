#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "parameter.h"
#include "initialize.h"
#include "mathlib.h"


void derivX(int nv, real *arrMain, real *arrHalo0, 
     real *arrHalo1, real *arrDeriv, int iVarStartMain,
     int iVarStartHalo)
{
     int npex = npe;
     // iVarStartMain: we sometimes do not need to calculate derivative for all variables
     for(int ivar=0; ivar<nv;ivar++)
     {
          for(size_t j=0; j<jSize;j++)
          {  
               for(size_t i=0;i<iSize;i++)
               {
                    cfdArr[i] = arrMain[IDXIJ(ivar+iVarStartMain,i,j)];
               }

               for(size_t i=0;i<NHALO;i++)
               {
                    cfdHalo0[i] = arrHalo0[IDX(ivar+iVarStartHalo,i,j,NHALO,jSize)];
                    cfdHalo1[i] = arrHalo1[IDX(ivar+iVarStartHalo,i,j,NHALO,jSize)];
               }

               derivCFD1st(iSize, dx, cfdArr, cfdHalo0, cfdHalo1, cfdRhs, cfdSol,
                         boundaryTypeX, iRank, npex);

               for(size_t i=0;i<iSize;i++)
               {
                    arrDeriv[IDXIJ(ivar,i,j)] = cfdSol[i];
               }
          }          
     }
}

void derivXX(int nv, real *arrMain, real *arrHalo0, 
     real *arrHalo1, real *arrDeriv, int iVarStartMain, 
     int iVarStartHalo)
{
     int npex = npe;
     // iVarStartMain: we sometimes do not need to calculate derivative for all variables
     for(int ivar=0; ivar<nv;ivar++)
     {

          for(size_t j=0; j<jSize;j++)
          {  
               for(size_t i=0;i<iSize;i++)
               {
                    cfdArr[i] = arrMain[IDXIJ(ivar+iVarStartMain,i,j)];
               }

               for(size_t i=0;i<NHALO;i++)
               {
                    cfdHalo0[i] = arrHalo0[IDX(ivar+iVarStartHalo,i,j,NHALO,jSize)];
                    cfdHalo1[i] = arrHalo1[IDX(ivar+iVarStartHalo,i,j,NHALO,jSize)];
               }

               derivCFD2nd(iSize, dx, cfdArr, cfdHalo0, cfdHalo1, cfdRhs, cfdSol,
                    boundaryTypeX, iRank, npex);

               for(size_t i=0;i<iSize;i++)
               {
                    arrDeriv[IDXIJ(ivar,i,j)] = cfdSol[i];
               }
          }        
     }
}


void derivY(int nv, real *arrMain, real *arrDeriv, int iVarStartMain)
{
     real realPart,imagPart;

     for(int ivar=0; ivar<nv;ivar++)
     {
          for(size_t i=0;i<iSize;i++)
          {
               // copy the variable to fftwIn
               for(size_t j=0;j<ny;j++)
               {
                    fftwIn[j] = arrMain[IDXIJ(ivar+iVarStartMain,i,j)];
               }

               fftw_execute(fftwPlanR2C);

               for(size_t ik=0; ik<nky; ik++)
               {
                    realPart = fftwOut[ik][0]; 
                    imagPart = fftwOut[ik][1];

                    //normalization by dividing n
                    fftwOut[ik][0] = -imagPart * 2 * PI / Ly * ik  / ny;
                    fftwOut[ik][1] = realPart * 2 * PI / Ly * ik  / ny;
               }

               fftw_execute(fftwPlanC2R);

               for(size_t j=0;j<ny;j++)
               {
                    arrDeriv[IDXIJ(ivar,i,j)] = fftwIn[j];
               }
          }   
     }
}



// 2nd order derivative for ux,uy,uz
void derivYY(int nv, real *arrMain, real *arrDeriv, int iVarStartMain) 
{
     // note: iVarStartMain: sometimes we may need to calculate derivatives for certain
     // variables, e.g. ux,uy,uz, thus we may need to start from iVar = 1
     real realPart,imagPart;

     for(int ivar=0; ivar<nv;ivar++)
     {
          for(size_t i=0;i<iSize;i++)
          {
               // copy the variable to fftwIn
               for(size_t j=0;j<ny;j++)
               {
                    fftwIn[j] = arrMain[IDXIJ(iVarStartMain + ivar,i,j)]; // take ux,uy,uz
               }

               fftw_execute(fftwPlanR2C);

               for(size_t ik=0; ik<nky; ik++)
               {
                    realPart = fftwOut[ik][0]; 
                    imagPart = fftwOut[ik][1];

                    //normalization by dividing n
                    fftwOut[ik][0] = -realPart * pow(2 * PI / Ly * ik, 2)  / ny;
                    fftwOut[ik][1] = -imagPart * pow(2 * PI / Ly * ik, 2)  / ny;
               }

               fftw_execute(fftwPlanC2R);

               for(size_t j=0;j<ny;j++)
               {
                    arrDeriv[IDXIJ(ivar,i,j)] = fftwIn[j];
               }
          }   
     }
}