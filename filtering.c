#include "parameter.h"
#include "initialize.h"
#include "mathlib.h"
#include "parallel.h"
#include "macros.h"



void filteringX(int nv, real *arrMain, real *arrHalo0, real *arrHalo1)
{
    int npex = npe;
    for(int ivar=0; ivar<nv;ivar++)
    {

        for(size_t j=0; j<jSize;j++)
        {  
            for(size_t i=0;i<iSize;i++)
            {
                    cfdArr[i] = arrMain[IDXIJ(ivar,i,j)];
            }

            for(size_t i=0;i<NHALO;i++)
            {
                    cfdHalo0[i] = arrHalo0[IDX(ivar,i,j,NHALO,jSize)];
                    cfdHalo1[i] = arrHalo1[IDX(ivar,i,j,NHALO,jSize)];
            }

            filterCFD(iSize, cfdArr, cfdHalo0, cfdHalo1, cfdRhs, cfdSol,
                    boundaryTypeX, iRank, npex);

            for(size_t i=0;i<iSize;i++)
            {
                    arrMain[IDXIJ(ivar,i,j)] = cfdSol[i];
            }
        }        
    }
}


void filteringY()
{
    real realPart,imagPart;

    for(int ivar=0; ivar<NVAR;ivar++)
    {
        for(int i=0;i<iSize;i++)
        {
            // copy the variable to fftwIn
            for(int j=0;j<ny;j++)
            {
                fftwIn[j] = uu[IDXIJ(ivar,i,j)];
            }

            // forward FFT
            fftw_execute(fftwPlanR2C);


            for(int ik=0; ik<nky; ik++)
            {
                if (filteringOptionY==0)
                {
                
                    if(ik<=floor(ny/3))
                    {
                        realPart = fftwOut[ik][0]; 
                        imagPart = fftwOut[ik][1];

                        //normalization by dividing n
                        fftwOut[ik][0] = realPart / ny;
                        fftwOut[ik][1] = imagPart / ny;    
                    }
                    else
                    {
                        fftwOut[ik][0] = 0;
                        fftwOut[ik][1] = 0;
                    }
                }
                else if (filteringOptionY==1)
                {
                    fftwOut[ik][0] = fftwOut[ik][0] / ny * filtY[ik];
                    fftwOut[ik][1] = fftwOut[ik][1] / ny * filtY[ik];
                }
            }
            
            // inverse FFT
            fftw_execute(fftwPlanC2R);

            for(int j=0;j<ny;j++)
            {
                uu[IDXIJ(ivar,i,j)] = fftwIn[j];
            }
        }  
    }
}



void filtering()
{
    filteringX(NVAR, uu, uuHaloX0, uuHaloX1);
    filteringY();
}