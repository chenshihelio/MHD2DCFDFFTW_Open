// 3-step explicit RK

#include "parameter.h"
#include "initialize.h"
#include "macros.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "divBcleaning.h"

#define CC10 8.0/15.0
#define CC20 5.0/12.0
#define CC30 0.75

#define DD20 -17./60.
#define DD30 -5./12.

real cc1[3], dd1[3];

void rkInit(real dt)
{
    cc1[0] = CC10*dt ; dd1[0] = 0.0;
    cc1[1] = CC20*dt ; dd1[1] = DD20*dt;
    cc1[2] = CC30*dt ; dd1[2] = DD30*dt;
    
    memset(dudtRK, 0, arrSize * sizeof(real));
}

void rkExecute(int iRK)
{
    for(int i=0; i<arrSize; i++)
    {
        uu[i] = uu[i] + cc1[iRK] * dudt[i] + dd1[iRK] * dudtRK[i];
        dudtRK[i] = dudt[i]; 
    }
}

real calcDt(real dt0)
{
    real dt, tempMax;

    real csound2, calfven2, cms2;
    real calfvenx,calfveny,calfvenz,cfastx,cslowx,cnsx2,
        cfasty,cslowy,cnsy2,cfastz,cslowz,cnsz2,
        ux,cmaxx,uy,cmaxy,uz,cmaxz,cmaxHall;
    real dtx, dty, dtz, dtmin, dtmin_iter;

    dtmin = -1;

    for(int i=0; i<iSize;i++)
    {
        for(int j=0; j<jSize;j++)
        {
                //sound speed squared = kappa * p / rho
            csound2 = adiabaticIdx * uu[P(i,j)] / uu[RHO(i,j)];

            // alfven wave speeds
            calfvenx =  uu[BX(i,j)] / sqrt(uu[RHO(i,j)]);
            calfveny =  uu[BY(i,j)] / sqrt(uu[RHO(i,j)]);
            calfvenz =  uu[BZ(i,j)] / sqrt(uu[RHO(i,j)]);
            
            calfven2 = calfvenx * calfvenx + calfveny * calfveny + calfvenz * calfvenz;
            cms2 = csound2 + calfven2;

            cnsx2 = sqrt(MAX( (cms2*cms2 - 4 * csound2 * calfvenx*calfvenx) , 0.));
            cnsy2 = sqrt(MAX( (cms2*cms2 - 4 * csound2 * calfveny*calfveny) , 0.));
            // cnsz2 = sqrt(MAX( (cms2*cms2 - 4 * csound2 * calfvenz*calfvenz) , 0.));

            cfastx = sqrt(cms2 + cnsx2) / sqrt(2.);
            cfasty = sqrt(cms2 + cnsy2) / sqrt(2.);
            // cfastz = sqrt(cms2 + cnsz2) / sqrt(2.);

            cslowx = sqrt( MAX( (cms2 - cnsx2) ,0.) ) / sqrt(2.);
            cslowy = sqrt( MAX( (cms2 - cnsy2) ,0.) ) / sqrt(2.);
            // cslowz = sqrt( MAX( (cms2 - cnsz2) ,0.) ) / sqrt(2.);

            ux = uu[UX(i,j)];
            uy = uu[UY(i,j)];
            // uz = uu[UZ(i,j)];

            // get the maximum speed along each direction
            tempMax = MAX( fabs(ux), fabs(ux - calfvenx) );
            tempMax = MAX( fabs(ux - cslowx), tempMax);
            tempMax = MAX( fabs(ux - cfastx), tempMax);
            tempMax = MAX( fabs(ux + calfvenx), tempMax);
            tempMax = MAX( fabs(ux + cslowx), tempMax);
            cmaxx = MAX( fabs(ux + cfastx), tempMax);

            tempMax = MAX( fabs(uy), fabs(uy - calfveny) );
            tempMax = MAX( fabs(uy - cslowy), tempMax);
            tempMax = MAX( fabs(uy - cfasty), tempMax);
            tempMax = MAX( fabs(uy + calfveny), tempMax);
            tempMax = MAX( fabs(uy + cslowy), tempMax);
            cmaxy = MAX( fabs(uy + cfasty), tempMax);


            if (Hall==1)
            {
                cmaxHall = di / (MIN(dx,dy)) / uu[RHO(i,j)]  * 
                    ( MAX(fabs(uu[BX(i,j)]) ,
                        ( MAX( fabs(uu[BY(i,j)]) , fabs(uu[BZ(i,j)]) ) ) ) );

                cmaxx = MAX(cmaxx,cmaxHall);
                cmaxy = MAX(cmaxy,cmaxHall);
            }

            if(divBcleaning==2)
            {
                // Ch = max_speed
                Ch_divBcleaning = MAX(cmaxx,cmaxy);
                calc_Cp_fix_Cr();
            }

            dtx = dx / cmaxx ;
            dty = dy / cmaxy ;

            if(resistivity > 0)
            {
                dtx = MIN(dtx, dx*dx/resistivity);
                dty = MIN(dty, dy*dy/resistivity);
            }

            if(viscosity > 0)
            {
                dtx = MIN(dtx, dx*dx/viscosity);
                dty = MIN(dty, dy*dy/viscosity);
            }

            dtmin_iter = MIN(dtx,dty);

            if ( dtmin < 0 || dtmin_iter < dtmin )
            {
                dtmin = dtmin_iter;
            }
        }
    }

    dtmin_iter = dtmin;
    MPI_Allreduce(&dtmin_iter, &dtmin, 1,mpi_real,MPI_MIN,
                MPI_COMM_WORLD);

    dtmin = dtmin * cfl;

    if (dt0<0.98*dtmin || dt0>1.02*dtmin)
    {
        dt = dtmin;
    }
    else
    {
        dt = dt0;
    }
    return dt;
}