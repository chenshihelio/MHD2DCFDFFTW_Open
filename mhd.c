#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "macros.h"
#include "parameter.h"
#include "initialize.h"
#include "writeOutput.h"
#include "usrInitialize.h"
#include "timeAdvance.h"
#include "parallel.h"
#include "calcDerivs.h"
#include "calcRHS.h"
#include "filtering.h"


void initialLog();
void writeLog(real timeSim, real dt, real diffTime, long iStep);
real checkDivB();

int main(int argc, char** argv)
{
    real timeSim = 0, dt = 0, tOut = 0, dtLog = 0, tLog = 0;
    int iOutput = 0;
    long iStep = 0;
    time_t tStart, tRunTime;
    real diffTime;
    real maxDivB;
    
    time(&tStart);

    // initialize the code: allocate arrays and FFTW plans
    initialize();

    // write the log
    initialLog();
    
    // Output the grid file
    writeGrid();

    // User-defined initialization function
    if (myRank==0)
    {
        printf("usrInitialize......\n");
    }
    usrInitialize();

    time(&tRunTime);
    diffTime = difftime(tRunTime, tStart);
    if(myRank==0)
    {
        printf("Initializing takes %.2f seconds......\n", diffTime);
    }
    // Initialize tStart again after all initialization modules
    time(&tStart); 

    // Calculate dt
    dt = calcDt(dt);
    if(myRank==0)
    {
        printf("  dt = %e\n", dt);
    }
    rkInit(dt); // initialize the Runge-Kutta method


    dtLog = dtOut / 10.0;

    // initialize the times
    timeSim = 0;
    iOutput = 0;
    tOut = timeSim ;
    tLog = timeSim ;

    // Output the main array at t = 0
    maxDivB = checkDivB();
    if(myRank==0)
        printf("  Output at time = %10.4f, max(div B) = %12.4e, dt = %12.4e\n", 
            timeSim, maxDivB, dt);

    writeArray(iOutput, timeSim);  
    if(outDivB==1)
    {
        writeDivB(iOutput, timeSim);
    }
    iOutput++;
    tOut += dtOut;
    
    // Output the log file at t=0
    if (myRank==0)
    {
        time(&tRunTime);
        diffTime = difftime(tRunTime, tStart);
        writeLog(timeSim, dt, diffTime, iStep);
    }
    tLog += dtLog;
    

    // Main loop -------------------------------
    while (timeSim <= Tmax)
    {
        // output main array
        if (timeSim >= tOut)
        {
            maxDivB = checkDivB();
            if (myRank==0)
                printf("  Output at time = %10.4f, max(div B) = %12.4e, dt = %12.4e\n", 
                    timeSim, maxDivB, dt);

            writeArray(iOutput, timeSim);
            if(outDivB==1)
            {
                writeDivB(iOutput, timeSim);
            }
            iOutput++;
            tOut += dtOut;
        }

        //write log file
        if (timeSim >= tLog)
        {
            time(&tRunTime);
            diffTime = difftime(tRunTime, tStart);
            if(myRank==0)
            {
                writeLog(timeSim, dt, diffTime, iStep);
            }
            tLog += dtLog;
        }

        // dt too small, exit...
        if (dt < 1e-8)
        {
            if(myRank==0)
            {
                printf("****************************************************\n");
                printf("\n");
                printf("  STOP: dt < 1e-8\n");
                printf("  time = %10.4f, dt = %12.4e\n", timeSim, dt);
                printf("  iteration = %14ld\n", iStep);
                printf("\n");
                printf("****************************************************\n");
                printf("  Output at time = %10.4f, max(div B) = %12.4e, dt = %12.4e\n", 
                    timeSim, maxDivB, dt);
            }
            writeArray(iOutput, timeSim);
            if(outDivB==1)
            {
                writeDivB(iOutput, timeSim);
            }
            break;
        }

        // evolution module ------------
        for (int iRK = 0; iRK<3; iRK ++)
        {
            //exchange halo points of the main array
            sendrecvHaloPointsX();

            //calculate derivative of uu
            derivX(NVAR, uu, uuHaloX0, uuHaloX1, dudx, 0, 0);
            derivY(NVAR, uu, dudy, 0);

            if(viscosity>0)
            {
                derivXX(3, uu, uuHaloX0, uuHaloX1, d2udx2, 1, 1);
                derivYY(3, uu, d2udy2, 1);
            }

            if((Hall==1) || (resistivity>0))
            {
                derivXX(3, uu, uuHaloX0, uuHaloX1, d2Bdx2, 4, 4); 
                derivYY(3, uu, d2Bdy2, 4);
            }

            if(divBcleaning == 1)
            {
                calcDivB();

                sendrecvHaloPointsdivBX();
                derivX(1, divB, divBHaloX0, divBHaloX1, ddivBdx,0,0);
                derivY(1, divB, ddivBdy, 0);
            }


            // calculate Efield and derivative (for magnetic field evolution)
            calcEfield();
            sendrecvHaloPointsEfieldX();
            derivX(3, Efield, EfieldHaloX0, EfieldHaloX1, dEdx,0,0);
            derivY(3, Efield, dEdy, 0);


            calcRHS(); // calculate RHS

            rkExecute(iRK); //Runge-Kutta method 

            //filtering
            sendrecvHaloPointsX();
            filtering(); 
        }

        timeSim+= dt;
        iStep++;

        // update dt
        dt = calcDt(dt);
        rkInit(dt);
        // -----------------------------
    }

    
    // after Tmax is reached, output one file and write Log
    maxDivB = checkDivB();
    if (myRank==0)
        printf("  Output at time = %10.4f, max(div B) = %12.4e, dt = %12.4e\n", 
            timeSim, maxDivB, dt);

    writeArray(iOutput, timeSim);
    if(outDivB==1)
    {
        writeDivB(iOutput, timeSim);
    }
    iOutput++;
    tOut += dtOut;

    if (myRank==0)
    {
        time(&tRunTime);
        diffTime = difftime(tRunTime, tStart);
        writeLog(timeSim, dt, diffTime, iStep);
    }
    
    tLog += dtLog;
    // -------------------------------------------------------------------


    // free the allocated memory and destroy FFTW plans
    finalize();

    if (myRank==0)
    {
        //tRunTime = time(NULL);
        time(&tRunTime);

        printf("Simulation end... Total time = %.2f sec.\n", 
            difftime(tRunTime,tStart));
    }
    
    return 0;
}


//------------------------------------------------------------

void initialLog()
{
    if (myRank==0)
    {
        printf("****************************************************\n");

        if (boundaryTypeX==0)
        {
            printf("  Boundary Condition X: Periodic.\n");
        }
        else if(boundaryTypeX==1)
        {
            printf("  Boundary Condition X: Open.\n");
        }
        else
        {
            printf("  Boundary Condition X Unknown!!!\nExisting......\n");
            exit(0);
        }

        printf("  Boundary Condition Y: Periodic.\n");

        
        printf("  npe = %5d\n",npe);
        printf("  cfl = %f\n", cfl);
        if(filteringOptionY==0)
        {
            printf("  Use cut-off method for dealiasing along Y.\n");
        }
        else if(filteringOptionY==1)
        {
            printf("  Use smooth function for dealiasing along Y. afy = %.4f.\n",
                afy);
        }
        printf("  Tmax = %f, dtOut = %f\n", Tmax, dtOut);
        printf("  (nx, ny) = (%ld, %ld)\n", nx, ny);
        printf("  Lx = %f, Ly = %f, dx = %e, dy = %e\n", Lx, Ly, dx, dy);
        printf("  adiabatic Index = %f\n", adiabaticIdx);
        if(Hall==0)
        {
            printf("  Hall term is OFF.\n");
        }
        else
        {
            printf("  Hall term is ON. di = %f\n", di);
        }
        printf("  resistivity = %e, viscosity = %e\n", resistivity, viscosity);

        if(divBcleaning==0)
        {
            printf("  DivB cleaning is off.\n");
        }
        else if (divBcleaning==1)
        {
            printf("  Use diffusion method, Cd = %.2e to clean DivB.\n", constDiffDivB);
        }
        else if (divBcleaning==2)
        {
            printf("  Use hyperbolic method to clean divB.\n");
        }
        
        printf("****************************************************\n");
    }
}

void writeLog(real timeSim, real dt, real diffTime, long iStep)
{
    long seconds = (long)diffTime;
    int hour, minute, second;

    hour = seconds/3600;
    minute = (seconds % 3600) / 60 ;
    second = seconds % 60;

    FILE *fpWrite = fopen("log","w");

    fprintf(fpWrite, "  Simulation time = %f\n", timeSim);
    fprintf(fpWrite, "  dt = %12.4e\n", dt);
    fprintf(fpWrite, "  Real time used (sec) = %12ld\n", seconds);
    fprintf(fpWrite, "  Real time used (hh:mm:ss) = (%03d:%02d:%02d)\n",
        hour,minute,second);
    fprintf(fpWrite, "  Iteration = %14ld\n", iStep);
    fclose(fpWrite);
    
}



real checkDivB()
{
    real maxDivB = 0;
    real tmp;

    for(size_t i=0;i<iSize;i++)
    {
        for(size_t j=0;j<jSize;j++)
        {
            divB[IDXIJ(0,i,j)] = dudx[BX(i,j)] + dudy[BY(i,j)];

            tmp = fabs(divB[IDXIJ(0,i,j)]);

            if(tmp > maxDivB )
            {
                maxDivB = tmp;
            }
        }
    }

    tmp = maxDivB;
    MPI_Allreduce(&tmp, &maxDivB, 1,mpi_real,MPI_MAX,
                MPI_COMM_WORLD);

    return maxDivB;
}