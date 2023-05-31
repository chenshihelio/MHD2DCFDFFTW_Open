// This should be the only .c file that is modified
// by the user to define the initial conditions
#include "parameter.h"
#include "initialize.h"
#include "macros.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>



// test cases------------------------------------------------
void usrInitialize()
{
    int testCase = 5;

    if(iRank==0)
    {
        printf("Running test case #%d......\n", testCase);
    }


    real rho0 = 1;
    real P0 = 1;
    real Bx0 = 1;
    real By0 = 1;

    int Nmax;
    real db,du;
    real xcent,a0,xlow,xup;
    real *randPhase;
    int iMode;

    real deltaP; 
    real Uy_shear;



    switch (testCase)
    {
    case 1: //Add circular Alfvenic wave along x
        Nmax = 1;
        db = 0.5;
        iMode = 16;


        // generate the random phase at Process #0
        randPhase = malloc(sizeof(real) * Nmax);
        memset(randPhase,0,sizeof(real) * Nmax);

        if (iRank==0)
        {
            for(int i=0;i<Nmax;i++)
            {
                randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
            }
        }
        // broadcast the random phase to all processors
        MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        for(size_t i=0;i<iSize;i++)
        {
            for(size_t j=0;j<jSize;j++)
            {
                uu[RHO(i,j)] = rho0;
                uu[UX(i,j)] = 0.0;

                uu[UY(i,j)] = uu[UY(i,j)] + db * cos( 
                    (double)iMode*2*PI/Lx*xgrid[i] + randPhase[iMode]);
                
                uu[UZ(i,j)] = uu[UZ(i,j)] + db * sin( 
                    (double)iMode*2*PI/Lx*xgrid[i] + randPhase[iMode]);

                uu[BY(i,j)] = -uu[UY(i,j)];
                uu[BZ(i,j)] = -uu[UZ(i,j)];
                
                uu[BX(i,j)] = Bx0;
                uu[P(i,j)] = P0;
            }
        }

        free(randPhase);

        break;
    
    case 2: // Add linear Alfvenic wave band along x
        Nmax = 4;
        db = 1e-1;

        randPhase = malloc(sizeof(real) * Nmax);
        memset(randPhase,0,sizeof(real) * Nmax);

        if (iRank==0)
        {
            for(int i=0;i<Nmax;i++)
            {
                randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
            }
        }

        MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);


        for(size_t i=0;i<iSize;i++)
        {
            for(size_t j=0;j<jSize;j++)
            {
                uu[RHO(i,j)] = rho0;
                uu[UX(i,j)] = 0.0;

                for(int M=0; M<Nmax;M++)
                {
                    uu[UY(i,j)] = uu[UY(i,j)] + db * cos( 
                        (double)(M+1)*2*PI/Lx*xgrid[i] + randPhase[M]);
                    
                    // uu[UZ(i,j)] = uu[UZ(i,j)] + db * sin( 
                    //     (double)(M+1)*2*PI/Lx*xgrid[i] + randPhase[M]);
                    uu[UZ(i,j)] = 0;
                }

                uu[BY(i,j)] = -uu[UY(i,j)];
                uu[BZ(i,j)] = -uu[UZ(i,j)];
                
                uu[BX(i,j)] = Bx0;
                uu[P(i,j)] = P0;
            }
        }

        free(randPhase);
        break;

    case 3: // add a pressure peak
        deltaP = 2.0;
        By0 = 3.0;
        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                uu[RHO(i,j)] = rho0;
                uu[UX(i,j)] = 0;
                uu[UY(i,j)] = 0;
                uu[UZ(i,j)] = 0;
                uu[BX(i,j)] = 0;
                uu[BY(i,j)] = By0;
                uu[BZ(i,j)] = 0;

                uu[P(i,j)] = P0 + deltaP * exp( - 
                    pow((xgrid[i] - Lx/2) / (Lx/8), 2) - 
                    pow((ygrid[j] - Ly/2) / (Ly/8), 2));
            }
        }
        break;

    case 4: // add shear flow - needs open boundary condition
        By0 = 0;
        Nmax = 32;
        Uy_shear = 1.0;
        du = 1e-3;

        xcent = Lx/2;
        a0 = Lx/40;

        // generate random phases 
        randPhase = malloc(sizeof(real) * Nmax);
        memset(randPhase,0,sizeof(real)*Nmax);

        if (myRank==0)
        {
            for(int i=0;i<Nmax;i++)
            {
                randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
            }
        }
        MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);


        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                // background field
                uu[RHO(i,j)] = rho0;
                uu[UX(i,j)] = 0;
                uu[UY(i,j)] = Uy_shear * tanh((xgrid[i]-xcent)/a0);
                uu[UZ(i,j)] = 0;
                uu[BX(i,j)] = 0;
                uu[BY(i,j)] = By0;
                uu[BZ(i,j)] = 0;
                uu[P(i,j)] = P0;

                // add fluctuations in U
                for(int ik=0;ik<Nmax;ik++)
                {
                    uu[UY(i,j)] = uu[UY(i,j)] + du * 
                        cos((double)(ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xcent)/a0 ,2))
                        * (-2 * (xgrid[i]-xcent)/a0/a0 );
                    uu[UX(i,j)] = uu[UX(i,j)] + du * ((ik+1)*2*PI/Ly)
                        * cos((double)(ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xcent)/a0 ,2));
                }

            }
        }
        
        free(randPhase);
        break;

    case 41: // add shear flow - double layer -- periodic BC
        By0 = 1.0;
        Nmax = 32;
        Uy_shear = 1.0;
        du = 1e-3;

        xlow = 0.25 * Lx;
        xup = 0.75 * Lx;
        xcent = 0.5 * Lx;
        a0 = Lx/60;


        // add background field
        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                
                uu[RHO(i,j)] = rho0;
                uu[UX(i,j)] = 0;
                if(xgrid[i]<xcent)
                {
                    uu[UY(i,j)] = Uy_shear * tanh((xgrid[i]-xlow)/a0);
                }
                else
                {
                    uu[UY(i,j)] = -Uy_shear * tanh((xgrid[i]-xup)/a0);
                }
                uu[UZ(i,j)] = 0;
                uu[BX(i,j)] = 0;
                uu[BY(i,j)] = By0;
                uu[BZ(i,j)] = 0;
                uu[P(i,j)] = P0;
            }
        }



        // add fluctuations in U
        
        randPhase = malloc(sizeof(real) * Nmax);
        memset(randPhase,0,sizeof(real)*Nmax);

        // upper layer----
        // generate random phases 
        if (myRank==0)
        {
            for(int i=0;i<Nmax;i++)
            {
                randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
            }
        }
        MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);

        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                for(int ik=0;ik<Nmax;ik++)
                {
                    uu[UY(i,j)] = uu[UY(i,j)] + du * 
                        cos((double)(ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xlow)/a0 ,2))
                        * (-2 * (xgrid[i]-xlow)/a0/a0 );
                    uu[UX(i,j)] = uu[UX(i,j)] + du * ((ik+1)*2*PI/Ly)
                        * cos((double)(ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xlow)/a0 ,2));
                }
            }
        }



        // lower layer----
        // generate random phases 
        if (myRank==0)
        {
            for(int i=0;i<Nmax;i++)
            {
                randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
            }
        }
        MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);

        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                for(int ik=0;ik<Nmax;ik++)
                {
                    uu[UY(i,j)] = uu[UY(i,j)] + du * 
                        cos((double)(ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xup)/a0 ,2))
                        * (-2 * (xgrid[i]-xup)/a0/a0 );
                    uu[UX(i,j)] = uu[UX(i,j)] + du * ((ik+1)*2*PI/Ly)
                        * cos((double)(ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xup)/a0 ,2));
                }
            }
        }

        
        free(randPhase);
        break;

    case 5: // a current sheet - needs open boundary 
        rho0 = 1.0;
        By0 = 1.0;
        P0 = 1.0;

        xcent = Lx/2;
        a0 = Lx/80;

        // add the background structure
        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {

                uu[BY(i,j)] = By0 * tanh((xgrid[i]-xcent)/a0);
                // use density to balance
                uu[RHO(i,j)] = rho0 * (1 + By0*By0/2/P0/ pow( 
                    cosh( (xgrid[i]-xcent)/a0  ) ,2) );

                uu[UX(i,j)] = 0; 
                uu[UY(i,j)] = 0; 
                uu[UZ(i,j)] = 0;
                uu[BX(i,j)] = 0; 
                uu[BZ(i,j)] = 0; 
                uu[P(i,j)] = P0;
            }
        }
        

        // add perturbations
        Nmax = 32;
        db = 1e-3;

        randPhase = malloc(sizeof(real) * Nmax);
        memset(randPhase,0,sizeof(real) * Nmax);
        if (myRank==0)
        {
            for(int i=0;i<Nmax;i++)
            {
                randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
            }
        }

        MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);

        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                for(int ik=0;ik<Nmax;ik++)
                {
                    uu[BY(i,j)] = uu[BY(i,j)] + db * 
                        cos((ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xcent)/a0 ,2))
                        * (-2 * (xgrid[i]-xcent)/a0/a0 );
                    uu[BX(i,j)] = uu[BX(i,j)] + db * ((ik+1)*2*PI/Ly)
                        * cos((ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
                        * 0.5 * exp(- pow( (xgrid[i]-xcent)/a0 ,2));
                }
            }
        }

        free(randPhase);
        break;


    default:
        for(int i=0;i<iSize;i++)
        {
            for(int j=0;j<jSize;j++)
            {
                uu[RHO(i,j)] = 1;
                uu[UX(i,j)] = 0;
                uu[UY(i,j)] = 0;
                uu[UZ(i,j)] = 0;
                uu[BX(i,j)] = 0;
                uu[BY(i,j)] = 0;
                uu[BZ(i,j)] = 0;
                uu[P(i,j)] = 1;
            }
        }


        break;
    }
}



// // ////////////////////////////////////////////////////////////////////////
// // add a jet along y-direction with Alfven waves on top
// void usrInitialize()
// {
//     real jet_width = 1.0;

//     for(int i=0;i<iSize;i++)
//     {
//         for(int j=0;j<jSize;j++)
//         {
//             uu[RHO(i,j)] = 1.0;
//             uu[UX(i,j)] = 0;
//             uu[UZ(i,j)] = 0;
//             uu[UY(i,j)] = 5.0 * exp(- 0.5 *  pow((xgrid[i] - 0.5*Lx)/(jet_width),2));
//             uu[BY(i,j)] = 1.0;
//             uu[P(i,j)] = 5.0;
//         }
//     }    


//     // add alfven waves 
//     int Nmax = 1;
//     real db = 0.5;

//     real randPhase[Nmax];


//     if (myRank==0)
//     {
//         for(int i=0;i<Nmax;i++)
//         {
//             randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//         }
//     }

//     MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);


//     for(int i=0;i<iSize;i++)
//     {
//         for(int j=0;j<jSize;j++)
//         {
//             for(int M=0; M<Nmax;M++)
//             {
//                 uu[UX(i,j)] = uu[UX(i,j)] + db * cos( 
//                     (M+1)*2*PI/Ly*ygrid[j] + randPhase[M]);
                
//                 uu[UZ(i,j)] = uu[UZ(i,j)] + db * sin( 
//                     (M+1)*2*PI/Ly*ygrid[j] + randPhase[M]);
//             }

//             uu[BX(i,j)] = -uu[UX(i,j)];
//             uu[BZ(i,j)] = -uu[UZ(i,j)];
//         }
//     }


//     // // add ux fluctuations
    
//     // int Nmax = 16;
//     // real dux = 0.1;

//     // real randPhase[Nmax];


//     // if (myRank==0)
//     // {
//     //     for(int i=0;i<Nmax;i++)
//     //     {
//     //         randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//     //     }
//     // }

//     // MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);


//     // for(int i=0;i<iSize;i++)
//     // {
//     //     for(int j=0;j<jSize;j++)
//     //     {
//     //         for(int M=0; M<Nmax;M++)
//     //         {
//     //             uu[UX(i,j)] = uu[UX(i,j)] + dux * cos( 
//     //                 (M+1)*2*PI/Ly*ygrid[j] + randPhase[M]) *
//     //                 0.5 * (exp(-0.5 * pow((xgrid[i]-0.5*Lx-jet_width)/(0.2*jet_width),2)) 
//     //                 + exp(-0.5 * pow((xgrid[i]-0.5*Lx+jet_width)/(0.2*jet_width),2)) );
                
//     //         }

//     //         uu[BX(i,j)] = 0;
//     //         uu[BY(i,j)] = 0;
//     //         uu[BZ(i,j)] = 0;
//     //     }
//     // }

// }
// ////////////////////////////////////////////////////////////////////////




// ////////////////////////////////////////////////////////////////////////
// // add a 1D pressure peak
// void usrInitialize()
// {
//     for(int i=0;i<iSize;i++)
//     {
//         for(int j=0;j<jSize;j++)
//         {
//             uu[RHO(i,j)] = 1.0;
//             uu[UX(i,j)] = 0;
//             uu[UY(i,j)] = 0;
//             uu[UZ(i,j)] = 0;
//             uu[BX(i,j)] = 0.0;
//             uu[BY(i,j)] = 0.0;
//             uu[BZ(i,j)] = 0;

//             uu[P(i,j)] = 1.0 + 3.0 * exp( - 
//                 pow((xgrid[i] - Lx/2) / (Lx/8), 2));
//             // uu[P(i,j)] = 1.0 + 1.0 * exp( - 
//             //     pow((ygrid[j] - Ly/2) / (Ly/8), 2));
//         }
//     }
// }






// // add double-Harris current sheet 
// void usrInitialize()
// {
//     real ninf,B0,T0,a0;
//     real ylow,yup,ycent;

//     ninf = 1.0;
//     B0 = 1.0;
//     a0 = 1.0;
//     T0 = 1.0;

//     ylow = 0.25 * Ly;
//     yup = 0.75 * Ly;
//     ycent = 0.5 * Ly ;


//     // add the background structure
//     for(int i=0;i<iSize;i++)
//     {
//         for(int j=0;j<jSize;j++)
//         {
//             if(ygrid[j]<ycent)
//             {
//                 uu[BX(i,j)] = B0 * tanh((ygrid[j]-ylow)/a0);
//                 uu[BZ(i,j)] = B0 / cosh((ygrid[j]-ylow)/a0);


//                 // use density to balance
//                 // uu[RHO(i,j)] = ninf + B0*B0/2/T0/ pow( 
//                 //     cosh( (ygrid[j]-ylow)/a0  ) ,2);
//             }
//             else
//             {
//                 uu[BX(i,j)] = -B0 * tanh((ygrid[j]-yup)/a0);
//                 uu[BZ(i,j)] = B0 / cosh((ygrid[j]-yup)/a0);

//                 // uu[RHO(i,j)] = ninf + B0*B0/2/T0/ pow( 
//                 //     cosh( (ygrid[j]-yup)/a0  ) ,2);
//             }
//             uu[UX(i,j)] = 0; uu[UY(i,j)] = 0; uu[UZ(i,j)] = 0;
//             uu[BY(i,j)] = 0; uu[RHO(i,j)] = ninf;
//             uu[P(i,j)] = uu[RHO(i,j)] * T0;
//         }
//     }
    

//     // add perturbations
//     int nk = 32;
//     real db = 1e-3;

//     real randPhase[nk];


//     if (myRank==0)
//     {
//         for(int i=0;i<nk;i++)
//         {
//             randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//         }
//     }

//     MPI_Bcast(randPhase, nk, mpi_real, 0, MPI_COMM_WORLD);

//     for(int ik=0;ik<nk;ik++)
//     {
//         for(int i=0;i<iSize;i++)
//         {
//             for(int j=0;j<jSize;j++)
//             {
//                 uu[BX(i,j)] = uu[BX(i,j)] + db * 
//                     cos((ik+1)*2*PI/Lx*xgrid[i] + randPhase[ik])
//                     * 0.5 * exp(- pow( (ygrid[j]-ylow)/a0 ,2))
//                     * (-2 * (ygrid[j]-ylow)/a0/a0 );
//                 uu[BY(i,j)] = uu[BY(i,j)] + db * ((ik+1)*2*PI/Lx)
//                     * cos((ik+1)*2*PI/Lx*xgrid[i] + randPhase[ik])
//                     * 0.5 * exp(- pow( (ygrid[j]-ylow)/a0 ,2));
//             }
//         }
//     }


//     if (myRank==0)
//     {
//         for(int i=0;i<nk;i++)
//         {
//             randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//         }
//     }

//     MPI_Bcast(randPhase, nk, mpi_real, 0, MPI_COMM_WORLD);

//     for(int ik=0;ik<nk;ik++)
//     {
//         for(int i=0;i<iSize;i++)
//         {
//             for(int j=0;j<jSize;j++)
//             {
//                 uu[BX(i,j)] = uu[BX(i,j)] + db * 
//                     cos((ik+1)*2*PI/Lx*xgrid[i] + randPhase[ik])
//                     * 0.5 * exp(- pow( (ygrid[j]-yup)/a0 ,2))
//                     * (-2 * (ygrid[j]-yup)/a0/a0 );
//                 uu[BY(i,j)] = uu[BY(i,j)] + db * ((ik+1)*2*PI/Lx)
//                     * cos((ik+1)*2*PI/Lx*xgrid[i] + randPhase[ik])
//                     * 0.5 * exp(- pow( (ygrid[j]-yup)/a0 ,2));
//             }
//         }
//     }
// }








// ////////////////////////////////////////////////////////////////////////
// // Add Alfvenic wave band along y
// void usrInitialize()
// {
//     int Nmax = 12;
//     real db = 1e-2;

//     real randPhase[Nmax];


//     if (myRank==0)
//     {
//         for(int i=0;i<Nmax;i++)
//         {
//             randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//         }
//     }

//     MPI_Bcast(randPhase, Nmax, mpi_real, 0, MPI_COMM_WORLD);


//     for(int i=0;i<iSize;i++)
//     {
//         for(int j=0;j<jSize;j++)
//         {
//             uu[RHO(i,j)] = 1.0;


//             uu[UY(i,j)] = 0.0;

//             for(int M=0; M<Nmax;M++)
//             {
//                 // uu[UX(i,j)] = uu[UX(i,j)] + db * cos( 
//                 //     (M+1)*2*PI/Ly*ygrid[j] + randPhase[M]);
                
//                 uu[UZ(i,j)] = uu[UZ(i,j)] + db * sin( 
//                     (M+1)*2*PI/Ly*ygrid[j] + randPhase[M]);
//             }

//             uu[BX(i,j)] = -uu[UX(i,j)];
//             uu[BZ(i,j)] = -uu[UZ(i,j)];
            
//             uu[BY(i,j)] = 1.0;
//             uu[P(i,j)] = 1.0;
//         }
//     }
// }
// ////////////////////////////////////////////////////////////////////////



// // add double-Harris current sheet along y direction
// void usrInitialize()
// {
//     real ninf,B0,T0,a0;
//     real xlow,xup,xcent;

//     ninf = 1.0;
//     B0 = 1.0;
//     a0 = 0.02;
//     T0 = 1.0;

//     xlow = 0.25 * Lx;
//     xup = 0.75 * Lx;
//     xcent = 0.5 * Lx;


//     // add the background structure
//     for(int i=0;i<iSize;i++)
//     {
//         for(int j=0;j<jSize;j++)
//         {
//             if(xgrid[i]<xcent)
//             {
//                 uu[BY(i,j)] = B0 * tanh((xgrid[i]-xlow)/a0);
//                 // uu[BZ(i,j)] = B0 / cosh((xgrid[i]-xlow)/a0);


//                 // use density to balance
//                 uu[RHO(i,j)] = ninf + B0*B0/2/T0/ pow( 
//                     cosh( (xgrid[i]-xlow)/a0  ) ,2);
//             }
//             else
//             {
//                 uu[BY(i,j)] = -B0 * tanh((xgrid[i]-xup)/a0);
//                 // uu[BZ(i,j)] = B0 / cosh((xgrid[i]-xup)/a0);

//                 uu[RHO(i,j)] = ninf + B0*B0/2/T0/ pow( 
//                     cosh( (xgrid[i]-xup)/a0  ) ,2);
//             }
//             uu[UX(i,j)] = 0; uu[UY(i,j)] = 0; uu[UZ(i,j)] = 0;
//             uu[BX(i,j)] = 0; uu[BZ(i,j)] = 0;  // uu[RHO(i,j)] = ninf;
//             uu[P(i,j)] = uu[RHO(i,j)] * T0;
//         }
//     }
    

//     // add perturbations
//     int nk = 128;
//     real db = 1e-4;

//     real randPhase[nk];


//     if (myRank==0)
//     {
//         for(int i=0;i<nk;i++)
//         {
//             randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//         }
//     }

//     MPI_Bcast(randPhase, nk, mpi_real, 0, MPI_COMM_WORLD);

//     for(int ik=0;ik<nk;ik++)
//     {
//         for(int i=0;i<iSize;i++)
//         {
//             for(int j=0;j<jSize;j++)
//             {
//                 uu[BY(i,j)] = uu[BY(i,j)] + db * 
//                     cos((ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
//                     * 0.5 * exp(- pow( (xgrid[i]-xlow)/a0 ,2))
//                     * (-2 * (xgrid[i]-xlow)/a0/a0 );
//                 uu[BX(i,j)] = uu[BX(i,j)] + db * ((ik+1)*2*PI/Ly)
//                     * cos((ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
//                     * 0.5 * exp(- pow( (xgrid[i]-xlow)/a0 ,2));
//             }
//         }
//     }


//     if (myRank==0)
//     {
//         for(int i=0;i<nk;i++)
//         {
//             randPhase[i] = 2* PI * (real)rand() / ( RAND_MAX + 1.0 );
//         }
//     }

//     MPI_Bcast(randPhase, nk, mpi_real, 0, MPI_COMM_WORLD);

//     for(int ik=0;ik<nk;ik++)
//     {
//         for(int i=0;i<iSize;i++)
//         {
//             for(int j=0;j<jSize;j++)
//             {
//                 uu[BY(i,j)] = uu[BY(i,j)] + db * 
//                     cos((ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
//                     * 0.5 * exp(- pow( (xgrid[i]-xup)/a0 ,2))
//                     * (-2 * (xgrid[i]-xup)/a0/a0 );
//                 uu[BX(i,j)] = uu[BX(i,j)] + db * ((ik+1)*2*PI/Ly)
//                     * cos((ik+1)*2*PI/Ly*ygrid[j] + randPhase[ik])
//                     * 0.5 * exp(- pow( (xgrid[i]-xup)/a0 ,2));
//             }
//         }
//     }
// }
