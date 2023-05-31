#include "initialize.h"
#include "macros.h"

//TODO: Implement hyperbolic divB cleaning to open boundary condition,
//     right now it works for periodic boundary condition,
//     At open boundary, we still treat psi as it is inner point.

void calcCharacteristicsX(size_t i, size_t j, int boundary,real *dudtChar);
void calcRHSInnerPoint(size_t i0, size_t i1, size_t j0, size_t j1);
void calcRHSOpenBoundaryX(size_t i, size_t j0, size_t j1);
void addSourceTerm();

void calcDivB()
{
    for(size_t i=0;i<iSize;i++)
    {
        for(size_t j=0;j<jSize;j++)
        {
            divB[IDXIJ(0,i,j)] = dudx[BX(i,j)] + dudy[BY(i,j)];
        }
    }
}

void calcEfield()
{
    for(size_t i=0;i<iSize;i++)
    {
        for(size_t j=0;j<jSize;j++)
        {
            // Ex = uz * By - uy * bz
            Efield[IDXIJ(0,i,j)] = uu[UZ(i,j)] * uu[BY(i,j)]
                - uu[UY(i,j)] * uu[BZ(i,j)];

            // Ey = ux * Bz - uz * Bx
            Efield[IDXIJ(1,i,j)] = uu[UX(i,j)] * uu[BZ(i,j)]
                - uu[UZ(i,j)] * uu[BX(i,j)];

            // Ez = uy * Bx - ux * By
            Efield[IDXIJ(2,i,j)] = uu[UY(i,j)] * uu[BX(i,j)]
                - uu[UX(i,j)] * uu[BY(i,j)];
        }
    }
}


void calcRHS()
{
    real divU, divB, tmpJx, tmpJy, tmpJz;
    int npex;
    npex = npe;

    if (boundaryTypeX==0)
    {
        calcRHSInnerPoint(0,iSize,0,jSize);
    }
    else if (boundaryTypeX==1)
    {
        if ( (iRank != 0) && (iRank != npex-1))
        {
            calcRHSInnerPoint(0,iSize,0,jSize);
        }
        else if (iRank == 0)
        {
            calcRHSInnerPoint(1,iSize,0,jSize);
            calcRHSOpenBoundaryX(0,0,jSize);    
        }
        else 
        {
            calcRHSInnerPoint(0,iSize-1,0,jSize);
            calcRHSOpenBoundaryX(iSize-1,0,jSize);
        }
    }

    addSourceTerm();
}

void calcRHSInnerPoint(size_t i0, size_t i1, size_t j0, size_t j1)
{
    real divB_point, divU, currentX, currentY, currentZ;

    for(size_t i=i0; i<i1; i++)
    {
        for(size_t j=j0; j<j1; j++)
        {
            divU = dudx[UX(i,j)] + dudy[UY(i,j)];
            divB_point = dudx[BX(i,j)] + dudy[BY(i,j)];
            currentX = dudy[BZ(i,j)];
            currentY = -dudx[BZ(i,j)];
            currentZ = dudx[BY(i,j)] - dudy[BX(i,j)];

            dudt[RHO(i,j)] = - uu[RHO(i,j)] * divU
                - uu[UX(i,j)] * dudx[RHO(i,j)]
                - uu[UY(i,j)] * dudy[RHO(i,j)];

            dudt[UX(i,j)] = -uu[UX(i,j)] * dudx[UX(i,j)]
                -uu[UY(i,j)] * dudy[UX(i,j)] 
                - dudx[P(i,j)]/uu[RHO(i,j)]
                + currentY * uu[BZ(i,j)] / uu[RHO(i,j)]
                - currentZ * uu[BY(i,j)] / uu[RHO(i,j)];

            dudt[UY(i,j)] = -uu[UX(i,j)] * dudx[UY(i,j)]
                -uu[UY(i,j)] * dudy[UY(i,j)] 
                - dudy[P(i,j)]/uu[RHO(i,j)]
                + currentZ * uu[BX(i,j)] / uu[RHO(i,j)]
                - currentX * uu[BZ(i,j)] / uu[RHO(i,j)];

            dudt[UZ(i,j)] = -uu[UX(i,j)] * dudx[UZ(i,j)]
                -uu[UY(i,j)] * dudy[UZ(i,j)]
                + currentX * uu[BY(i,j)] / uu[RHO(i,j)]
                - currentY * uu[BX(i,j)] / uu[RHO(i,j)];

            
            // // expand all terms : do not need to calculate Efield and its gradients
            // dudt[BX(i,j)] = -uu[UX(i,j)] * dudx[BX(i,j)]
            //     -uu[UY(i,j)] * dudy[BX(i,j)] 
            //     +uu[BX(i,j)] * dudx[UX(i,j)] 
            //     +uu[BY(i,j)] * dudy[UX(i,j)]
            //     -divU * uu[BX(i,j)] + divB_point * uu[UX(i,j)]; 

            // dudt[BY(i,j)] = -uu[UX(i,j)] * dudx[BY(i,j)]
            //     -uu[UY(i,j)] * dudy[BY(i,j)] 
            //     +uu[BX(i,j)] * dudx[UY(i,j)] 
            //     +uu[BY(i,j)] * dudy[UY(i,j)]
            //     -divU * uu[BY(i,j)] + divB_point * uu[UY(i,j)];
            
            // dudt[BZ(i,j)] = -uu[UX(i,j)] * dudx[BZ(i,j)]
            //     -uu[UY(i,j)] * dudy[BZ(i,j)] 
            //     +uu[BX(i,j)] * dudx[UZ(i,j)] 
            //     +uu[BY(i,j)] * dudy[UZ(i,j)]
            //     -divU * uu[BZ(i,j)] + divB_point * uu[UZ(i,j)];
            


            // dB/dt = - curl(E) : smaller div(B) errors
           dudt[BX(i,j)] = -dEdy[EZ(i,j)];
           dudt[BY(i,j)] = dEdx[EZ(i,j)];
           dudt[BZ(i,j)] = dEdy[EX(i,j)] - dEdx[EY(i,j)];
           

            dudt[P(i,j)] = -uu[UX(i,j)] * dudx[P(i,j)]
                -uu[UY(i,j)] * dudy[P(i,j)]
                - adiabaticIdx * uu[P(i,j)] * divU;

            if (divBcleaning==2) 
            {
                // hyperbolic divB cleaning (Dedner+2002) (eq 24) 
                // TODO: Try (eq 36).
                dudt[UX(i,j)] = dudt[UX(i,j)] -  divB_point * uu[BX(i,j)];
                dudt[UY(i,j)] = dudt[UY(i,j)] -  divB_point * uu[BY(i,j)];
                dudt[UZ(i,j)] = dudt[UZ(i,j)] -  divB_point * uu[BZ(i,j)];


                dudt[BX(i,j)] = dudt[BX(i,j)] - dudx[PSI(i,j)];
                dudt[BY(i,j)] = dudt[BY(i,j)] - dudy[PSI(i,j)];

                dudt[PSI(i,j)] = - Ch_divBcleaning*Ch_divBcleaning * (divB_point 
                    + uu[PSI(i,j)]/Cp_divBcleaning/Cp_divBcleaning );
            }
        }
    }
}

void calcRHSOpenBoundaryX(size_t i, size_t j0, size_t j1)
{
    // terms for open boundary, i.e. no explicit d/dx terms
    real divB, divU, tmpJx, tmpJy, tmpJz;
    real dudtChar[NVAR],divBpoint;

    for(size_t j=j0;j<j1;j++)
    {
        // dBx/dx needs not to be decomposed ///////////////
        divU = dudx[UX(i,j)] + dudy[UY(i,j)];
        divB = dudx[BX(i,j)] + dudy[BY(i,j)];

        dudt[BX(i,j)] = -uu[UX(i,j)] * dudx[BX(i,j)]
                -uu[UY(i,j)] * dudy[BX(i,j)]
                +uu[BX(i,j)] * dudx[UX(i,j)] 
                +uu[BY(i,j)] * dudy[UX(i,j)]
                -divU * uu[BX(i,j)];
        ////////////////////////////////////////////////////

        // do not include any d/dx terms
        tmpJx = dudy[BZ(i,j)];
        tmpJy = 0;
        tmpJz = -dudy[BX(i,j)];

        divU = dudy[UY(i,j)];
        divB = dudy[BY(i,j)];

        dudt[P(i,j)] = -uu[UY(i,j)] * dudy[P(i,j)]
            - adiabaticIdx * uu[P(i,j)] * divU;

        dudt[RHO(i,j)] = - uu[RHO(i,j)] * divU
            - uu[UY(i,j)] * dudy[RHO(i,j)];

        dudt[UX(i,j)] = -uu[UY(i,j)] * dudy[UX(i,j)] 
            + tmpJy * uu[BZ(i,j)] / uu[RHO(i,j)]
            - tmpJz * uu[BY(i,j)] / uu[RHO(i,j)];

        dudt[UY(i,j)] = -uu[UY(i,j)] * dudy[UY(i,j)] 
            - dudy[P(i,j)]/uu[RHO(i,j)]
            + tmpJz * uu[BX(i,j)] / uu[RHO(i,j)]
            - tmpJx * uu[BZ(i,j)] / uu[RHO(i,j)];

        dudt[UZ(i,j)] = -uu[UY(i,j)] * dudy[UZ(i,j)]
            + tmpJx * uu[BY(i,j)] / uu[RHO(i,j)]
            - tmpJy * uu[BX(i,j)] / uu[RHO(i,j)];

        dudt[BY(i,j)] = -uu[UY(i,j)] * dudy[BY(i,j)] 
            +uu[BY(i,j)] * dudy[UY(i,j)]
            -divU * uu[BY(i,j)];
        
        dudt[BZ(i,j)] = -uu[UY(i,j)] * dudy[BZ(i,j)] 
            +uu[BY(i,j)] * dudy[UZ(i,j)]
            -divU * uu[BZ(i,j)];

        if(i==0)
        {
            calcCharacteristicsX(i, j, 0, dudtChar);
            dudt[RHO(i,j)] = dudt[RHO(i,j)] + dudtChar[0];
            dudt[UX(i,j)] = dudt[UX(i,j)] + dudtChar[1];
            dudt[UY(i,j)] = dudt[UY(i,j)] + dudtChar[2];
            dudt[UZ(i,j)] = dudt[UZ(i,j)] + dudtChar[3];
            dudt[BY(i,j)] = dudt[BY(i,j)] + dudtChar[4];
            dudt[BZ(i,j)] = dudt[BZ(i,j)] + dudtChar[5];
            dudt[P(i,j)] = dudt[P(i,j)] + dudtChar[6];
        }
        else if (i==iSize-1)
        {
            calcCharacteristicsX(i, j, 1, dudtChar);
            dudt[RHO(i,j)] = dudt[RHO(i,j)] + dudtChar[0];
            dudt[UX(i,j)] = dudt[UX(i,j)] + dudtChar[1];
            dudt[UY(i,j)] = dudt[UY(i,j)] + dudtChar[2];
            dudt[UZ(i,j)] = dudt[UZ(i,j)] + dudtChar[3];
            dudt[BY(i,j)] = dudt[BY(i,j)] + dudtChar[4];
            dudt[BZ(i,j)] = dudt[BZ(i,j)] + dudtChar[5];
            dudt[P(i,j)] = dudt[P(i,j)] + dudtChar[6];
        }
        else
        {
            printf("Value of i error. Inside calcRHSOpenBoundaryX(...): at %s Line %d",
                __FILE__, __LINE__);
        }

        if(divBcleaning==2)
        {
            divBpoint = dudx[BX(i,j)] + dudy[BY(i,j)];
            dudt[BX(i,j)] = dudt[BX(i,j)] - dudx[PSI(i,j)];
            dudt[BY(i,j)] = dudt[BY(i,j)] - dudy[PSI(i,j)];

            dudt[PSI(i,j)] = - Ch_divBcleaning*Ch_divBcleaning * (divBpoint
                    + uu[PSI(i,j)]/Cp_divBcleaning/Cp_divBcleaning);
        }
    }
}


void calcCharacteristicsX(size_t i, size_t j, int boundary, real *dudtChar)
{
    real charSlowM, charFastM, charAlfM, charEntro,
        charSlowP, charFastP, charAlfP;
    real charByzP, charUyzP, charByzM, charUyzM;
    real cSound, cSound2, cAlfx, cAlfx2, cAlfy, cAlfz ,cAlf2, 
        cAlfperp,cM2, cN2, cFastx, cSlowx,cFlow;
    real charAlpha1, charAlpha2;
    real charSign, angleY, angleZ;

    real ux,uy,uz,bx,by,bz,rho,p;
    real dux, duy, duz, dbx, dby, dbz, drho, dp;
    
    ux = uu[UX(i,j)];
    uy = uu[UY(i,j)];
    uz = uu[UZ(i,j)];
    bx = uu[BX(i,j)];
    by = uu[BY(i,j)];
    bz = uu[BZ(i,j)];
    rho = uu[RHO(i,j)];
    p = uu[P(i,j)];

    dux = dudx[UX(i,j)];
    duy = dudx[UY(i,j)];
    duz = dudx[UZ(i,j)];
    dbx = dudx[BX(i,j)];
    dby = dudx[BY(i,j)];
    dbz = dudx[BZ(i,j)];
    drho = dudx[RHO(i,j)];
    dp = dudx[P(i,j)];

    // calculate various wave speeds
    cFlow = ux;

    cAlfx = bx / sqrt(rho);
    cAlfy = by / sqrt(rho);
    cAlfz = bz / sqrt(rho);
    cSound2 = adiabaticIdx * p / rho;
    cSound = sqrt(cSound2);
    cAlfx2 = cAlfx * cAlfx;

    if (cAlfx>=0)
    {
        charSign = 1.0;
    }
    else
    {
        charSign = -1.0;
    }


    cAlf2 = cAlfx2 + cAlfy*cAlfy + cAlfz * cAlfz;
    cAlfperp = sqrt( cAlfy*cAlfy + cAlfz * cAlfz );

    if (cAlfperp < 1e-7)
    {
        // angleY = 1/sqrt(2.0);
        // angleZ = 1/sqrt(2.0);

        // Chen Shi: Need to be cautious here.
        angleY = 0.0;
        angleZ = 0.0;;
    }
    else
    {
        angleY = cAlfy/cAlfperp;
        angleZ = cAlfz/cAlfperp;
    }
    

    cM2 = cSound2 + cAlf2;
    cN2 = sqrt( MAX( (cM2*cM2 - 4 * cSound2 * cAlfx2) , 0. ) );

    cFastx = sqrt( MAX((cM2 + cN2) , 0.) ) / sqrt(2.0);
    cSlowx = sqrt( MAX((cM2 - cN2) , 0.) ) / sqrt(2.0);

    // need abs(): sometimes we get a -0.0 inside
    charAlpha1 = sqrt( fabs((cFastx*cFastx - cAlfx2) / 
        (cFastx*cFastx - cSlowx * cSlowx)) );  
    charAlpha2 = sqrt( fabs((cFastx*cFastx - cSound2) / 
        (cFastx*cFastx - cSlowx * cSlowx)) );


    // Some quantities used in calc the final characteristics
    charByzP = angleZ * dbz + angleY * dby;
    charUyzP = angleZ * duz + angleY * duy;
    charByzM = angleY * dbz - angleZ * dby;
    charUyzM = angleY * duz - angleZ * duy;


    // calculate the 7 characteristics
    charFastP = (cFlow + cFastx) * ( charAlpha1 * ( dp / rho / cFastx + dux) 
        + charAlpha2 * ( charByzP / sqrt(rho) - cAlfx/cFastx * charUyzP ) );

    charAlfP = (cFlow + charSign * cAlfx) * ( charUyzM - charSign / sqrt(rho) * charByzM);

    charSlowP = (cFlow + cSlowx) * (-charAlpha2 * ( dp / rho / cSound + 
        charSign * cAlfx / cFastx * dux) + charAlpha1 * ( cSound/cFastx/sqrt(rho) * 
        charByzP - charSign * charUyzP ));

    charEntro = ux * cSound * (dp/adiabaticIdx/p - drho/rho);

    charSlowM = (cFlow - cSlowx) * (-charAlpha2 * ( dp / rho / cSound - 
        charSign * cAlfx / cFastx * dux) + charAlpha1 * ( cSound/cFastx/sqrt(rho) 
        * charByzP + charSign * charUyzP ));

    charAlfM = (cFlow - charSign * cAlfx) * ( charUyzM + charSign / sqrt(rho) * charByzM);

    charFastM = (cFlow - cFastx) * ( charAlpha1 * ( dp / rho / cFastx - dux) 
        + charAlpha2 * ( charByzP / sqrt(rho) + cAlfx/cFastx * charUyzP ) );



    // calculate the dudt according to the boundary type
    if (boundary==0)
    {
        // left boundary
        if (cFlow + cFastx >= 0)
        {
            charFastP = 0;
        }
        
        if (cFlow + charSign * cAlfx >= 0)
        {
            charAlfP = 0;
        }

        if (cFlow + cSlowx >= 0)
        {
            charSlowP = 0;
        }

        if (cFlow >= 0)
        {
            charEntro = 0;
        }

        if (cFlow - cSlowx >= 0)
        {
            charSlowM = 0;
        }

        if (cFlow - charSign * cAlfx >=0)
        {
            charAlfM = 0;
        }

        if (cFlow - cFastx >=0)
        {
            charFastM = 0;
        }


        // inject waves inside the simulation domain
        
    }
    else if (boundary==1)
    {
        // right boundary
        if (cFlow + cFastx <= 0)
        {
            charFastP = 0;
        }
        
        if (cFlow + charSign * cAlfx <= 0)
        {
            charAlfP = 0;
        }

        if (cFlow + cSlowx <= 0)
        {
            charSlowP = 0;
        }

        if (cFlow <= 0)
        {
            charEntro = 0;
        }

        if (cFlow - cSlowx <= 0)
        {
            charSlowM = 0;
        }

        if (cFlow - charSign * cAlfx <=0)
        {
            charAlfM = 0;
        }

        if (cFlow - cFastx <=0)
        {
            charFastM = 0;
        }
    }


    // drho/dt
    dudtChar[0] = rho / cSound * charEntro 
        + rho * (charAlpha2/2/cSound 
        * (charSlowP + charSlowM) - charAlpha1/2/cFastx
        * (charFastP + charFastM));
    
    // dux/dt
    dudtChar[1] = cSlowx * charAlpha2 /2/cSound *
        (charSlowP - charSlowM) - charAlpha1/2 *
        (charFastP - charFastM);

    // duy/dt
    dudtChar[2] = angleZ/2 * (charAlfP + charAlfM)
        + angleY * (charSign * charAlpha1 / 2 *
        (charSlowP - charSlowM) + cAlfx*charAlpha2/2/cFastx * 
        (charFastP - charFastM));

    // duz/dt
    dudtChar[3] = -angleY/2 * (charAlfP + charAlfM)
        + angleZ * (charSign * charAlpha1 / 2 *
        (charSlowP - charSlowM) + cAlfx*charAlpha2/2/cFastx * 
        (charFastP - charFastM));

    // dBy/dt
    dudtChar[4] = -charSign * sqrt(rho) * angleZ/2
        * (charAlfP - charAlfM) - sqrt(rho) 
        * angleY * (cSound * charAlpha1 /2/cFastx 
        * (charSlowP + charSlowM) + charAlpha2 /2 
        * (charFastP + charFastM) );
    
    // dBz/dt
    dudtChar[5] = charSign * sqrt(rho) * angleY/2
        * (charAlfP - charAlfM) - sqrt(rho)
        * angleZ * (cSound * charAlpha1 /2/cFastx 
        * (charSlowP + charSlowM) + charAlpha2 /2 
        * (charFastP + charFastM) );
    
    // dp/dt
    dudtChar[6] = adiabaticIdx * p * (charAlpha2/2/cSound 
        * (charSlowP + charSlowM) - charAlpha1/2/cFastx
        * (charFastP + charFastM));
}


void addSourceTerm()
{
    // add additional source terms to the RHS
    if(viscosity>0)
    {
        for(size_t i=0;i<iSize;i++)
        {
            for(size_t j=0;j<jSize;j++)
            {
                dudt[UX(i,j)] = dudt[UX(i,j)] + viscosity/uu[RHO(i,j)] 
                    * (d2udx2[IDXIJ(0,i,j)] + d2udy2[IDXIJ(0,i,j)]);
                dudt[UY(i,j)] = dudt[UY(i,j)] + viscosity/uu[RHO(i,j)] 
                    * (d2udx2[IDXIJ(1,i,j)] + d2udy2[IDXIJ(1,i,j)]);
                dudt[UZ(i,j)] = dudt[UZ(i,j)] + viscosity/uu[RHO(i,j)] 
                    * (d2udx2[IDXIJ(2,i,j)] + d2udy2[IDXIJ(2,i,j)]);
            }
        }
    }

    if(resistivity>0)
    {
        for(size_t i=0;i<iSize;i++)
        {
            for(size_t j=0;j<jSize;j++)
            {
                dudt[BX(i,j)] = dudt[BX(i,j)] + resistivity * 
                    (d2Bdx2[IDXIJ(0,i,j)] + d2Bdy2[IDXIJ(0,i,j)]);
                dudt[BY(i,j)] = dudt[BY(i,j)] + resistivity * 
                    (d2Bdx2[IDXIJ(1,i,j)] + d2Bdy2[IDXIJ(1,i,j)]);
                dudt[BZ(i,j)] = dudt[BZ(i,j)] + resistivity * 
                    (d2Bdx2[IDXIJ(2,i,j)] + d2Bdy2[IDXIJ(2,i,j)]);
            }
        }
    }

    if(divBcleaning == 1)
    {
        for(size_t i=0;i<iSize;i++)
        {
            for(size_t j=0;j<jSize;j++)
            {
                dudt[BX(i,j)] = dudt[BX(i,j)] + constDiffDivB * ddivBdx[IDXIJ(0,i,j)];
                dudt[BY(i,j)] = dudt[BY(i,j)] + constDiffDivB * ddivBdy[IDXIJ(0,i,j)];
            }
        }
    }
}