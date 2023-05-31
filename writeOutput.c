#include <stdio.h>
#include <mpi.h>
#include "initialize.h"

void writeGrid()
{
    MPI_Datatype typeContX;
    MPI_File fpWrite;

    if (myRank==0)
    {
        FILE *fpWrite0;

        // write ygrid
        fpWrite0 = fopen("ygrid.dat","wb");
        fwrite(&ny, sizeof(size_t), 1, fpWrite0);
        fwrite(ygrid, sizeof(real), ny, fpWrite0);
        fclose(fpWrite0);

        // write open xgrid file and write nx
        fpWrite0 = fopen("xgrid.dat","wb");
        fwrite(&nx, sizeof(size_t), 1, fpWrite0);
        fclose(fpWrite0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Write xgrid
    MPI_Type_contiguous(iSize,mpi_real,&typeContX);
    MPI_Type_commit(&typeContX);

    MPI_Offset displacement0 = sizeof(size_t);
    MPI_Offset displacement;
    displacement = displacement0 + iOffset * sizeof(real);

    MPI_File_open(MPI_COMM_WORLD, "xgrid.dat",
        MPI_MODE_WRONLY, MPI_INFO_NULL, &fpWrite);

    MPI_File_set_view(fpWrite, displacement, mpi_real, 
        typeContX,"native",MPI_INFO_NULL);

    MPI_File_write_all(fpWrite, xgrid, 1, typeContX, MPI_STATUS_IGNORE);

    MPI_File_close(&fpWrite);
}

void writeArray(int iOutput, real timeSim)
{
    int bufsz = snprintf(NULL, 0, "out%04d.dat", iOutput);
    char fileName[bufsz + 1]; 
    int snprintfReturn = snprintf(fileName,bufsz + 1, "out%04d.dat",iOutput);    


    if (myRank==0)
    {
        FILE *fpWrite = fopen(fileName,"wb");
        fwrite(&timeSim, sizeof(real), 1, fpWrite);
        fclose(fpWrite);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File fpWrite;
    MPI_Offset displacement = sizeof(real);

    MPI_File_open(MPI_COMM_WORLD, fileName,
        MPI_MODE_WRONLY, MPI_INFO_NULL, &fpWrite);

    MPI_File_set_view(fpWrite, displacement, mpi_real, 
        typeMainArr,"native",MPI_INFO_NULL);

    MPI_File_write_all(fpWrite, uu, arrSize, 
        mpi_real, MPI_STATUS_IGNORE);

    MPI_File_close(&fpWrite);
}

void writeDivB(int iOutput, real timeSim)
{
    int bufsz = snprintf(NULL, 0, "divB%04d.dat", iOutput);
    char fileName[bufsz + 1]; 
    int snprintfReturn = snprintf(fileName,bufsz + 1, "divB%04d.dat",iOutput);    


    if (myRank==0)
    {
        FILE *fpWrite = fopen(fileName,"wb");
        fwrite(&timeSim, sizeof(real), 1, fpWrite);
        fclose(fpWrite);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File fpWrite;
    MPI_Offset displacement = sizeof(real);

    MPI_File_open(MPI_COMM_WORLD, fileName,
        MPI_MODE_WRONLY, MPI_INFO_NULL, &fpWrite);

    MPI_File_set_view(fpWrite, displacement, mpi_real, 
        typeDivBArr,"native",MPI_INFO_NULL);

    MPI_File_write_all(fpWrite, divB, iSize * jSize, 
        mpi_real, MPI_STATUS_IGNORE);

    MPI_File_close(&fpWrite);
}


// for debug 
void writeDerivX(int iOutput, real timeSim)
{
    int bufsz = snprintf(NULL, 0, "derivX%04d.dat", iOutput);
    char fileName[bufsz + 1]; 
    int snprintfReturn = snprintf(fileName,bufsz + 1, "derivX%04d.dat",iOutput);    


    if (myRank==0)
    {
        FILE *fpWrite = fopen(fileName,"wb");
        fwrite(&timeSim, sizeof(real), 1, fpWrite);
        fclose(fpWrite);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File fpWrite;
    MPI_Offset displacement = sizeof(real);

    MPI_File_open(MPI_COMM_WORLD, fileName,
        MPI_MODE_WRONLY, MPI_INFO_NULL, &fpWrite);

    MPI_File_set_view(fpWrite, displacement, mpi_real, 
        typeMainArr,"native",MPI_INFO_NULL);

    MPI_File_write_all(fpWrite, dudx, arrSize, 
        mpi_real, MPI_STATUS_IGNORE);

    MPI_File_close(&fpWrite);
}




void writeDerivY(int iOutput, real timeSim)
{
    int bufsz = snprintf(NULL, 0, "derivY%04d.dat", iOutput);
    char fileName[bufsz + 1]; 
    int snprintfReturn = snprintf(fileName,bufsz + 1, "derivY%04d.dat",iOutput);    


    if (myRank==0)
    {
        FILE *fpWrite = fopen(fileName,"wb");
        fwrite(&timeSim, sizeof(real), 1, fpWrite);
        fclose(fpWrite);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_File fpWrite;
    MPI_Offset displacement = sizeof(real);

    MPI_File_open(MPI_COMM_WORLD, fileName,
        MPI_MODE_WRONLY, MPI_INFO_NULL, &fpWrite);

    MPI_File_set_view(fpWrite, displacement, mpi_real, 
        typeMainArr,"native",MPI_INFO_NULL);

    MPI_File_write_all(fpWrite, dudy, arrSize, 
        mpi_real, MPI_STATUS_IGNORE);

    MPI_File_close(&fpWrite);
}

// void writeDerivXX(int iOutput, real timeSim)
// {
//     int bufsz = snprintf(NULL, 0, "derivXX%04d.dat", iOutput);
//     char fileName[bufsz + 1]; 
//     int snprintfReturn = snprintf(fileName,bufsz + 1, "derivXX%04d.dat",iOutput);    


//     if (myRank==0)
//     {
//         FILE *fpWrite = fopen(fileName,"wb");
//         fwrite(&timeSim, sizeof(real), 1, fpWrite);
//         fclose(fpWrite);
//     }

//     MPI_Barrier(MPI_COMM_WORLD);

//     MPI_File fpWrite;
//     MPI_Offset displacement0 = sizeof(real);
//     MPI_Offset displacement;
//     displacement = displacement0 + iOffset * sizeof(real);

//     MPI_File_open(MPI_COMM_WORLD, fileName,
//         MPI_MODE_WRONLY, MPI_INFO_NULL, &fpWrite);

//     MPI_File_set_view(fpWrite, displacement, mpi_real, 
//         typeMainArr,"native",MPI_INFO_NULL);

//     // MPI_File_write_all(fpWrite, uu, NVAR*iSize, mpi_real, MPI_STATUS_IGNORE);
//     MPI_File_write_all(fpWrite, d2udx2, NVAR, typeBlockMain, MPI_STATUS_IGNORE);
//     MPI_File_close(&fpWrite);
// }