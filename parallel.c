#include "parameter.h"
#include "initialize.h"

void sendrecvHaloPointsX()
{
    int npex;
    npex = npe;

    int rightRank = (iRank + 1) % npex;
    int leftRank = (iRank -1 + npex) % npex;
    int tag = 0;

    if (boundaryTypeX==0) // periodic BC
    {
        // send right-boundary points out and receive Halo points on the left-boundary
        MPI_Sendrecv( uu , 1, typeHaloSendX1, rightRank, tag, 
            uuHaloX0, NVAR * NHALO * jSize , mpi_real, leftRank, tag, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // vice-versa
        MPI_Sendrecv(uu, 1, typeHaloSendX0, leftRank, tag, 
            uuHaloX1, NVAR * NHALO * jSize, mpi_real, rightRank, tag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (boundaryTypeX==1) // open boundary
    {
        // send right-boundary points out and receive Halo points on the left-boundary
        if ( (iRank != 0) && (iRank != (npex - 1) )  )
        {    
            MPI_Sendrecv( uu , 1, typeHaloSendX1, rightRank, tag, 
                uuHaloX0, NVAR * NHALO * jSize, mpi_real, leftRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == 0) // left-most processor
        {
            MPI_Send( uu , 1, typeHaloSendX1, rightRank, tag, 
                MPI_COMM_WORLD);
        }
        else if (iRank == (npex-1))
        {
            MPI_Recv( uuHaloX0, NVAR * NHALO * jSize, mpi_real, leftRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // vice versa
        if ( (iRank != 0) && (iRank != (npex - 1) )  )
        {    
            MPI_Sendrecv(uu, 1, typeHaloSendX0, leftRank, tag, 
                uuHaloX1, NVAR * NHALO * jSize, mpi_real, rightRank, tag,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == 0) // left-most processor
        {
            MPI_Recv( uuHaloX1, NVAR * NHALO * jSize, mpi_real, rightRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == (npex-1))
        {
            MPI_Send( uu, 1, typeHaloSendX0, leftRank, tag, 
                MPI_COMM_WORLD);
        }
    }
}


void sendrecvHaloPointsdivBX()
{
    int npex;
    npex = npe;

    int rightRank = (iRank + 1) % npex;
    int leftRank = (iRank -1 + npex) % npex;
    int tag = 0;

    if (boundaryTypeX==0) // periodic BC
    {
        // send right-boundary points out and receive Halo points on the left-boundary
        MPI_Sendrecv(divB , 1, typeHalodivBSendX1, rightRank, tag, 
            divBHaloX0, NHALO * jSize , mpi_real, leftRank, tag, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // vice-versa
        MPI_Sendrecv(divB, 1, typeHalodivBSendX0, leftRank, tag, 
            divBHaloX1, NHALO * jSize, mpi_real, rightRank, tag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (boundaryTypeX==1) // open boundary
    {
        // send right-boundary points out and receive Halo points on the left-boundary
        if ( (iRank != 0) && (iRank != (npex - 1) )  )
        {    
            MPI_Sendrecv( divB , 1, typeHalodivBSendX1, rightRank, tag, 
                divBHaloX0, NHALO * jSize, mpi_real, leftRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == 0) // left-most processor
        {
            MPI_Send( divB , 1, typeHalodivBSendX1, rightRank, tag, 
                MPI_COMM_WORLD);
        }
        else if (iRank == (npex-1))
        {
            MPI_Recv( divBHaloX0, NHALO * jSize, mpi_real, leftRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // vice versa
        if ( (iRank != 0) && (iRank != (npex - 1) )  )
        {    
            MPI_Sendrecv(divB, 1, typeHalodivBSendX0, leftRank, tag, 
                divBHaloX1, NHALO * jSize, mpi_real, rightRank, tag,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == 0) // left-most processor
        {
            MPI_Recv( divBHaloX1, NHALO * jSize, mpi_real, rightRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == (npex-1))
        {
            MPI_Send( divB, 1, typeHalodivBSendX0, leftRank, tag, 
                MPI_COMM_WORLD);
        }
    }
}




void sendrecvHaloPointsEfieldX()
{
    int npex;
    npex = npe;

    int rightRank = (iRank + 1) % npex;
    int leftRank = (iRank -1 + npex) % npex;
    int tag = 0;

    if (boundaryTypeX==0) // periodic BC
    {
        // send right-boundary points out and receive Halo points on the left-boundary
        MPI_Sendrecv(Efield , 1, typeHaloEfieldSendX1, rightRank, tag, 
            EfieldHaloX0, 3 * NHALO * jSize , mpi_real, leftRank, tag, 
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // vice-versa
        MPI_Sendrecv(Efield, 1, typeHaloEfieldSendX0, leftRank, tag, 
            EfieldHaloX1, 3 * NHALO * jSize, mpi_real, rightRank, tag,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (boundaryTypeX==1) // open boundary
    {
        // send right-boundary points out and receive Halo points on the left-boundary
        if ( (iRank != 0) && (iRank != (npex - 1) )  )
        {    
            MPI_Sendrecv( Efield , 1, typeHaloEfieldSendX1, rightRank, tag, 
                EfieldHaloX0, 3 * NHALO * jSize, mpi_real, leftRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == 0) // left-most processor
        {
            MPI_Send( Efield , 1, typeHaloEfieldSendX1, rightRank, tag, 
                MPI_COMM_WORLD);
        }
        else if (iRank == (npex-1))
        {
            MPI_Recv( EfieldHaloX0, 3 * NHALO * jSize, mpi_real, leftRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // vice versa
        if ( (iRank != 0) && (iRank != (npex - 1) )  )
        {    
            MPI_Sendrecv(Efield, 1, typeHaloEfieldSendX0, leftRank, tag, 
                EfieldHaloX1, 3 * NHALO * jSize, mpi_real, rightRank, tag,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == 0) // left-most processor
        {
            MPI_Recv( EfieldHaloX1,3 * NHALO * jSize, mpi_real, rightRank, tag, 
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else if (iRank == (npex-1))
        {
            MPI_Send( Efield, 1, typeHaloEfieldSendX0, leftRank, tag, 
                MPI_COMM_WORLD);
        }
    }
}