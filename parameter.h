// Define the basic parameters, such as the type of real,
// Should not be touched unless necessary

#ifndef _PARAMETER_
#define _PARAMETER_

#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <limits.h>

#define PI acos(-1)

// number of variables. For MHD nvar = 8. (7 for 1DMHD)
// #define NVAR 8 // this is a variable now

#define NHALO 2

#define NHALOFILTER 2

typedef double real;

// CFD default constants
#define ALPHAINNER1 0.375
#define ALPHABOUND1 3.0
#define ALPHASUB1 0.0
#define ALPHAINNER2 0.375
#define ALPHABOUND2 5.0
#define ALPHASUB2 0.0
#define ALPHAINNERF 0.49
#define ALPHASUBF 0.95
#define ALPHABOUNDF0 1.05
#define ALPHABOUNDF1 0.49
#define ALPHABOUNDF2 0.49

// CFD inner scheme for points near boundary 
// (with less points needed)
// for points i=2 (P2)--------
// 1st deriv
#define ALPHA1P2 1./3.
#define A1P2 14./9.
#define B1P2 1./9.
// 2nd deriv
#define ALPHA2P2 2./11.
#define A2P2 12./11.
#define B2P2 3./11.
//-----------------------------

// for points i=1 (P1)---------
// 1st deriv
#define ALPHA1P1 1./4.
#define A1P1 1.5
// 2nd deriv
#define ALPHA2P1 0.1
#define A2P1 1.2
//-----------------------------

// MPI variables
#define mpi_real MPI_DOUBLE  // do not use MPI_REAL (already defined)
#if SIZE_MAX == UCHAR_MAX
   #define mpiSIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define mpiSIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define mpiSIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define mpiSIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define mpiSIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif


#endif