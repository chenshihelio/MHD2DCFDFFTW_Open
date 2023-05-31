#ifndef _MACROS_
#define _MACROS_

#define RHO(i,j) ((i) * jSize + (j))
#define UX(i,j) ((iSize + (i)) * jSize + (j))
#define UY(i,j) ((2 * iSize + (i)) * jSize + (j))
#define UZ(i,j) ((3 * iSize + (i)) * jSize + (j))
#define BX(i,j) ((4 * iSize + (i)) * jSize + (j))
#define BY(i,j) ((5 * iSize + (i)) * jSize + (j))
#define BZ(i,j) ((6 * iSize + (i)) * jSize + (j))
#define P(i,j) ((7 * iSize + (i)) * jSize + (j))
#define PSI(i,j) ((8 * iSize + (i)) * jSize + (j))

#define JX(i,j) ((i) * jSize + (j))
#define JY(i,j) ((iSize + (i)) * jSize + (j))
#define JZ(i,j) ((2 * iSize + (i)) * jSize + (j))

#define EX(i,j) ((i) * jSize + (j))
#define EY(i,j) ((iSize + (i)) * jSize + (j))
#define EZ(i,j) ((2 * iSize + (i)) * jSize + (j))


#define IDX(v,i,j,n,m) (((v) * (n) + (i)) * (m) + (j))
#define IDXIJ(v,i,j) (((v) * iSize + (i)) * jSize + (j))

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

#endif