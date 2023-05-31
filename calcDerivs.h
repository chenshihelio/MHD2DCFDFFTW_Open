#ifndef _CALCDERIVS_
#define _CALCDERIVS_

#include "parameter.h"

void derivX(int nv, real *arrMain, real *arrHalo0, 
     real *arrHalo1, real *arrDeriv, int iVarStartMain,
     int iVarStartHalo);

void derivXX(int nv, real *arrMain, real *arrHalo0, 
     real *arrHalo1, real *arrDeriv, int iVarStartMain, 
     int iVarStartHalo);

void derivY(int nv, real *arrMain, real *arrDeriv, int iVarStartMain);
void derivYY(int nv, real *arrMain, real *arrDeriv, int iVarStartMain);

#endif