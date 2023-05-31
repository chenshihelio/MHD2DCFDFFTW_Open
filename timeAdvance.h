#ifndef _TIMEADVANCE_
#define _TIMEADVANCE_

#include "parameter.h"

void rkInit(real dt);
void rkExecute(int iRK);
real calcDt(real dt0);

#endif