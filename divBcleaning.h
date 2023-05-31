#ifndef _DIVBCLEANING_
#define _DIVBCLEANING_

#include "parameter.h"

// ratio between the diffusion and convection
// times of divB (Dedner+2002)
// smaller: stronger diffusion

extern real Cr_divBCleaning;
void calc_Cp_fix_Cr();

#endif