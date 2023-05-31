#include "initialize.h"
#include <math.h>

real Cr_divBCleaning = 0.18; 


void calc_Cp_fix_Cr()
{
    Cp_divBcleaning = sqrt(Cr_divBCleaning * Ch_divBcleaning);
}