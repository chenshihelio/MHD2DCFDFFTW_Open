#ifndef _WRITEOUTPUT_
#define _WRITEOUTPUT_

#include "parameter.h"

void writeGrid();
void writeArray(int iOutput, real time);
void writeDivB(int iOutput, real time);
void writeDerivX(int iOutput, real time);
void writeDerivY(int iOutput, real time);
void writeDerivXX(int iOutput, real time);
#endif