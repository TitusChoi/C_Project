//-----------------------------------------------------------------------
// Created  : 2019.12.09
// Modified : 2019.12.16
// Version  : Microsoft Visual Studio 2019
// Title    : myComplex.cpp
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Header file declaration
//-----------------------------------------------------------------------
#ifndef		_MY_COMPLEX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_COMPLEX_H

#include <iostream>
#include <string>
#include <fstream>

typedef struct {
	double* re;
	double* im;
	int cols;
}Complex;

using namespace std;

extern Complex createCmp(int _cols);

extern void freeCmp(Complex _A);

extern void printCmp(Complex _A);

extern Complex onesCmp(int _cols);

extern Complex zerosCmp(int _cols);

extern Complex polar2Cmp(double radius, double theta, int _cols);
#endif