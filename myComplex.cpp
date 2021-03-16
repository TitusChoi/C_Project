//-----------------------------------------------------------------------
// Created  : 2019.12.09
// Modified : 2019.12.16
// Version  : Microsoft Visual Studio 2019
// Title    : myComplex.cpp
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------

#include "myComplex.h"

//-----------------------------------------------------------------------
// Function definition
//-----------------------------------------------------------------------
Complex createCmp(int _cols) {
	// check vector dimension
	if (_cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createCmp' function");
		printf("\n****************************************************\n");
		return createCmp(0);
	}

	Complex Output;

	for (int i = 0; i < _cols; i++) {
		Output.re = (double*)malloc(sizeof(double) * _cols);
		Output.im = (double*)malloc(sizeof(double) * _cols);
	}

	Output.cols = _cols;

	return Output;
}//allocating Complex vector

Complex onesCmp(int _cols) {
	if (_cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'onesCmp' function");
		printf("\n****************************************************\n");
		return createCmp(0);
	}

	Complex Result = createCmp(_cols);

	for (int i = 0; i < Result.cols; i++) {
		Result.re[i] = 1;
		Result.im[i] = 1;
	}
	return Result;
}

Complex zerosCmp(int _cols) {
	if (_cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'zerosCmp' function");
		printf("\n****************************************************\n");
		return createCmp(0);
	}

	Complex Result = createCmp(_cols);

	for (int i = 0; i < Result.cols; i++) {
		Result.re[i] = 0;
		Result.im[i] = 0;
	}
	return Result;
}

void printCmp(Complex _A) {
	for (int i = 0; i < _A.cols; i++) {
		printf("%f	+	%fi\n", _A.re[i], _A.im[i]);
	}
}

void freeCmp(Complex _A){
	free(_A.re);
	free(_A.im);
}

Complex polar2Cmp(double radius, double theta, int _cols) {
	Complex Result = createCmp(_cols);
	for (int i = 0; i < Result.cols; i++) {
		Result.re[i] = radius * cos(theta);
		Result.im[i] = radius * sin(theta);
	}
	return Result;
}