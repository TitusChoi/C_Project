//-----------------------------------------------------------------------
// Created  : 2018.03.26
// Modified : 2018.12.31
// Version  : Microsoft Visual Studio 2019
// Title    : myMatrix.h
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Header file declaration
//-----------------------------------------------------------------------
#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

typedef struct { 
	double** at;
	int rows, cols;
}Matrix;

typedef struct {
	double* at;
	int cols;
}Vector;

using namespace std;

// Create Matrix with specified size
extern Matrix createMat(int _rows, int _cols);

// Create Vector with specified size
extern Vector createVec(int _cols);

// Free a memory allocated vector
extern void	freeVec(Vector _A);

extern void initVec(Vector _A, double _initVal);

//// Print vector
extern void printVec(Vector _A);

// Free a memory allocated matrix
extern void	freeMat(Matrix _A);

// Create a matrix from a text file
extern Matrix txt2Mat(string _filePath, string _fileName);

extern void	initMat(Matrix _A, double _val);

// Copy matrix Elements from A to B
extern void	copyVal(Matrix _A, Matrix _B);

extern Matrix arr2Mat(double* _array, int _rows, int _cols);

//// Print matrix
extern void printMat(Matrix _A);

extern void printMat2(Matrix _A, Matrix _B);

// Create matrix of all zeros
extern Matrix zeros(int _rows, int _cols);

extern Matrix zeros(Matrix _A);

// Create matrix of all ones
extern Matrix ones(int _rows, int _cols);

// Create Transpose matrix
extern Matrix transMat(Matrix _A);

// Create identity 
extern Matrix eye(int _rows, int _cols);

// Create Determinant 
double det(Matrix _A, int dim);

double norm(Matrix _A);
#endif