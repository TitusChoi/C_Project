//-----------------------------------------------------------------------
// Created  : 2018.03.26
// Modified : 2019.12.09
// Version  : Microsoft Visual Studio 2019
// Title    : myMatrix.cpp
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------

#include "myMatrix.h"

// Create Matrix with specified size
Matrix createMat(int _rows, int _cols){
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Output;
	// 1. Allocate row array first
	Output.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Output.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Output.rows = _rows;
	Output.cols = _cols;

	return Output;
}

Vector createVec(int _cols) {
	// check vector dimension
	if (_cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createVec' function");
		printf("\n****************************************************\n");
		return createVec(0);
	}

	Vector Output;

	// 1. Then, allocate column 
	for (int i = 0; i < _cols; i++)
		Output.at = (double*)malloc(sizeof(double) * _cols);
	// 2. Initialize row & column values of a vector
	Output.cols = _cols;

	return Output;
}//allocating vector

void freeVec(Vector _A) {
	// 1. Free allocated column memory
	free(_A.at);
}//memory free vector

void printVec(Vector _A) {
	for (int i = 0; i < _A.cols; i++) {
		printf("%f\n", _A.at[i]);
	}
}//printing vector

void initVec(Vector _A, double _initVal) {
	for (int i = 0; i < _A.cols; i++) {
		_A.at[i] = _initVal;
	}
}//initiating vector

// Free a memory allocated matrix
void freeMat(Matrix _A){
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix txt2Mat(string _filePath, string _fileName){
	ifstream file;
	string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Output = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Output.at[i][j] = stof(temp_string);
		}
	file.close();

	return Output;
}

void initMat(Matrix _A, double _initVal){
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _initVal;
}

// Copy Matrix Elements
void copyVal(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n**************************************************");
		printf("\n| ERROR!!: dimension error at 'copyVal' function |");
		printf("\n**************************************************\n");
		return;
	}

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_B.at[i][j] = _A.at[i][j];
}

// 1D-array to Matrix
Matrix arr2Mat(double* _array, int _rows, int _cols){
	Matrix Result = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Result.at[i][j] = _array[i*_cols + j];

	return Result;
}

// Print Matrix
void printMat(Matrix _A){
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%f\t", _A.at[i][j]);
		printf("\n");
	}
}

// Print 2 Matrices
void printMat2(Matrix _A, Matrix _B){
	if (_A.rows != _B.rows) {
		printf("\n****************************************************");
		printf("\n| ERROR!!: dimension error at 'printMat2' function |");
		printf("\n****************************************************\n");
		return;
	}

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%.2f\t", _A.at[i][j]);
		printf("|\t");

		for (int j = 0; j < _B.cols; j++)
			printf("%.2f\t", _B.at[i][j]);
		printf("\n");
	}
}

Matrix zeros(int _rows, int _cols){
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'zeros' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Result = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Result.at[i][j] = 0;

	return Result;
}

Matrix ones(int _rows, int _cols) {
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'ones' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Result = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Result.at[i][j] = 1;

	return Result;
}

Matrix zeros(Matrix _A) {
	if (_A.rows =! _A.cols) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'zeros' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Result = createMat(_A.cols, _A.rows);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Result.at[i][j] = 0;

	return Result;
}

// Transpose
Matrix transMat(Matrix _A) {
	if (_A.rows < 0 || _A.cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'transMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Result = createMat(_A.cols, _A.rows);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Result.at[j][i] = _A.at[i][j];

	return Result;
}

// Identity Matrix
Matrix eye(int _rows, int _cols) {
	if (_rows != _cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'eye' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Output = createMat(_rows, _cols);
	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			if (i == j) {
				Output.at[i][j] = 1.0;
			}
			else {
				Output.at[i][j] = 0.0;
			}

	return Output;
}

// Determinant
double det(Matrix _A, int dim){
	// check matrix dimension
	if (_A.rows =!_A.cols) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'det' function");
		printf("\n****************************************************\n");
		return 0;
	}
	
	else{
		Matrix _D = createMat(dim, dim);
		int m, n;
		double Output;
		double temp = 1;

		if (dim <= 1) {
			printf("\n****************************************************");
			printf("\n  ERROR!!: Please dimension more than 2 ");
			printf("\n****************************************************\n");
			return 0;
		}

		else if (dim == 2){
			Output = 0;
			Output = _A.at[0][0] * _A.at[1][1] - _A.at[0][1] * _A.at[1][0]; // definition of  2 by 2 matrix determinant
			return Output;
		}

		else{
			Output = 0;
			for (int i = 0; i < dim; i++){
				m = 0, n = 0;
				for (int j = 0; j < dim; j++){
					for (int k = 0; k < dim; k++){
						if (j != 0 && k != i){ // storing determinant factor of A
							 _D.at[m][n] = _A.at[j][k];
							n = n + 1;
							if (n > dim - 2){
								m = m + 1;
								n = 0;
							}
						}
					}
				}
				Output = Output + temp * (_A.at[0][i] * det(_D, dim - 1)); // call itself
				temp = -1 * temp;
			}
			return Output;
		}
		
	}
}

double norm(Matrix _A) {
	double output = 0;
		
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {
			output = output + _A.at[i][j] * _A.at[i][j];
		}
	}
	return sqrt(output);
}
