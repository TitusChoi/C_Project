//-----------------------------------------------------------------------
// Created  : 2018.03.26
// Modified : 2019.12.16
// Version  : Microsoft Visual Studio 2019
// Title    : myMath.cpp
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Header file declaration
//-----------------------------------------------------------------------

#include "myMath.h"
#define Max 10000
#define PI (double)3.141592653589793238460
#define log2(x) log(x)/log(2)

//-----------------------------------------------------------------------
// Function definition
//-----------------------------------------------------------------------

Matrix addMat(Matrix _A, Matrix _B){
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Output = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Output.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Output;
}

Matrix subtractMat(Matrix _A, Matrix _B){
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Output = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Output.at[i][j] = _A.at[i][j] - _B.at[i][j];

	return Output;
}

Matrix multiplyMat(Matrix _A, Matrix _B){
	if (_A.cols != _B.rows) {
		printf("\n*************************************************************");
		printf("\n  ERROR!!: dimension error at 'multiplication Mat' function");
		printf("\n*************************************************************\n");
		return createMat(0, 0);
	}

	double sum;
	Matrix Output = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _B.cols; j++) {
			sum = 0.0;
			for (int k = 0; k<_B.rows; k++)
				sum += _A.at[i][k] * _B.at[k][j];
			Output.at[i][j] = sum;
		}
	}

	return Output;
}

Vector addVec(Vector _A, Vector _B) {
	if (_A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addVec' function");
		printf("\n*************************************************\n");
		return createVec (0);
	}

	Vector Output = createVec(_A.cols);
	for (int j = 0; j < _A.cols; j++)
		Output.at[j] = _A.at[j] + _B.at[j];

	return Output;
}

Vector subtractVec(Vector _A, Vector _B) {
	if (_A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'subtractVec' function");
		printf("\n*************************************************\n");
		return createVec(0);
	}

	Vector Output = createVec(_A.cols);
	for (int j = 0; j < _A.cols; j++)
		Output.at[j] = _A.at[j] - _B.at[j];

	return Output;
}

double dot(Vector _A, Vector _B) {
	if (_A.cols != _B.cols) {
		printf("\n*************************************************************");
		printf("\n			 ERROR!!: dimension error at 'dot' function         ");
		printf("\n*************************************************************\n");
		return 0;
	}

	else {
		double sum = 0;
		for (int i = 0; i < _A.cols; i++) {
			sum += _A.at[i] * _B.at[i];
		}
		return sum;
	}
}

Matrix multiplysMat(Matrix _A, double w) {
	Matrix Output = createMat(_A.rows, _A.cols);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {
				Output.at[i][j] = w * _A.at[i][j];
		}
	}
	return Output;
}

Complex addCmp(Complex _A, Complex _B) {
	if (_A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addCmp' function");
		printf("\n*************************************************\n");
		return createCmp(0);
	}

	Complex Output = createCmp(_A.cols);
	for (int i = 0; i < _A.cols; i++) {
		Output.re[i] = _A.re[i] + _B.re[i];
		Output.im[i] = _A.im[i] + _B.im[i];
	}

	return Output;
}

Complex subtractCmp(Complex _A, Complex _B) {
	if (_A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'subtractCmp' function");
		printf("\n*************************************************\n");
		return createCmp(0);
	}

	Complex Output = createCmp(_A.cols);
	for (int i = 0; i < _A.cols; i++) {
		Output.re[i] = _A.re[i] - _B.re[i];
		Output.im[i] = _A.im[i] - _B.im[i];
	}

	return Output;
}

Complex multiplyCmp(Complex _A, Complex _B) {
	if (_A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'multiplyCmp' function");
		printf("\n*************************************************\n");
		return createCmp(0);
	}

	Complex Output = createCmp(_A.cols);
	for (int i = 0; i < _A.cols; i++) {
		Output.re[i] = _A.re[i] * _B.re[i] - _A.im[i] * _B.im[i];
		Output.im[i] = _A.im[i] * _B.re[i] + _A.re[i] * _B.im[i];
	}

	return Output;
}

Complex multiplysMat(Complex _A, double w) {
	Complex Output = createCmp(_A.cols);

	for (int i = 0; i < _A.cols; i++) {
		Output.re[i] = w * _A.re[i];
		Output.im[i] = w * _A.im[i];
	}
	return Output;
}

double myfunc(double w) {
	double f;           //Initial myfunc
	double m = 2500;    //Initial mass
	double k = 300000;  //Initial spring constant
	double c = 36000;   //Initial damper
	f = w * w * w - w * w - 1;
	return f;
}

double myfuncDerv(double w) {
	double Df;          //Initial myfunc
	Df = 3 * w * w - 2 * w;
	return Df;
}

void bisecRoot(double a, double b, int maxlimit, double Xn, double tolf) {
	if (myfunc(a) * myfunc(b) >= 0)
		printf("Please retry right a and b.");
	else {
		printf("\n iteration    a           b           Xn          tolerance          using bisection method\n");
		for (int iteration = 0; iteration < maxlimit; iteration++) {
			Xn = (a + b) / 2;                 //The definition of bisection method
			double tol = (b - a) / 2;          //The tolerance of this function
			printf("%3d           %5.5f    %5.5f    %5.5f    %5.5f\n", iteration, a, b, Xn, tol);
			if (myfunc(Xn) == 0) {
				break;
			}
			if (tol < tolf) {
				break;
			}
			if (iteration == maxlimit) {
				break;
			}
			if (myfunc(a) * myfunc(Xn) < 0) {
				b = Xn;
			}
			else {
				a = Xn;
			}
		}
	}
}

void newtonRoot(double Xi, int maxlimit, double Xn, double tolf) {
	printf("\n iteration    Xn              F(Xn)          dF/dx(Xn)         abs(Xm - Xn)          using Newton-Raphson method\n");
	double Xm;               //Xn are pass through Xm
	int iteration;          //iteration is a variable without forloop
	Xn = Xi;       //The estimated solution
	for (iteration = 0; iteration < maxlimit; iteration++) {
		if (myfuncDerv(Xn) == 0) {
			printf("%3d        %5.5f       %5.5f         %5.5f          %5.5f\n", iteration, Xn ,myfunc(Xn), myfuncDerv(Xn), INFINITY);
			printf("The derivative is zero and value is infinity.\n");
		}
		else {
			Xm = Xn - myfunc(Xn) / myfuncDerv(Xn);
			printf("%3d         %5.5f         %5.5f           %5.5f           %5.5f\n", iteration, Xn, myfunc(Xn), myfuncDerv(Xn), fabs(myfunc(Xn) / myfuncDerv(Xn)));
			if (fabs(myfunc(Xn) / myfuncDerv(Xn)) < tolf) {
				Xn = Xm;
				break;
			}
			Xn = Xm;
		}
	}
	if (iteration == maxlimit) {
		printf("Please enter a longer interation.");
	}
}

void gaussjElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn) {
	initMat(_U, 0);
	initMat(_bn, 0);
	Matrix _S = createMat(_b.rows, _b.cols);
	if (_A.cols != _A.rows || _A.cols != _b.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'gaussjElim' function");
		printf("\n*************************************************\n");
		system("pause");
	}
	else {
		for (int i = 0; i < _A.rows; i++) { //making U
			for (int j = 0; j < _A.cols; j++) {
				_U.at[i][j] = _A.at[i][j];
			}
		}

		for (int i = 0; i < _A.cols; i++) { //making bn
			_bn.at[i][0] = _b.at[i][0];
		}

		scalefactor(_U, _bn, _S);
		partialpivot(_U, _bn, _S);
		for (int k = 0; k < _U.rows; k++) {
			for (int i = 0; i < _U.rows; i++) {
				if (i != k) { // diagonal matrix of gauss-jordan elimination
					if (_U.at[k][k] != 0) { // preventing divided by zero
						double mol = _U.at[i][k] / _U.at[k][k];
						for (int j = 0; j < _U.cols; j++) {
							_U.at[i][j] = _U.at[i][j] - mol * _U.at[k][j];
						}
						_bn.at[i][0] = _bn.at[i][0] - mol * _bn.at[k][0];
					}
					else {
						printf("\n we can`t continue gauss elimination because we can`t divide by zero.");
						system("pause");
					}
				}
			}
		}

		printf("\n\n[ Method : Gauss-Jordan Elimination ]\n\n");
		printf("[prob_matU] =\n");
		for (int i = 0; i < _A.rows; i++) { //printing U
			for (int j = 0; j < _A.cols; j++) {
				printf("%f\t", _U.at[i][j]);
			}
			printf("\n");
		}

		printf("\n[prob_vecbn] =\n");
		for (int i = 0; i < _A.cols; i++) { //printing bn
			printf("%f\n", _bn.at[i][0]);
		}
	}

	freeMat(_S);
}

void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn){
	initMat(_U, 0);
	initMat(_bn, 0);
	Matrix _S = createMat(_b.rows, _b.cols);
	if (_A.cols != _A.rows || _A.cols != _b.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'gaussElim' function");
		printf("\n*************************************************\n");
		system("pause");
	}
	else {
		for (int i = 0; i < _A.rows; i++) { //making U
			for (int j = 0; j < _A.cols; j++) {
				_U.at[i][j] = _A.at[i][j];
			}
		}

		for (int i = 0; i < _A.cols; i++) { //making bn
			_bn.at[i][0] = _b.at[i][0];
		}

		scalefactor(_U, _bn, _S);
		partialpivot(_U, _bn, _S);
		for (int k = 0; k < _U.rows - 1; k++) {
			for (int i = 0; i < _U.rows; i++) {
				if (i > k) { // upper triangular matrix of gauss elimination
					if (_U.at[k][k] != 0) { // preventing divided by zero
						double mol = _U.at[i][k] / _U.at[k][k];
						for (int j = 0; j < _U.cols; j++) {
							_U.at[i][j] = _U.at[i][j] - mol * _U.at[k][j];
						}
						_bn.at[i][0] = _bn.at[i][0] - mol * _bn.at[k][0];
					}
					else {
						printf("\n we can`t continue gauss elimination because we can`t divide by zero.");
						system("pause");
					}
				}
			}
		}

		printf("\n\n[ Method : Gauss Elimination ]\n\n");
		printf("[prob_matU] =\n");
		for (int i = 0; i < _A.rows; i++) { //printing U
			for (int j = 0; j < _A.cols; j++) {
				printf("%f\t", _U.at[i][j]);
			}
			printf("\n");
		}

		printf("\n[prob_vecbn] =\n");
		for (int i = 0; i < _A.cols; i++) { //printing bn
			printf("%f\n", _bn.at[i][0]);
		}
		printf("\n");
	}

	freeMat(_S);
}

void backsubj(Matrix _U, Matrix _bn, Matrix _x) {
	for (int i = 0; i < _U.rows; i++) {
		if (_U.at[i][i] != 0) { // preventing divided by zero
			_x.at[i][0] = _bn.at[i][0] / _U.at[i][i];
		}
		else {
			printf("\n we can`t continue gauss elimination because we can`t divide by zero.");
			system("pause");
		}
	}
}

void scalefactor(Matrix _A, Matrix _b, Matrix _S){
	_S = createMat(_b.rows, _b.cols);

	for (int i = 0; i < _A.cols; i++){
		_S.at[i][0] = _A.at[i][0];
		for (int j = 1; j < _A.cols; ++j){//	Comparing between each scale
			if (_S.at[i][0] < fabs(_A.at[i][j]))
				_S.at[i][0] = fabs(_A.at[i][j]);
		}

	}
}

void partialpivot(Matrix _A, Matrix _b, Matrix _S){
	Matrix _Aswap = createMat(_A.rows, _b.cols);
	double _bswap;

	for (int j = 0; j < _A.rows; j++){
		for (int i = j ; i < _A.rows; i++){
			if(_S.at[i][0] == 0){
				if (_b.at[i][0] == 0){
					printf("\nSystem have infinity solution.");
				}
				else {
					printf("\nSystem doesn`t have a solution.");
					system("pause");
				}
			}
			
			if (fabs(_A.at[i][j] / _S.at[i][0]) > fabs(_A.at[j][j] / _S.at[j][0])){//	Scaled Partial Pivoting
				for (int k = 0; k < _A.cols; k++){//	swap each eliments
					_Aswap.at[i][0] = _A.at[i][k];
					_A.at[i][k] = _A.at[j][k];
					_A.at[j][k] = _Aswap.at[i][0];
				}
				_bswap = _b.at[i][0];
				_b.at[i][0] = _b.at[j][0];
				_b.at[j][0] = _bswap;
			}
		}

		if (_A.at[j][j] == 0)//	Singular solution is detected.
		{
			printf("\nSystem is singular.\n");
			system("pause");
		}
	}
	freeMat(_Aswap);
}

void backsub(Matrix _U, Matrix _bn, Matrix _x) {
	for (int i = _U.rows - 1; i >= 0; i--) {
		double sum = 0;
		for (int j = i + 1; j < _U.cols; j++) {
			sum = sum + _U.at[i][j] * _x.at[j][0];
		}
		if (_U.at[i][i] != 0) { // preventing divided by zero
			_x.at[i][0] = (_bn.at[i][0] - sum) / _U.at[i][i];
		}
		else {
			printf("\n we can`t continue gauss elimination because we can`t divide by zero.");
			system("pause");
		}
	}
}

void forwsub(Matrix _L, Matrix _bn, Matrix _y){
	for (int i = 0; i < _L.cols; i++){
		_y.at[i][0] = _bn.at[i][0];
		for (int j = 0; j < i; j++){
			_y.at[i][0] = _y.at[i][0] - _L.at[i][j] * _y.at[j][0];
		}
	}

	printf("\n[prob_vecy] =\n");
	for (int i = 0; i < _L.cols; i++)
	{
		printf("%.6f\n", _y.at[i][0]);
	}
}

void backsubl(Matrix _U, Matrix _y, Matrix _x) {
	for (int i = _U.rows - 1; i >= 0; i--) {
		_x.at[i][0] = _y.at[i][0];
		for (int j = i + 1; j < _U.cols; j++) {
			_x.at[i][0] = _x.at[i][0] - _U.at[i][j] * _x.at[j][0];
		}
		if (_U.at[i][i] != 0) { // preventing divided by zero
			_x.at[i][0] = _x.at[i][0] / _U.at[i][i];
		}
		else {
			printf("\n we can`t continue LU decomposition because we can`t divide by zero.");
			system("pause");
		}
	}
}

void LUdcmp(Matrix _A, Matrix _L, Matrix _U) {
	initMat(_L, 0);
	initMat(_U, 0);
	if (_A.cols != _A.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'LUdcmp' function");
		printf("\n*************************************************\n");
		system("pause");
	}
	else {//	LU decomposition is derived from Crout's algorithm

		for (int k = 0; k < _A.cols; k++) { 
			for (int i = 0; i < _A.cols; i++) {//	Upper triangular matrix
				if (k > i)
					_U.at[k][k] = 0.0;
				else {
					double sum = 0;
					for (int j = 0; j < k; j++)
						sum = sum + _L.at[k][j] * _U.at[j][i];
					_U.at[k][i] = _A.at[k][i] - sum;
				}
			}

			for (int i = 0; i < _A.cols; i++) {//	Lower triangular matrix
				if (k == i)
					_L.at[k][k] = 1.0;
				else {
					double sum = 0;
					for (int j = 0; j < k; j++)
						sum = sum + (_L.at[i][j] * _U.at[j][k]);
					_L.at[i][k] = (_A.at[i][k] - sum) / _U.at[k][k];
				}
			}
		}

		printf("\n\n[ Method : LU decomposition ]\n\n");
		printf("[prob_matL] =\n");
		for (int i = 0; i < _A.cols; i++){
			for (int j = 0; j < _A.cols; j++)
				printf("%f\t", _L.at[i][j]);
			printf("\n");
		}

		printf("\n[prob_matU] =\n");
		for (int i = 0; i < _A.cols; i++){
			for (int j = 0; j < _A.cols; j++) {
				printf("%f\t", _U.at[i][j]);
			}
			printf("\n");
		}
	}
}

Matrix solveLinear(Matrix _A, Matrix _b, char method) {
	//-----------------------------------------------------------------------
	// Initialization
	//-----------------------------------------------------------------------
	Matrix _L = createMat(_A.rows, _A.cols);
	Matrix _U = createMat(_A.rows, _A.cols);
	Matrix _bn = createMat(_b.rows, _b.cols);
	Matrix _x = createMat(_b.rows, _b.cols);
	Matrix _y = createMat(_b.rows, _b.rows);
	Matrix _S = createMat(_A.cols, _b.cols);
	//-----------------------------------------------------------------------
	// Using each functions 
	//-----------------------------------------------------------------------

	switch (method) {
	case 'g':	//Gauss Elimination with scaled pivoting
		gaussElim(_A, _b, _U, _bn);
		backsub(_U, _bn, _x);
		return _x;
		freeMat(_U);
		freeMat(_bn);
		break;

	case 'j':	//Gauss-Jordan Elimination with scaled pivoting
		gaussjElim(_A, _b, _U, _bn);
		backsubj(_U, _bn, _x);
		return _x;
		freeMat(_U);
		freeMat(_bn);
		break;

	case 'L':	//LU decomposition with scaled pivoting
		scalefactor(_A, _b, _S);
		partialpivot(_A, _b, _S);
		LUdcmp(_A, _L, _U);
		forwsub(_L, _b, _y);
		backsubl(_U, _y, _x);
		return _x;
		freeMat(_L);
		freeMat(_U);
		freeMat(_y);
		break;

	default:
		printf("\n Please retry input g, j, or L to solve the equations.");
		return createMat(0, 0);
	}
}

//Inverse Matrix using gaussElim
Matrix inverseMat(Matrix _A) {
	if (_A.cols != _A.rows) {
		printf("\n***************************************************");
		printf("\n  ERROR!!: dimension error at 'inverseMat' function");
		printf("\n***************************************************\n");
		return createMat(0, 0);
	}

	else if (det(_A, _A.cols) == 0) {
		printf("\n****************************************************************");
		printf("\n  The matrix is a singular solution because determinant is zero.");
		printf("\n**************************************************************\n");
		return createMat(0, 0);
	}

	else {
		Matrix _Ai = createMat(_A.rows, _A.cols);
		copyVal(_A, _Ai);
		Matrix _I = eye(_A.rows, _A.cols);
		for (int k = 0; k < _Ai.rows; k++) {
			for (int i = 0; i < _Ai.rows; i++) {
				if (i != k) {
					double mol = _Ai.at[i][k] / _Ai.at[k][k];
					for (int j = 0; j < _Ai.cols; j++) {
						_Ai.at[i][j] = _Ai.at[i][j] - mol * _Ai.at[k][j];
						_I.at[i][j] = _I.at[i][j] - mol * _I.at[k][j];//	It is different to original code.
					}
				}
			}
		}

		for (int k = 0; k < _Ai.rows; k++) {
			for (int j = 0; j < _Ai.rows; j++) {
				_I.at[k][j] = _I.at[k][j] / _Ai.at[k][k];
			}
			_Ai.at[k][k] = _Ai.at[k][k] / _Ai.at[k][k];//	as 1
		}
		freeMat(_Ai);
		return _I;
	}
}

Matrix QRdcmph(Matrix _A, Matrix _Q, Matrix _R) {
	if (_A.cols != _A.rows) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: Input matrix is not square matrix at 'QRdcmph' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}
	else {//   QR decomposition householder method
		Matrix _E = createMat(_A.rows, _A.cols);

		Matrix _I = eye(_A.rows, _A.cols);
		Matrix _H = createMat(_A.rows, _A.cols);
		Matrix _c = createMat(_A.rows, _A.cols);
		Matrix _e = createMat(_A.rows, _A.cols);
		Matrix _v = createMat(_A.rows, _A.cols);
		Matrix _vt = createMat(_A.cols, _A.rows);

		initMat(_c, 0);
		initMat(_e, 0);
		initMat(_v, 0);
		initMat(_vt, 0);
		initMat(_H, 0);

		printf("\n[prob1_A] =\n");
		for (int i = 0; i < _A.rows; i++) {
			for (int j = 0; j < _A.cols; j++) {
				printf("%f\t\t", _A.at[i][j]);
				_E.at[i][j] = _A.at[i][j];
			}
			printf("\n");
		}

		for (int m = 0; m < Max; m++) {// iteration

			_Q = eye(_A.rows, _A.cols); //making Q
			for (int i = 0; i < _A.rows; i++) { //making R
				for (int j = 0; j < _A.cols; j++) {
					_R.at[i][j] = _E.at[i][j];
				}
			}

			for (int k = 0; k < _A.rows - 1; k++) {
				initMat(_c, 0);
				initMat(_e, 0);
				initMat(_v, 0);
				initMat(_vt, 0);
				initMat(_H, 0);

				for (int i = k; i < _A.rows; i++) {
					_c.at[i][0] = _R.at[i][k];
				}

				if (_c.at[k][0] >= 0) {
					_e.at[k][0] = 1;
				}
				else {
					_e.at[k][0] = -1;
				}

				double cnorm = 0;

				for (int i = 0; i < _A.rows; i++) {
					cnorm = cnorm + _c.at[i][0] * _c.at[i][0];
				}

				for (int j = 0; j < _A.rows; j++)
					_v.at[j][0] = _c.at[j][0] + (sqrt(cnorm) * _e.at[j][0]);

				for (int j = 0; j < _A.rows; j++)
				{
					_vt.at[0][j] = _v.at[j][0];
				}

				double vnorm = 0;

				for (int j = 0; j < _A.rows; j++) {
					vnorm = vnorm + _vt.at[0][j] * _v.at[j][0];
				}


				for (int i = 0; i < _A.rows; i++) {
					for (int j = 0; j < _A.cols; j++) {
						_H.at[i][j] = _I.at[i][j] - ((2 / vnorm)* _v.at[i][0] * _vt.at[0][j]);
					}
				}

				_R = multiplyMat(_H, _R);
				_Q = multiplyMat(_Q, _H);
			}
			_E = multiplyMat(_R, _Q);
		}

		printf("\n[prob1_finalQ] =\n");
		for (int i = 0; i < _A.rows; i++) {
			for (int j = 0; j < _A.cols; j++)
				printf("%f\t\t", _Q.at[i][j]);
			printf("\n");
		}

		printf("\n[prob1_finalR] =\n");
		for (int i = 0; i < _A.rows; i++) {
			for (int j = 0; j < _A.cols; j++)
				printf("%f\t\t", _R.at[i][j]);
			printf("\n");
		}

		printf("\n[prob1_eigenvalue matrix] =\n");
		for (int i = 0; i < _A.rows; i++) {
			for (int j = 0; j < _A.cols; j++)
				printf("%f\t\t", _E.at[i][j]);
			printf("\n");
		}

		return _E;
	}
}

void eigVec(Matrix _A, Matrix _E, Matrix _V) {
	Matrix _B = createMat(_A.rows, _A.cols);
	Matrix _L = createMat(_A.rows, _A.cols);
	Matrix _I = eye(_A.rows, _A.cols);
	_V = eye(_A.rows, _A.cols);

	Matrix temp1 = createMat(_A.rows - 1, 1);
	Matrix temp2 = createMat(_A.rows - 1, _A.cols - 1);

	Matrix _v = createMat(_A.rows - 1, 1);

	Matrix VDVt = createMat(_A.rows, _A.cols);

	int p = 0;
	int q = 0;

	Matrix _EigVal = createMat(_E.rows, 1);
	for (int i = 0; i < _E.rows; i++) {
		_EigVal.at[i][0] = _E.at[i][i];
	}

	for (int k = 0; k < _I.rows; k++) {
		initMat(_L, 0);
		initMat(_B, 0);
		initMat(_v, 0);

		_L = multiplysMat(_I, _EigVal.at[k][0]);
		_B = subtractMat(_A, _L);

		for (int i = 0; i < _B.rows; i++) {
			if (k != i) {
				for (int j = 0; j < _B.cols; j++) {
					if (k != j) {
						temp2.at[p][q] = _B.at[i][j];
						q = q + 1;
					}
					else {
						temp1.at[p][0] = -1 * _B.at[i][j];
					}
				}
				q = 0;
				p = p + 1;
			}
		}
		p = 0;

		Matrix temp2i = inverseMat(temp2);//	Inverse matrix using GaussElim

		for (int i = 0; i < temp2i.rows; ++i) {
			for (int j = 0; j < temp1.cols; ++j) {
				for (int m = 0; m < temp2i.cols; ++m) {
					_v.at[i][j] = _v.at[i][j] + temp2i.at[i][m] * temp1.at[m][j];
				}
			}
		}

		for (int i = 0; i < _V.rows; i++) {
			if (k != i) {
				_V.at[i][k] = _v.at[p][0];
				p = p + 1;
			}
		}
		p = 0;

		for (int i = 0; i < _V.rows; i++) {
			double Vnorm = 0;
			for (int j = 0; j < _V.cols; j++) {
				Vnorm = Vnorm + _V.at[j][i] * _V.at[j][i];
			}

			for (int j = 0; j < _V.rows; j++) {
				_V.at[j][i] = _V.at[j][i] / sqrt(Vnorm);
			}
		}
	}

	VDVt = multiplyMat(multiplyMat(_V, _E), transMat(_V));

	printf("\n[eigenvector matirx] = \n");
	for (int i = 0; i < _V.rows; i++)
	{
		for (int j = 0; j < _V.cols; j++)
		{
			printf("%lf\t", _V.at[i][j]);
		}
		printf("\n");
	}

	printf("\n[V*D*V'] = \n");
	for (int i = 0; i < _V.rows; i++)
	{
		for (int j = 0; j < _V.cols; j++)
		{
			printf("%lf\t", VDVt.at[i][j]);
		}
		printf("\n");
	}

	freeMat(_B);
	freeMat(_v);
	freeMat(_I);
	freeMat(temp1);
	freeMat(temp2);
	freeMat(VDVt);
}

Matrix normalizeVec(Matrix _V) {
	Matrix output = createMat(_V.rows, 1);
	initMat(output, 0);

	for (int i = 0; i < _V.rows; i++) {
		for (int j = 0; j < _V.cols; j++) {
			output.at[j][0] = _V.at[j][i];
		}

		double vnorm = 0;

		for (int l = 0; l < _V.rows; l++) {
			vnorm = vnorm + output.at[0][l] * output.at[l][0];
		}

		copyVal(multiplysMat(output, 1.0 / sqrt(vnorm)), output);
		for (int j = 0; j < _V.cols; j++) {
			_V.at[j][i] = output.at[j][0];
		}
	}
	return _V;
}

//Inverse Matrix with cofactor
Matrix inverseMatwc(Matrix _A) {
	if (_A.cols != _A.rows) {
		printf("\n***************************************************");
		printf("\n  ERROR!!: dimension error at 'inverseMatwc' function");
		printf("\n***************************************************\n");
		return createMat(0, 0);
	}

	else if (det(_A, _A.cols) == 0) {
		printf("\n****************************************************************");
		printf("\n  The matrix is a singular solution because determinant is zero.");
		printf("\n**************************************************************\n");
		return createMat(0, 0);
	}

	else {
		Matrix Output = createMat(_A.rows, _A.cols);
		for (int i = 0; i < _A.rows; i++){
			for (int j = 0; j < _A.cols; j++){
				Output.at[i][j] = multiplysMat(transMat(cofactor(_A)), 1 / det(_A, _A.rows)).at[i][j];
			}
		}
		return Output;
	}
}

Matrix cofactor(Matrix _A) {
	Matrix Output = createMat(_A.rows, _A.cols);
	Matrix temp = createMat(_A.rows, _A.cols);

	for (int q = 0; q < _A.rows; q++){
		for (int p = 0; p < _A.rows; p++){
			int m = 0;
			int n = 0;
			for (int i = 0; i < _A.rows; i++){
				for (int j = 0; j < _A.rows; j++){
					if (i != q && j != p){
						temp.at[m][n] = _A.at[i][j];
						if (n < (_A.rows - 2)){
							n++;
						}
						else{
							n = 0;
							m++;
						}
					}
				}
			}
			Output.at[q][p] = pow(-1, q + p) * det(temp, _A.rows - 1);
		}
	}
	return Output;
}

// Rank
int rankMat(Matrix _A) {

	Matrix _Q = createMat(_A.rows, _A.cols);
	Matrix _R = createMat(_A.rows, _A.cols);

	initMat(_Q, 0);
	initMat(_R, 0);

	Matrix factor = QRdcmph(_A, _Q, _R); //we find eigenvalue

	Matrix _EigVal = createMat(factor.rows, 1);
	initMat(_EigVal, 0);


	int Output = 0;
	for (int i = 0; i < factor.rows; i++) {
		_EigVal.at[i][0] = factor.at[i][i];
		if (_EigVal.at[i][0] != 0) {
			Output = Output + 1; // we count each 
		}
		else {
			printf("Rank is smaller than given matrix dimension.");
			break;
		}
	}

	return Output;

}

double factorial(int n) {
	double result = 1;
	if (n == 0) {
		return 1;
	}
	else {
		return n * factorial(n - 1);
	}
}

void polyfit(Matrix _X, Matrix _Y, Matrix _z, int n) {
	if (_X.rows != _Y.rows || _X.cols != _Y.cols) {
		printf("\n******************************************************************");
		printf("\n  ERROR!!: Input matrix is not same dimension at 'polyfit' function");
		printf("\n****************************************************************\n");
		system("pause");
	}

	else if (_X.rows <= n) {
		printf("\n******************************************************************");
		printf("\n  ERROR!!: You take values like n < m.");
		printf("\n****************************************************************\n");
		system("pause");
	}

	else {
		int m = _X.rows;

		Matrix Xsum = createMat(2 * n, 1);
		initMat(Xsum, 0);

		Matrix _S = createMat(n + 1, n + 1);
		initMat(_S, 0);

		Matrix _b = createMat(n + 1, 1);
		initMat(_b, 0);

		initMat(_z, 0);


		for (int i = 0; i < 2 * n; i++) {
			for (int j = 0; j < m; j++) {
				Xsum.at[i][0] = Xsum.at[i][0] + pow((_X.at[j][0]), i + 1);
			}
		}

		_S.at[0][0] = m;

		for (int j = 1; j < n + 1; j++) {
			_S.at[0][j] = Xsum.at[j - 1][0];
		}

		for (int j = 0; j < m; j++) {
			_b.at[0][0] = _b.at[0][0] + _Y.at[j][0];
		}

		for (int i = 1; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				_S.at[i][j] = Xsum.at[i + j - 1][0];
			}
			for (int j = 0; j < m; j++) {
				_b.at[i][0] = _b.at[i][0] + pow((_X.at[j][0]), i) * _Y.at[j][0];
			}
		}

		Matrix _U = createMat(_S.rows, _S.cols);
		initMat(_U, 0);

		Matrix _bn = createMat(_b.rows, _b.cols);
		initMat(_bn, 0);

		gaussjElim(_S, _b, _U, _bn);
		backsubj(_U, _bn, _z);

		printf("\n************************ Matrix S %d X %d **********************\n", n + 1, n + 1);
		printMat(_S);

		printf("\n************************ Vector b %d X %d **********************\n", n + 1, 1);
		printMat(_b);

		printf("\n************************ Vector z %d X %d **********************\n", n + 1, 1);
		printMat(_z);

		freeMat(Xsum);
		freeMat(_S);
		freeMat(_b);
		freeMat(_U);
		freeMat(_bn);
	}
}

Matrix deriv(Matrix _X, Matrix _Y, double h, char method) {
	int n = _X.rows;
	//-----------------------------------------------------------------------
	// Using each functions 
	//-----------------------------------------------------------------------
	switch (method) {
	case 'f':
		Matrix _Vf = createMat(n - 1, 1);
		initMat(_Vf, 0);
		for (int i = 0; i < n - 1; i++) {
			_Vf.at[i][0] = (_Y.at[i + 1][0] - _Y.at[i][0]) / (_X.at[i + 1][0] - _X.at[i][0]);
		}
		return _Vf;
	case 'b':
		Matrix _Vb = createMat(n - 1, 1);
		initMat(_Vb, 0);
		for (int i = 1; i < n; i++) {
			_Vb.at[i][0] = (_Y.at[i][0] - _Y.at[i - 1][0]) / (_X.at[i][0] - _X.at[i - 1][0]);
		}
		return _Vb;
	case 'c':
		Matrix _Vc = createMat(n - 2, 1);
		initMat(_Vc, 0);
		for (int i = 1; i < n - 1; i++) {
			_Vc.at[i][0] = (_Y.at[i + 1][0] - _Y.at[i - 1][0]) / (_X.at[i + 1][0] - _X.at[i - 1][0]);
		}
		return _Vc;
	case 'F':
		Matrix _VF = createMat(n - 2, 1);
		initMat(_VF, 0);
		for (int i = 0; i < n - 2; i++) {
			_VF.at[i][0] = (-3 * _Y.at[i][0] + 4 * _Y.at[i + 1][0] - _Y.at[i + 2][0]) / (_X.at[i + 2][0] - _X.at[i][0]);
		}
		return _VF;
	case 'B':
		Matrix _VB = createMat(n - 2, 1);
		initMat(_VB, 0);
		for (int i = 2; i < n; i++) {
			_VB.at[i][0] = (_Y.at[i - 2][0] - 4 * _Y.at[i - 1][0] + 3 * _Y.at[i][0]) / (_X.at[i][0] - _X.at[i - 2][0]);
		}
		return _VB;
	default:
		printf("\n Please retry choose a method.");
		system("pause");
		return createMat(0, 0);
	}
}

double deriv(double x, double myfunc(double x), double h, char method) {
	double v;
	
	//-----------------------------------------------------------------------
	// Using each functions 
	//-----------------------------------------------------------------------
	switch (method) {
	case 'f':
		v = ( myfunc(x + h) - myfunc(x) ) / h;
		return v;
	case 'b':
		v = (myfunc(x) - myfunc(x - h)) / h;
		return v;
	case 'c':
		v = (myfunc(x + h) - myfunc(x - h)) /(2 * h);
		return v;
	case 'F':
		v = (- 3 * myfunc(x) + 4 * myfunc(x + h) - myfunc(x + 2 * h)) / (2 * h);
		return v;
	case 'B':
		v = (myfunc(x - 2 * h) - 4 * myfunc(x - h) + 3 * myfunc(x)) / (2 * h);
		return v;
	default:
		printf("\n Please retry choose a method.");
		system("pause");
		return 0;
	}
}

Matrix deriv2(Matrix _X, Matrix _Y, double h, char method) {
	int n = _X.rows;
	//-----------------------------------------------------------------------
	// Using each functions 
	//-----------------------------------------------------------------------
	switch (method) {
	case 'f':
		Matrix _Af = createMat(n - 2, 1);
		initMat(_Af, 0);
		for (int i = 0; i < n - 2; i++) {
			_Af.at[i][0] = (_Y.at[i + 2][0] - 2 * _Y.at[i + 1][0] + _Y.at[i][0]) / ((_X.at[i + 1][0] - _X.at[i][0])*(_X.at[i][0] - _X.at[i - 1][0]));
		}
		return _Af;
	case 'b':
		Matrix _Ab = createMat(n - 2, 1);
		initMat(_Ab, 0);
		for (int i = 2; i < n; i++) {
			_Ab.at[i][0] = (_Y.at[i - 2][0] - 2 * _Y.at[i - 1][0] + _Y.at[i][0]) / ((_X.at[i + 1][0] - _X.at[i][0])*(_X.at[i][0] - _X.at[i - 1][0]));
		}
		return _Ab;
	case 'c':
		Matrix _Ac = createMat(n - 2, 1);
		initMat(_Ac, 0);
		for (int i = 1; i < n - 1; i++) {
			_Ac.at[i][0] = (_Y.at[i + 1][0]   + _Y.at[i - 1][0]) / ((_X.at[i + 1][0] - _X.at[i][0])*(_X.at[i][0] - _X.at[i - 1][0]));
		}
		return _Ac;

	default:
		printf("\n Please retry choose a method.");
		system("pause");
		return createMat(0, 0);
	}
}

double deriv2(double x, double myfunc(double x), double h, char method) {
	double a;

	//-----------------------------------------------------------------------
	// Using each functions 
	//-----------------------------------------------------------------------
	switch (method) {
	case 'f':
		a = (myfunc(x + 2 * h) - 2 * myfunc(x + h) + myfunc(x)) / (h * h);
		return a;
	case 'b':
		a = (myfunc(x - 2 * h) - 2 * myfunc(x - h) + myfunc(x)) / (h * h);
		return a;
	case 'c':
		a = (myfunc(x + h) - 2 * myfunc(x) + myfunc(x - h)) / (h * h);
		return a;
	default:
		printf("\n Please retry choose a method.");
		system("pause");
		return 0;
	}
}

double myfunc1(double x) {
	return pow(sin(x), 2);
}

double myfunc2(double x) {
	return cos(2 * x) * exp(-1 * x);
}

double trapz(Matrix _X, Matrix _Y){
	int N = _X.rows - 1;
	double sum = 0;
	double h = (_X.at[N][0] - _X.at[0][0]) / N;

	for (int i = 0; i < N + 1; i++) {
		if (i == 0 || i == N) {
			sum = sum + _Y.at[i][0] / 2;
		}
		else {
			sum = sum + _Y.at[i][0];
		}
	}

	return sum * h;
}

double trapz(double myfunc(double x), double a, double b, int N) {
	double sum = 0;
	double h = (b - a) / N;
	Matrix _X = createMat(N + 1, 1);
	Matrix _Y = createMat(N + 1, 1);

	for (int i = 0; i < N + 1; i++) {
		_X.at[i][0] = a + i * h;
		_Y.at[i][0] = myfunc(_X.at[i][0]);
	}

	for (int i = 0; i < N + 1; i++) {
		if (i == 0 || i == N) {
			sum = sum + _Y.at[i][0] / 2;
		}
		else {
			sum = sum + _Y.at[i][0];
		}
	}

	return sum * h;
}

double simpson13(Matrix _X, Matrix _Y) {
	int N = _X.rows - 1;
	if (N % 2 == 1) {	// N is must be even number.
		N = N + 1;
	}

	double sum_odd = 0;
	double sum_even = 0;
	double h = (_X.at[N][0] - _X.at[0][0]) / N;

	for (int i = 1; i < N; i++) {
		if (i % 2 == 0) {
			sum_even = sum_even + _Y.at[i][0];
		}
		else {
			sum_odd = sum_odd + _Y.at[i][0];
		}
	}

	return h / 3 * (_Y.at[0][0] + _Y.at[N][0] + 4 * sum_odd + 2 * sum_even);
}

double simpson13(double myfunc(double x), double a, double b, int N) {
	if (N % 2 == 1){	// N is must be even number.
		N = N + 1;
	}

	Matrix _X = createMat(N + 1, 1);
	Matrix _Y = createMat(N + 1, 1);
	
	double sum_odd = 0;
	double sum_even = 0;
	double h = (b - a) / N;
	
	for (int i = 0; i < N + 1; i++){
		_X.at[i][0] = a + i * h;
		_Y.at[i][0] = myfunc(_X.at[i][0]);
	}

	for (int i = 1; i < N; i++){
		if (i % 2 == 0){
			sum_even = sum_even + _Y.at[i][0];
		}
		else{
			sum_odd = sum_odd + _Y.at[i][0];
		}
	}
	
	return h / 3 * (_Y.at[0][0] + _Y.at[N][0] + 4 * sum_odd + 2 * sum_even);
}

double dy_test(double x, double y) {
	double m = 10;
	double k = 1000;
	double c = 140;
	double A = 100;
	double f = 0.5;
	double F = A * cos(2 * PI * f * x);
	double dy = - k * y / m + F / m;
	return dy;
}

double dy2_test(double x, double y, double dy) {
	double m = 10;
	double k = 1000;
	double c = 140;
	double A = 100;
	double f = 0.5;
	double F = A * cos(2 * PI * f * x);
	double dy2 = -c * dy / m - k * y / m + F / m;
	return dy2;
}

void ODEeuler(Matrix _X, Matrix _YE, double x0, double xf, double y0, double h) {
	_X.at[0][0] = x0;
	_YE.at[0][0] = y0;
	int N = (int)((xf - x0) / h);

	for (int i = 0; i < N; i++){
		_X.at[i + 1][0] = _X.at[i][0] + h;
		_YE.at[i + 1][0] = _YE.at[i][0] + h * dy_test(_X.at[i][0], _YE.at[i][0]);
	}
}

void ODEeulerM(Matrix _X, Matrix _YEM, double x0, double xf, double y0, double h) {
	_X.at[0][0] = x0;
	_YEM.at[0][0] = y0;
	double slope1 = 0;
	double _YEMu = 0;
	double slope2 = 0;
	int N = (int)((xf - x0) / h);

	for (int i = 0; i < N; i++) {
		_X.at[i + 1][0] = _X.at[i][0] + h;
		slope1 = dy_test(_X.at[i][0], _YEM.at[i][0]);
		_YEMu = _YEM.at[i][0] + h * slope1;
		slope2 = dy_test(_X.at[i + 1][0], _YEMu);
		_YEM.at[i + 1][0] = _YEM.at[i][0] + (slope1 + slope2) * h / 2;
	}
}

void ODEmid(Matrix _X, Matrix _YM, double x0, double xf, double y0, double h) {
	_X.at[0][0] = x0;
	_YM.at[0][0] = y0;
	double _Xm = 0;
	double _Ym = 0;
	double slopemid = 0;

	int N = (int)((xf - x0) / h);
	
	for (int i = 0; i < N; i++) {
		_Xm = _X.at[i][0] + h / 2;
		_Ym = _YM.at[i][0] + h * dy_test(_X.at[i][0], _YM.at[i][0]) / 2;
		slopemid = dy_test(_Xm, _Ym);
		_YM.at[i + 1][0] = _YM.at[i][0] + slopemid * h;
		_X.at[i + 1][0] = _X.at[i][0] + h;
	}
}

void ODERK2sys1(Matrix _X, Matrix _YRK, double x0, double xf, double y0, double h) {
	_X.at[0][0] = x0;
	_YRK.at[0][0] = y0;

	double slope1 = 0;
	double slope2 = 0;
	double slope = 0;

	int N = (int)((xf - x0) / h);

	for (int i = 0; i < N; i++) {
		slope1 = dy_test(_X.at[i][0], _YRK.at[i][0]);
		slope2 = dy_test(_X.at[i][0] + h, _YRK.at[i][0] + slope1 * h);
		slope = (slope1 + slope2) / 2;
		_YRK.at[i + 1][0] = _YRK.at[i][0] + slope * h;
		_X.at[i + 1][0] = _X.at[i][0] + h;
	}
}

void ODERK4sys1(Matrix _X, Matrix _YRK, double x0, double xf, double y0, double h) {
	_X.at[0][0] = x0;
	_YRK.at[0][0] = y0;

	double slope1 = 0;
	double slope2 = 0;
	double slope3 = 0;
	double slope4 = 0;
	double slope = 0;

	int N = (int)((xf - x0) / h);

	for (int i = 0; i < N; i++) {
		slope1 = dy_test(_X.at[i][0], _YRK.at[i][0]);
		slope2 = dy_test(_X.at[i][0] + h / 2, _YRK.at[i][0] + slope1 * h / 2);
		slope3 = dy_test(_X.at[i][0] + h / 2, _YRK.at[i][0] + slope2 * h / 2);
		slope4 = dy_test(_X.at[i][0] + h, _YRK.at[i][0] + slope3 * h);
		slope = (slope1 + 2 * slope2 + 2 * slope3 + slope4) / 6;
		_YRK.at[i + 1][0] = _YRK.at[i][0] + slope * h;
		_X.at[i + 1][0] = _X.at[i][0] + h;
	}
}

void ODERK2sys2(Matrix _X, Matrix _YRK, Matrix _dYRK, double x0, double xf, double y0, double dy0, double h) {
	_X.at[0][0] = x0;
	_YRK.at[0][0] = y0;
	_dYRK.at[0][0] = dy0;

	int N = (int)((xf - x0) / h);

	double temp11 = 0;
	double temp12 = 0;

	double temp21 = 0;
	double temp22 = 0;

	for (int i = 0; i < N; i++) {
		temp11 = h * _dYRK.at[i][0];
		temp12 = h * dy2_test(_X.at[i][0], _YRK.at[i][0], _dYRK.at[i][0]);
		temp21 = h * (_dYRK.at[i][0] + temp12);
		temp22 = h * dy2_test(_X.at[i][0] + h, _YRK.at[i][0] + temp11, _dYRK.at[i][0] + temp12);

		_YRK.at[i + 1][0] = _YRK.at[i][0] + (temp11 + temp21) / 2;
		_dYRK.at[i + 1][0] = _dYRK.at[i][0] + (temp12 + temp22) / 2;
		_X.at[i + 1][0] = _X.at[i][0] + h;
	}
}

void ODERK4sys2(Matrix _X, Matrix _YRK, Matrix _dYRK, double x0, double xf, double y0, double dy0, double h) {
	_X.at[0][0] = x0;
	_YRK.at[0][0] = y0;
	_dYRK.at[0][0] = dy0;

	int N = (int)((xf - x0) / h);

	double temp11 = 0;
	double temp12 = 0;
	double temp13 = 0;
	double temp14 = 0;

	double temp21 = 0;
	double temp22 = 0;
	double temp23 = 0;
	double temp24 = 0;

	for (int i = 0; i < N; i++) {
		temp11 = h * _dYRK.at[i][0];
		temp21 = h * dy2_test(_X.at[i][0], _YRK.at[i][0], _dYRK.at[i][0]);
		temp12 = h * (_dYRK.at[i][0] + temp21 / 2);
		temp22 = h * dy2_test(_X.at[i][0] + h / 2, _YRK.at[i][0] + temp11 / 2, _dYRK.at[i][0] + temp21 / 2);
		temp13 = h * (_dYRK.at[i][0] + temp22 / 2);
		temp23 = h * dy2_test(_X.at[i][0] + h / 2, _YRK.at[i][0] + temp12 / 2, _dYRK.at[i][0] + temp22 / 2);
		temp14 = h * (_dYRK.at[i][0] + temp23);
		temp24 = h * dy2_test(_X.at[i][0] + h, _YRK.at[i][0] + temp13, _dYRK.at[i][0] + temp23);

		_YRK.at[i + 1][0] = _YRK.at[i][0] + (temp11 + 2 * temp12 + 2 * temp13 + temp14) / 6;
		_dYRK.at[i + 1][0] = _dYRK.at[i][0] + (temp21 + 2 * temp22 + 2 * temp23 + temp24) / 6;
		_X.at[i + 1][0] = _X.at[i][0] + h;
	}
}

//Single-Sided DFT
Complex dft(Complex _Y) {
	int n = _Y.cols;

	Complex _W = createCmp(n);

	for (int i = 0; i < n; i++) {
		_W.re[i] = 0;
		_W.im[i] = 0;
		for (int k = 0; k < n; k++) {
			double theta = - 2 * PI * k * i / n;
			_W.re[i] += _Y.re[k] * cos(theta) - _Y.im[k] * sin(theta);
			_W.im[i] += _Y.re[k] * sin(theta) + _Y.im[k] * cos(theta);
		}
	}

	return _W;
}

//Cooley-Tukey Algorithm : must use power of 2
Complex fft(Complex _Y) {
	int n = _Y.cols;
	Complex _W = createCmp(n);

	if (n == 1) {
		_W.re[0] = _Y.re[0];
		_W.im[0] = _Y.im[0];
		return _W;
	}
	
	else {
		Complex _e = createCmp(n / 2);
		Complex _E = createCmp(n / 2);
		Complex _o = createCmp(n / 2);
		Complex _O = createCmp(n / 2);

		for (int k = 0; k < n / 2; k++) {
			_e.re[k] = _Y.re[2 * k];
			_e.im[k] = _Y.im[2 * k];
			_o.re[k] = _Y.re[2 * k + 1];
			_o.im[k] = _Y.im[2 * k + 1];
		}

		_E = fft(_e);
		_O = fft(_o);

		Complex _P = createCmp(n / 2);
		for (int k = 0; k < _P.cols; k++) {
			_P.re[k] = cos(-2 * PI * k / n);
			_P.im[k] = sin(-2 * PI * k / n);
		}

		_O = multiplyCmp(_O, _P);

		for (int k = 0; k < n / 2; k++) {
			_W.re[k] = _E.re[k] + _O.re[k];
			_W.im[k] = _E.im[k] + _O.im[k];
			_W.re[k + n / 2] = _E.re[k] - _O.re[k];
			_W.im[k + n / 2] = _E.im[k] - _O.im[k];
		}

		freeCmp(_o);
		freeCmp(_e);
		freeCmp(_O);
		freeCmp(_E);
		freeCmp(_P);

		return _W;
	}
}

Complex idft(Complex _Y) {
	int n = _Y.cols;

	Complex _W = createCmp(n);

	for (int i = 0; i < n; i++) {
		_W.re[i] = 0;
		_W.im[i] = 0;
		for (int k = 0; k < n; k++) {
			double theta = 2 * PI * k * i / n;
			_W.re[i] += _Y.re[k] * cos(theta) - _Y.im[k] * sin(theta);
			_W.im[i] += _Y.re[k] * sin(theta) + _Y.im[k] * cos(theta);
		}
	}

	for (int i = 0; i < n; i++) {
		_W.re[i] = _W.re[i] / n;
		_W.im[i] = _W.im[i] / n;
	}

	return _W;
}

Complex ifft(Complex _Y) {
	int n = _Y.cols;
	Complex _W = createCmp(n);

	if (n == 1) {
		_W.re[0] = _Y.re[0];
		_W.im[0] = _Y.im[0];
		return _W;
	}

	else {
		Complex _e = createCmp(n / 2);
		Complex _E = createCmp(n / 2);
		Complex _o = createCmp(n / 2);
		Complex _O = createCmp(n / 2);

		for (int k = 0; k < n / 2; k++) {
			_e.re[k] = _Y.re[2 * k];
			_e.im[k] = _Y.im[2 * k];
			_o.re[k] = _Y.re[2 * k + 1];
			_o.im[k] = _Y.im[2 * k + 1];
		}

		_E = fft(_e);
		_O = fft(_o);

		Complex _P = createCmp(n / 2);
		for (int k = 0; k < _P.cols; k++) {
			_P.re[k] = cos(-2 * PI * k / n);
			_P.im[k] = sin(-2 * PI * k / n);
		}

		_O = multiplyCmp(_O, _P);

		for (int k = 0; k < n / 2; k++) {
			_W.re[k] = _E.re[k] + _O.re[k];
			_W.im[k] = _E.im[k] + _O.im[k];
			_W.re[k + n / 2] = _E.re[k] - _O.re[k];
			_W.im[k + n / 2] = _E.im[k] - _O.im[k];

			_W.re[k] = _W.re[k] / n;
			_W.im[k] = _W.im[k] / n;
			_W.re[k + n / 2] = _W.re[k + n / 2] / n;
			_W.im[k + n / 2] = _W.im[k + n / 2] / n;
		}

		for (int k = 0; k < n / 2; k++) {
			_W.re[k] = _E.re[k] + _O.re[k];
			_W.im[k] = _E.im[k] + _O.im[k];
			_W.re[k + n / 2] = _E.re[k] - _O.re[k];
			_W.im[k + n / 2] = _E.im[k] - _O.im[k];

			_W.re[k] = _W.re[k] / n;
			_W.im[k] = _W.im[k] / n;
			_W.re[k + n / 2] = _W.re[k + n / 2] / n;
			_W.im[k + n / 2] = _W.im[k + n / 2] / n;

		}

		double temp_re;
		double temp_im;

		for (int i = 1; i < n / 2; i++) {	// Index Reverse
			temp_re = _W.re[i];
			_W.re[i] = _W.re[n - i];
			_W.re[n - i] = temp_re;

			temp_im = _W.im[i];
			_W.im[i] = _W.im[n - i];
			_W.im[n - i] = temp_im;
		}

		freeCmp(_o);
		freeCmp(_e);
		freeCmp(_O);
		freeCmp(_E);
		freeCmp(_P);


		return _W;
	}
}

/* Fourier Transform example

	clock_t start1, start2, end1, end2;
	float res1, res2;

	int L = 4096;												// Length of signal -> must be power of 2
	Vector _T = createVec(L);									// Time vector //(0:L - 1) * T;
	Complex _Y = createCmp(L);									// Signal vector //Y = 0.7 * sin(2 * pi * 50 * t) + sin(2 * pi * 120 * t);

	for (int i = 0; i < L; i++) {
		_T.at[i] = 0.001 * i;
		_Y.re[i] = 0.7 * sin(2 * PI * 50 * _T.at[i]) + sin(2 * PI * 120 * _T.at[i]);
		_Y.im[i] = 0;
	}

	start1 = clock();
	Complex result1 = dft(_Y);
	end1 = clock();
	res1 = (float)(end1 - start1) / CLOCKS_PER_SEC;
	
	start2 = clock();
	Complex result2 = fft(_Y);
	end2 = clock();
	res2 = (float)(end2 - start2) / CLOCKS_PER_SEC;

	printf(" dft : %.3fs \n", res1);
	printf(" fft : %.3fs \n", res2);

	freeCmp(result1);
	freeCmp(result2);
*/

/* Inverse Fourier Transform example
int L = 1000;												// Length of signal -> must be power of 2
Vector _T = createVec(L);									// Time vector //(0:L - 1) * T;
Complex _Y = createCmp(L);									// Signal vector //Y = 0.7 * sin(2 * pi * 50 * t) + sin(2 * pi * 120 * t);

for (int i = 0; i < L; i++) {
	_T.at[i] = 0.001 * i;
	_Y.re[i] = 0.7 * sin(2 * PI * 50 * _T.at[i]) + sin(2 * PI * 120 * _T.at[i]);
	_Y.im[i] = 0;
}

printf("\n======Original Singal======\n");
printCmp(_Y);

Complex result1 = dft(_Y);

printf("\n======dft Singal======\n");
printCmp(result1);

Complex result2 = idft(result1);

printf("\n======Recovery Singal======\n");
printCmp(result2);

freeCmp(result1);
freeCmp(result2);
*/