//-----------------------------------------------------------------------
// Created  : 2018.03.26
// Modified : 2019.12.16
// Version  : Microsoft Visual Studio 2019
// Title    : myMath.h
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------

#ifndef		_MY_MATH_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATH_H
//-----------------------------------------------------------------------
// Header file declaration
//-----------------------------------------------------------------------
#define _CRT_SECURE_NO_WARNINGS // A way to avoid a vaccine program

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <time.h>
#include <windows.h>
#include "myMatrix.h"
#include "myComplex.h"

extern Matrix addMat(Matrix _A, Matrix _B);

extern Matrix subtractMat(Matrix _A, Matrix _B);

extern Matrix multiplyMat(Matrix _A, Matrix _B);

extern Matrix multiplysMat(Matrix _A, double c);

extern Vector addVec(Vector _A, Vector _B);

extern Vector subtractVec(Vector _A, Vector _B);

extern Complex addCmp(Complex _A, Complex _B);

extern Complex subtractCmp(Complex _A, Complex _B);

extern Complex multiplyCmp(Complex _A, Complex _B);

extern Complex multiplysMat(Complex _A, double w);

extern double dot(Vector _A, Vector _B);

extern double myfunc(double w);

extern double myfuncDerv(double w);

extern void bisecRoot(double a, double b, int maxlimit, double Xn, double tolf);

extern void newtonRoot(double Xi, int maxlimit, double Xn, double tolf);

extern void gaussjElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn);

extern void backsubj(Matrix _U, Matrix _bn, Matrix _x);

extern void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _bn);

extern void backsub(Matrix _U, Matrix _bn, Matrix _x);

extern void forwsub(Matrix _L, Matrix _bn, Matrix _y);

extern void backsubl(Matrix _U, Matrix _y, Matrix _x);

extern void partialpivot(Matrix _A, Matrix _b, Matrix _S);

extern void scalefactor(Matrix _A, Matrix _b, Matrix _S);

extern void LUdcmp(Matrix _A, Matrix _L, Matrix _U);

extern Matrix solveLinear(Matrix _A, Matrix _b, char method);

extern Matrix inverseMat(Matrix _A);

extern Matrix QRdcmph(Matrix _A, Matrix _Q, Matrix _R);

extern void eigVec(Matrix _A, Matrix _E, Matrix _V);

extern Matrix normalizeVec(Matrix _V);

extern Matrix inverseMatwc(Matrix _A);

extern Matrix cofactor(Matrix _A);

extern int rankMat(Matrix _A);

extern double factorial(int n);

extern void polyfit(Matrix _X, Matrix _Y, Matrix _z, int n);

extern Matrix deriv(Matrix _X, Matrix _Y, double h, char method);

extern double deriv(double x, double myfunc(double x), double h, char method);

extern Matrix deriv2(Matrix _X, Matrix _Y, double h, char method);

extern double deriv2(double x, double myfunc(double x), double h, char method);

extern double myfunc1(double x);

extern double myfunc2(double x);

extern double trapz(Matrix _X, Matrix _Y);

extern double trapz(double myfunc(double x), double a, double b, int N);

extern double simpson13(Matrix _X, Matrix _Y);

extern double simpson13(double myfunc(double x), double a, double b, int N);

extern void ODEeuler(Matrix _X, Matrix _YE, double x0, double xf, double y0, double h);

extern void ODEeulerM(Matrix _X, Matrix _YEM, double x0, double xf, double y0, double h);

extern void ODEmid(Matrix _X, Matrix _YM, double x0, double xf, double y0, double h);

extern void ODERK2sys1(Matrix _X, Matrix _YRK, double x0, double xf, double y0, double h);

extern void ODERK4sys1(Matrix _X, Matrix _YRK, double x0, double xf, double y0, double h);

extern void ODERK2sys2(Matrix _X, Matrix _YRK, Matrix _dYRK, double x0, double xf, double y0, double dy0, double h);

extern void ODERK4sys2(Matrix _X, Matrix _YRK, Matrix _dYRK, double x0, double xf, double y0, double dy0, double h);

extern Complex dft(Complex _Y);

//must use power of 2
extern Complex fft(Complex _Y);

extern Complex idft(Complex _Y);

extern Complex ifft(Complex _Y);
#endif