//-----------------------------------------------------------------------
// Created  : 2018.03.26
// Modified : 2019.12.16
// Version  : Microsoft Visual Studio 2019
// Title    : myMain.cpp
//-----------------------------------------------------------------------
// Copyright by Titus Choi
//-----------------------------------------------------------------------


#define _CRT_SECURE_NO_WARNINGS // A way to avoid a vaccine program
#include "myMath.h"

#define PI (double)3.14159265359
#define N  256						// Multiple of 2

int main(){
	/*==========================================================================*/
	/*									Main Run								*/
	/*==========================================================================*/

	clock_t start1, start2, end1, end2;
	float res1, res2;

	int L =4096;												// Length of signal -> must be power of 2
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
	

	printf("Press any key to terminate your program.\n");
	return 0;
}
