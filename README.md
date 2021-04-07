# Numerical Method [![Build Status][travis-image]][travis-url]
The Numerical Method Library (NM) is a self-contained C++ library for engineering mathematics.

## Installation
Version : C++17<br>
IDE : [Microsoft Visual Studio Community 2019](https://visualstudio.microsoft.com/ko/vs/community/)

## Usage
It is possible to use libraries that directly implement numerical method.<br>
myMain.cpp : This is a main function file for packages.<br>
myComplex.cpp : Complex functions are implemented.<br>
myMath.cpp : Numerical Linear Algebra functions are implemented.<br>
myMatrix.cpp : Matrix functions are implemented.<br>
myComplex.h : It works as a header file of 'myComplex.cpp'.
myMath.h : It works as a header file of 'myMath.cpp'.
myMatrix.h : It works as a header file of 'myMatrix.cpp'.

## Contents
1. Bisection method
2. Newton-Raphson method
3. Gauss Elimination method
4. Gauss-Jordan Elimination method
5. LU decomposition method
6. Inverse matrix solving method
7. QR decomposition for Eigenvalues and Eigenvector
8. Linear polynomial regression method
9. Numerical differentiation method(forward divided method, backward divided method, center divided method)
10. Numerical integration method(Trapezoidal rule, Simpson method)
11. Ordinary differential equation method(Euler method, Runge–Kutta method)
12. Digital signal processing(DFT, FFT(Cooley-Tukey Algorithm))


<!-- Markdown link & img dfn's -->
[travis-image]: https://img.shields.io/travis/dbader/node-datadog-metrics/master.svg?style=flat-square
[travis-url]: https://travis-ci.org/dbader/node-datadog-metrics
