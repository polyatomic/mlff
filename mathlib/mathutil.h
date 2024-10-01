#pragma once

#include "mkl.h"

#define IDX(n,i,j) (2*n-i-1)*i/2+j-i-1

double sq_dist(const double *v1, const double *v2, int n);
double rp_3_6_kernel(double x1, double x2);
double rp_3_5_kernel(double x1, double x2);
double rp_3_4_kernel(double x1, double x2);
double rp_3_3_kernel(double x1, double x2);
double rp_3_2_kernel(double x1, double x2);
double rp_3_1_kernel(double x1, double x2);
double rp_3_0_kernel(double x1, double x2);
double symmetrized_kernel3_0(double *x1, double *x2, double (*kernel)(double, double));
double symmetrized_kernel3_1(double *x1, double *x2, double (*kernel)(double, double));
double symmetrized_kernel3_2(double *x1, double *x2, double (*kernel)(double, double));
double symmetrized_kernel4_0(double *x1, double *x2, double (*kernel)(double, double));
double symmetrized_kernel4_1(double *x1, double *x2, double (*kernel)(double, double));
double symmetrized_kernel4_2(double *x1, double *x2, double (*kernel)(double, double));
double symmetrized_kernel4_3(double *x1, double *x2, double (*kernel)(double, double));
void scale(const double *mat, int m, int n, double *cmat, double *avgs, double *stdevs);
void transpose(const double *mat_in, double *mat_out, int m, int n);
double dot_prod(const double *v1, const double *v2, int n);
double dot_prod_blas(const double *v1, const double *v2, int n);
double dot_prod_kahan(const double *v1, const double *v2, int n);
bool eigs(const double m[], double v[], double e[], int n);
void MatByVec(const double a[], const double b[], double c[], int m, int n);
bool svd(const double mt[], double u[], double sigma[], double vt[], lapack_int m, lapack_int n, bool transpose);
bool multiple_linear_regression_svd(const int M, const int N, const double *X, const double *Y, double *B, int& d, bool transpose);
double find_prediction_rss(const double *X, const double *u, const double *s, const double *vt, const int M, const int N, const double *Y, double *utb, double *y, double *B);
bool is_singular(const int M, const int N, const double *s);
