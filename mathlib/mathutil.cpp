#include <cmath>
#include <iostream>

#include "mathutil.h"

using std::cerr;
using std::endl;

static const double g_3_over_56 = 3.0/56.0;
static const double g_4_over_3 = 4.0/3.0;
static const double g_7_over_15 = 7.0/15.0;
static const double g_3_over_35 = 3.0/35.0;
static const double g_5_over_12 = 5.0/12.0;
static const double g_8_over_7 = 8.0/7.0;
static const double g_5_over_14 = 5.0/14.0;
static const double g_2_over_7 = 2.0/7.0;
static const double g_1_over_28 = 1.0/28.0;
static const double g_28_over_55 = 28.0/55.0;

double sq_dist(const double *v1, const double *v2, int n) {
   int i;
   double tmp, sum;
   for (i = 0, sum = 0.0; i < n; i++) {
      tmp = v1[i] - v2[i];
      tmp *= tmp;
      sum += tmp;
   }
   return sum;
}

double rp_3_6_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = xl*xl*xl;
   t1 = t1*t1;
   t1 = t1*xl;
   t1 = 1.0/t1;
   t2 = xs/xl;
   return g_1_over_28*t1*(1.0 - t2*(1.4 - g_28_over_55*t2));
}

double rp_3_5_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = xl*xl;
   t1 = t1*t1*t1;
   t1 = 1.0/t1;
   t2 = xs/xl;
   return g_3_over_56*t1*(1.0 - t2*(g_4_over_3 - g_7_over_15*t2));
}

double rp_3_4_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = xl*xl;
   t1 = t1*t1*xl;
   t1 = 1.0/t1;
   t2 = xs/xl;
   return g_3_over_35*t1*(1.0 - t2*(1.25 - g_5_over_12*t2));
}

double rp_3_3_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = xl*xl;
   t1 = t1*t1;
   t1 = 1.0/t1;
   t2 = xs/xl;
   return 0.15*t1*(1.0 - t2*(g_8_over_7 - g_5_over_14*t2));
}

double rp_3_2_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = xl*xl*xl;
   t1 = 1.0/t1;
   t2 = xs/xl;
   return 0.3*t1*(1.0 - t2*(1.0 - g_2_over_7*t2));
}

double rp_3_1_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = xl*xl;
   t1 = 1.0/t1;
   t2 = xs/xl;
   return 0.75*t1*(1.0 - t2*(0.8 - 0.2*t2));
}

double rp_3_0_kernel(double x1, double x2) {
   double xl, xs, t1, t2;
   if (x1 < x2) {
      xs = x1;
      xl = x2;
   } else {
      xs = x2;
      xl = x1;
   }
   t1 = 1.0/xl;
   t2 = xs/xl;
   return 3.0*t1*(1.0 - t2*(0.5 - 0.1*t2));
}

// Symmetrized 3D kernel for AAA system.
double symmetrized_kernel3_0(double *x1, double *x2, double (*kernel)(double, double)) {
   double ret, k00, k11, k22, k01, k12, k21, k10, k20, k02;
   k00 = kernel(x1[0], x2[0]);
   k11 = kernel(x1[1], x2[1]);
   k22 = kernel(x1[2], x2[2]);
   k01 = kernel(x1[0], x2[1]);
   k12 = kernel(x1[1], x2[2]);
   k21 = kernel(x1[2], x2[1]);
   k10 = kernel(x1[1], x2[0]);
   k20 = kernel(x1[2], x2[0]);
   k02 = kernel(x1[0], x2[2]);
   ret = k00*(k11*k22 + k12*k21) +
         k01*(k10*k22 + k12*k20) +
         k02*(k10*k21 + k11*k20);
   return ret;
}

// Symmetrized 3D kernel for ABB system.
double symmetrized_kernel3_1(double *x1, double *x2, double (*kernel)(double, double)) {
   double ret;
   ret = (kernel(x1[0], x2[0])*kernel(x1[1], x2[1]) +
          kernel(x1[0], x2[1])*kernel(x1[1], x2[0]))*
          kernel(x1[2], x2[2]);
   return ret;
}

// Symmetrized 4D kernel for AAAA system.
double symmetrized_kernel4_0(double *x1, double *x2, double (*kernel)(double, double)) {
   double ret, k00, k11, k22, k33, k44, k55, k12, k21, k34, k43, k01, k10, k45, k54, k20, k35, k02, k53,
          k13, k24, k31, k42, k14, k23, k32, k41, k03, k52, k04, k51, k25, k30, k15, k40, k05, k50;
   k00 = kernel(x1[0], x2[0]);
   k11 = kernel(x1[1], x2[1]);
   k22 = kernel(x1[2], x2[2]);
   k33 = kernel(x1[3], x2[3]);
   k44 = kernel(x1[4], x2[4]);
   k55 = kernel(x1[5], x2[5]);
   k12 = kernel(x1[1], x2[2]);
   k21 = kernel(x1[2], x2[1]);
   k34 = kernel(x1[3], x2[4]);
   k43 = kernel(x1[4], x2[3]);
   k01 = kernel(x1[0], x2[1]);
   k10 = kernel(x1[1], x2[0]);
   k45 = kernel(x1[4], x2[5]);
   k54 = kernel(x1[5], x2[4]);
   k20 = kernel(x1[2], x2[0]);
   k35 = kernel(x1[3], x2[5]);
   k02 = kernel(x1[0], x2[2]);
   k53 = kernel(x1[5], x2[3]);
   k13 = kernel(x1[1], x2[3]);
   k24 = kernel(x1[2], x2[4]);
   k31 = kernel(x1[3], x2[1]);
   k42 = kernel(x1[4], x2[2]);
   k14 = kernel(x1[1], x2[4]);
   k23 = kernel(x1[2], x2[3]);
   k32 = kernel(x1[3], x2[2]);
   k41 = kernel(x1[4], x2[1]);
   k03 = kernel(x1[0], x2[3]);
   k52 = kernel(x1[5], x2[2]);
   k04 = kernel(x1[0], x2[4]);
   k51 = kernel(x1[5], x2[1]);
   k25 = kernel(x1[2], x2[5]);
   k30 = kernel(x1[3], x2[0]);
   k15 = kernel(x1[1], x2[5]);
   k40 = kernel(x1[4], x2[0]);
   k05 = kernel(x1[0], x2[5]);
   k50 = kernel(x1[5], x2[0]);
   ret = k00*k55*(k11*k22*k33*k44 +
                  k12*k21*k34*k43 +
                  k13*k24*k31*k42 +
                  k14*k23*k32*k41) +
         k01*k54*(k10*k22*k33*k45 +
                  k12*k20*k35*k43 +
                  k13*k25*k30*k42 +
                  k15*k23*k32*k40) +
         k02*k53*(k10*k21*k34*k45 +
                  k11*k20*k35*k44 +
                  k14*k25*k30*k41 +
                  k15*k24*k31*k40) +
         k03*k52*(k10*k24*k31*k45 +
                  k14*k20*k35*k41 +
                  k11*k25*k30*k44 +
                  k15*k21*k34*k40) +
         k04*k51*(k10*k23*k32*k45 +
                  k13*k20*k35*k42 +
                  k12*k25*k30*k43 +
                  k15*k22*k33*k40) +
         k05*k50*(k11*k23*k32*k44 +
                  k13*k21*k34*k42 +
                  k12*k24*k31*k43 +
                  k14*k22*k33*k41);
   return ret;
}

// Symmetrized 4D kernel for ABBB system.
double symmetrized_kernel4_1(double *x1, double *x2, double (*kernel)(double, double)) {
   double ret, k00, k11, k22, k33, k44, k55, k12, k21, k34, k43, k01, k10, k45, k54, k20, k35, k02, k53;
   k00 = kernel(x1[0], x2[0]);
   k11 = kernel(x1[1], x2[1]);
   k22 = kernel(x1[2], x2[2]);
   k33 = kernel(x1[3], x2[3]);
   k44 = kernel(x1[4], x2[4]);
   k55 = kernel(x1[5], x2[5]);
   k12 = kernel(x1[1], x2[2]);
   k21 = kernel(x1[2], x2[1]);
   k34 = kernel(x1[3], x2[4]);
   k43 = kernel(x1[4], x2[3]);
   k01 = kernel(x1[0], x2[1]);
   k10 = kernel(x1[1], x2[0]);
   k45 = kernel(x1[4], x2[5]);
   k54 = kernel(x1[5], x2[4]);
   k20 = kernel(x1[2], x2[0]);
   k35 = kernel(x1[3], x2[5]);
   k02 = kernel(x1[0], x2[2]);
   k53 = kernel(x1[5], x2[3]);
   ret = k00*k55*(k11*k22*k33*k44 +
                  k12*k21*k34*k43) +
         k01*k54*(k10*k22*k33*k45 +
                  k12*k20*k35*k43) +
         k02*k53*(k10*k21*k34*k45 +
                  k11*k20*k35*k44);
   return ret;
}

// Symmetrized 4D kernel for AABB system.
double symmetrized_kernel4_2(double *x1, double *x2, double (*kernel)(double, double)) {
   double ret;
   ret = kernel(x1[0], x2[0])*kernel(x1[5], x2[5])*
         (kernel(x1[1], x2[1])*kernel(x1[2], x2[2])*kernel(x1[3], x2[3])*kernel(x1[4], x2[4]) +
          kernel(x1[1], x2[2])*kernel(x1[2], x2[1])*kernel(x1[3], x2[4])*kernel(x1[4], x2[3]) +
          kernel(x1[1], x2[3])*kernel(x1[2], x2[4])*kernel(x1[3], x2[1])*kernel(x1[4], x2[2]) +
          kernel(x1[1], x2[4])*kernel(x1[2], x2[3])*kernel(x1[3], x2[2])*kernel(x1[4], x2[1]));
   return ret;
}

void scale(const double *mat, int m, int n, double *cmat, double *avgs, double *stdevs) {
   int i, j;
   double avg, stdev;
   // m - number of cols
   // n - number of rows
   for (i = 0; i < m; i++) {
      for (j = 0, avg = 0.0; j < n; j++) avg += mat[j*m+i];
      if (!avg) {
         avgs[i] = 0.0;
         stdevs[i] = 0.0;
         for (j = 0; j < n; j++) cmat[j*m+i] = 0.0;
         continue;
      }
      avg /= n;
      for (j = 0; j < n; j++) cmat[j*m+i] = mat[j*m+i] - avg;
      avgs[i] = avg;
      for (j = 0, stdev = 0.0; j < n; j++) stdev += cmat[j*m+i]*cmat[j*m+i];
      stdev = sqrt(stdev/n);
      stdevs[i] = stdev;
      for (j = 0; j < n; j++) cmat[j*m+i] /= stdev;
   }
}

void transpose(const double *mat_in, double *mat_out, int m, int n) {
   int i, j, in;
   for (i = 0; i < m; i++) {
      in = i*n;
      for (j = 0; j < n; j++) {
         mat_out[j*m+i] = mat_in[in+j];
      }
   }
}

// Dot product of n-dim vectors.
double dot_prod(const double *v1, const double *v2, int n) {
   int i;
   double ret = 0.0;
   for (i = 0; i < n; i++) ret += v1[i]*v2[i];
   return ret;
}

// Dot product of n-dim vectors with Kahan summation
double dot_prod_kahan(const double *v1, const double *v2, int n) {
   int i;
   double y, f, t;
   double ret = 0.0;
   double c = 0.0;
   for (i = 0; i < n; i++) {
      f = v1[i]*v2[i];
      y = f - c;
      t = ret + y;
      c = (t - ret) - y;
      ret = t;
   }
   return ret;
}

bool eigs(const double m[], double v[], double e[], int n) {
   int i, info;
   for (i = 0; i < n*n; i++) v[i] = m[i];
   info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, v, n, e);
   if (info > 0) return false;
   return true;
}

void MatByVec(const double a[], const double b[], double c[], int m, int n) {
   int i, j;
   const double *v;
   for (i = 0; i < m; i++) {
      v = a + i*n;
      c[i] = 0.0;
      for (j = 0; j < n; j++)
         c[i] += v[j]*b[j];
   }
}

bool svd(const double mt[], double u[], double sigma[], double vt[], lapack_int m, lapack_int n, bool transpose) {
   int i, j;
   long long il;
   bool ret = true;
   lapack_int info;
   lapack_int minmn = m;
   if (n < m) minmn = n;
   double *mtc = new double[m*((long long)n)];
   double *superb = new double[minmn-1];
   if (transpose) {
      for (i = 0; i < m; i++) {
         for (j = 0; j < n; j++)
            mtc[j*((long long)m)+i] = mt[i*((long long)n)+j];
      }
   } else {
      for (il = 0; il < ((long long)m)*n; il++) mtc[il] = mt[il];
   }
   info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, mtc, m, sigma, u, m, vt, minmn, superb);
   if (info > 0)
      ret = false;
   delete [] mtc;
   delete [] superb;
   return ret;
}

bool multiple_linear_regression_svd(const int M, const int N, const double *X, const double *Y, double *B, int& d, bool transpose) {
   // M - number of cases
   // N - number of variables
   int i, rank;
   bool ret = true;
   int maxdim, mindim;
   double tol;
   if (N > M) {
      maxdim = N;
      mindim = M;
      cerr << "Underdefined problem\n";
   } else {
      maxdim = M;
      mindim = N;
   }
   double *u = new double[mindim*M];
   double *vt = new double[mindim*N];
   double *s = new double[mindim];
   double *utb = new double[mindim];
   double *y = new double[mindim];
   bool r = svd(X, u, s, vt, M, N, transpose);
   if (!r) {
      ret = false;
      goto end;
   }
   tol = 1.0e-15*maxdim*s[0];
   for (i = 0; i < mindim; i++)
      utb[i] = dot_prod(u+i*M, Y, M);
   rank = mindim;
   for (i = 0; i < mindim; i++) {
      if (s[i] > tol) y[i] = utb[i]/s[i];
      else {
         y[i] = 0.0;
         rank--;
      }
   }
   d = rank;
   if (rank < mindim)
      cerr << "Singular matrix, rank = " << rank << endl;
   for (i = 0; i < N; i++)
      B[i] = dot_prod(vt+i*mindim, y, rank);
end:
   delete [] u;
   delete [] vt;
   delete [] s;
   delete [] utb;
   delete [] y;
   return ret;
}

double find_prediction_rss(const double *X, const double *u, const double *s, const double *vt, const int M, const int N, const double *Y, double *utb, double *y, double *B) {
   int i, maxdim, mindim, rank, d;
   double tol, rss, yp, tmp;
   const double *dp;
   if (N > M) {
      return 0.0;
   }
   maxdim = M;
   mindim = N;
   tol = 1.0e-15*maxdim*s[0];
   for (i = 0; i < mindim; i++)
      utb[i] = dot_prod(u+i*((long long)M), Y, M);
   rank = mindim;
   for (i = 0; i < mindim; i++) {
      if (s[i] > tol) y[i] = utb[i]/s[i];
      else {
         y[i] = 0.0;
         rank--;
      }
   }
   d = rank;
   if (rank < mindim)
      cerr << "Singular matrix, rank = " << rank << endl;
   for (i = 0; i < N; i++)
      B[i] = dot_prod(vt+i*mindim, y, rank);
   for (i = 0, rss = 0.0; i < M; i++) {
      dp = X + i*((long long)N);
      yp = dot_prod(dp, B, N);
      tmp = Y[i] - yp;
      rss += tmp*tmp;
   }
   return rss;
   // return sqrt(rss/M);
}

bool is_singular(const int M, const int N, const double *s) {
   double tol;
   int i, rank;
   tol = 1.0e-15*M*s[0];
   rank = N;
   for (i = 0; i < N; i++) {
      if (s[i] <= tol) rank--;
   }
   return rank < N;
}
