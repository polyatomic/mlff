#pragma once

#define MAX_THREADS 16

class FunctorDaDaI {
public:
   virtual void operator()(double *darg1, double *darg2, int iarg) = 0;
   virtual void Call(double *darg1, double *darg2, int iarg) = 0;
};

template<class FClass>
class TFunctorDaDaI : public FunctorDaDaI {
public:
   TFunctorDaDaI() {
      pt2Object = 0;
      fpt = 0;
   }

   TFunctorDaDaI(FClass *_pt2Object, void(FClass::*_fpt)(double *darg1, double *darg2, int iarg)) {
      pt2Object = _pt2Object;
      fpt = _fpt;
   }

   void init(FClass *_pt2Object, void(FClass::*_fpt)(double *darg1, double *darg2, int iarg)) {
      pt2Object = _pt2Object;
      fpt = _fpt;
   }

   virtual void operator()(double *darg1, double *darg2, int iarg) {
      (*pt2Object.*fpt)(darg1, darg2, iarg);
   }

   virtual void Call(double *darg1, double *darg2, int iarg) {
      (*pt2Object.*fpt)(darg1, darg2, iarg);
   }

private:
   void (FClass::*fpt)(double *darg1, double *darg2, int iarg);
   FClass *pt2Object;
};

bool calculate_descriptors(int N, int M, int olength, double *r, FunctorDaDaI *dcalc, double *stdevs, double *avgs, /*double *wmatc,*/ int nthreads, int maxrows);
bool order_by_leverages(int N, int M, double *wmatc, int *ids_selected, double *leverages);
bool order_by_leverages(int N, int M[], int nsub, const char *wmatc_name, const char *wmatct_name, int *ids_selected[], double *leverages[], int ranks[], int bufsize, int nthreads);
void sketch_matrices(int N, int M[], int nsub, const char *wmatc_name, int *ids_selected[], int ranks[], int n_test_rows, int nthreads, int *selected_variables[]);
bool calculate_linear_regression(int N, int M, double *wmatc, int *ids, int nt, double *Yc);
bool calculate_linear_regression(int N, int M, int Mtot, const char *wmatc_name, int *ids, int nt, double *Yc, double& rmse_train, double& rmse_test, double& rmse_total, double& mae_train, double& mae_test, double& mae_total, int& ndim, bool calculate_test_stats = true, bool print_data = true);
bool build_regression(int N, int M[], int Mtot, int substart, int nsub, const char *wmatc_name, double *Yc, int *ids_selected[], int n_incr, int n_steps, int& n_training, int *ids_training);
