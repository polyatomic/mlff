#pragma once

#include "functor.h"

#define MAX_THREADS 16

bool calculate_descriptors(int N, int M, int olength, double *r, FunctorDaDaI *dcalc, double *stdevs, double *avgs, /*double *wmatc,*/ int nthreads, int maxrows, const char *outp_file, const char *stdevs_file = nullptr, const char *avgs_file = nullptr);
