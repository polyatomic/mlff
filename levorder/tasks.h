#pragma once

#define MAX_THREADS 16

bool order_by_leverages(int N, int M, const char *wmatc_name, const char *wmatct_name, int *ids_selected, double *leverages, int &rank, int bufsize, int nthreads);
