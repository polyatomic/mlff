#pragma once

#include "functor.h"

#define MAX_THREADS 16

bool calculate_descriptors(int olength, double *r, FunctorDaDaI *dcalc, double *dval, int batch);
