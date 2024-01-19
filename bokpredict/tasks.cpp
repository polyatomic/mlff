#include "tasks.h"

bool calculate_descriptors(int olength, double *r, FunctorDaDaI *dcalc, double *dval) {
   (*dcalc)(r, dval, 0);
   return true;
}
