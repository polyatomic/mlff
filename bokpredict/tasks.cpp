#include "tasks.h"

bool calculate_descriptors(int olength, double *r, FunctorDaDaI *dcalc, double *dval, int batch) {
   (*dcalc)(r, dval, batch);
   return true;
}
