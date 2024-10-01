#include <iostream>
#include <iomanip>
#include <thread>

#include "RunTasks.h"
#include "tasks.h"

extern PTypeExtended g_params;

using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;
using std::thread;

void calculate_descriptors_for_rows(int n, int ndesc, double *v, double *B, double *avgs, double avg, double *stdevs, FunctorDaDaI* dcalc, double *dval, double *Yp, int istart, int iend, int batch) {
   int i, j;
   double sum;
   double *r;
   for (i = istart; i < iend; i++) {
      r = v + i*n;
      calculate_descriptors(n, r, dcalc, dval, batch);
      for (sum = 0.0, j = 0; j < ndesc; j++) {
         dval[j] = (dval[j] - avgs[j])/stdevs[j];
         sum += dval[j]*B[j];
      }
      Yp[i] = sum + avg;
   }
}

void RunTasks() {
   int i, nstr, n, ndesc, nr, nc, istart, iend, nthreads, tc;
   double *v, *B, *stdevs, *avgs, *Y, *Yp;
   double *dvals[MAX_THREADS];
   int *tsizes;
   int nblocks[3];
   char *blocks;
   double avg, rss, dif, rss2;
   thread threads[MAX_THREADS];
   Iterator<Mol> mols;
   stdevs = avgs = B = Y = Yp = v = 0;
   tsizes = 0;
   nthreads = 0;
   Descriptors descs;
   TFunctorDaDaI<Descriptors> dcalc;
   string s = g_params.linearCoeffs;
   if (!calculation_prepare(&g_params, mols, descs, nstr, n, v, Y, nblocks, &blocks)) goto end;
   ndesc = descs.GetNDescriptors();
   cout << "Number of descriptors: " << ndesc << endl;
   dcalc.init(&descs, &Descriptors::Calculate);
   B = new double[ndesc+1];
   if (s.size() > 0) {
      File2Matrix(s.c_str(), nr, nc, B);
   } else {
      cerr << "Descriptor coeffs undefined" << endl;
      goto end;
   }
   stdevs = new double[ndesc];
   avgs = new double[ndesc];
   File2Matrix("avgs.bin", nr, nc, avgs);
   File2Matrix("stdevs.bin", nr, nc, stdevs);
   avg = g_params.bias;
   Yp = new double[nstr];
   nthreads = g_params.nthreads;
   cout << "Using " << nthreads << " threads" << endl;
   tsizes = new int[nthreads];
   for (i=0; i < nthreads; i++) dvals[i] = new double[ndesc];
   divide(tsizes, nstr, nthreads);
   istart = 0;
   iend = tsizes[0];
   for (tc = 0; tc < nthreads; tc++) {
      threads[tc] = thread(calculate_descriptors_for_rows, n, ndesc, v, B, avgs, avg, stdevs, &dcalc, dvals[tc], Yp, istart, iend, tc);
      istart = iend;
      if (tc < nthreads-1) iend = istart + tsizes[tc+1];
   }
   for (tc = 0; tc < nthreads; tc++)
      threads[tc].join();
   for (i=0, rss=0.0, rss2=0.0; i < nstr; i++) {
      cout << Yp[i] << endl;
      dif = Y[i] - Yp[i];
      rss += dif*dif;
      rss2 += fabs(dif);
   }
   cout << "rmse: " << sqrt(rss/nstr) << endl;
   cout << "mae: " << rss2/nstr << endl;
end:
   for (i = 0; i < nthreads; i++) delete [] dvals[i];
   delete [] tsizes;
   delete [] Yp;
   delete [] avgs;
   delete [] stdevs;
   delete [] B;
   delete [] blocks;
   delete [] Y;
   delete [] v;
}
