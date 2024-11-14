#include <iostream>
#include <thread>
#include <string>
#include <fstream>

#include "util.h"
#include "tasks.h"

using std::cout;
using std::endl;
using std::thread;
using std::string;
using std::ofstream;
using std::ios;

void runFunctorDaDaI(int ivars[], double *dvars[], FunctorDaDaI *dcalc) {
   int j;
   double *rcoords, *wp;
   int istart = ivars[0];
   int iend = ivars[1];
   int olength = ivars[2];
   int nidi = ivars[3];
   int wsi = ivars[4];
   int ii = ivars[5];
   double *rp = dvars[0];
   double *wmat = dvars[1];
   for (j = istart; j < iend; j++,ii++) {
      rcoords = rp + j*olength;
      wp = wmat + ii*nidi;
      (*dcalc)(rcoords, wp, wsi);
   }
}

bool calculate_descriptors(int N, int M, int olength, double *r, FunctorDaDaI *dcalc, double *stdevs, double *avgs, /*double* wmatc,*/ int nthreads, int maxrows, const char *outp_file, const char *stdevs_file, const char *avgs_file) {
   /*int i;
   double *wp, *wmat, *rcoords;
   cout << "#####   Starting calculate_descriptors   #####\n";
   if (nthreads > 1) cout << "Using " << nthreads << " threads" << endl;
   wmat = new double[M*N];
   for (i = 0; i < N; i++) {
      wp = wmat + i*M;
      rcoords = r + i*olength;
      (*dcalc)(rcoords, wp, 0);
   }
   scale(wmat, M, N, wmatc, avgs, stdevs);
   //exit(0);
   cout << "#####   Finished calculate_descriptors   #####\n";
   delete [] wmat;*/
   int tc, jstart, jend, cur_rows, rowsum, ii, jj, matcount, i, j, nr, nc, cur_cols, colsum;
   int *tsizes;
   double s;
   double *wmat, *wmatc, *wmatct;
   long long il;
   string fn;
   thread threads[MAX_THREADS];
   int ivars[MAX_THREADS][6];
   double *dvars[MAX_THREADS][2];
   cout << "#####   Starting calculate_descriptors   #####\n";
   if (nthreads > 1) cout << "Using " << nthreads << " threads" << endl;
   cout << "Row buffer size is: " << maxrows << endl;
   tsizes = new int[nthreads];
   wmat = new double[M*maxrows];
   wmatc = new double[M*maxrows];
   wmatct = new double[((long long)N)*maxrows];
   if (stdevs_file) {
      File2Matrix(avgs_file, i, j, avgs);
      File2Matrix(stdevs_file, i, j, stdevs);
   } else {
      for (i=0; i < M; i++) {
         avgs[i] = 0.0;
         stdevs[i] = 0.0;
      }
   }
   rowsum = 0;
   matcount = 0;
   for (;;) {
      if (N - rowsum < maxrows)
         cur_rows = N - rowsum;
      else
         cur_rows = maxrows;
      divide(tsizes, cur_rows, nthreads);
      jstart = 0;
      jend = tsizes[0];
      for (tc = 0, ii = 0; tc < nthreads; tc++) {
         ivars[tc][0] = rowsum + jstart;
         ivars[tc][1] = rowsum + jend;
         ivars[tc][2] = olength;
         ivars[tc][3] = M;
         ivars[tc][4] = tc;
         ivars[tc][5] = ii;
         dvars[tc][0] = r;
         dvars[tc][1] = wmat;
         threads[tc] = thread(runFunctorDaDaI, ivars[tc], dvars[tc], dcalc);
         ii += tsizes[tc];
         jstart = jend;
         if (tc < nthreads-1) jend = jstart + tsizes[tc+1];
      }
      for (tc = 0; tc < nthreads; tc++)
         threads[tc].join();
      fn = "_mat." + std::to_string(matcount) + ".bin";
      Matrix2File(wmat, cur_rows, M, fn.c_str());
      if (!stdevs_file) {
         for (i = 0; i < cur_rows; i++) {
            for (j = 0; j < M; j++) {
               avgs[j] += wmat[i*M+j];
            }
         }
      }
      rowsum += cur_rows;
      if (rowsum == N) break;
      matcount++;
   }
   if (!stdevs_file) {
      for (i=0; i < M; i++) avgs[i] = avgs[i]/N;
   }
   rowsum = 0;
   matcount = 0;
   for (;;) {
      if (N - rowsum < maxrows)
         cur_rows = N - rowsum;
      else
         cur_rows = maxrows;
      fn = "_mat." + std::to_string(matcount) + ".bin";
      File2Matrix(fn.c_str(), nr, nc, wmat);
      if (!stdevs_file) {
         for (i = 0; i < cur_rows; i++) {
            for (j = 0; j < M; j++) {
               s = wmat[i*M+j] - avgs[j];
               stdevs[j] += s*s;
            }
         }
      }
      rowsum += cur_rows;
      if (rowsum == N) break;
      matcount++;
   }
   if (!stdevs_file) {
      for (i = 0; i < M; i++) stdevs[i] = sqrt(stdevs[i]/N);
   }
   rowsum = 0;
   matcount = 0;
   for (;;) {
      if (N - rowsum < maxrows)
         cur_rows = N - rowsum;
      else
         cur_rows = maxrows;
      fn = "_mat." + std::to_string(matcount) + ".bin";
      File2Matrix(fn.c_str(), nr, nc, wmat);
      remove(fn.c_str());
      for (i = 0; i < cur_rows; i++) {
         for (j = 0; j < M; j++) {
            wmatc[i*M+j] = (wmat[i*M+j] - avgs[j])/stdevs[j];
         }
      }
      fn = "_matc." + std::to_string(matcount) + ".bin";
      Matrix2File(wmatc, cur_rows, M, fn.c_str());
      rowsum += cur_rows;
      if (rowsum == N) break;
      matcount++;
   }
   fn = string(outp_file) + ".bin";
   ofstream ofile(fn.c_str(), ios::out | ios::binary);
   nr = N;
   nc = M;
   ofile.write((char*)&nr, sizeof(int));
   ofile.write((char*)&nc, sizeof(int));
   rowsum = 0;
   matcount = 0;
   for (;;) {
      if (N - rowsum < maxrows)
         cur_rows = N - rowsum;
      else
         cur_rows = maxrows;
      fn = "_matc." + std::to_string(matcount) + ".bin";
      File2Matrix(fn.c_str(), nr, nc, wmatc);
//      remove(fn.c_str());
      for (i = 0; i < cur_rows; i++) {
         for (j = 0; j < M; j++)
            ofile.write((char*)&wmatc[i*M+j], sizeof(double));
      }
      rowsum += cur_rows;
      if (rowsum == N) break;
      matcount++;
   }
   ofile.close();
   // cerr << "Press enter to continue:";
   // std::cin.get();
   fn = string(outp_file) + "t.bin";
   ofile.open(fn.c_str(), ios::out | ios::binary);
   nr = M;
   nc = N;
   ofile.write((char*)&nr, sizeof(int));
   ofile.write((char*)&nc, sizeof(int));
   colsum = 0;
   for (;;) {
      rowsum = 0;
      matcount = 0;
      if (M - colsum < maxrows)
         cur_cols = M - colsum;
      else
         cur_cols = maxrows;
      for (;;) {
         if (N - rowsum < maxrows)
            cur_rows = N - rowsum;
         else
            cur_rows = maxrows;
         fn = "_matc." + std::to_string(matcount) + ".bin";
         // Read cur_rows*M matrix into wmatc
         File2Matrix(fn.c_str(), nr, nc, wmatc);
         // Fill columns rowsum until rowsum + cur_rows of wmatct
         for (i=0, ii = colsum; i < cur_cols; i++,ii++) {
            il = i*((long long)N);
            for (j=0, jj = rowsum; j < cur_rows; j++,jj++) {
               wmatct[il+jj] = wmatc[j*M+ii];
            }
         }
         rowsum += cur_rows;
         if (rowsum == N) break;
         matcount++;
      }
      for (il=0; il < ((long long)N)*cur_cols; il++) ofile.write((char*)&wmatct[il], sizeof(double));
      colsum += cur_cols;
      if (colsum == M) break;
   }
   ofile.close();
   rowsum = 0;
   matcount = 0;
   for (;;) {
      if (N - rowsum < maxrows)
         cur_rows = N - rowsum;
      else
         cur_rows = maxrows;
      fn = "_matc." + std::to_string(matcount) + ".bin";
      remove(fn.c_str());
      rowsum += cur_rows;
      if (rowsum == N) break;
      matcount++;
   }
   cout << "#####   Finished calculate_descriptors   #####\n";
   delete [] wmatct;
   delete [] wmatc;
   delete [] wmat;
   delete [] tsizes;
   return true;
}
