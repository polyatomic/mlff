#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include <map>
#include <vector>
#include <thread>

#include "util.h"
#include "tasks.h"
#include "mathutil.h"

using std::cout;
using std::cerr;
using std::endl;
using std::sort;
using std::ifstream;
using std::ios;
using std::set;
using std::map;
using std::vector;
using std::thread;

// Calculates portion of Gram matrix
// X1 - nr1*N matrix slice
// X2 - nr2*N matrix slice
// K - nr*nr Gram matrix
// r1, r2 - starting indexes of X1 and X2 slices
void runKCalc(double *X1, double *X2, double *K, int r1, int r2, int nr1, int nr2, int N, int nr) {
   int i, j, ii, jj;
   double *dp1, *dp2;
   for (i = 0, ii = r1; i < nr1; i++, ii++) {
      dp1 = X1 + i*N;
      for (j = 0, jj = r2; j < nr2; j++, jj++) {
         dp2 = X2 + j*N;
         //K[ii*nr+jj] = dot_prod(dp1, dp2, N);
         K[ii*nr+jj] = dot_prod_kahan(dp1, dp2, N);
      }
   }
}

void runLeverageCalc(int jstart, int jend, int rowsum, int rank, int Ms, int ntot, double *levmat, double *wmatc, double *wwp, double *leverages) {
   int i, j;
   for (i = rowsum + jstart, j = jstart; j < jend; i++, j++) {
      MatByVec(levmat, wmatc + j*ntot, wwp, rank, Ms);
      leverages[i] = dot_prod(wwp, wwp, rank);
   }
}

bool order_by_leverages(int N, int M, const char *wmatc_name, const char *wmatct_name, int *ids_selected, double *leverages, int &rank, int bufsize, int nthreads) {
   int i, j, ntotsub, rowsum, cur_rows, nr, nc, matid, next_rowsum, rowsum2, cur_rows2,
       rslast, crlast, npairs, id, nr1, nr2, iM, ii, bufsize2, jstart, jend;
   double tmp, tol, ss;
   bool res = true;
   double *dp1, *dp2;
   double *Kds = 0;
   double **wmatct = 0;
   double *wmatc = 0;
   double *evec = 0;
   double *eval = 0;
   double *sevec = 0;
   double *seval = 0;
   double *levmat = 0;
   double **wwp = 0;
   int *tsizes = 0;
   bool next_init = false;
   set<int> ids;
   map<int, int> matm, matmr;
   vector<int> ids1, ids2;
   thread threads[MAX_THREADS];
   cout << "#####   Starting order_by_leverages   #####\n";
   if (nthreads > 1) cout << "Using " << nthreads << " threads" << endl;
   cout << "Buffer size: " << bufsize << endl;
   cout << "Number of rows: " << N << endl;
   cout << "Number of columns: " << M << endl;
   ifstream ifile;
   ifstream ifilet(wmatct_name, ios::in | ios::binary);
   if (!ifilet) {
      cerr << "Cannot open " << wmatct_name << endl;
      return false;
   }
   Kds = new double[M*M];
   wmatct = new double*[nthreads+1];
   for (j=0; j < M*M; j++) Kds[j] = 0.0;
   bufsize2 = nthreads*bufsize;
   wmatc = new double[M*bufsize2];
   for (i=0; i <= nthreads; i++) wmatct[i] = new double[N*bufsize];
   ntotsub = M;
   rowsum = 0;
   npairs = 0;
   id = 0;
   matid = 0;
   rslast = -1;
   crlast = 0;
   for (;;) {
      ifilet.seekg(0, ifilet.beg);
      if (ntotsub - rowsum < bufsize)
         cur_rows = ntotsub - rowsum;
      else
         cur_rows = bufsize;
      ifilet.read((char*)&nr, sizeof(int));
      ifilet.read((char*)&nc, sizeof(int));
      if (nr != M || nc != N) {
         cerr << "Incorrect number of rows or columns" << endl;
         res = false;
         goto end;
      }
      ifilet.seekg(((long long)rowsum)*N*sizeof(double), ifilet.cur);
      dp1 = wmatct[matid];
      for (i = 0; i < cur_rows; i++) {
         for (j = 0; j < N; j++) ifilet.read((char*)&dp1[i*N+j], sizeof(double));
      }
      if (matid != 0) {
         next_init = true;
         next_rowsum = rowsum;
      }
      matid++;
      rowsum2 = rowsum;
      for (;;) {
         if (ntotsub - rowsum2 < bufsize) {
            cur_rows2 = ntotsub - rowsum2;
            rslast = rowsum2;
            crlast = cur_rows2;
         } else {
            cur_rows2 = bufsize;
         }
         if (rowsum2 > rowsum) {
            dp1 = wmatct[matid];
            for (i = 0; i < cur_rows2; i++) {
               for (j = 0; j < N; j++) ifilet.read((char*)&dp1[i*N+j], sizeof(double));
            }
            matid++;
         }
         npairs++;
         auto r = ids.insert(rowsum);
         if (r.second) {
            matm[id] = rowsum;
            matmr[rowsum] = id;
            ids1.push_back(id);
            id++;
         } else {
            ids1.push_back(matmr[*r.first]);
         }
         r = ids.insert(rowsum2);
         if (r.second) {
            matm[id] = rowsum2;
            matmr[rowsum2] = id;
            ids2.push_back(id);
            id++;
         } else {
            ids2.push_back(matmr[*r.first]);
         }
         if (npairs == nthreads || (rowsum+cur_rows == ntotsub && rowsum2+cur_rows2 == ntotsub)) {
            // Computing
            for (i = 0; i < npairs; i++) {
               nr1 = nr2 = bufsize;
               if (matm[ids1[i]] == rslast) nr1 = crlast;
               if (matm[ids2[i]] == rslast) nr2 = crlast;
               threads[i] = thread(runKCalc, wmatct[ids1[i]], wmatct[ids2[i]], Kds, matm[ids1[i]], matm[ids2[i]], nr1, nr2, N, ntotsub);
            }
            for (i=0; i < npairs; i++) threads[i].join();
            npairs = 0;
            if (rowsum2+cur_rows2 == ntotsub) matid = 0;
            else matid = 1;
            if (next_init) {
               next_init = false;
               if (matid != 0) {
                  dp1 = wmatct[0];
                  dp2 = wmatct[matmr[next_rowsum]];
                  nr2 = bufsize;
                  if (matm[matmr[next_rowsum]] == rslast) nr2 = crlast;
                  for (i=0; i < nr2*N; i++) dp1[i] = dp2[i];
               }
            }
            ids.clear();
            matm.clear();
            matmr.clear();
            ids1.clear();
            ids2.clear();
            id = 0;
         }
         rowsum2 += cur_rows2;
         if (rowsum2 == ntotsub) break;
      }
      rowsum += cur_rows;
      if (rowsum == ntotsub) break;
   }
   ifilet.close();
   for (i = 1; i < M; i++) {
      iM = i*M;
      for (j = 0; j < i; j++) Kds[iM+j] = Kds[j*M+i];
   }
   evec = new double[M*M];
   eval = new double[M];
   sevec = new double[M*M];
   seval = new double[M];
   levmat = new double[M*M];
   wwp = new double*[nthreads];
   for (i = 0; i < nthreads; i++) wwp[i] = new double[M];
   tsizes = new int[nthreads];
   ifile.open(wmatc_name, ios::in | ios::binary);
   if (!ifile) {
      cerr << "Cannot open " << wmatc_name << endl;
      res = false;
      goto end;
   }
   cout << "Starting eigendecomposition" << endl;
   cout << "eigs returned " << eigs(Kds, evec, eval, M) << endl;
   tol = 1.0e-15*M*sqrt(eval[M-1]);
   rank = 0;
   for (i = 0; i < M; i++) {
      if (eval[i] > 0.0) {
         tmp = sqrt(eval[i]);
         if (tmp > tol) rank++;
      }
   }
   cout << "rank = " << rank << endl;
   cout << "Finished eigendecomposition" << endl;
   for (ii = 0; ii < M; ii++) {
      i = M-ii-1;
      for (j = 0; j < M; j++) sevec[ii*M+j] = evec[i*M+j];
      seval[ii] = eval[i];
   }
   for (i = 0; i < M*M; i++) levmat[i] = 0.0;
   for (i = 0; i < rank; i++) {
      ss = sqrt(seval[i]);
      dp1 = levmat + i*M;
      dp2 = sevec + i*M;
      for (j = 0; j < M; j++) dp1[j] = dp2[j]/ss;
   }
   ifile.read((char*)&nr, sizeof(int));
   ifile.read((char*)&nc, sizeof(int));
   if (nr != N || nc != M) {
      cerr << "Incorrect number of rows or columns" << endl;
      res = false;
      goto end;
   }
   rowsum = 0;
   for (;;) {
      if (N - rowsum < bufsize2)
         cur_rows = N - rowsum;
      else
         cur_rows = bufsize2;
      divide(tsizes, cur_rows, nthreads);
      for (i = 0; i < cur_rows; i++) {
         for (j = 0; j < M; j++) ifile.read((char*)&wmatc[i*M+j], sizeof(double));
      }
      jstart = 0;
      jend = tsizes[0];
      for (i = 0; i < nthreads; i++) {
         threads[i] = thread(runLeverageCalc, jstart, jend, rowsum, rank, M, M, levmat, wmatc, wwp[i], leverages);
         jstart = jend;
         if (i < nthreads-1) jend = jstart + tsizes[i+1];
      }
      for (i=0; i < nthreads; i++) threads[i].join();
//      for (i = rowsum, j=0; j < cur_rows; i++, j++) {
//         MatByVec(levmat, wmatc + j*ntot + start, wwp[0], rank, Ms);
//         leverages[k][i] = dot_prod(wwp[0], wwp[0], rank);
//      }
      rowsum += cur_rows;
      if (rowsum == N) break;
   }
   for (i = 0; i < N; i++) ids_selected[i] = i;
   sort(ids_selected, ids_selected+N, [leverages](int i, int j) -> bool { return leverages[i] > leverages[j]; });
   cout << "#####   Finished order_by_leverages   #####\n";
   ifile.close();
end:
   delete [] tsizes;
   if (wwp) {
      for (i = 0; i < nthreads; i++) delete [] wwp[i];
   }
   delete [] wwp;
   delete [] levmat;
   delete [] seval;
   delete [] sevec;
   delete [] eval;
   delete [] evec;
   for (i=0; i <= nthreads; i++) delete [] wmatct[i];
   delete [] wmatc;
   delete [] wmatct;
   delete [] Kds;
   return res;
}
