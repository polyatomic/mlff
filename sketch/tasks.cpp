#include <iostream>
#include <fstream>
#include <algorithm>
#include <thread>
#include <set>
#include <mutex>

#include "tasks.h"
#include "mathutil.h"
#include "util.h"

using std::cout;
using std::endl;
using std::sort;
using std::ifstream;
using std::ios;
using std::set;
using std::thread;
using std::mutex;

mutex mtx;

void sketch_matrix(int nrows, int ncols, int ncols2, const double *X, const double *X2, int *selected_variables, int nsel_max, int& nsel_real, int thr_id) {
   int i, j, jj, ntot, pcount, nsel, imax;
   long long il;
   int *varids, *priority_queue;
   double rssmax, rss;
   double *X_selected_t, *X_selected, *u, *s, *vt, *priority_queue_weights, *dp1;
   double *utb;
   double *yy;
   double *b;
   double *sr;
   mtx.lock();
   cout << thr_id << ": sketch_matrix:" << endl;
   cout << thr_id << ":    nrows: " << nrows << endl;
   cout << thr_id << ":    ncols: " << ncols << endl;
   cout << thr_id << ":    ncols2: " << ncols2 << endl;
   mtx.unlock();
   X_selected_t = new double[((long long)nrows)*(nsel_max+ncols2)];
   X_selected = new double[((long long)nrows)*(nsel_max+ncols2)];
   u = new double[((long long)nrows)*(nsel_max+ncols2)];
   s = new double[nsel_max+ncols2];
   vt = new double[(nsel_max+ncols2)*(nsel_max+ncols2)];
   priority_queue_weights = new double[ncols];
   priority_queue = new int[ncols];
   varids = new int[ncols];
   sr = new double[nrows];
   utb = new double[nsel_max+ncols2];
   yy = new double[nsel_max+ncols2];
   b = new double[nsel_max+ncols2];
   for (il=0; il < ((long long)nrows)*ncols2; il++) X_selected_t[il] = X2[il];
   for (i=0; i < nrows; i++) {
      for (j = 0; j < ncols2; j++) {
         X_selected[j*((long long)nrows)+i] = X2[i*((long long)ncols2)+j];
      }
   }
   // Order variables of X according to correlation with X2 variables
   // Apply lazy greedy
   ntot = ncols2;
   nsel = 0;
   svd(X_selected_t, u, s, vt, nrows, ntot, true);
   for (i = 0; i < ncols; i++) priority_queue_weights[i] = 0.0;
   for (i = 0, pcount = 0; i < ncols; i++) {
      priority_queue[pcount++] = i;
      for (j = 0; j < nrows; j++) sr[j] = X[j*((long long)ncols)+i];
      rss = find_prediction_rss(X_selected_t, u, s, vt, nrows, ntot, sr, utb, yy, b);
      priority_queue_weights[i] = rss;
   }
   sort(priority_queue, priority_queue + pcount, [priority_queue_weights](int i, int j) -> bool { return priority_queue_weights[i] > priority_queue_weights[j]; });
   imax = priority_queue[0];
   rssmax = priority_queue_weights[imax];
   pcount--;
   for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
   mtx.lock();
   cout << thr_id << ": Variable(" << nsel + 1 << "): " << imax << " " << rssmax << endl;
   mtx.unlock();
   varids[nsel] = imax;
   dp1 = X_selected + ((long long)ntot)*nrows;
   for (i = 0; i < nrows; i++) dp1[i] = X[i*((long long)ncols)+imax];
   nsel++;
   ntot++;
   for (i = 0; i < ntot; i++) {
      for (j = 0; j < nrows; j++) X_selected_t[j*((long long)ntot)+i] = X_selected[i*((long long)nrows)+j];
   }
   svd(X_selected_t, u, s, vt, nrows, ntot, true);
   for (;;) {
      jj = priority_queue[0];
      for (i = 0; i < nrows; i++) sr[i] = X[i*((long long)ncols)+jj];
      rss = find_prediction_rss(X_selected_t, u, s, vt, nrows, ntot, sr, utb, yy, b);
      priority_queue_weights[jj] = rss;
      sort(priority_queue, priority_queue + pcount, [priority_queue_weights](int i, int j) -> bool { return priority_queue_weights[i] > priority_queue_weights[j]; });
      if (priority_queue[0] != jj) continue;
      imax = priority_queue[0];
      rssmax = priority_queue_weights[imax];
      pcount--;
      for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
      mtx.lock();
      cout << thr_id << ": Variable(" << nsel + 1 << "): " << imax << " " << rssmax << endl;
      mtx.unlock();
      varids[nsel] = imax;
      dp1 = X_selected + ((long long)ntot)*nrows;
      for (i = 0; i < nrows; i++) dp1[i] = X[i*((long long)ncols)+imax];
      nsel++;
      ntot++;
      //sort(varids_sorted, varids_sorted + nsel);
      if (nsel == nsel_max || rssmax/nrows < 1.0e-15) {
         nsel_real = nsel;
         for (i = 0; i < nsel; i++) selected_variables[i] = varids[i];
         break;
      }
      for (i = 0; i < ntot; i++) {
         for (j = 0; j < nrows; j++) X_selected_t[j*((long long)ntot)+i] = X_selected[i*((long long)nrows)+j];
      }
      svd(X_selected_t, u, s, vt, nrows, ntot, true);
   }
   delete [] b;
   delete [] yy;
   delete [] utb;
   delete [] sr;
   delete [] varids;
   delete [] priority_queue;
   delete [] priority_queue_weights;
   delete [] vt;
   delete [] s;
   delete [] u;
   delete [] X_selected;
   delete [] X_selected_t;
}

void sketch_matrices(int N, vector<string>& mat_names, int *ids_selected_total, int ranks[], int n_test_rows, int nthreads, int *selected_variables[], int max_cols, bool select_from_leading_matrix) {
   int i, id, j, k, nr, nc, nmat, rs, ncols, mid, max_rank, md, ms, is, i1, i2;
   double *Xprev, *dp, *Xs, *avgs, *stdevs, *avgstosave, *stdevstosave;
   double **X;
   int *mat_ids, *mat_sizes, *nsels, *svsorted;
   ifstream ifile;
   thread threads[MAX_THREADS];
   cout << "#####   Starting sketch_matrices   #####\n";
   cout << "Number of rows: " << N << endl;
   cout << "Number of rows to use: " << n_test_rows << endl;
   cout << "Number of threads: " << nthreads << endl;
   nmat = mat_names.size();
   cout << "Number of matrices: " << nmat << endl;
   ifile.open("descs.bin", ios::in | ios::binary);
   ifile.read((char*)&nr, sizeof(int));
   ifile.read((char*)&nc, sizeof(int));
   X = new double*[nthreads];
   for (i = 0; i < nthreads; i++) {
      X[i] = new double[((long long)n_test_rows)*max_cols];
   }
   max_rank = 0;
   for (i=0; i < nmat; i++) {
      if (max_rank < ranks[i]) max_rank = ranks[i];
   }
   Xs = new double[n_test_rows*max_rank];
   Xprev = new double[((long long)n_test_rows)*nc];
   avgs = new double[max_cols];
   stdevs = new double[max_cols];
   avgstosave = new double[max_rank];
   stdevstosave = new double[max_rank];
   if (select_from_leading_matrix) {
      for (i = 0, j = 0; i < N; i++) {
         if (j < n_test_rows && ids_selected_total[j] == i) {
            for (k=0; k < nc; k++) ifile.read((char*)&Xprev[((long long)j)*nc+k], sizeof(double));
            j++;
            if (j == n_test_rows) break;
         } else {
            ifile.seekg(nc*sizeof(double), ifile.cur);
         }
      }
      Matrix2File(Xprev, n_test_rows, nc, "descssel.bin");
   } else {
      for (i = 0; i < n_test_rows; i++) {
         for (k=0; k < nc; k++) ifile.read((char*)&Xprev[((long long)i)*nc+k], sizeof(double));
      }
   }
   ifile.close();
   mat_ids = new int[nmat];
   mat_sizes = new int[nmat];
   nsels = new int[nmat];
   svsorted = new int[max_rank];
   rs = 0;
   for (id = 0; id < nmat; id++) {
      ifile.open("mat" + mat_names[id] + ".bin", ios::in | ios::binary);
      ifile.read((char*)&nr, sizeof(int));
      ifile.read((char*)&ncols, sizeof(int));
      dp = X[rs];
      for (i = 0, j = 0; i < N; i++) {
         if (j < n_test_rows && ids_selected_total[j] == i) {
            for (k=0; k < ncols; k++) ifile.read((char*)&dp[((long long)j)*ncols+k], sizeof(double));
            j++;
            if (j == n_test_rows) break;
         } else {
            ifile.seekg(ncols*sizeof(double), ifile.cur);
         }
      }
      mat_ids[rs] = id;
      mat_sizes[rs] = ncols;
      rs++;
      ifile.close();
      if (rs == nthreads || id == nmat-1) {
         for (i = 0; i < rs; i++) {
            threads[i] = thread(sketch_matrix, n_test_rows, mat_sizes[i], nc, X[i], Xprev, selected_variables[mat_ids[i]], ranks[mat_ids[i]], std::ref(nsels[i]), i);
         }
         for (i=0; i < rs; i++) threads[i].join();
         for (mid=0; mid < rs; mid++) {
            ranks[mat_ids[mid]] = nsels[mid];
            md = nsels[mid];
            std::copy(selected_variables[mat_ids[mid]], selected_variables[mat_ids[mid]]+md, svsorted);
            sort(svsorted, svsorted+md);
            // X[i] -> Xs
            dp = X[mid];
            ms = mat_sizes[mid];
            File2Matrix(("avgs" + mat_names[mat_ids[mid]] + ".bin").c_str(), i1, i2, avgs);
            File2Matrix(("stdevs" + mat_names[mat_ids[mid]] + ".bin").c_str(), i1, i2, stdevs);
            for (i=0, is=0; i < ms; i++) {
               if (is < md && svsorted[is] == i) {
                  for (j = 0; j < n_test_rows; j++) {
                     Xs[j*md+is] = dp[j*ms+i];
                  }
                  avgstosave[is] = avgs[i];
                  stdevstosave[is] = stdevs[i];
                  is++;
                  if (is == md) break;
               }
            }
            Matrix2File(Xs, n_test_rows, md, ("mat" + mat_names[mat_ids[mid]] + "sk.bin").c_str());
            Matrix2File(avgstosave, md, 1, ("avgs" + mat_names[mat_ids[mid]] + "sk.bin").c_str());
            Matrix2File(stdevstosave, md, 1, ("stdevs" + mat_names[mat_ids[mid]] + "sk.bin").c_str());
         }
         rs = 0;
      }
   }
   cout << "#####   Finished sketch_matrices   #####\n";
   delete [] svsorted;
   delete [] nsels;
   delete [] mat_sizes;
   delete [] mat_ids;
   delete [] stdevstosave;
   delete [] avgstosave;
   delete [] stdevs;
   delete [] avgs;
   delete [] Xprev;
   delete [] Xs;
   for (i = 0; i < nthreads; i++) {
      delete [] X[i];
   }
   delete [] X;
}
