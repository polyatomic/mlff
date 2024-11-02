#include <iostream>
#include <algorithm>
#include <thread>
#include <string>
#include <fstream>
#include <set>
#include <map>

#include "util.h"
#include "tasks.h"
#include "mathutil.h"
#include "ThreadPool.h"

using std::cout;
using std::cerr;
using std::endl;
using std::sort;
using std::set_difference;
using std::set_union;
using std::lower_bound;
using std::thread;
using std::string;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::set;
using std::map;

void find_slot_and_pos(int n_slots, int slot_sums[], int id, int& slot, int& pos) {
   int *pt;
   pt = lower_bound(slot_sums, slot_sums+n_slots, id+1);
   slot = pt - slot_sums;
   if (slot > 0)
      pos = id - slot_sums[slot-1];
   else
      pos = id;
}

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

// wmatct - matrix with cur_rows columns (cases) and M rows (variables)
// Calculates dot products for subset of Ms variables from istart until istart+Ms 
//void runKCalc(double *wmatct, double *Kd, int cur_rows, int istart, int Ms) {
//   int i, j, istarti, istartj, iM;
//   istarti = istart;
//   for (i = 0; i < Ms; i++) {
//      iM = i*Ms;
//      istartj = istarti;
//      for (j = i; j < Ms; j++) {
//         Kd[iM+j] += dot_prod(wmatct+istarti*cur_rows, wmatct+istartj*cur_rows, cur_rows);
//         istartj++;
//      }
//      istarti++;
//   }
//}

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

void runLeverageCalc(int jstart, int jend, int rowsum, int rank, int Ms, int ntot, int start, double *levmat, double *wmatc, double *wwp, double *leverages) {
   int i, j;
   for (i = rowsum + jstart, j = jstart; j < jend; i++, j++) {
      MatByVec(levmat, wmatc + j*ntot + start, wwp, rank, Ms);
      leverages[i] = dot_prod(wwp, wwp, rank);
   }
}

bool calculate_descriptors(int N, int M, int olength, double *r, FunctorDaDaI *dcalc, double *stdevs, double *avgs, /*double* wmatc,*/ int nthreads, int maxrows, const char *stdevs_file, const char *avgs_file) {
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
   fn = "descs.bin";
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
   fn = "descst.bin";
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

bool order_by_leverages(int N, int M, double *wmatc, int *ids_selected, double *leverages) {
   int i, j, iM, rank, ii;
   double *wmatct, *Kd, *evec, *eval, *sevec, *seval, *levmat, *dp1, *dp2, *wwp;
   double tol, tmp, ss;
   cout << "#####   Starting order_by_leverages   #####\n";
   cout << "Number of rows: " << N << endl;
   cout << "Number of basis functions: " << M << endl;
   wmatct = new double[((long long)M)*N];
   Kd = new double[M*M];
   evec = new double[M*M];
   eval = new double[M];
   sevec = new double[M*M];
   seval = new double[M];
   levmat = new double[M*M];
   wwp = new double[M];
   transpose(wmatc, wmatct, N, M);
   for (i = 0; i < M*M; i++) Kd[i] = 0.0;
   for (i = 0; i < M; i++) {
      iM = i*M;
      for (j = i; j < M; j++) {
         Kd[iM+j] += dot_prod(wmatct + ((long long)i)*N, wmatct + ((long long)j)*N, N);
      }
   }
   for (i = 1; i < M; i++) {
      iM = i*M;
      for (j = 0; j < i; j++) Kd[iM+j] = Kd[j*M+i];
   }
   cout << "Starting eigendecomposition" << endl;
   cout << "eigs returned " << eigs(Kd, evec, eval, M) << endl;
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
   for (ii=0; ii < M; ii++) {
      i = M - ii - 1;
      for (j=0; j < M; j++) sevec[ii*M+j] = evec[i*M+j];
      seval[ii] = eval[i];
   }
   for (i=0; i < M*M; i++) levmat[i] = 0.0;
   for (i = 0; i < rank; i++) {
      ss = sqrt(seval[i]);
      dp1 = levmat + i*M;
      dp2 = sevec + i*M;
      for (j=0; j < M; j++) dp1[j] = dp2[j]/ss;
   }
   for (i = 0; i < N; i++) {
      MatByVec(levmat, wmatc+((long long)i)*M, wwp, rank, M);
      leverages[i] = dot_prod(wwp, wwp, rank);
   }
   auto rule = [leverages](int i, int j) -> bool { return leverages[i] > leverages[j]; };
   for (i=0; i < N; i++) ids_selected[i] = i;
   sort(ids_selected, ids_selected+N, rule);
   cout << "#####   Finished order_by_leverages   #####\n";
   delete [] wwp;
   delete [] levmat;
   delete [] seval;
   delete [] sevec;
   delete [] eval;
   delete [] evec;
   delete [] Kd;
   delete [] wmatct;
   return true;
}

// Performs leverage ordering of submatrices of a large matrix saved as file
// N - number of matrix rows
// nsub - number of submatrices
// M[] - array of submatrix widths
// wmatc_name - name of matrix file
// ids_selected - row numbers ordered according to leverage values for each submatrix (out)
// leverages - calculated leverage values for each submatrix (out)
// bufsize - row buffer size
bool order_by_leverages(int N, int M[], int nsub, const char *wmatc_name, const char *wmatct_name, int *ids_selected[], double *leverages[], int ranks[], int bufsize, int nthreads) {
   int i, j, k, ntot, Ms, rowsum, cur_rows, nr, nc, iM, Ms_max, rank, ii, start, jstart, jend, bufsize2, ntotsub, npairs,
       id, matid, rslast, crlast, next_rowsum, rowsum2, cur_rows2, nr1, nr2;
   int *tsizes = 0;
   double tol, tmp, ss;
   double *Kds, *dp1, *dp2, *lp;
   double *wmatc = 0;
   double *evec = 0;
   double *eval = 0;
   double *sevec = 0;
   double *seval = 0;
   double *levmat = 0;
   double **wwp = 0;
   double **wmatct = 0;
   double **Kd = 0;
   bool next_init = false;
   bool res = true;
   set<int> ids;
   map<int, int> matm, matmr;
   vector<int> ids1, ids2;
   thread threads[MAX_THREADS];
   cout << "#####   Starting order_by_leverages   #####\n";
   if (nthreads > 1) cout << "Using " << nthreads << " threads" << endl;
   cout << "Buffer size: " << bufsize << endl;
   cout << "Number of rows: " << N << endl;
   cout << "Number of submatrices: " << nsub << endl;
   ifstream ifilet(wmatct_name, ios::in | ios::binary);
   if (!ifilet) return false;
   Kd = new double*[nsub];
   wmatct = new double*[nthreads+1];
   ntot = 0;
   Ms_max = 0;
   Kds = 0;
   for (i=0; i < nsub; i++) {
      Ms = M[i];
      if (Ms > Ms_max) Ms_max = Ms;
      if (Ms) Kds = Kd[i] = new double[Ms*Ms];
      else Kd[i] = 0;
      for (j=0; j < Ms*Ms; j++) Kds[j] = 0.0;
      ntot += Ms;
   }
   cout << "Total number of basis functions: " << ntot << endl;
   bufsize2 = nthreads*bufsize;
   wmatc = new double[ntot*bufsize2];
   for (i=0; i <= nthreads; i++) wmatct[i] = new double[N*bufsize];
   for (k=0; k < nsub; k++) {
      for (i=0,start=0; i < k; i++) start += M[i];
      ntotsub = M[k];
      if (!ntotsub) continue;
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
         ifilet.seekg(((long long)(rowsum+start))*N*sizeof(double), ifilet.cur);
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
                  threads[i] = thread(runKCalc, wmatct[ids1[i]], wmatct[ids2[i]], Kd[k], matm[ids1[i]], matm[ids2[i]], nr1, nr2, N, ntotsub);
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
   }
   ifilet.close();
   for (k = 0; k < nsub; k++) {
      Ms = M[k];
      Kds = Kd[k];
      for (i = 1; i < Ms; i++) {
         iM = i*Ms;
         for (j = 0; j < i; j++) Kds[iM+j] = Kds[j*Ms+i];
      }
   }
   evec = new double[Ms_max*Ms_max];
   eval = new double[Ms_max];
   sevec = new double[Ms_max*Ms_max];
   seval = new double[Ms_max];
   levmat = new double[Ms_max*Ms_max];
   wwp = new double*[nthreads];
   for (i = 0; i < nthreads; i++) wwp[i] = new double[Ms_max];
   tsizes = new int[nthreads];
   ifstream ifile(wmatc_name, ios::in | ios::binary);
   if (!ifile) {
      res = false;
      goto end;
   }
   start = 0;
   for (k=0; k < nsub; k++) ranks[k] = 0;
   for (k=0; k < nsub; k++) {
      Ms = M[k];
      if (!Ms) continue;
      Kds = Kd[k];
      cout << "Starting eigendecomposition(" << k << ")" << endl;
      cout << "eigs returned " << eigs(Kds, evec, eval, Ms) << endl;
      tol = 1.0e-15*Ms*sqrt(eval[Ms-1]);
      rank = 0;
      for (i = 0; i < Ms; i++) {
         if (eval[i] > 0.0) {
            tmp = sqrt(eval[i]);
            if (tmp > tol) rank++;
         }
      }
      cout << "rank = " << rank << endl;
      ranks[k] = rank;
      cout << "Finished eigendecomposition" << endl;
      for (ii = 0; ii < Ms; ii++) {
         i = Ms-ii-1;
         for (j = 0; j < Ms; j++) sevec[ii*Ms+j] = evec[i*Ms+j];
         seval[ii] = eval[i];
      }
      for (i = 0; i < Ms*Ms; i++) levmat[i] = 0.0;
      for (i = 0; i < rank; i++) {
         ss = sqrt(seval[i]);
         dp1 = levmat + i*Ms;
         dp2 = sevec + i*Ms;
         for (j = 0; j < Ms; j++) dp1[j] = dp2[j]/ss;
      }
      ifile.seekg(0, ifile.beg);
      ifile.read((char*)&nr, sizeof(int));
      ifile.read((char*)&nc, sizeof(int));
      rowsum = 0;
      for (;;) {
         if (N - rowsum < bufsize2)
            cur_rows = N - rowsum;
         else
            cur_rows = bufsize2;
         divide(tsizes, cur_rows, nthreads);
         for (i = 0; i < cur_rows; i++) {
            for (j = 0; j < ntot; j++) ifile.read((char*)&wmatc[i*ntot+j], sizeof(double));
         }
         jstart = 0;
         jend = tsizes[0];
         for (i = 0; i < nthreads; i++) {
            threads[i] = thread(runLeverageCalc, jstart, jend, rowsum, rank, Ms, ntot, start, levmat, wmatc, wwp[i], leverages[k]);
            jstart = jend;
            if (i < nthreads-1) jend = jstart + tsizes[i+1];
         }
         for (i=0; i < nthreads; i++) threads[i].join();
//         for (i = rowsum, j=0; j < cur_rows; i++, j++) {
//            MatByVec(levmat, wmatc + j*ntot + start, wwp[0], rank, Ms);
//            leverages[k][i] = dot_prod(wwp[0], wwp[0], rank);
//         }
         rowsum += cur_rows;
         if (rowsum == N) break;
      }
      lp = leverages[k];
      auto rule = [lp](int i, int j) -> bool { return lp[i] > lp[j]; };
      for (i = 0; i < N; i++) ids_selected[k][i] = i;
      sort(ids_selected[k], ids_selected[k]+N, rule);
      start += Ms;
   }
   cout << "#####   Finished order_by_leverages   #####\n";
   ifile.close();
end:
   delete [] tsizes;
   for (i=0; i < nthreads; i++) delete [] wwp[i];
   delete [] wwp;
   delete [] levmat;
   delete [] seval;
   delete [] sevec;
   delete [] eval;
   delete [] evec;
   for (i=0; i <= nthreads; i++) delete [] wmatct[i];
   delete [] wmatc;
   for (i=0; i < nsub; i++) delete [] Kd[i];
   delete [] wmatct;
   delete [] Kd;
   return res;
}

// Select columns of nrows x ncols matrix X that are least correlated with each other and with the columns of nrows x ncols2 matrix X2
// nrows - number of rows of matrices X and X2
// ncols - number of columns of matrix X
// ncols2 - number of columns of matrix X2
// selected_variables - selected columns of matrix X
// nsel_max - maximum number of columns to select
// nsel_real - number of columns actually selected
// nthreads - number of threads to use
void sketch_matrix(int nrows, int ncols, int ncols2, const double *X, const double *X2, int *selected_variables, int nsel_max, int& nsel_real, int nthreads) {
   int i, j, jj, k, kk, ntot, pcount, nsel, imax;
   int js[MAX_THREADS];
   long long il;
   int *varids_sorted, *priority_queue;
   double rssmax;
   double rss[MAX_THREADS];
   double *X_selected_t, *X_selected, *u, *s, *vt, *priority_queue_weights, *dp1;
   double *utb[MAX_THREADS] = {};
   double *yy[MAX_THREADS] = {};
   double *b[MAX_THREADS] = {};
   double *sr[MAX_THREADS] = {};
   cout << "sketch_matrix:" << endl;
   cout << "   nrows: " << nrows << endl;
   cout << "   ncols: " << ncols << endl;
   cout << "   ncols2: " << ncols2 << endl;
   cout << "   nthreads: " << nthreads << endl;
   X_selected_t = new double[((long long)nrows)*(nsel_max+ncols2)];
   X_selected = new double[((long long)nrows)*(nsel_max+ncols2)];
   u = new double[((long long)nrows)*(nsel_max+ncols2)];
   s = new double[nsel_max+ncols2];
   vt = new double[(nsel_max+ncols2)*(nsel_max+ncols2)];
   priority_queue_weights = new double[ncols];
   priority_queue = new int[ncols];
   auto rule2 = [priority_queue_weights](int i, int j) -> bool { return priority_queue_weights[i] > priority_queue_weights[j]; };
   varids_sorted = new int[ncols];
   for (i=0; i < MAX_THREADS; i++) {
      sr[i] = new double[nrows];
      utb[i] = new double[nsel_max+ncols2];
      yy[i] = new double[nsel_max+ncols2];
      b[i] = new double[nsel_max+ncols2];
   }
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
      for (j = 0; j < nrows; j++) sr[0][j] = X[j*((long long)ncols)+i];
      rss[0] = find_prediction_rss(X_selected_t, u, s, vt, nrows, ntot, sr[0], utb[0], yy[0], b[0]);
      priority_queue_weights[i] = rss[0];
   }
   sort(priority_queue, priority_queue + pcount, rule2);
   imax = priority_queue[0];
   rssmax = priority_queue_weights[imax];
   pcount--;
   for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
   cout << "Variable(" << nsel + 1 << "): " << imax << " " << rssmax << endl;
   varids_sorted[nsel] = imax;
   dp1 = X_selected + ((long long)ntot)*nrows;
   for (i = 0; i < nrows; i++) dp1[i] = X[i*((long long)ncols)+imax];
   nsel++;
   ntot++;
   for (i = 0; i < ntot; i++) {
      for (j = 0; j < nrows; j++) X_selected_t[j*((long long)ntot)+i] = X_selected[i*((long long)nrows)+j];
   }
   svd(X_selected_t, u, s, vt, nrows, ntot, true);
   ThreadPool tp(nthreads);
   tp.Start();
   for (;;) {
      jj = priority_queue[0];
      for (k=0,kk=0;;) {
         j = priority_queue[k];
         k++;
         for (i = 0; i < nrows; i++) sr[kk][i] = X[i*((long long)ncols)+j];
         auto job = [&](int jid){ rss[jid] = find_prediction_rss(X_selected_t, u, s, vt, nrows, ntot, sr[jid], utb[jid], yy[jid], b[jid]); };
         js[kk] = j;
         tp.QueueJob(job, kk);
         // job(kk);
         kk++;
         if (kk == nthreads || k == pcount) break;
      }
      while (tp.busy());
      for (i = 0; i < kk; i++) {
         priority_queue_weights[js[i]] = rss[i];
      }
      sort(priority_queue, priority_queue + pcount, rule2);
      if (priority_queue[0] != jj) continue;
      imax = priority_queue[0];
      rssmax = priority_queue_weights[imax];
      pcount--;
      for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
      cout << "Variable(" << nsel + 1 << "): " << imax << " " << rssmax << endl;
      varids_sorted[nsel] = imax;
      dp1 = X_selected + ((long long)ntot)*nrows;
      for (i = 0; i < nrows; i++) dp1[i] = X[i*((long long)ncols)+imax];
      nsel++;
      ntot++;
      sort(varids_sorted, varids_sorted + nsel);
      if (nsel == nsel_max || rssmax/nrows < 1.0e-15) {
         nsel_real = nsel;
         for (i = 0; i < nsel; i++) selected_variables[i] = varids_sorted[i];
         break;
      }
      for (i = 0; i < ntot; i++) {
         for (j = 0; j < nrows; j++) X_selected_t[j*((long long)ntot)+i] = X_selected[i*((long long)nrows)+j];
      }
      svd(X_selected_t, u, s, vt, nrows, ntot, true);
   }
   tp.Stop();
   for (i=0; i < MAX_THREADS; i++) {
      delete [] b[i];
      delete [] yy[i];
      delete [] utb[i];
      delete [] sr[i];
   }
   delete [] varids_sorted;
   delete [] priority_queue;
   delete [] priority_queue_weights;
   delete [] vt;
   delete [] s;
   delete [] u;
   delete [] X_selected;
   delete [] X_selected_t;
}

void sketch_matrix2(int nrows, int ncols, int ncols2, const double* X, const double* X2, int* selected_variables, int nsel_max, int& nsel_real, int thr_id) {
   int i, j, jj, ntot, pcount, nsel, imax, imin;
   long long il;
   int *varids_sorted, *priority_queue;
   double rssmax, rssmin, rss;
   double *X_selected_t, *X_selected, *u, *s, *vt, *priority_queue_weights, *dp1;
   double *utb;
   double *yy;
   double *b;
   double *sr;
   cout << thr_id << ": sketch_matrix:" << endl;
   cout << thr_id << ":    nrows: " << nrows << endl;
   cout << thr_id << ":    ncols: " << ncols << endl;
   cout << thr_id << ":    ncols2: " << ncols2 << endl;
   X_selected_t = new double[((long long)nrows)*(nsel_max+ncols2)];
   X_selected = new double[((long long)nrows)*(nsel_max+ncols2)];
   u = new double[((long long)nrows)*(nsel_max+ncols2)];
   s = new double[nsel_max+ncols2];
   vt = new double[(nsel_max+ncols2)*(nsel_max+ncols2)];
   priority_queue_weights = new double[ncols];
   priority_queue = new int[ncols];
   auto rule2 = [priority_queue_weights](int i, int j) -> bool { return priority_queue_weights[i] > priority_queue_weights[j]; };
   varids_sorted = new int[ncols];
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
   sort(priority_queue, priority_queue + pcount, rule2);
   imax = priority_queue[0];
   rssmax = priority_queue_weights[imax];
   pcount--;
   for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
   cout << thr_id << ": Variable(" << nsel + 1 << "): " << imax << " " << rssmax << endl;
   varids_sorted[nsel] = imax;
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
      sort(priority_queue, priority_queue + pcount, rule2);
      if (priority_queue[0] != jj) continue;
      imax = priority_queue[0];
      imin = priority_queue[pcount-1];
      rssmax = priority_queue_weights[imax];
      rssmin = priority_queue_weights[imin];
      pcount--;
      for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
      cout << thr_id << ": Variable(" << nsel + 1 << "): " << imax << " " << rssmax << " " << rssmin << endl;
      varids_sorted[nsel] = imax;
      dp1 = X_selected + ((long long)ntot)*nrows;
      for (i = 0; i < nrows; i++) dp1[i] = X[i*((long long)ncols)+imax];
      nsel++;
      ntot++;
      sort(varids_sorted, varids_sorted + nsel);
      if (nsel == nsel_max || rssmax/nrows < 1.0e-15) {
         nsel_real = nsel;
         for (i = 0; i < nsel; i++) selected_variables[i] = varids_sorted[i];
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
   delete [] varids_sorted;
   delete [] priority_queue;
   delete [] priority_queue_weights;
   delete [] vt;
   delete [] s;
   delete [] u;
   delete [] X_selected;
   delete [] X_selected_t;
}

void sketch_matrix(int nrows, int ncols, int ncols_retain, const double *X, const double *y, int *selected_variables) {
   double betamax, beta, betamin, rss, rssmax;
   double *dvec, *X_selected, *X_selected_t, *dp1, *u, *vt, *s, *priority_queue_weights, *sr, *utb, *yy, *b;
   int i, j, imax, imin, nsel, scount, pcount;
   int *varids_sorted, *priority_queue;
   // Let's find the variable most correlated with target
   betamax = 0.0;
   imax = -1;
   dvec = new double[nrows];
   X_selected = new double[ncols_retain*nrows];
   X_selected_t = new double[ncols_retain*nrows];
   varids_sorted = new int[ncols_retain];
   u = new double[ncols_retain*nrows];
   vt = new double[ncols_retain*ncols_retain];
   s = new double[ncols_retain];
   priority_queue_weights = new double[ncols_retain];
   priority_queue = new int[ncols_retain];
   auto rule2 = [priority_queue_weights](int i, int j) -> bool { return priority_queue_weights[i] > priority_queue_weights[j]; };
   sr = new double[nrows];
   utb = new double[ncols_retain];
   yy = new double[ncols_retain];
   b = new double[ncols_retain];
   for (i = 0; i < ncols; i++) {
      for (j = 0; j < nrows; j++) dvec[j] = X[j*ncols+i];
      beta = fabs(dot_prod(dvec, y, nrows)/dot_prod(dvec, dvec, nrows));
      if (beta > betamax) {
         betamax = beta;
         imax = i;
      }
   }
   cout << "First variable: " << imax << endl;
   for (i = 0; i < nrows; i++) X_selected[i] = X[i*ncols+imax];
   // Let's find another variable that is least correlated with the first variable
   betamin = 1.0e10;
   imin = -1;
   for (i = 0; i < ncols; i++) {
      if (i == imax) continue;
      for (j = 0; j < nrows; j++) dvec[j] = X[j*ncols+i];
      beta = fabs(dot_prod(dvec, X_selected, nrows)/dot_prod(dvec, dvec, nrows));
      if (beta < betamin) {
         betamin = beta;
         imin = i;
      }
   }
   cout << "Second variable: " << imin << endl;
   dp1 = X_selected + nrows;
   for (i = 0; i < nrows; i++) dp1[i] = X[i*ncols+imin];
   nsel = 2;
   for (i = 0; i < nsel; i++) {
      for (j = 0; j < nrows; j++) X_selected_t[j*nsel+i] = X_selected[i*nrows+j];
   }
   varids_sorted[0] = imax;
   varids_sorted[1] = imin;
   sort(varids_sorted, varids_sorted + 2);
   svd(X_selected_t, u, s, vt, nrows, nsel, true);
   for (i = 0; i < ncols; i++) priority_queue_weights[i] = 0.0;
   for (i = 0, scount = 0, pcount = 0; i < ncols; i++) {
      if (scount < nsel && i == varids_sorted[scount]) {
         scount++;
         continue;
      }
      priority_queue[pcount++] = i;
      for (j = 0; j < nrows; j++) sr[j] = X[j*ncols+i];
      rss = find_prediction_rss(X_selected_t, u, s, vt, nrows, nsel, sr, utb, yy, b);
      priority_queue_weights[i] = rss;
   }
   sort(priority_queue, priority_queue + pcount, rule2);
   imax = priority_queue[0];
   rssmax = priority_queue_weights[imax];
   pcount--;
   for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
   cout << "Variable(" << nsel+1 << "): " << imax << " " << rssmax << endl;
   varids_sorted[nsel] = imax;
   dp1 = X_selected + nsel*nrows;
   for (i = 0; i < nrows; i++) dp1[i] = X[i*ncols+imax];
   nsel++;
   for (i = 0; i < nsel; i++) {
      for (j = 0; j < nrows; j++) X_selected_t[j*nsel+i] = X_selected[i*nrows+j];
   }
   svd(X_selected_t, u, s, vt, nrows, nsel, true);
   sort(varids_sorted, varids_sorted + nsel);
   for (;;) {
      // Use lazy greedy algorithm
      j = priority_queue[0];
      for (i = 0; i < nrows; i++) sr[i] = X[i*ncols+j];
      rss = find_prediction_rss(X_selected_t, u, s, vt, nrows, nsel, sr, utb, yy, b);
      priority_queue_weights[j] = rss;
      sort(priority_queue, priority_queue + pcount, rule2);
      if (priority_queue[0] != j) continue;
      imax = priority_queue[0];
      rssmax = priority_queue_weights[imax];
      pcount--;
      for (i = 0; i < pcount; i++) priority_queue[i] = priority_queue[i+1];
      cout << "Variable(" << nsel+1 << "): " << imax << " " << rssmax << endl;
      varids_sorted[nsel] = imax;
      dp1 = X_selected + nsel*nrows;
      for (i = 0; i < nrows; i++) dp1[i] = X[i*ncols+imax];
      nsel++;
      for (i = 0; i < nsel; i++) {
         for (j = 0; j < nrows; j++) X_selected_t[j*nsel+i] = X_selected[i*nrows+j];
      }
      svd(X_selected_t, u, s, vt, nrows, nsel, true);
      sort(varids_sorted, varids_sorted + nsel);
      if (nsel == ncols_retain) {
         for (i = 0; i < nsel; i++) selected_variables[i] = varids_sorted[i];
         break;
      }
   }
   delete [] b;
   delete [] yy;
   delete [] utb;
   delete [] sr;
   delete [] priority_queue;
   delete [] priority_queue_weights;
   delete [] s;
   delete [] vt;
   delete [] u;
   delete [] varids_sorted;
   delete [] X_selected_t;
   delete [] X_selected;
   delete [] dvec;
}

// N - number of rows
// M - number of columns for each block
// nsub - number of blocks
// wmatc_name - name of the matrix file
// ids_selected - leverage ordered ids for each block
// ranks - rank for each block
// n_test_rows - number of rows to select from each block
// nthreads - number of threads to use
// selected_variables - selected column ids for each block
void sketch_matrices(int N, int M[], int nblocks[], int nsub, const char *wmatc_name, int *ids_selected[], int ranks[], int n_test_rows, int nthreads, int *selected_variables[]) {
   int i, j, k, Mtot, ranktot, maxcols, mat_id, nr, nc, r_skip1, nrows, n_selected_total, r_skip2, ncols, ncolsprev, ncolsnew, ncolsprevupd, nsel;
   long long il;
   int *ids_selected_total_cp, *idit, *ids_selectedo, *ids_selected_total;
   double *X, *Xprev, *Xtmp, *Xprevupd;
   cout << "#####   Starting sketch_matrices   #####\n";
   cout << "Number of rows: " << N << endl;
   cout << "Number of rows to use: " << n_test_rows << endl;
   cout << "Numbers of columns and corresponding ranks" << endl;
   for (i=0, Mtot=0, ranktot=0, maxcols=0; i < nsub; i++) {
      if (M[i] > maxcols) maxcols = M[i];
      cout << M[i] << " " << ranks[i] << endl;
      Mtot += M[i];
      ranktot += ranks[i];
   }
   ifstream ifile(wmatc_name, ios::in | ios::binary);
   if (!ifile) return;
   ids_selectedo = new int[n_test_rows];
   ids_selected_total = new int[nsub*n_test_rows];
   ids_selected_total_cp = new int[nsub*n_test_rows];
   ifile.read((char*)&nr, sizeof(int));
   ifile.read((char*)&nc, sizeof(int));
   nrows = n_test_rows;
   for (i = 0; i < nrows; i++) ids_selected_total_cp[i] = ids_selected[0][i];
   sort(ids_selected_total_cp, ids_selected_total_cp + nrows);
   n_selected_total = nrows;
   for (mat_id = 1; mat_id < nsub; mat_id++) {
      for (i = 0; i < nrows; i++) ids_selectedo[i] = ids_selected[mat_id][i];
      sort(ids_selectedo, ids_selectedo + nrows);
      idit = set_union(ids_selectedo, ids_selectedo + nrows, ids_selected_total_cp, ids_selected_total_cp + n_selected_total, ids_selected_total);
      n_selected_total = idit - ids_selected_total;
      for (i=0; i < n_selected_total; i++) ids_selected_total_cp[i] = ids_selected_total[i];
   }
   nrows = n_selected_total;
   cout << "Total number of rows: " << nrows << endl;
   X = new double[((long long)nrows)*maxcols];
   Xprev = new double[((long long)nrows)*ranktot];
   Xprevupd = new double[((long long)nrows)*ranktot];
   Xtmp = new double[((long long)nrows)*ranktot];
   r_skip1 = 0;
   for (mat_id = 0; mat_id < nsub; mat_id++) {
      ncols = M[mat_id];
      if (mat_id > 0) r_skip1 += M[mat_id-1];
      for (i = mat_id + 1, r_skip2 = 0; i < nsub; i++) r_skip2 += M[i];
      for (i = 0, j = 0; i < N; i++) {
         if (j < nrows && ids_selected_total[j] == i) {
            ifile.seekg(r_skip1*sizeof(double), ifile.cur);
            for (k=0; k < ncols; k++) ifile.read((char*)&X[((long long)j)*ncols+k], sizeof(double));
            j++;
            if (j == nrows) break;
            ifile.seekg(r_skip2*sizeof(double), ifile.cur);
         } else {
            ifile.seekg(Mtot*sizeof(double), ifile.cur);
         }
      }
      ifile.seekg(0, ifile.beg);
      ifile.read((char*)&nr, sizeof(int));
      ifile.read((char*)&nc, sizeof(int));
      if (mat_id == 0) {
         for (il=0; il < ((long long)ncols)*nrows; il++) Xprevupd[il] = X[il];
         ncolsprevupd = ncols;
      } else if (mat_id < nblocks[0]) {
         ncolsnew = ncolsprevupd + ncols;
         for (il = 0; il < ((long long)nrows)*ncolsprevupd; il++) Xtmp[il] = Xprevupd[il];
         for (i = 0; i < nrows; i++) {
            for (j = 0; j < ncolsprevupd; j++) {
               Xprevupd[i*((long long)ncolsnew)+j] = Xtmp[i*((long long)ncolsprevupd)+j];
            }
            for (j = 0; j < ncols; j++) {
               Xprevupd[i*((long long)ncolsnew)+ncolsprevupd+j] = X[((long long)i)*ncols+j];
            }
         }
         if (mat_id == nblocks[0]-1) {
            for (il=0; il < ((long long)ncolsnew)*nrows; il++) Xprev[il] = Xprevupd[il];
            ncolsprev = ncolsnew;
         }
         ncolsprevupd = ncolsnew;
      } else {
         // In case there is a need to save memory:
         // save Xprevupd and Xtmp to file
         // delete [] Xprevupd
         // delete [] Xtmp
         sketch_matrix(nrows, ncols, ncolsprev, X, Xprev, selected_variables[mat_id], ranks[mat_id], nsel, nthreads);
         // Xprevupd = new double[((long long)nrows)*ranktot];
         // Xtmp = new double[((long long)nrows)*ranktot];
         // read Xprevupd and Xtmp from file
         ranks[mat_id] = nsel;
         // Add X selected cols to Xprevupd
         ncolsnew = ncolsprevupd + ranks[mat_id];
         for (il = 0; il < ((long long)nrows)*ncolsprevupd; il++) Xtmp[il] = Xprevupd[il];
         for (i = 0; i < nrows; i++) {
            for (j = 0; j < ncolsprevupd; j++) {
               Xprevupd[i*((long long)ncolsnew)+j] = Xtmp[i*((long long)ncolsprevupd)+j];
            }
            for (j = 0; j < ranks[mat_id]; j++) {
               Xprevupd[i*((long long)ncolsnew)+ncolsprevupd+j] = X[i*((long long)ncols)+selected_variables[mat_id][j]];
            }
         }
         if (mat_id == nblocks[0]+nblocks[1]-1) {
            for (il = 0; il < ncolsnew*((long long)nrows); il++) Xprev[il] = Xprevupd[il];
            ncolsprev = ncolsnew;
         }
         ncolsprevupd = ncolsnew;
      }
   }
   cout << "#####   Finished sketch_matrices   #####\n";
   ifile.close();
   delete [] Xtmp;
   delete [] Xprevupd;
   delete [] Xprev;
   delete [] X;
   delete [] ids_selected_total_cp;
   delete [] ids_selected_total;
   delete [] ids_selectedo;
}

void sketch_matrices2(int N, int M[], int nblocks[], int nsub, const char *wmatc_name, int *ids_selected[], int ranks[], int n_test_rows, int nthreads, int *selected_variables[]) {
   int i, j, k, nr, nc, nrows, n_selected_total, mat_id, ncols, r_skip1, r_skip2, Mtot, rs, ncolsprevupd, ranktot, ncolsnew, ncolsprev, id, mid;
   long long il;
   int *ids_selectedo, *ids_selected_total, *ids_selected_total_cp, *idit, *sm_order, *mat_ids, *nsels;
   double *Xprevupd, *Xtmp, *Xprev, *dp;
   double **X;
   thread threads[MAX_THREADS];
   cout << "#####   Starting sketch_matrices   #####\n";
   cout << "Number of rows: " << N << endl;
   cout << "Number of rows to use: " << n_test_rows << endl;
   cout << "Number of threads: " << nthreads << endl;
   cout << "Numbers of columns and corresponding ranks" << endl;
   sm_order = new int[nblocks[2]];
   for (i=0, Mtot=0, ranktot=0; i < nsub; i++) {
      if (i >= nblocks[0] + nblocks[1]) sm_order[i-nblocks[0]-nblocks[1]] = i;
      else ranktot += ranks[i];
      cout << M[i] << " " << ranks[i] << endl;
      Mtot += M[i];
   }
   sort(sm_order, sm_order+nblocks[2], [M](int i, int j) -> bool { return M[i] > M[j]; });
   ifstream ifile(wmatc_name, ios::in | ios::binary);
   if (!ifile) {
      delete [] sm_order;
      return;
   }
   ids_selectedo = new int[n_test_rows];
   ids_selected_total = new int[nsub*n_test_rows];
   ids_selected_total_cp = new int[nsub*n_test_rows];
   ifile.read((char*)&nr, sizeof(int));
   ifile.read((char*)&nc, sizeof(int));
   nrows = n_test_rows;
   for (i = 0; i < nrows; i++) ids_selected_total_cp[i] = ids_selected[0][i];
   sort(ids_selected_total_cp, ids_selected_total_cp + nrows);
   n_selected_total = nrows;
   for (mat_id = 1; mat_id < nsub; mat_id++) {
      for (i = 0; i < nrows; i++) ids_selectedo[i] = ids_selected[mat_id][i];
      sort(ids_selectedo, ids_selectedo + nrows);
      idit = set_union(ids_selectedo, ids_selectedo + nrows, ids_selected_total_cp, ids_selected_total_cp + n_selected_total, ids_selected_total);
      n_selected_total = idit - ids_selected_total;
      for (i=0; i < n_selected_total; i++) ids_selected_total_cp[i] = ids_selected_total[i];
   }
   nrows = n_selected_total;
   cout << "Total number of rows: " << nrows << endl;
   X = new double*[nthreads];
   for (i = 0; i < nthreads; i++) {
      X[i] = new double[((long long)nrows)*M[sm_order[i]]];
   }
   Xprev = new double[((long long)nrows)*ranktot];
   Xprevupd = new double[((long long)nrows)*ranktot];
   Xtmp = new double[((long long)nrows)*ranktot];
   mat_ids = new int[nthreads];
   nsels = new int[nthreads];
   rs = 0;
   for (id = 0; id < nsub; id++) {
      mat_id = id;
      if (id >= nblocks[0]+nblocks[1]) mat_id = sm_order[id-nblocks[0]-nblocks[1]];
      ncols = M[mat_id];
      for (i=0, r_skip1=0; i < mat_id; i++) r_skip1 += M[i];
      for (i = mat_id + 1, r_skip2 = 0; i < nsub; i++) r_skip2 += M[i];
      dp = X[rs];
      for (i = 0, j = 0; i < N; i++) {
         if (j < nrows && ids_selected_total[j] == i) {
            ifile.seekg(r_skip1*sizeof(double), ifile.cur);
            for (k=0; k < ncols; k++) ifile.read((char*)&dp[((long long)j)*ncols+k], sizeof(double));
            j++;
            if (j == nrows) break;
            ifile.seekg(r_skip2*sizeof(double), ifile.cur);
         } else {
            ifile.seekg(Mtot*sizeof(double), ifile.cur);
         }
      }
      ifile.seekg(0, ifile.beg);
      ifile.read((char*)&nr, sizeof(int));
      ifile.read((char*)&nc, sizeof(int));
      if (mat_id == 0) {
         for (il=0; il < ((long long)ncols)*nrows; il++) Xprevupd[il] = dp[il];
         ncolsprevupd = ncols;
      } else if (mat_id < nblocks[0]) {
         ncolsnew = ncolsprevupd + ncols;
         for (il = 0; il < ((long long)nrows)*ncolsprevupd; il++) Xtmp[il] = Xprevupd[il];
         for (i = 0; i < nrows; i++) {
            for (j = 0; j < ncolsprevupd; j++) {
               Xprevupd[i*((long long)ncolsnew)+j] = Xtmp[i*((long long)ncolsprevupd)+j];
            }
            for (j = 0; j < ncols; j++) {
               Xprevupd[i*((long long)ncolsnew)+ncolsprevupd+j] = dp[((long long)i)*ncols+j];
            }
         }
         if (mat_id == nblocks[0]-1) {
            for (il=0; il < ((long long)ncolsnew)*nrows; il++) Xprev[il] = Xprevupd[il];
            ncolsprev = ncolsnew;
         }
         ncolsprevupd = ncolsnew;
      } else {
         mat_ids[rs] = mat_id;
         rs++;
         if (rs == nthreads || id == nblocks[0]+nblocks[1]-1 || id == nsub-1) {
            // X[0] until X[rs-1] now contain matrices
            // Start computations here
            for (i = 0; i < rs; i++) {
               threads[i] = thread(sketch_matrix2, nrows, M[mat_ids[i]], ncolsprev, X[i], Xprev, selected_variables[mat_ids[i]], ranks[mat_ids[i]], std::ref(nsels[i]), i);
            }
            for (i=0; i < rs; i++) threads[i].join();
            for (mid=0; mid < rs; mid++) {
               ranks[mat_ids[mid]] = nsels[mid];
               ncolsnew = ncolsprevupd + ranks[mat_ids[mid]];
               for (il = 0; il < ((long long)nrows)*ncolsprevupd; il++) Xtmp[il] = Xprevupd[il];
               dp = X[mid];
               for (i = 0; i < nrows; i++) {
                  for (j = 0; j < ncolsprevupd; j++) {
                     Xprevupd[i*((long long)ncolsnew)+j] = Xtmp[i*((long long)ncolsprevupd)+j];
                  }
                  for (j = 0; j < ranks[mat_ids[mid]]; j++) {
                     Xprevupd[i*((long long)ncolsnew)+ncolsprevupd+j] = dp[i*((long long)ncols)+selected_variables[mat_ids[mid]][j]];
                  }
               }
               if (mid == rs-1 && id == nblocks[0]+nblocks[1]-1) {
                  for (il = 0; il < ncolsnew*((long long)nrows); il++) Xprev[il] = Xprevupd[il];
                  ncolsprev = ncolsnew;
               }
               ncolsprevupd = ncolsnew;
            }
            rs = 0;
         }
      }
   }
   cout << "#####   Finished sketch_matrices   #####\n";
   ifile.close();
   delete [] nsels;
   delete [] mat_ids;
   delete [] Xtmp;
   delete [] Xprevupd;
   delete [] Xprev;
   for (i = 0; i < nthreads; i++) {
      delete [] X[i];
   }
   delete [] X;
   delete [] ids_selected_total_cp;
   delete [] ids_selected_total;
   delete [] ids_selectedo;
   delete [] sm_order;
}

bool calculate_linear_regression(int N, int M, double *wmatc, int *ids, int nt, double *Yc) {
   int i, j, ii, ndim;
   double rss, rss1, rss2, tmp, rs1, rs2, rs;
   double *wmatt, *dp, *dp2, *b, *y;
   cout << "#####   Starting calculate_linear_regression   #####\n";
   cout << "Number of rows: " << N << endl;
   cout << "Number of basis functions: " << M << endl;
   cout << "Training set size: " << nt << endl;
   wmatt = new double[M*nt];
   b = new double[M];
   y = new double[nt];
   for (i = 0; i < nt; i++) {
      ii = ids[i];
      dp = wmatc + ii*M;
      dp2 = wmatt + i*M;
      for (j=0; j < M; j++) dp2[j] = dp[j];
      y[i] = Yc[ii];
   }
   multiple_linear_regression_svd(nt, M, wmatt, y, b, ndim, true);
   for (i = 0,rss1 = 0.0,rs1 = 0.0; i < nt; i++) {
      dp = wmatt + i*M;
      tmp = y[i] - dot_prod(b, dp, M);
      rs1 += fabs(tmp);
      rss1 += tmp*tmp;
   }
   cout << "RMSE(train) = " << sqrt(rss1/nt) << endl;
   cout << "MAE(train) = " << rs1/nt << endl;
   for (i = nt,rss2 = 0.0,rs2 = 0.0; i < N; i++) {
      ii = ids[i];
      dp = wmatc + ii*M;
      tmp = Yc[ii] - dot_prod(b, dp, M);
      rs2 += fabs(tmp);
      rss2 += tmp*tmp;
   }
   cout << "RMSE(test) = " << sqrt(rss2/(N-nt)) << endl;
   cout << "MAE(test) = " << rs2/(N-nt) << endl;
   rss = rss1 + rss2;
   rs = rs1 + rs2;
   cout << "RMSE(total) = " << sqrt(rss/N) << endl;
   cout << "MAE(total) = " << rs/N << endl;
   cout << "#####   Finished calculate_linear_regression   #####\n";
   delete [] y;
   delete [] b;
   delete [] wmatt;
   return true;
}

bool calculate_linear_regression(int N, int M, int Mtot, const char *wmatc_name, int *ids, int nt, double *Yc, double& rmse_train, double& rmse_test, double& rmse_total, double& mae_train, double& mae_test, double& mae_total, int& ndim, bool calculate_test_stats, bool verbose, bool save_coeffs) {
   int i, j, k, nr, nc, mtmm;
   double rss1, rs1, tmp, rs2, rss2, rss, rs;
   double *wmatt, *b, *y, *dp, *dpv;
   if (verbose) {
      cout << "#####   Starting calculate_linear_regression   #####\n";
      cout << "Number of rows: " << N << endl;
      cout << "Number of basis functions: " << M << endl;
      cout << "Training set size: " << nt << endl;
   }
   ifstream ifile(wmatc_name, ios::in | ios::binary);
   if (!ifile) return false;
   ifile.read((char*)&nr, sizeof(int));
   ifile.read((char*)&nc, sizeof(int));
   wmatt = new double[M*nt];
   b = new double[M];
   y = new double[nt];
   dpv = new double[M];
   mtmm = Mtot - M;
   for (i = 0, k = 0; i < N; i++) {
      if (k < nt && ids[k] == i) {
         for (j = 0; j < M; j++) ifile.read((char*)&wmatt[k*M+j], sizeof(double));
         ifile.seekg(mtmm*sizeof(double), ifile.cur);
         y[k] = Yc[i];
         k++;
         if (k == nt) break;
      } else {
         ifile.seekg(Mtot*sizeof(double), ifile.cur);
      }
   }
   ifile.seekg(0, ifile.beg);
   ifile.read((char*)&nr, sizeof(int));
   ifile.read((char*)&nc, sizeof(int));
   multiple_linear_regression_svd(nt, M, wmatt, y, b, ndim, true);
   for (i = 0, rss1 = 0.0, rs1 = 0.0; i < nt; i++) {
      dp = wmatt + i*M;
      tmp = y[i] - dot_prod(b, dp, M);
      rs1 += fabs(tmp);
      rss1 += tmp*tmp;
   }
   rmse_train = sqrt(rss1/nt);
   mae_train = rs1/nt;
   rmse_test = 0.0;
   mae_test = 0.0;
   rmse_total = 0.0;
   mae_total = 0.0;
   if (verbose) {
      cout << "RMSE(train) = " << rmse_train << endl;
      cout << "MAE(train) = " << mae_train << endl;
   }
   if (!calculate_test_stats) goto end;
   for (i = 0, k = 0, rs2 = 0.0, rss2 = 0.0; i < N; i++) {
      if (k < nt && ids[k] == i) {
         ifile.seekg(Mtot*sizeof(double), ifile.cur);
         k++;
      } else {
         for (j = 0; j < M; j++) ifile.read((char*)&dpv[j], sizeof(double));
         ifile.seekg(mtmm*sizeof(double), ifile.cur);
         tmp = Yc[i] - dot_prod(b, dpv, M);
         rs2 += fabs(tmp);
         rss2 += tmp*tmp;
      }
   }
   rmse_test = sqrt(rss2/(N-nt));
   mae_test = rs2/(N-nt);
   if (verbose) {
      cout << "RMSE(test) = " << rmse_test << endl;
      cout << "MAE(test) = " << mae_test << endl;
   }
   rss = rss1 + rss2;
   rs = rs1 + rs2;
   rmse_total = sqrt(rss/N);
   mae_total = rs/N;
   if (verbose) {
      cout << "RMSE(total) = " << rmse_total << endl;
      cout << "MAE(total) = " << mae_total << endl;
   }
end:
   if (verbose)
      cout << "#####   Finished calculate_linear_regression   #####\n";
   if (save_coeffs) {
      Matrix2File(b, M, 1, "coeffs.bin");
   }
   ifile.close();
   delete [] dpv;
   delete [] y;
   delete [] b;
   delete [] wmatt;
   return true;
}

bool build_regression(int N, int M[], int Mtot, int nsub, const char *wmatc_name, double *Yc, int *ids_selected[], int n_training, int *ids_training) {
   int i, ndim, nr, nc;
   int *isel;
   double rmse_train, rmse_test, rmse_total, mae_train, mae_test, mae_total;
   double *lev, *wmatc;
   bool res = true;
   set<int> tr_ids;
   cout << "#####   Starting build_regression   #####\n";
   cout << "Number of rows: " << N << endl;
   cout << "Total number of basis functions: " << Mtot << endl;
   cout << "Number of submatrices: " << nsub << endl;
   cout << "Size of training set: " << n_training << endl;
   isel = new int[N];
   lev = new double[N];
   wmatc = new double[((long long)Mtot)*N];
   File2Matrix(wmatc_name, nr, nc, wmatc);
   order_by_leverages(N, Mtot, wmatc, isel, lev);
   for (i=0; i < n_training; i++) ids_training[i] = isel[i];
   sort(ids_training, ids_training+n_training);
   delete [] wmatc;
   delete [] lev;
   delete [] isel;
   calculate_linear_regression(N, Mtot, Mtot, wmatc_name, ids_training, n_training, Yc, rmse_train, rmse_test, rmse_total, mae_train, mae_test, mae_total, ndim, true, false, true);
   cout << "rmse(train): " << rmse_train << endl;
   cout << "mae(train): " << mae_train << endl;
   cout << "rmse(test): " << rmse_test << endl;
   cout << "mae(test): " << mae_test << endl;
   cout << "rmse(total): " << rmse_total << endl;
   cout << "mae(total): " << mae_total << endl;
   cout << "#####   Finished build_regression   #####\n";
   return res;
}
