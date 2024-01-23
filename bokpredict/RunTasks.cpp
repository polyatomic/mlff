#include <iostream>
#include <iomanip>
#include <thread>

#include "chemutil.h"
#include "RunTasks.h"
#include "util.h"
#include "chemutil.h"
#include "functor.h"
#include "Descriptors.h"
#include "tasks.h"

extern PType g_params;

using std::cerr;
using std::cout;
using std::endl;
using std::setprecision;
using std::numeric_limits;
using std::streamsize;
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
   int i, j, k, l, o, nstr, na, na3, n, nt1, nt2, nt3, m, tp, nt4, na4, ii, jj, kk, ndesc, nr, nc, istart, iend, nthreads, tc;
   double *v, *dmin, *dmax, *B, *stdevs, *avgs, *Y, *Yp;
   double *dvals[MAX_THREADS];
   int *types2, *types3, *stypes3, *types4, *stypes4, *ngp, *tsizes;
   double avg, rss, dif, rss2;
   int pids[6] = {};
   thread threads[MAX_THREADS];
   Iterator<Mol> mols;
   vector<string> tkns, tkns2;
   vector<int> t1;
   vector<Pair> t2;
   vector<Triplet> t3;
   vector<Quadruplet> t4;
   Pair p;
   Triplet t;
   Quadruplet q;
   vector<Pair>::iterator t2pos;
   vector<Triplet>::iterator t3pos;
   vector<Quadruplet>::iterator t4pos;
   stdevs = avgs = B = Y = Yp = 0;
   tsizes = 0;
   nthreads = 0;
   int prec = numeric_limits<double>::max_digits10;
   streamsize oldprec = cout.precision();
   string s = g_params.molFile;
   mols = ReadXYZFile(s.c_str());
   if (!mols) {
      cerr << "Error opening " << s << endl;
      return;
   }
   nstr = mols.Size();
   cout << "Number of configurations = " << nstr << endl;
   na = mols[0].NAtoms();
   cout << "Number of atoms in configuration = " << na << endl;
   n = 3*na;
   j = nstr*n;
   v = new double[j];
   for (i = 0; i < nstr; i++) mols[i].GetCoordinates(v + i*n);
   mols.Resize(1);
   Descriptors descs;
   TFunctorDaDaI<Descriptors> dcalc;
   tkns = Tokenize(g_params.t1, "\r\n");
   nt1 = tkns.size();
   cout << "Number of atom types = " << nt1 << endl;
   cout << "Atom types:" << endl;
   for (i = 0; i < nt1; i++) {
      j = atoi(tkns[i].c_str());
      t1.push_back(j);
      cout << Atom::GetAtomicSymbol(j) << endl;
   }
   tkns = Tokenize(g_params.t2, "\r\n");
   nt2 = tkns.size();
   cout << "Number of atom type pairs = " << nt2 << endl;
   for (i = 0; i < nt2; i++) {
      tkns2 = Tokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      p.v[0] = j;
      j = atoi(tkns2[1].c_str());
      p.v[1] = j;
      t2.push_back(p);
   }
   dmin = new double[nt2];
   dmax = new double[nt2];
   types2 = new int[na*(na-1)/2];
   for (i = 0, m = 0; i < na-1; i++) {
      k = mols[0].GetAtomicNumber(i);
      for (j = i+1; j < na; j++) {
         p.v[0] = k;
         p.v[1] = mols[0].GetAtomicNumber(j);
         if (p.v[0] > p.v[1]) { l = p.v[0]; p.v[0] = p.v[1]; p.v[1] = l; }
         t2pos = lower_bound(t2.begin(), t2.end(), p, PairCompare());
         tp = t2pos - t2.begin();
         types2[m++] = tp;
      }
   }
   tkns = Tokenize(g_params.t3, "\r\n");
   nt3 = tkns.size();
   cout << "Number of atom type triplets = " << nt3 << endl;
   for (i = 0; i < nt3; i++) {
      tkns2 = Tokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      t.v[0] = j;
      j = atoi(tkns2[1].c_str());
      t.v[1] = j;
      j = atoi(tkns2[2].c_str());
      t.v[2] = j;
      t3.push_back(t);
   }
   na3 = na*(na-1)*(na-2)/6;
   types3 = new int[4*na3];
   stypes3 = new int[nt3];
   for (i = 0; i < nt3; i++) stypes3[i] = GetSubType(t3[i]);
   for (i = 0, o = 0; i < na-2; i++) {
      l = mols[0].GetAtomicNumber(i);
      for (j = i+1; j < na-1; j++) {
         m = mols[0].GetAtomicNumber(j);
         for (k = j+1; k < na; k++) {
            t.v[0] = l;
            t.v[1] = m;
            t.v[2] = mols[0].GetAtomicNumber(k);
            analyze_triplet(t, na, i, j, k, pids);
            t3pos = find(t3.begin(), t3.end(), t);
            tp = t3pos - t3.begin();
            types3[o] = tp;
            types3[na3+3*o] = pids[0]; types3[na3+3*o+1] = pids[1]; types3[na3+3*o+2] = pids[2];
            o++;
         }
      }
   }
   tkns = Tokenize(g_params.t4, "\r\n");
   nt4 = tkns.size();
   cout << "Number of atom type quadruplets = " << nt4 << endl;
   for (i = 0; i < nt4; i++) {
      tkns2 = Tokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      q.v[0] = j;
      j = atoi(tkns2[1].c_str());
      q.v[1] = j;
      j = atoi(tkns2[2].c_str());
      q.v[2] = j;
      j = atoi(tkns2[3].c_str());
      q.v[3] = j;
      t4.push_back(q);
   }
   na4 = na*(na-1)*(na-2)*(na-3)/24;
   types4 = new int[7*na4];
   stypes4 = new int[nt4];
   for (i = 0; i < nt4; i++) stypes4[i] = GetSubType(t4[i]);
   for (i = 0, o = 0; i < na-3; i++) {
      ii = mols[0].GetAtomicNumber(i);
      for (j = i+1; j < na-2; j++) {
         jj = mols[0].GetAtomicNumber(j);
         for (k = j+1; k < na-1; k++) {
            kk = mols[0].GetAtomicNumber(k);
            for (l = k+1; l < na; l++) {
               q.v[0] = ii;
               q.v[1] = jj;
               q.v[2] = kk;
               q.v[3] = mols[0].GetAtomicNumber(l);
               analyze_quadruplet(q, na, i, j, k, l, pids);
               t4pos = find(t4.begin(), t4.end(), q);
               tp = t4pos - t4.begin();
               types4[o] = tp;
               types4[na4+6*o] = pids[0]; types4[na4+6*o+1] = pids[1]; types4[na4+6*o+2] = pids[2]; types4[na4+6*o+3] = pids[3]; types4[na4+6*o+4] = pids[4]; types4[na4+6*o+5] = pids[5];
               o++;
            }
         }
      }
   }
   ngp = new int[nt2];
   tkns = Tokenize(g_params.ngp, "\r\n");
   if (tkns.size() != nt2) {
      cerr << "Input error reading NUMBER_OF_GRID_POINTS" << endl;
      goto end;
   }
   for (i = 0; i < nt2; i++)
      ngp[i] = atoi(tkns[i].c_str());
   descs.Init(na, nt2, types2, nt3, stypes3, types3, nt4, stypes4, types4, ngp, g_params.nthreads);
   tkns = Tokenize(g_params.dist, "\r\n");
   if (tkns.size() == nt2) {
      for (i = 0; i < nt2; i++) {
         tkns2 = Tokenize(tkns[i], " ");
         dmin[i] = atof(tkns2[0].c_str());
         dmax[i] = atof(tkns2[1].c_str());
      }
   } else {
      cerr << "Input error reading DISTANCES" << endl;
      goto end;
   }
   cout << "Atom type pairs, minimum and maximum distances, max-min and number of grid points:" << endl;
   for (i = 0; i < nt2; i++)
      cout << Atom::GetAtomicSymbol(t2[i].v[0]) << " " << Atom::GetAtomicSymbol(t2[i].v[1]) << " " << setprecision(prec) << dmin[i] << " " << dmax[i] << " " << setprecision(oldprec) << dmax[i] - dmin[i] << " " << ngp[i] << endl;
   descs.SetGrid(dmin, dmax);
   Descriptors::SetKernel2B(g_params.kernel_2b);
   Descriptors::SetKernel3B(g_params.kernel_3b);
   Descriptors::SetKernel4B(g_params.kernel_4b);
   ndesc = descs.GetNDescriptors();
   cout << "Number of descriptors: " << ndesc << endl;
   B = new double[ndesc+1];
   dcalc.init(&descs, &Descriptors::Calculate);
   s = g_params.linearCoeffs;
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
   s = g_params.eFile;
   Y = new double[nstr];
   Yp = new double[nstr];
   if (s.size() > 0) {
      if (!File2Array(s, Y)) {
         cerr << "Error reading " << s << endl;
         goto end;
      }
   } else {
      cerr << "Input error" << endl;
      goto end;
   }
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
   delete [] Y;
   delete [] avgs;
   delete [] stdevs;
   delete [] B;
   delete [] ngp;
   delete [] stypes4;
   delete [] types4;
   delete [] stypes3;
   delete [] types3;
   delete [] types2;
   delete [] dmax;
   delete [] dmin;
   delete [] v;
}
