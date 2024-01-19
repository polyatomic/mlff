#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "chemutil.h"
#include "RunTasks.h"
#include "util.h"
#include "Descriptors.h"
#include "tasks.h"
#include "mathutil.h"

extern PType g_params;

using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::streamsize;
using std::setprecision;
using std::set_difference;
using std::sort;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::ios;

void RunTasks() {
   int i, j, k, l, nstr, na, n, nt1, nt2, nt4, tp, m, ndesc, ndesc_sketch, ndesc20, ndesc21, ndesc22, ndesc4, ndesc40, ndesc41, nsub, ntrain,
       nt3, na3, o, ndesc2, na4, ii, jj, kk, ndesc3, ndesc30, ndesc31, ndesc32, ndesc33, ndesc42, ndesc43, ndesc44, max_case_count, nr, nc, np;
   int *types2, *ngp, *types3, *stypes3, *ids_train, *types4, *stypes4;
   double avg, x;
   double *Y, *v, *dmin, *dmax, *rcoords, *stdevs, *avgs, *wmatc, *Yc, *stdevs_to_save, *avgs_to_save;
   int pids[6] = {};
   int n_train3B[4] = {};
   int n_train4B[5] = {};
   int ndescs[12] = {};
   int ndescs_sketch[12] = {};
   int ranks[12] = {};
   int *ids_selected[12] = {};
   int *vars_selected[12] = {};
   double *levs[12] = {};
   const char *blocks[12] = {
      "g20",
      "g21",
      "g22",
      "g30",
      "g31",
      "g32",
      "g33",
      "g40",
      "g41",
      "g42",
      "g43",
      "g44"
   };
   Pair p;
   Triplet t;
   Quadruplet q;
   vector<int> t1;
   vector<Pair> t2;
   vector<Triplet> t3;
   vector<Quadruplet> t4;
   Iterator<Mol> mols;
   vector<string> tkns, tkns2;
   vector<Pair>::iterator t2pos;
   vector<Triplet>::iterator t3pos;
   vector<Quadruplet>::iterator t4pos;
   set<int> pss[3];
   set<Triplet, TripletCompare> tss[4];
   set<Triplet, TripletCompare> tss_ref[4];
   set<Triplet, TripletCompare>::iterator it, it2;
   set<Sixtuplet, SixtupletCompare> sss[5];
   set<Sixtuplet, SixtupletCompare> sss_ref[5];
   set<Sixtuplet, SixtupletCompare>::iterator its, its2;
   set<int>::iterator itp;
   vector<int> ivec;
   ofstream ofile;
   ifstream ifile;
   istringstream isstream;
   string line;
   int prec = numeric_limits<double>::max_digits10;
   streamsize oldprec = cout.precision();
   string s = g_params.molFile;
   mols = ReadXYZFile(s.c_str());
   if (!mols) {
      cerr << "Error opening " << s << endl;
      return;
   }
   nstr = mols.Size();
   cout << "Number of training configurations = " << nstr << endl;
   na = mols[0].NAtoms();
   cout << "Number of atoms in configuration = " << na << endl;
   s = g_params.eFile;
   Y = new double[nstr];
   if (s.size() > 0) {
      if (!File2Array(s, Y)) {
         cerr << "Error reading " << s << endl;
         delete [] Y;
         return;
      }
   } else {
      cerr << "Input error" << endl;
      delete [] Y;
      return;
   }
   n = 3*na;
   j = nstr*n;
   v = new double[j];
   for (i=0; i < nstr; i++) mols[i].GetCoordinates(v + i*n);
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
   for (i=0; i < nt3; i++) stypes3[i] = GetSubType(t3[i]);
   for (i = 0,o = 0; i < na-2; i++) {
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
   for (i=0; i < nt4; i++) stypes4[i] = GetSubType(t4[i]);
   for (i=0,o=0; i < na-3; i++) {
      ii = mols[0].GetAtomicNumber(i);
      for (j=i+1; j < na-2; j++) {
         jj = mols[0].GetAtomicNumber(j);
         for (k=j+1; k < na-1; k++) {
            kk = mols[0].GetAtomicNumber(k);
            for (l=k+1; l < na; l++) {
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
   stdevs = avgs = stdevs_to_save = avgs_to_save = wmatc = Yc = 0;
   ids_train = 0;
   if (tkns.size() != nt2) {
      cerr << "Input error reading NUMBER_OF_GRID_POINTS" << endl;
      goto end;
   }
   for (i = 0; i < nt2; i++)
      ngp[i] = atoi(tkns[i].c_str());
   descs.Init(na, nt2, types2, nt3, stypes3, types3, nt4, stypes4, types4, ngp, g_params.nthreads);
   if (g_params.dist.size() == 0) {
      for (i = 0; i < nt2; i++) {
         dmin[i] = 1000.0;
         dmax[i] = 0.0;
      }
      for (i = 0; i < nstr; i++) {
         rcoords = v + i*n;
         descs.GetDistanceRanges(rcoords, dmin, dmax);
      }
   }
   cout << "Atom type pairs, minimum and maximum distances, max-min and number of grid points:" << endl;
   for (i = 0; i < nt2; i++)
      cout << Atom::GetAtomicSymbol(t2[i].v[0]) << " " << Atom::GetAtomicSymbol(t2[i].v[1]) << " " << setprecision(prec) << dmin[i] << " " << dmax[i] << " " << setprecision(oldprec) << dmax[i] - dmin[i] << " " << ngp[i] << endl;
   descs.SetGrid(dmin, dmax);
   if (g_params.preselect_descriptors) {
      for (i = 0; i < nstr; i++) {
         rcoords = v + i*n;
         descs.FindClosestMultiplets(rcoords, pss, tss, sss);
      }
      cout << "2B descriptors selected according to closeness to training points" << endl;
      for (i = 0; i < 3; i++) {
         cout << "Descriptors of type " << i << endl;
         if (i == 0) ofile.open("g20.txt");
         else if (i == 1) ofile.open("g21.txt");
         else ofile.open("g22.txt");
         cout << pss[i].size() << endl;
         ofile << pss[i].size() << endl;
         for (itp = pss[i].begin(); itp != pss[i].end(); itp++) {
            cout << *itp << endl;
            ofile << *itp << endl;
         }
         ofile.close();
      }
      descs.Get3BGridPoints(tss_ref);
      cout << "3B descriptors selected according to closeness to training points" << endl;
      for (i = 0; i < 4; i++) {
         ivec.clear();
         cout << "Descriptors of type " << i << endl;
         cout << tss[i].size() << endl;
         j = 0;
         for (it = tss[i].begin(); it != tss[i].end(); it++) {
            it2 = tss_ref[i].find(*it);
            if (it2 != tss_ref[i].end()) {
               cout << (*it).v[0] << " " << (*it).v[1] << " " << (*it).v[2] << endl;
               ivec.push_back((*it).v[0]); ivec.push_back((*it).v[1]); ivec.push_back((*it).v[2]);
               j++;
            } else {
               cerr << (*it).v[0] << " " << (*it).v[1] << " " << (*it).v[2] << " not found for type " << i << endl;
            }
         }
         cout << j << endl;
         if (j > 0) {
            if (i == 0) ofile.open("g30.txt");
            else if (i == 1) ofile.open("g31.txt");
            else if (i == 2) ofile.open("g32.txt");
            else ofile.open("g33.txt");
            ofile << j << endl;
            for (k=0; k < j; k++) {
               ofile << ivec[k*3] << " " << ivec[k*3+1] << " " << ivec[k*3+2] << endl;
            }
            ofile.close();
         }
      }
      descs.Get4BGridPoints(sss_ref);
      cout << "4B descriptors selected according to closeness to training points" << endl;
      for (i = 0; i < 5; i++) {
         ivec.clear();
         cout << "Descriptors of type " << i << endl;
         cout << sss[i].size() << endl;
         j = 0;
         for (its = sss[i].begin(); its != sss[i].end(); its++) {
            its2 = sss_ref[i].find(*its);
            if (its2 != sss_ref[i].end()) {
               cout << (*its).v[0] << " " << (*its).v[1] << " " << (*its).v[2] << " " << (*its).v[3] << " " << (*its).v[4] << " " << (*its).v[5] << endl;
               ivec.push_back((*its).v[0]); ivec.push_back((*its).v[1]); ivec.push_back((*its).v[2]); ivec.push_back((*its).v[3]); ivec.push_back((*its).v[4]); ivec.push_back((*its).v[5]);
               j++;
            } else {
               cerr << (*its).v[0] << " " << (*its).v[1] << " " << (*its).v[2] << " " << (*its).v[3] << " " << (*its).v[4] << " " << (*its).v[5] << " not found for type " << i << endl;
            }
         }
         cout << j << endl;
         if (j > 0) {
            if (i == 0) ofile.open("g40.txt");
            else if (i == 1) ofile.open("g41.txt");
            else if (i == 2) ofile.open("g42.txt");
            else if (i == 3) ofile.open("g43.txt");
            else ofile.open("g44.txt");
            ofile << j << endl;
            for (k = 0; k < j; k++) {
               ofile << ivec[k*6] << " " << ivec[k*6+1] << " " << ivec[k*6+2] << " " << ivec[k*6+3] << " " << ivec[k*6+4] << " " << ivec[k*6+5] << endl;
            }
            ofile.close();
         }
      }
      ivec.clear();
      goto end;
   }
   Descriptors::SetKernel2B(g_params.kernel_2b);
   Descriptors::SetKernel3B(g_params.kernel_3b);
   Descriptors::SetKernel4B(g_params.kernel_4b);
   ndesc = descs.GetNDescriptors();
   stdevs = new double[ndesc];
   avgs = new double[ndesc];
   stdevs_to_save = new double[ndesc];
   avgs_to_save = new double[ndesc];
   for (i=0; i < 12; i++) ids_selected[i] = new int[nstr];
   for (i=0; i < 12; i++) levs[i] = new double[nstr];
   // wmatc = new double[nstr*ndesc];
   cout << "Number of descriptors: " << ndesc << endl;
   dcalc.init(&descs, &Descriptors::Calculate);
   if (!g_params.read_descriptors) {
      calculate_descriptors(nstr, ndesc, n, v, &dcalc, stdevs, avgs, g_params.nthreads, g_params.maxrows);
      // Matrix2File(wmatc, nstr, ndesc, "descs.bin");
      Matrix2File(stdevs, ndesc, 1, "stdevs.bin");
      Matrix2File(avgs, ndesc, 1, "avgs.bin");
      goto end;
   } else {
      // Add error handling here for the case bin files don't exist
      // File2Matrix("descs.bin", nstr, ndesc, wmatc);
      File2Matrix("stdevs.bin", ndesc, i, stdevs);
      File2Matrix("avgs.bin", ndesc, i, avgs);
   }
   ndesc2 = descs.GetN2BDescriptors();
   cout << "Number of 2B descriptors: " << ndesc2 << endl;
   ndesc20 = descs.GetN2BDescriptors(0);
   ndesc21 = descs.GetN2BDescriptors(1);
   ndesc22 = descs.GetN2BDescriptors(2);
   cout << "   type 0: " << ndesc20 << endl;
   cout << "   type 1: " << ndesc21 << endl;
   cout << "   type 2: " << ndesc22 << endl;
   Yc = new double[nstr];
   for (i = 0, avg = 0.0; i < nstr; i++) avg += Y[i];
   avg /= nstr;
   for (i = 0; i < nstr; i++) Yc[i] = Y[i] - avg;
   ndescs[0] = ndesc20;
   ndescs[1] = ndesc21;
   ndescs[2] = ndesc22;
   ndesc3 = descs.GetN3BDescriptors();
   cout << "Number of 3B descriptors: " << ndesc3 << endl;
   ndesc30 = descs.GetN3BDescriptors(0);
   ndesc31 = descs.GetN3BDescriptors(1);
   ndesc32 = descs.GetN3BDescriptors(2);
   ndesc33 = descs.GetN3BDescriptors(3);
   cout << "   type 0: " << ndesc30 << endl;
   cout << "   type 1: " << ndesc31 << endl;
   cout << "   type 2: " << ndesc32 << endl;
   cout << "   type 3: " << ndesc33 << endl;
   ndescs[3] = ndesc30;
   ndescs[4] = ndesc31;
   ndescs[5] = ndesc32;
   ndescs[6] = ndesc33;
   ndesc4 = descs.GetN4BDescriptors();
   cout << "Number of 4B descriptors: " << ndesc4 << endl;
   ndesc40 = descs.GetN4BDescriptors(0);
   ndesc41 = descs.GetN4BDescriptors(1);
   ndesc42 = descs.GetN4BDescriptors(2);
   ndesc43 = descs.GetN4BDescriptors(3);
   ndesc44 = descs.GetN4BDescriptors(4);
   cout << "   type 0: " << ndesc40 << endl;
   cout << "   type 1: " << ndesc41 << endl;
   cout << "   type 2: " << ndesc42 << endl;
   cout << "   type 3: " << ndesc43 << endl;
   cout << "   type 4: " << ndesc44 << endl;
   ndescs[7] = ndesc40;
   ndescs[8] = ndesc41;
   ndescs[9] = ndesc42;
   ndescs[10] = ndesc43;
   ndescs[11] = ndesc44;
   if (!ndesc3)
      nsub = 3;
   else if (!ndesc4)
      nsub = 7;
   else
      nsub = 12;
   if (!g_params.build_regression) {
      if (g_params.order_by_leverages) {
         order_by_leverages(nstr, ndescs, nsub, "descs.bin", "descst.bin", ids_selected, levs, ranks, g_params.maxrows, g_params.nthreads);
         ofile.open("g20s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[0][i] << endl;
         ofile.close();
         ofile.open("g21s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[1][i] << endl;
         ofile.close();
         ofile.open("g22s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[2][i] << endl;
         ofile.close();
         ofile.open("g30s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[3][i] << endl;
         ofile.close();
         ofile.open("g31s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[4][i] << endl;
         ofile.close();
         ofile.open("g32s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[5][i] << endl;
         ofile.close();
         ofile.open("g33s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[6][i] << endl;
         ofile.close();
         ofile.open("g40s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[7][i] << endl;
         ofile.close();
         ofile.open("g41s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[8][i] << endl;
         ofile.close();
         ofile.open("g42s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[9][i] << endl;
         ofile.close();
         ofile.open("g43s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[10][i] << endl;
         ofile.close();
         ofile.open("g44s.txt");
         for (i = 0; i < nstr; i++) ofile << ids_selected[11][i] << endl;
         ofile.close();
         ofile.open("ranks.txt");
         for (i = 0; i < nsub; i++) ofile << ranks[i] << endl;
         ofile.close();
      } else {
         if (!File2Array("g20s.txt", ids_selected[0])) goto end;
         if (!File2Array("g21s.txt", ids_selected[1])) goto end;
         if (!File2Array("g22s.txt", ids_selected[2])) goto end;
         if (!File2Array("g30s.txt", ids_selected[3])) goto end;
         if (!File2Array("g31s.txt", ids_selected[4])) goto end;
         if (!File2Array("g32s.txt", ids_selected[5])) goto end;
         if (!File2Array("g33s.txt", ids_selected[6])) goto end;
         if (!File2Array("g40s.txt", ids_selected[7])) goto end;
         if (!File2Array("g41s.txt", ids_selected[8])) goto end;
         if (!File2Array("g42s.txt", ids_selected[9])) goto end;
         if (!File2Array("g43s.txt", ids_selected[10])) goto end;
         if (!File2Array("g44s.txt", ids_selected[11])) goto end;
         if (!File2Array("ranks.txt", ranks)) goto end;
         for (i=0; i < nsub; i++) {
            vars_selected[i] = new int[ranks[i]];
            for (j = 0; j < ranks[i]; j++) vars_selected[i][j] = -1;
         }
         sketch_matrices(nstr, ndescs, nsub, "descs.bin", ids_selected, ranks, g_params.nsketch, g_params.nthreads, vars_selected);
         for (i = 0, l = 0, m = 0; i < nsub; i++) {
            if (vars_selected[i][0] == -1) {
               ndescs_sketch[i] = ndescs[i];
               copy_file((string(blocks[i]) + ".txt").c_str(), (string(blocks[i]) + "f.txt").c_str());
               for (j = 0; j < ndescs_sketch[i]; j++) {
                  avgs_to_save[l] = avgs[m];
                  stdevs_to_save[l++] = stdevs[m++];
               }
            } else {
               ndescs_sketch[i] = ranks[i];
               ifile.open((string(blocks[i]) + ".txt").c_str(), ios::in);
               ofile.open((string(blocks[i]) + "f.txt").c_str(), ios::out);
               getline(ifile, line, '\n');
               isstream.str(line);
               isstream >> np;
               isstream.clear();
               ofile << ndescs_sketch[i] << endl;
               for (j = 0, k = 0; j < np; j++) {
                  getline(ifile, line, '\n');
                  if (k < ndescs_sketch[i] && j == vars_selected[i][k]) {
                     ofile << line << endl;
                     k++;
                     avgs_to_save[l] = avgs[m];
                     stdevs_to_save[l++] = stdevs[m++];
                  } else {
                     m++;
                  }
               }
               ifile.close();
               ofile.close();
            }
         }
         for (i=0, ndesc_sketch=0; i < nsub; i++) {
            ndesc_sketch += ndescs_sketch[i];
         }
         if (l != ndesc_sketch) cerr << "Internal error l != ndesc_sketch" << endl;
         Matrix2File(avgs_to_save, ndesc_sketch, 1, "avgsf.bin");
         Matrix2File(stdevs_to_save, ndesc_sketch, 1, "stdevsf.bin");
         ofile.open("descs_sketch.bin", ios::out | ios::binary);
         ofile.write((char*)&nstr, sizeof(int));
         ofile.write((char*)&ndesc_sketch, sizeof(int));
         ifile.open("descs.bin", ios::in | ios::binary);
         ifile.read((char*)&nr, sizeof(int));
         ifile.read((char*)&nc, sizeof(int));
         for (i = 0; i < nstr; i++) {
            for (j=0; j < nsub; j++) {
               if (vars_selected[j][0] == -1) {
                  for (k = 0; k < ndescs[j]; k++) {
                     ifile.read((char*)&x, sizeof(double));
                     ofile.write((char*)&x, sizeof(double));
                  }
               } else {
                  for (k=0,kk=0; k < ndescs[j]; k++) {
                     ifile.read((char*)&x, sizeof(double));
                     if (kk < ndescs_sketch[j] && k == vars_selected[j][kk]) {
                        kk++;
                        ofile.write((char*)&x, sizeof(double));
                     }
                  }
               }
            }
         }
         ofile.close();
         ifile.close();
         ofile.open("ndescs_sketch.txt");
         for (i=0; i < nsub; i++) ofile << ndescs_sketch[i] << endl;
         ofile.close();
      }
   } else {
      max_case_count = (g_params.rbs_2b + g_params.rbs_3b + g_params.rbs_4b)*100;
      cout << "max_case_count: " << max_case_count << endl;
      if (!File2Array("g20s.txt", ids_selected[0])) goto end;
      if (!File2Array("g21s.txt", ids_selected[1])) goto end;
      if (!File2Array("g22s.txt", ids_selected[2])) goto end;
      if (!File2Array("g30s.txt", ids_selected[3])) goto end;
      if (!File2Array("g31s.txt", ids_selected[4])) goto end;
      if (!File2Array("g32s.txt", ids_selected[5])) goto end;
      if (!File2Array("g33s.txt", ids_selected[6])) goto end;
      if (!File2Array("g40s.txt", ids_selected[7])) goto end;
      if (!File2Array("g41s.txt", ids_selected[8])) goto end;
      if (!File2Array("g42s.txt", ids_selected[9])) goto end;
      if (!File2Array("g43s.txt", ids_selected[10])) goto end;
      if (!File2Array("g44s.txt", ids_selected[11])) goto end;
      if (!File2Array("ndescs_sketch.txt", ndescs)) goto end;
      ids_train = new int[max_case_count];
      for (i=0, ndesc=0; i < nsub; i++) ndesc += ndescs[i];
      ntrain = 0;
      build_regression(nstr, ndescs, ndesc, 0, 3, "descs_sketch.bin", Yc, ids_selected, 100, g_params.rbs_2b, ntrain, ids_train);
      if (nsub > 3) {
         build_regression(nstr, ndescs, ndesc, 3, 4, "descs_sketch.bin", Yc, ids_selected, 100, g_params.rbs_3b, ntrain, ids_train);
         if (nsub > 7)
            build_regression(nstr, ndescs, ndesc, 7, 5, "descs_sketch.bin", Yc, ids_selected, 100, g_params.rbs_4b, ntrain, ids_train, true);
      }
      cout << "Dependent variable average: " << setprecision(prec) << avg << setprecision(oldprec) << endl;
   }
end:
   delete [] ids_train;
   delete [] Yc;
   for (i=0; i < 12; i++) {
      delete [] levs[i];
      delete [] vars_selected[i];
      delete [] ids_selected[i];
   }
   delete [] avgs_to_save;
   delete [] stdevs_to_save;
   delete [] avgs;
   delete [] stdevs;
   delete [] ngp;
   delete [] stypes4;
   delete [] types4;
   delete [] stypes3;
   delete [] types3;
   delete [] types2;
   delete [] dmax;
   delete [] dmin;
   delete [] v;
   delete [] Y;
}
