#include <iostream>
#include <sstream>
#include <iomanip>

#include "RunTasks.h"
#include "tasks.h"

extern PTypeExtended g_params;

using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::streamsize;
using std::setprecision;
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::ios;

void RunTasks() {
   int i, j, k, l, nstr, n, m, ndesc, ndesc_sketch, ndesc20, ndesc21, ndesc22, ndesc4, ndesc40, ndesc41, nsub, ntrain,
       ndesc2, kk, ndesc3, ndesc30, ndesc31, ndesc32, ndesc33, ndesc42, ndesc43, ndesc44, max_case_count, nr, nc, np;
   int *ids_train;
   double avg, x;
   double *Y, *v, *rcoords, *stdevs, *avgs, *wmatc, *Yc, *stdevs_to_save, *avgs_to_save;
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
   Iterator<Mol> mols;
   Descriptors descs;
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
   Y = v = wmatc = Yc = stdevs = avgs = stdevs_to_save = avgs_to_save = 0;
   ids_train = 0;
   int prec = numeric_limits<double>::max_digits10;
   streamsize oldprec = cout.precision();
   TFunctorDaDaI<Descriptors> dcalc;
   if (!calculation_prepare(&g_params, mols, descs, nstr, n, v, Y)) goto end;
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
         for (j = 0; j < 12; j++) {
            ofile.open((string(blocks[j]) + "s.txt").c_str());
            for (i = 0; i < nstr; i++) ofile << ids_selected[j][i] << endl;
            ofile.close();
         }
         ofile.open("ranks.txt");
         for (i = 0; i < nsub; i++) ofile << ranks[i] << endl;
         ofile.close();
      } else {
         for (j = 0; j < 12; j++) {
            if (!File2Array(string(blocks[j]) + "s.txt", ids_selected[j])) goto end;
         }
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
      max_case_count = (g_params.rbs_2b + g_params.rbs_3b + g_params.rbs_4b)*1000;
      cout << "max_case_count: " << max_case_count << endl;
      for (j = 0; j < 12; j++) {
         if (!File2Array(string(blocks[j]) + "s.txt", ids_selected[j])) goto end;
      }
      if (!File2Array("ndescs_sketch.txt", ndescs)) goto end;
      ids_train = new int[max_case_count];
      for (i=0, ndesc=0; i < nsub; i++) ndesc += ndescs[i];
      ntrain = 0;
      build_regression(nstr, ndescs, ndesc, 0, 3, "descs_sketch.bin", Yc, ids_selected, 1000, g_params.rbs_2b, ntrain, ids_train);
      if (nsub > 3) {
         build_regression(nstr, ndescs, ndesc, 3, 4, "descs_sketch.bin", Yc, ids_selected, 1000, g_params.rbs_3b, ntrain, ids_train);
         if (nsub > 7)
            build_regression(nstr, ndescs, ndesc, 7, 5, "descs_sketch.bin", Yc, ids_selected, 1000, g_params.rbs_4b, ntrain, ids_train, true);
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
   delete [] v;
   delete [] Y;
}
