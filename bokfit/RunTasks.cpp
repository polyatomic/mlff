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
   int i, j, k, l, nstr, n, m, ndesc, ndesc_sketch, ndesc4, nsub, ntrain,
       ndesc2, kk, ndesc3, nr, nc, np, nb_total;
   int *ids_train, *ndescs, *ndescs_sketch, *ranks;
   int **ids_selected, **vars_selected;
   double avg, x, rmse;
   double *Y, *v, *rcoords, *stdevs, *avgs, *wmatc, *Yc, *stdevs_to_save, *avgs_to_save;
   double **levs;
   int nblocks[3];
   char *blocks;
   set<int> *pss = nullptr;
   set<Triplet, TripletCompare> *tss = nullptr;
   set<Triplet, TripletCompare> *tss_ref = nullptr;
   set<Sixtuplet, SixtupletCompare> *sss = nullptr;
   set<Sixtuplet, SixtupletCompare> *sss_ref = nullptr;
   set<Triplet, TripletCompare>::iterator it, it2;
   set<Sixtuplet, SixtupletCompare>::iterator its, its2;
   Iterator<Mol> mols;
   Descriptors descs;
   set<int>::iterator itp;
   vector<int> ivec;
   ofstream ofile;
   ifstream ifile;
   istringstream isstream;
   string line;
   Y = v = wmatc = Yc = stdevs = avgs = stdevs_to_save = avgs_to_save = 0;
   ids_train = ndescs = ndescs_sketch = ranks = 0;
   ids_selected = vars_selected = 0;
   blocks = 0;
   levs = 0;
   int prec = numeric_limits<double>::max_digits10;
   streamsize oldprec = cout.precision();
   TFunctorDaDaI<Descriptors> dcalc;
   descs.SetUseInverseSpace(g_params.uis);
   if (!calculation_prepare(&g_params, mols, descs, nstr, n, v, Y, nblocks, &blocks)) goto end;
   pss = new set<int>[nblocks[0]];
   tss = new set<Triplet, TripletCompare>[nblocks[1]];
   tss_ref = new set<Triplet, TripletCompare>[nblocks[1]];
   sss = new set<Sixtuplet, SixtupletCompare>[nblocks[2]];
   sss_ref = new set<Sixtuplet, SixtupletCompare>[nblocks[2]];
   nb_total = nblocks[0] + nblocks[1] + nblocks[2];
   ndescs = new int[nb_total];
   ndescs_sketch = new int[nb_total];
   ranks = new int[nb_total];
   ids_selected = new int*[nb_total];
   vars_selected = new int*[nb_total];
   levs = new double*[nb_total];
   for (i=0; i < nb_total; i++) {
      ndescs[i] = 0;
      ndescs_sketch[i] = 0;
      ranks[i] = 0;
      ids_selected[i] = 0;
      vars_selected[i] = 0;
      levs[i] = 0;
   }
   if (g_params.preselect_descriptors) {
      for (i = 0; i < nstr; i++) {
         rcoords = v + i*n;
         descs.FindClosestMultiplets(rcoords, pss, tss, sss);
      }
      cout << "2B descriptors selected according to closeness to training points" << endl;
      for (i = 0; i < nblocks[0]; i++) {
         cout << "Descriptors of type " << i << endl;
         ofile.open(string(blocks+5*i) + ".txt");
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
      for (i = 0; i < nblocks[1]; i++) {
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
            ofile.open(string(blocks+5*(nblocks[0]+i)) + ".txt");
            ofile << j << endl;
            for (k=0; k < j; k++) {
               ofile << ivec[k*3] << " " << ivec[k*3+1] << " " << ivec[k*3+2] << endl;
            }
            ofile.close();
         }
      }
      descs.Get4BGridPoints(sss_ref);
      cout << "4B descriptors selected according to closeness to training points" << endl;
      for (i = 0; i < nblocks[2]; i++) {
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
            ofile.open(string(blocks+5*(nblocks[0]+nblocks[1]+i)) + ".txt");
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
   for (i=0; i < nb_total; i++) ids_selected[i] = new int[nstr];
   for (i=0; i < nb_total; i++) levs[i] = new double[nstr];
   // wmatc = new double[nstr*ndesc];
   cout << "Number of descriptors: " << ndesc << endl;
   dcalc.init(&descs, &Descriptors::Calculate);
   if (!g_params.read_descriptors) {
      if (g_params.ues) {
         calculate_descriptors(nstr, ndesc, n, v, &dcalc, stdevs, avgs, g_params.nthreads, g_params.maxrows, "stdevs.bin", "avgs.bin");
      } else {
         calculate_descriptors(nstr, ndesc, n, v, &dcalc, stdevs, avgs, g_params.nthreads, g_params.maxrows);
         Matrix2File(stdevs, ndesc, 1, "stdevs.bin");
         Matrix2File(avgs, ndesc, 1, "avgs.bin");
      }
      // Matrix2File(wmatc, nstr, ndesc, "descs.bin");
      goto end;
   } else {
      // Add error handling here for the case bin files don't exist
      // File2Matrix("descs.bin", nstr, ndesc, wmatc);
      File2Matrix("stdevs.bin", ndesc, i, stdevs);
      File2Matrix("avgs.bin", ndesc, i, avgs);
   }
   ndesc2 = descs.GetN2BDescriptors();
   cout << "Number of 2B descriptors: " << ndesc2 << endl;
   for (i = 0; i < nblocks[0]; i++) {
      ndescs[i] = descs.GetN2BDescriptors(i);
      cout << "   type " << i << ": " << ndescs[i] << endl;
   }
   Yc = new double[nstr];
   for (i = 0, avg = 0.0; i < nstr; i++) avg += Y[i];
   avg /= nstr;
   for (i = 0, rmse = 0.0; i < nstr; i++) {
      Yc[i] = Y[i] - avg;
      rmse += Yc[i]*Yc[i];
   }
   rmse = sqrt(rmse/nstr);
   ndesc3 = descs.GetN3BDescriptors();
   cout << "Number of 3B descriptors: " << ndesc3 << endl;
   for (i = nblocks[0], j=0; i < nblocks[0] + nblocks[1]; i++, j++) {
      ndescs[i] = descs.GetN3BDescriptors(j);
      cout << "   type " << j << ": " << ndescs[i] << endl;
   }
   ndesc4 = descs.GetN4BDescriptors();
   cout << "Number of 4B descriptors: " << ndesc4 << endl;
   for (i = nblocks[0] + nblocks[1], j=0; i < nb_total; i++, j++) {
      ndescs[i] = descs.GetN4BDescriptors(j);
      cout << "   type " << j << ": " << ndescs[i] << endl;
   }
   if (!ndesc3 && !ndesc4)
      nsub = nblocks[0];
   else if (!ndesc4)
      nsub = nblocks[0] + nblocks[1];
   else
      nsub = nb_total;
   if (!g_params.build_regression) {
      if (g_params.order_by_leverages) {
         order_by_leverages(nstr, ndescs, nsub, "descs.bin", "descst.bin", ids_selected, levs, ranks, g_params.maxrows, g_params.nthreads);
         for (j = 0; j < nb_total; j++) {
            if (ndescs[j]) {
               ofile.open((string(blocks+5*j) + "s.txt").c_str());
               for (i = 0; i < nstr; i++) ofile << ids_selected[j][i] << endl;
               ofile.close();
            }
         }
         ofile.open("ranks.txt");
         for (i = 0; i < nsub; i++) ofile << ranks[i] << endl;
         ofile.close();
      } else {
         for (j = 0; j < nb_total; j++) {
            if (!File2Array(string(blocks+5*j) + "s.txt", ids_selected[j])) goto end;
         }
         if (!File2Array("ranks.txt", ranks)) goto end;
         for (i=0; i < nsub; i++) {
            vars_selected[i] = new int[ranks[i]];
            for (j = 0; j < ranks[i]; j++) vars_selected[i][j] = -1;
         }
         //sketch_matrices(nstr, ndescs, nblocks, nsub, "descs.bin", ids_selected, ranks, g_params.nsketch, g_params.nthreads, vars_selected);
         sketch_matrices2(nstr, ndescs, nblocks, nsub, "descs.bin", ids_selected, ranks, g_params.nsketch, g_params.nthreads, vars_selected);
         for (i = 0, l = 0, m = 0; i < nsub; i++) {
            if (vars_selected[i][0] == -1) {
               ndescs_sketch[i] = ndescs[i];
               copy_file((string(blocks+5*i) + ".txt").c_str(), (string(blocks+5*i) + "f.txt").c_str());
               for (j = 0; j < ndescs_sketch[i]; j++) {
                  avgs_to_save[l] = avgs[m];
                  stdevs_to_save[l++] = stdevs[m++];
               }
            } else {
               ndescs_sketch[i] = ranks[i];
               ifile.open((string(blocks+5*i) + ".txt").c_str(), ios::in);
               ofile.open((string(blocks+5*i) + "f.txt").c_str(), ios::out);
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
      for (j = 0; j < nb_total; j++) {
         if (!File2Array(string(blocks+5*j) + "s.txt", ids_selected[j])) goto end;
      }
      if (!File2Array("ndescs_sketch.txt", ndescs)) goto end;
      ntrain = g_params.tss;
      ids_train = new int[ntrain];
      for (i=0, ndesc=0; i < nsub; i++) ndesc += ndescs[i];
      build_regression(nstr, ndescs, ndesc, nsub, "descs_sketch.bin", Yc, ids_selected, ntrain, ids_train);
      cout << "Dependent variable average: " << setprecision(prec) << avg << endl;
      cout << "Dependent variable rmse: " << rmse << setprecision(oldprec) << endl;
   }
end:
   delete [] ids_train;
   delete [] Yc;
   for (i=0; i < nb_total; i++) {
      delete [] levs[i];
      delete [] vars_selected[i];
      delete [] ids_selected[i];
   }
   delete [] avgs_to_save;
   delete [] stdevs_to_save;
   delete [] avgs;
   delete [] stdevs;
   delete [] levs;
   delete [] vars_selected;
   delete [] ids_selected;
   delete [] ranks;
   delete [] ndescs_sketch;
   delete [] ndescs;
   delete [] sss_ref;
   delete [] sss;
   delete [] tss_ref;
   delete [] tss;
   delete [] pss;
   delete [] blocks;
   delete [] v;
   delete [] Y;
}
