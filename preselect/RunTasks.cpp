#include <iostream>

#include "RunTasks.h"

extern PTypeExtended g_params;

using std::cout;
using std::endl;
using std::ofstream;

void RunTasks() {
   int i, j, k, nstr, n;
   double *Y, *v, *rcoords;
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
   vector<int> ivec;
   set<int>::iterator itp;
   ofstream ofile;
   Y = v = 0;
   blocks = 0;
   descs.SetUseInverseSpace(g_params.uis);
   if (!calculation_prepare(&g_params, mols, descs, nstr, n, v, Y, nblocks, &blocks)) goto end;
   pss = new set<int>[nblocks[0]];
   tss = new set<Triplet, TripletCompare>[nblocks[1]];
   tss_ref = new set<Triplet, TripletCompare>[nblocks[1]];
   sss = new set<Sixtuplet, SixtupletCompare>[nblocks[2]];
   sss_ref = new set<Sixtuplet, SixtupletCompare>[nblocks[2]];
   for (i = 0; i < nstr; i++) {
      rcoords = v + i*n;
      descs.FindClosestMultiplets(rcoords, pss, tss, sss);
   }
   cout << "2B descriptors selected according to closeness to training points" << endl;
   for (i = 0; i < nblocks[0]; i++) {
      cout << "Descriptors of type " << i << ": " << pss[i].size() << endl;
      ofile.open(string(blocks+5*i) + ".txt");
      ofile << pss[i].size() << endl;
      for (itp = pss[i].begin(); itp != pss[i].end(); itp++) {
         ofile << *itp << endl;
      }
      ofile.close();
   }
   descs.Get3BGridPoints(tss_ref);
   cout << "3B descriptors selected according to closeness to training points" << endl;
   for (i = 0; i < nblocks[1]; i++) {
      ivec.clear();
      j = 0;
      for (it = tss[i].begin(); it != tss[i].end(); it++) {
         it2 = tss_ref[i].find(*it);
         if (it2 != tss_ref[i].end()) {
            ivec.push_back((*it).v[0]); ivec.push_back((*it).v[1]); ivec.push_back((*it).v[2]);
            j++;
         }
      }
      cout << "Descriptors of type " << i << ": " << j << endl;
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
      j = 0;
      for (its = sss[i].begin(); its != sss[i].end(); its++) {
         its2 = sss_ref[i].find(*its);
         if (its2 != sss_ref[i].end()) {
            ivec.push_back((*its).v[0]); ivec.push_back((*its).v[1]); ivec.push_back((*its).v[2]); ivec.push_back((*its).v[3]); ivec.push_back((*its).v[4]); ivec.push_back((*its).v[5]);
            j++;
         }
      }
      cout << "Descriptors of type " << i << ": " << j << endl;
      if (j > 0) {
         ofile.open(string(blocks+5*(nblocks[0]+nblocks[1]+i)) + ".txt");
         ofile << j << endl;
         for (k = 0; k < j; k++) {
            ofile << ivec[k*6] << " " << ivec[k*6+1] << " " << ivec[k*6+2] << " " << ivec[k*6+3] << " " << ivec[k*6+4] << " " << ivec[k*6+5] << endl;
         }
         ofile.close();
      }
   }
end:
   delete [] sss_ref;
   delete [] sss;
   delete [] tss_ref;
   delete [] tss;
   delete [] pss;
   delete [] blocks;
   delete [] v;
   delete [] Y;
}
