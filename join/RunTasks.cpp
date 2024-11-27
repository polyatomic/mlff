#include <iostream>
#include <fstream>

#include "RunTasks.h"

extern PType g_params;

using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std:: cerr;
using std::endl;

void RunTasks() {
   int i, j, k, nrows, nc, pos, ir, ic;
   double x;
   double *avgs, *stdevs;
   vector<string> tkns;
   ifstream ifile[20];
   ofstream ofile;
   int ncols[20] = {};
   nrows = 0;
   tkns = Tokenize(g_params.matrix_names, "\r\n");
   ofile.open("descs.bin", ios::out | ios::binary);
   for (i = 0, nc = 0; i < tkns.size(); i++) {
      ifile[i].open(("mat" + tkns[i] + ".bin").c_str(), ios::in | ios::binary);
      ifile[i].read((char*)&nrows, sizeof(int));
      ifile[i].read((char*)&ncols[i], sizeof(int));
      nc += ncols[i];
   }
   avgs = new double[nc];
   stdevs = new double[nc];
   ofile.write((char*)&nrows, sizeof(int));
   ofile.write((char*)&nc, sizeof(int));
   for (j = 0, pos = 0; j < tkns.size(); j++) {
      if (!File2Matrix(("avgs" + tkns[j] + ".bin").c_str(), ir, ic, avgs+pos)) {
         cerr << "Couldn't open avgs" + tkns[j] + ".bin" << endl;
         goto end;
      }
      if (!File2Matrix(("stdevs" + tkns[j] + ".bin").c_str(), ir, ic, stdevs+pos)) {
         cerr << "Couldn't open stdevs" + tkns[j] + ".bin";
         goto end;
      }
      pos += ncols[j];
   }
   for (i = 0; i < nrows; i++) {
      for (j=0; j < tkns.size(); j++) {
         for (k = 0; k < ncols[j]; k++) {
            ifile[j].read((char*)&x, sizeof(double));
            ofile.write((char*)&x, sizeof(double));
         }
      }
   }
   for (i = 0; i < tkns.size(); i++) {
      ifile[i].close();
   }
   ofile.close();
   Matrix2File(avgs, nc, 1, "avgs.bin");
   Matrix2File(stdevs, nc, 1, "stdevs.bin");
   cout << "Saved " << nrows << "x" << nc << " matrix" << endl;
end:
   delete [] stdevs;
   delete [] avgs;
}
