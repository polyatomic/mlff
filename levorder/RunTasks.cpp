#include <iostream>
#include <fstream>

#include "RunTasks.h"
#include "tasks.h"
#include "util.h"

extern PType g_params;

using std::cout;
using std::endl;
using std::ofstream;

void RunTasks() {
   int i, j, rank;
   int *ids_selected;
   double *levs;
   ids_selected = 0;
   levs = 0;
   ids_selected = new int[g_params.nrows];
   levs = new double[g_params.nrows];
   vector<string> tkns, tkns2;
   ofstream ofile;
   tkns = Tokenize(g_params.matrices, "\r\n");
   tkns2 = Tokenize(g_params.ncols, "\r\n");
   for (i=0; i < tkns.size(); i++) {
      if (!order_by_leverages(g_params.nrows, atoi(tkns2[i].c_str()), ("mat" + tkns[i] + ".bin").c_str(), ("mat" + tkns[i] + "t.bin").c_str(), ids_selected, levs, rank, g_params.maxrows, g_params.nthreads)) break;
      ofile.open((tkns[i] + "s.txt").c_str());
      for (j = 0; j < g_params.nrows; j++) ofile << ids_selected[j] << endl;
      ofile.close();
   }
   delete [] levs;
   delete [] ids_selected;
}
