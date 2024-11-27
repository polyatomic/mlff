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
   int i, j, nsub;
   int **vars_selected;
   int *ranks, *ids_selected, *ip;
   ofstream ofile;
   vars_selected = 0;
   auto names = Tokenize(g_params.matrices, "\r\n");
   nsub = names.size();
   vars_selected = new int*[nsub];
   for (i = 0; i < nsub; i++) {
      vars_selected[i] = 0;
   }
   ranks = new int[nsub];
   auto rs = Tokenize(g_params.ranks, "\r\n");
   for (i = 0; i < nsub; i++) {
      ranks[i] = atoi(rs[i].c_str());
   }
   for (i=0; i < nsub; i++) {
      vars_selected[i] = new int[ranks[i]];
   }
   ids_selected = new int[g_params.nsketch];
   if (File2Array("sorted_ids.txt", ids_selected)) {
      sketch_matrices(g_params.nrows, names, ids_selected, ranks, g_params.nsketch, g_params.nthreads, vars_selected, g_params.maxcols, g_params.select_from_leading_matrix);
      for (i = 0; i < nsub; i++) {
         ip = vars_selected[i];
         ofile.open(names[i] + "sk.txt");
         for (j = 0; j < ranks[i]; j++) {
            ofile << ip[j] << endl;
         }
         ofile.close();
      }
   }
   delete [] ids_selected;
   for (i = 0; i < nsub; i++) {
      delete [] vars_selected[i];
   }
   delete [] ranks;
   delete [] vars_selected;
}
