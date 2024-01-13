#pragma once

#include <string>

using std::string;

struct PType {
   string molFile;
   string eFile;
   string t1;
   string t2;
   string t3;
   string t4;
   string dist;
   string ngp;
   int rbs_2b;
   int rbs_3b;
   int rbs_4b;
   int nthreads;
   int maxrows;
   int kernel_2b;
   int kernel_3b;
   int kernel_4b;
   int nsketch = 0;
   bool preselect_descriptors;
   bool read_descriptors;
   bool build_regression;
   bool order_by_leverages;
};

void RunTasks();
