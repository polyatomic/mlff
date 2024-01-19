#pragma once

#include <string>

using std::string;

struct PType {
   string molFile;
   string t1;
   string t2;
   string t3;
   string t4;
   string ngp;
   string dist;
   string linearCoeffs;
   string eFile;
   int nthreads;
   int kernel_2b;
   int kernel_3b;
   int kernel_4b;
   double bias;
};

void RunTasks();
