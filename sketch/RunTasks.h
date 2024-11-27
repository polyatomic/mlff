#pragma once

#include <string>

using std::string;

struct PType {
   int nrows;
   int nsketch;
   int nthreads;
   int maxcols;
   bool select_from_leading_matrix;
   string matrices;
   string ranks;
};

void RunTasks();
