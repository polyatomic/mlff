#pragma once

#include <string>

using std::string;

struct PType {
   int nrows;
   string ncols;
   int maxrows;
   int nthreads;
   string matrices;
};

void RunTasks();
