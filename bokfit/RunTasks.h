#pragma once

#include <string>
#include "chemutil.h"

using std::string;

struct PTypeExtended : public PType {
   int rbs_2b;
   int rbs_3b;
   int rbs_4b;
   int maxrows;
   int nsketch;
   bool preselect_descriptors;
   bool read_descriptors;
   bool build_regression;
   bool order_by_leverages;
};

void RunTasks();
