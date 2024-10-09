#pragma once

#include <string>
#include "chemutil.h"

using std::string;

struct PTypeExtended : public PType {
   int tss;
   int maxrows;
   int nsketch;
   bool preselect_descriptors;
   bool read_descriptors;
   bool build_regression;
   bool order_by_leverages;
   bool uis;
   bool ues;
};

void RunTasks();
