#pragma once

#include <string>
#include "chemutil.h"

using std::string;

struct PTypeExtended : public PType {
   string linearCoeffs;
   double bias;
};

void RunTasks();
