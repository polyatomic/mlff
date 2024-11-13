#pragma once

#include <string>
#include "chemutil.h"

using std::string;

struct PTypeExtended : public PType {
   bool uis;
};

void RunTasks();
