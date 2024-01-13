#pragma once

#include <fstream>

#include "Iterator.h"
#include "Mol.h"

using std::istream;

Iterator<Mol> ReadXYZFile(const char* fileName);
void ReadXYZFile(Iterator<Mol>& molecules, istream& InputStream);
