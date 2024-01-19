#pragma once

#include <fstream>

#include "Iterator.h"
#include "Mol.h"
#include "util.h"

using std::istream;

Iterator<Mol> ReadXYZFile(const char* fileName);
void ReadXYZFile(Iterator<Mol>& molecules, istream& InputStream);
int GetSubType(Triplet& t);
void analyze_triplet(Triplet& t, int na, int i, int j, int k, int pids[]);
void analyze_quadruplet(Quadruplet& q, int na, int i, int j, int k, int l, int pids[]);
int GetSubType(Quadruplet& q);
