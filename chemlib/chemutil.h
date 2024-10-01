#pragma once

#include <fstream>

#include "Iterator.h"
#include "Mol.h"
#include "util.h"
#include "Descriptors.h"

using std::istream;

struct PType {
   string molFile;
   string t1;
   string t2;
   string t3;
   string t4;
   string ngp;
   string gsp;
   string eFile;
   string dist;
   int kernel_2b;
   int kernel_3b;
   int kernel_4b;
   int nthreads;
   double gss;
};

Iterator<Mol> ReadXYZFile(const char* fileName);
void ReadXYZFile(Iterator<Mol>& molecules, istream& InputStream);
int GetSubType(Triplet& t);
void analyze_triplet(Triplet& t, int na, int i, int j, int k, int pids[]);
void analyze_quadruplet(Quadruplet& q, int na, int i, int j, int k, int l, int pids[]);
int GetSubType(Quadruplet& q);
bool calculation_prepare(PType *params, Iterator<Mol>& mols, Descriptors& descs, int& nstr, int& n, double*& v, double*& Y, int nblocks[], char **blocks);
