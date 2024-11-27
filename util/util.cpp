#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "util.h"

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;
using std::ios;

bool File2Array(const string& fn, double arr[]) {
   int i;
   size_t len;
   double fval;
   string line;
   ifstream efile(fn.c_str());
   i = 0;
   if (efile) {
      for (;;) {
         if (efile.eof()) break;
         getline(efile, line, '\n');
         len = line.size();
         if (len == 0) continue;
         fval = atof(line.c_str());
         arr[i++] = fval;
      }
      efile.close();
   } else {
      cerr << "Error opening " << fn << endl;
      return false;
   }
   return true;
}

bool File2Array(const string& fn, int arr[]) {
   int i, ival;
   size_t len;
   string line;
   ifstream efile(fn.c_str());
   i = 0;
   if (efile) {
      for (;;) {
         if (efile.eof()) break;
         getline(efile, line, '\n');
         len = line.size();
         if (len == 0) continue;
         ival = atoi(line.c_str());
         arr[i++] = ival;
      }
      efile.close();
   } else {
      cerr << "Error opening " << fn << endl;
      return false;
   }
   return true;
}

vector<string> Tokenize(const string& str, const string& delimiters) {
   vector<string> tokens;
   string::size_type lastPos = str.find_first_not_of(delimiters, 0);
   string::size_type pos = str.find_first_of(delimiters, lastPos);
   while (string::npos != pos || string::npos != lastPos) {
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      lastPos = str.find_first_not_of(delimiters, pos);
      pos = str.find_first_of(delimiters, lastPos);
   }
   return tokens;
}

void Matrix2File(double *X, int nrows, int ncols, const char *fn) {
   long long i;
   ofstream ofile(fn, ios::out | ios::binary);
   ofile.write((char*)&nrows, sizeof(int));
   ofile.write((char*)&ncols, sizeof(int));
   for (i = 0; i < ((long long)nrows)*ncols; i++) ofile.write((char*)&X[i], sizeof(double));
   ofile.close();
}

bool File2Matrix(const char *fn, int& nrows, int& ncols, double *X) {
   long long i;
   ifstream ifile(fn, ios::in | ios::binary);
   if (!ifile) return false;
   ifile.read((char*)&nrows, sizeof(int));
   ifile.read((char*)&ncols, sizeof(int));
   for (i = 0; i < ((long long)nrows)*ncols; i++) ifile.read((char*)&X[i], sizeof(double));
   ifile.close();
   return true;
}

void copy_file(const char *srce_file, const char *dest_file) {
   ifstream srce(srce_file, ios::binary);
   ofstream dest(dest_file, ios::binary);
   dest << srce.rdbuf();
}

void divide(int *tsizes, int total, int nthreads) {
   int i, rem, ts;
   ts = total / nthreads;
   rem = total % nthreads;
   for (i = 0; i < rem; i++) tsizes[i] = ts + 1;
   for (i = rem; i < nthreads; i++) tsizes[i] = ts;
}
