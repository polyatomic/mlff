#include <iostream>
#include <fstream>
#include <iomanip>

#ifdef _WIN32
   #include <windows.h>
#endif

#include "RunTasks.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;

PType g_params;

bool ReadParams(const char *fn) {
   int ival;
   size_t len, startpos;
   bool res, bt;
   string line, par, val;
   res = true;
   ifstream ifile(fn);
   if (!ifile) {
      res = false;
      goto end;
   }
   g_params.nrows = 0;
   g_params.nsketch = 0;
   g_params.nthreads = 1;
   g_params.maxcols = 0;
   g_params.select_from_leading_matrix = true;
   for (;;) {
      if (ifile.eof()) break;
      getline(ifile, line, '\n');
      len = line.size();
      if (len == 0) continue;
      if (line[0] == '#') continue;
      startpos = line.find("=");
      if (startpos != string::npos) {
         par = line.substr(0, startpos);
         while (par.size() > 0 && par[par.size()-1] == ' ') par.erase(par.size()-1, 1);
         startpos++;
         val = line.substr(startpos, len-startpos);
         while (val.size() > 0 && val[0] == ' ') val.erase(0, 1);
         bt = false;
         while (val.size() > 0 && val[val.size()-1] == '\\') {
            val.erase(val.size()-1, 1);
            bt = true;
         }
         if (bt) continue;
      } else if (line[len-1] == '\\') {
         while (line.size() > 0 && line[line.size()-1] == '\\') line.erase(line.size()-1, 1);
         val += '\n';
         val += line;
         continue;
      } else if (par.size() > 0) {
         val += '\n';
         val += line;
      }
      if (par == "N_ROWS") {
         ival = atoi(val.c_str());
         g_params.nrows = ival;
      } else if (par == "N_SKETCHING_CASES") {
         ival = atoi(val.c_str());
         g_params.nsketch = ival;
      } else if (par == "N_THREADS") {
         ival = atoi(val.c_str());
         g_params.nthreads = ival;
      } else if (par == "MATRIX_NAMES") {
         g_params.matrices = val;
      } else if (par == "RANKS") {
         g_params.ranks = val;
      } else if (par == "MAX_COLS") {
         ival = atoi(val.c_str());
         g_params.maxcols = ival;
      } else if (par == "SELECT_FROM_LEADING_MATRIX") {
         ival = atoi(val.c_str());
         if (ival == 0)
            g_params.select_from_leading_matrix = false;
         else
            g_params.select_from_leading_matrix = true;
      }
      par = "";
   }
end:
   ifile.close();
   return res;
}

int main(int argc, char *argv[]) {
   int i;
   double sUser = 0;
#ifdef _WIN32
   HANDLE hp;
   FILETIME ftimec, ftimee, fsys, fuser;
   __int64 i64User, i64UserPrev;
#endif
   string pfile;
   if (argc < 3) {
      cerr << "sketch version 1.0.0.1" << endl;
      cerr << "Usage: " << argv[0] << " -p <parameter file>" << endl;
      return 1;
   }
#ifdef _WIN32
   hp = GetCurrentProcess();
   i64UserPrev = 0;
#endif
   for (i = 1; i < 3; i++) {
      if (string(argv[i]) == "-p") {
         if (i+1 == argc) {
            cerr << "Errors in command line" << endl;
            return 1;
         }
         pfile = argv[i+1];
      }
   }
   const char *pfname = pfile.c_str();
   if (!ReadParams(pfname)) {
      cerr << "Error reading parameter file" << endl;
      return 1;
   }
   RunTasks();
#ifdef _WIN32
   GetProcessTimes(hp, &ftimec, &ftimee, &fsys, &fuser);
   i64User = *((__int64*)&fuser);
   sUser = ((double)(i64User - i64UserPrev)) / 10000000;
#endif
   cout << std::setiosflags(std::ios::fixed) << std::setprecision(1);
   cout << "CPU time: " << sUser << " seconds" << endl;
   return 0;
}
