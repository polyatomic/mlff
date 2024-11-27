#include <iostream>
#include <fstream>

#include "RunTasks.h"

using std::cerr;
using std::endl;
using std::ifstream;

PType g_params;

bool ReadParams(const char *fn) {
   size_t len, startpos;
   bool res, bt;
   string line, par, val;
   res = true;
   ifstream ifile(fn);
   if (!ifile) {
      res = false;
      goto end;
   }
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
      if (par == "MATRIX_NAMES") {
         g_params.matrix_names = val;
      }
      par = "";
   }
end:
   ifile.close();
   return res;
}

int main(int argc, char *argv[]) {
   int i;
   string pfile;
   if (argc < 3) {
      cerr << "join version 1.0.0.0" << endl;
      cerr << "Usage: " << argv[0] << " -p <parameter file>" << endl;
      return 1;
   }
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
   return 0;
}
