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

PTypeExtended g_params;

bool ReadParams(const char *fn) {
   int ival;
   double dval;
   size_t len, startpos;
   bool res, bt;
   string line, par, val;
   res = true;
   ifstream ifile(fn);
   if (!ifile) {
      res = false;
      goto end;
   }
   g_params.preselect_descriptors = false;
   g_params.read_descriptors = false;
   g_params.nthreads = 1;
   g_params.maxrows = 1000;
   g_params.build_regression = false;
   g_params.order_by_leverages = false;
   g_params.tss = 1000;
   g_params.kernel_2b = 0;
   g_params.kernel_3b = 0;
   g_params.kernel_4b = 0;
   g_params.nsketch = 0;
   g_params.uis = false;
   g_params.gss = 0.0;
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
      if (par == "GEOMETRIES") {
         g_params.molFile = val;
      } else if (par == "ENERGIES") {
         g_params.eFile = val;
      } else if (par == "ATOM_TYPES") {
         g_params.t1 = val;
      } else if (par == "ATOM_TYPE_PAIRS") {
         g_params.t2 = val;
      } else if (par == "ATOM_TYPE_TRIPLETS") {
         g_params.t3 = val;
      } else if (par == "ATOM_TYPE_QUADRUPLETS") {
         g_params.t4 = val;
      } else if (par == "TRAINING_SET_SIZE") {
         ival = atoi(val.c_str());
         g_params.tss = ival;
      } else if (par == "DISTANCES") {
         g_params.dist = val;
      } else if (par == "NUMBER_OF_GRID_POINTS") {
         g_params.ngp = val;
      } else if (par == "GRID_START_POINTS") {
         g_params.gsp = val;
      } else if (par == "GRID_STEP_SIZE") {
         dval =  atof(val.c_str());
         g_params.gss = dval;
      } else if (par == "PRESELECT_DESCRIPTORS") {
         ival = atoi(val.c_str());
         if (ival == 1) g_params.preselect_descriptors = true;
      } else if (par == "READ_DESCRIPTOR_VALUES") {
         ival = atoi(val.c_str());
         if (ival == 1) g_params.read_descriptors = true;
      } else if (par == "N_THREADS") {
         ival = atoi(val.c_str());
         g_params.nthreads = ival;
      } else if (par == "MAX_ROWS") {
         ival = atoi(val.c_str());
         g_params.maxrows = ival;
      } else if (par == "BUILD_REGRESSION") {
         ival = atoi(val.c_str());
         if (ival == 1) g_params.build_regression = true;
      } else if (par == "KERNEL_2B") {
         ival = atoi(val.c_str());
         g_params.kernel_2b = ival;
      } else if (par == "KERNEL_3B") {
         ival = atoi(val.c_str());
         g_params.kernel_3b = ival;
      } else if (par == "KERNEL_4B") {
         ival = atoi(val.c_str());
         g_params.kernel_4b = ival;
      } else if (par == "N_SKETCHING_CASES") {
         ival = atoi(val.c_str());
         g_params.nsketch = ival;
      } else if (par == "ORDER_BY_LEVERAGES") {
         ival = atoi(val.c_str());
         if (ival == 1) g_params.order_by_leverages = true;
      } else if (par == "USE_INVERSE_SPACE") {
         ival = atoi(val.c_str());
         if (ival == 1) g_params.uis = true;
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
      cerr << "bokfit version 1.0.0.1" << endl;
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
