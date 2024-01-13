#include <sstream>
//#include <iostream>

#include "Iterator.h"
#include "chemutil.h"

using std::ifstream;
using std::istringstream;
using std::ios;

Iterator<Mol> ReadXYZFile(const char *fileName) {
   Iterator<Mol> ret(0);
   ifstream ifile(fileName);
   if (!ifile) {
      goto end;
   }
   ReadXYZFile(ret, ifile);
end:
   ifile.close();
   return ret;
}

void ReadXYZFile(Iterator<Mol>& molecules, istream& InputStream) {
   int i, j, na, nmols;
   string line, ele;
   double coords[3];
   istringstream isstream;
   getline(InputStream, line, '\n');
   na = atoi(line.c_str());
   for (nmols=0;;) {
      if (InputStream.eof()) break;
      getline(InputStream, line, '\n');
      nmols++;
   }
   nmols /= (na+2);
   InputStream.clear();
   InputStream.seekg(0, ios::beg);
   molecules.Resize(nmols);
   for (j=0; j < nmols; j++) {
      getline(InputStream, line, '\n');
      getline(InputStream, line, '\n');
      //std::cout << line << '\n';
      molecules[j].Init(na, 0);
      for (i=0; i < na; i++) {
         getline(InputStream, line, '\n');
         isstream.str(line);
         isstream >> ele >> coords[0] >> coords[1] >> coords[2];
         isstream.clear();
         molecules[j].SetCoords(i, coords);
         molecules[j].SetAtom(i, ele);
      }
   }
}
