#pragma once

using std::string;
using std::vector;

struct Pair {
   int v[2];
};

struct PairCompare {
   bool operator()(const Pair& x, const Pair& y) const {
      if (x.v[0] == y.v[0]) {
         return x.v[1] < y.v[1];
      }
      return x.v[0] < y.v[0];
   }
};

struct Triplet {
   int v[3];

   bool operator==(const Triplet& t) {
      return (t.v[0] == v[0] && t.v[1] == v[1] && t.v[2] == v[2]);
   }
};

struct TripletCompare {
   bool operator()(const Triplet& x, const Triplet& y) const {
      if (x.v[0] == y.v[0]) {
         if (x.v[1] == y.v[1])
            return x.v[2] < y.v[2];
         else
            return x.v[1] < y.v[1];
      }
      return x.v[0] < y.v[0];
   }
};

struct Quadruplet {
   int v[4];

   bool operator==(const Quadruplet& q) {
      return (q.v[0] == v[0] && q.v[1] == v[1] && q.v[2] == v[2] && q.v[3] == v[3]);
   }
};

struct Sixtuplet {
   int v[6];
};

struct SixtupletCompare {
   bool operator()(const Sixtuplet& x, const Sixtuplet& y) const {
      if (x.v[0] == y.v[0]) {
         if (x.v[1] == y.v[1]) {
            if (x.v[2] == y.v[2]) {
               if (x.v[3] == y.v[3]) {
                  if (x.v[4] == y.v[4])
                     return x.v[5] < y.v[5];
                  else
                     return x.v[4] < y.v[4];
               } else
                  return x.v[3] < y.v[3];
            } else
               return x.v[2] < y.v[2];
         } else
            return x.v[1] < y.v[1];
      }
      return x.v[0] < y.v[0];
   }
};

bool File2Array(string& fn, double arr[]);
bool File2Array(const string& fn, int arr[]);
vector<string> Tokenize(const string& str, const string& delimiters);
void Matrix2File(double *X, int nrows, int ncols, const char *fn);
void File2Matrix(const char *fn, int& nrows, int& ncols, double *X);
