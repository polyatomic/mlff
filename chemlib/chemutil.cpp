#include <iostream>
#include <iomanip>

#include "chemutil.h"

#include "mathutil.h"

using std::numeric_limits;
using std::streamsize;
using std::setprecision;
using std::cout;
using std::cerr;
using std::endl;

int GetSubType(Triplet& t) {
   if (t.v[0] == t.v[1]) {
      if (t.v[1] == t.v[2]) {
         return 0;
      }
   } else if (t.v[1] != t.v[2])
      return 2;
   return 1;
}

void analyze_triplet(Triplet& t, int na, int i, int j, int k, int pids[]) {
   int itmp;
   int trp[3];
   if (t.v[0] == t.v[1]) {
      if (t.v[1] == t.v[2]) {
         // AAA
         pids[0] = IDX(na, i, j);
         pids[1] = IDX(na, i, k);
         pids[2] = IDX(na, j, k);
      } else {
         // AAB -> kij
         itmp = t.v[0];
         t.v[0] = t.v[2];
         t.v[2] = itmp;
         pids[0] = IDX(na, i, k);
         pids[1] = IDX(na, j, k);
         pids[2] = IDX(na, i, j);
      }
   } else {
      if (t.v[1] == t.v[2]) {
         // ABB
         pids[0] = IDX(na, i, j);
         pids[1] = IDX(na, i, k);
         pids[2] = IDX(na, j, k);
      } else if (t.v[0] == t.v[2]) {
         // ABA -> jik
         itmp = t.v[0];
         t.v[0] = t.v[1];
         t.v[1] = itmp;
         pids[0] = IDX(na, i, j);
         pids[1] = IDX(na, j, k);
         pids[2] = IDX(na, i, k);
      } else {
         if (t.v[0] < t.v[1]) {
            if (t.v[1] < t.v[2]) {
               // ABC
               pids[0] = IDX(na, i, j);
               pids[1] = IDX(na, i, k);
               pids[2] = IDX(na, j, k);
            } else if (t.v[0] < t.v[2]) {
               // ACB -> ikj
               itmp = t.v[1];
               t.v[1] = t.v[2];
               t.v[2] = itmp;
               pids[0] = IDX(na, i, k);
               pids[1] = IDX(na, i, j);
               pids[2] = IDX(na, j, k);
            } else {
               // BCA -> kij
               trp[0] = t.v[0];
               trp[1] = t.v[1];
               trp[2] = t.v[2];
               t.v[0] = trp[2];
               t.v[1] = trp[0];
               t.v[2] = trp[1];
               pids[0] = IDX(na, i, k);
               pids[1] = IDX(na, j, k);
               pids[2] = IDX(na, i, j);
            }
         } else {
            if (t.v[0] > t.v[1] && t.v[0] > t.v[2]) {
               if (t.v[1] < t.v[2]) {
                  // CAB -> jki
                  trp[0] = t.v[0];
                  trp[1] = t.v[1];
                  trp[2] = t.v[2];
                  t.v[0] = trp[1];
                  t.v[1] = trp[2];
                  t.v[2] = trp[0];
                  pids[0] = IDX(na, j, k);
                  pids[1] = IDX(na, i, j);
                  pids[2] = IDX(na, i, k);
               } else {
                  // CBA -> kji
                  itmp = t.v[0];
                  t.v[0] = t.v[2];
                  t.v[2] = itmp;
                  pids[0] = IDX(na, j, k);
                  pids[1] = IDX(na, i, k);
                  pids[2] = IDX(na, i, j);
               }
            } else {
               // BAC -> jik
               itmp = t.v[0];
               t.v[0] = t.v[1];
               t.v[1] = itmp;
               pids[0] = IDX(na, i, j);
               pids[1] = IDX(na, j, k);
               pids[2] = IDX(na, i, k);
            }
         }
      }
   }
}

void analyze_quadruplet(Quadruplet& q, int na, int i, int j, int k, int l, int pids[]) {
   int itmp;
   int qrp[4];
   if (q.v[0] == q.v[1]) {
      if (q.v[1] == q.v[2]) {
         if (q.v[2] == q.v[3]) {
            // AAAA
            pids[0] = IDX(na, i, j);
            pids[1] = IDX(na, i, k);
            pids[2] = IDX(na, i, l);
            pids[3] = IDX(na, j, k);
            pids[4] = IDX(na, j, l);
            pids[5] = IDX(na, k, l);
         } else {
            // AAAB -> lijk
            itmp = q.v[0];
            q.v[0] = q.v[3];
            q.v[3] = itmp;
            pids[0] = IDX(na, i, l);
            pids[1] = IDX(na, j, l);
            pids[2] = IDX(na, k, l);
            pids[3] = IDX(na, i, j);
            pids[4] = IDX(na, i, k);
            pids[5] = IDX(na, j, k);
         }
      } else {
         if (q.v[2] == q.v[3]) {
            if (q.v[1] > q.v[2]) {
               // BBAA -> klij
               itmp = q.v[0];
               q.v[0] = q.v[2];
               q.v[2] = itmp;
               itmp = q.v[1];
               q.v[1] = q.v[3];
               q.v[3] = itmp;
               pids[0] = IDX(na, k, l);
               pids[1] = IDX(na, i, k);
               pids[2] = IDX(na, j, k);
               pids[3] = IDX(na, i, l);
               pids[4] = IDX(na, j, l);
               pids[5] = IDX(na, i, j);
            } else {
               // AABB
               pids[0] = IDX(na, i, j);
               pids[1] = IDX(na, i, k);
               pids[2] = IDX(na, i, l);
               pids[3] = IDX(na, j, k);
               pids[4] = IDX(na, j, l);
               pids[5] = IDX(na, k, l);
            }
         } else {
            if (q.v[0] != q.v[3]) {
               if (q.v[2] < q.v[3]) {
                  // CCAB -> klij
                  qrp[0] = q.v[0];
                  qrp[1] = q.v[1];
                  qrp[2] = q.v[2];
                  qrp[3] = q.v[3];
                  q.v[0] = qrp[2];
                  q.v[1] = qrp[3];
                  q.v[2] = qrp[0];
                  q.v[3] = qrp[1];
                  pids[0] = IDX(na, k, l);
                  pids[1] = IDX(na, i, k);
                  pids[2] = IDX(na, j, k);
                  pids[3] = IDX(na, i, l);
                  pids[4] = IDX(na, j, l);
                  pids[5] = IDX(na, i, j);
               } else {
                  // CCBA -> lkij
                  qrp[0] = q.v[0];
                  qrp[1] = q.v[1];
                  qrp[2] = q.v[2];
                  qrp[3] = q.v[3];
                  q.v[0] = qrp[3];
                  q.v[1] = qrp[2];
                  q.v[2] = qrp[0];
                  q.v[3] = qrp[1];
                  pids[0] = IDX(na, k, l);
                  pids[1] = IDX(na, i, l);
                  pids[2] = IDX(na, j, l);
                  pids[3] = IDX(na, i, k);
                  pids[4] = IDX(na, j, k);
                  pids[5] = IDX(na, i, j);
               }
            } else {
               // AABA -> kjil
               itmp = q.v[0];
               q.v[0] = q.v[2];
               q.v[2] = itmp;
               pids[0] = IDX(na, j, k);
               pids[1] = IDX(na, i, k);
               pids[2] = IDX(na, k, l);
               pids[3] = IDX(na, i, j);
               pids[4] = IDX(na, j, l);
               pids[5] = IDX(na, i, l);
            }
         }
      }
   } else {
      if (q.v[1] == q.v[2]) {
         if (q.v[2] == q.v[3]) {
            // ABBB
            pids[0] = IDX(na, i, j);
            pids[1] = IDX(na, i, k);
            pids[2] = IDX(na, i, l);
            pids[3] = IDX(na, j, k);
            pids[4] = IDX(na, j, l);
            pids[5] = IDX(na, k, l);
         } else {
            if (q.v[0] != q.v[3]) {
               if (q.v[0] < q.v[3]) {
                  // ACCB -> iljk
                  itmp = q.v[1];
                  q.v[1] = q.v[3];
                  q.v[3] = itmp;
                  pids[0] = IDX(na, i, l);
                  pids[1] = IDX(na, i, j);
                  pids[2] = IDX(na, i, k);
                  pids[3] = IDX(na, j, l);
                  pids[4] = IDX(na, k, l);
                  pids[5] = IDX(na, j, k);
               } else {
                  // BCCA -> lijk
                  qrp[0] = q.v[0];
                  qrp[1] = q.v[1];
                  qrp[3] = q.v[3];
                  q.v[0] = qrp[3];
                  q.v[1] = qrp[0];
                  q.v[3] = qrp[1];
                  pids[0] = IDX(na, i, l);
                  pids[1] = IDX(na, j, l);
                  pids[2] = IDX(na, k, l);
                  pids[3] = IDX(na, i, j);
                  pids[4] = IDX(na, i, k);
                  pids[5] = IDX(na, j, k);
               }
            } else {
               if (q.v[0] < q.v[1]) {
                  // ABBA -> iljk 
                  itmp = q.v[1];
                  q.v[1] = q.v[3];
                  q.v[3] = itmp;
                  pids[0] = IDX(na, i, l);
                  pids[1] = IDX(na, i, j);
                  pids[2] = IDX(na, i, k);
                  pids[3] = IDX(na, j, l);
                  pids[4] = IDX(na, k, l);
                  pids[5] = IDX(na, j, k);
               } else {
                  // BAAB -> kjil
                  itmp = q.v[0];
                  q.v[0] = q.v[2];
                  q.v[2] = itmp;
                  pids[0] = IDX(na, j, k);
                  pids[1] = IDX(na, i, k);
                  pids[2] = IDX(na, k, l);
                  pids[3] = IDX(na, i, j);
                  pids[4] = IDX(na, j, l);
                  pids[5] = IDX(na, i, l);
               }
            }
         }
      } else {
         if (q.v[2] == q.v[3]) {
            if (q.v[0] != q.v[2]) {
               if (q.v[0] < q.v[1]) {
                  // ABCC
                  pids[0] = IDX(na, i, j);
                  pids[1] = IDX(na, i, k);
                  pids[2] = IDX(na, i, l);
                  pids[3] = IDX(na, j, k);
                  pids[4] = IDX(na, j, l);
                  pids[5] = IDX(na, k, l);
               } else {
                  // BACC -> jikl
                  itmp = q.v[1];
                  q.v[1] = q.v[0];
                  q.v[0] = itmp;
                  pids[0] = IDX(na, i, j);
                  pids[1] = IDX(na, j, k);
                  pids[2] = IDX(na, j, l);
                  pids[3] = IDX(na, i, k);
                  pids[4] = IDX(na, i, l);
                  pids[5] = IDX(na, k, l);
               }
            } else {
               // ABAA -> jikl
               itmp = q.v[1];
               q.v[1] = q.v[0];
               q.v[0] = itmp;
               pids[0] = IDX(na, i, j);
               pids[1] = IDX(na, j, k);
               pids[2] = IDX(na, j, l);
               pids[3] = IDX(na, i, k);
               pids[4] = IDX(na, i, l);
               pids[5] = IDX(na, k, l);
            }
         } else {
            if (q.v[0] != q.v[2]) {
               if (q.v[0] != q.v[3]) {
                  if (q.v[0] < q.v[2]) {
                     // ACBC -> ikjl
                     itmp = q.v[1];
                     q.v[1] = q.v[2];
                     q.v[2] = itmp;
                     pids[0] = IDX(na, i, k);
                     pids[1] = IDX(na, i, j);
                     pids[2] = IDX(na, i, l);
                     pids[3] = IDX(na, j, k);
                     pids[4] = IDX(na, k, l);
                     pids[5] = IDX(na, j, l);
                  } else {
                     // BCAC -> kijl
                     qrp[0] = q.v[0];
                     qrp[1] = q.v[1];
                     qrp[2] = q.v[2];
                     q.v[0] = qrp[2];
                     q.v[1] = qrp[0];
                     q.v[2] = qrp[1];
                     pids[0] = IDX(na, i, k);
                     pids[1] = IDX(na, j, k);
                     pids[2] = IDX(na, k, l);
                     pids[3] = IDX(na, i, j);
                     pids[4] = IDX(na, i, l);
                     pids[5] = IDX(na, j, l);
                  }
               } else {
                  if (q.v[1] < q.v[2]) {
                     // CABC -> jkil
                     qrp[0] = q.v[0];
                     qrp[1] = q.v[1];
                     qrp[2] = q.v[2];
                     q.v[0] = qrp[1];
                     q.v[1] = qrp[2];
                     q.v[2] = qrp[0];
                     pids[0] = IDX(na, j, k);
                     pids[1] = IDX(na, i, j);
                     pids[2] = IDX(na, j, l);
                     pids[3] = IDX(na, i, k);
                     pids[4] = IDX(na, k, l);
                     pids[5] = IDX(na, i, l);
                  } else {
                     // CBAC -> kjil
                     itmp = q.v[0];
                     q.v[0] = q.v[2];
                     q.v[2] = itmp;
                     pids[0] = IDX(na, j, k);
                     pids[1] = IDX(na, i, k);
                     pids[2] = IDX(na, k, l);
                     pids[3] = IDX(na, i, j);
                     pids[4] = IDX(na, j, l);
                     pids[5] = IDX(na, i, l);
                  }
               }
            } else {
               if (q.v[1] == q.v[3]) {
                  if (q.v[1] < q.v[2]) {
                     // BABA -> ljki
                     itmp = q.v[0];
                     q.v[0] = q.v[3];
                     q.v[3] = itmp;
                     pids[0] = IDX(na, j, l);
                     pids[1] = IDX(na, k, l);
                     pids[2] = IDX(na, i, l);
                     pids[3] = IDX(na, j, k);
                     pids[4] = IDX(na, i, j);
                     pids[5] = IDX(na, i, k);
                  } else {
                     // ABAB -> ikjl
                     itmp = q.v[1];
                     q.v[1] = q.v[2];
                     q.v[2] = itmp;
                     pids[0] = IDX(na, i, k);
                     pids[1] = IDX(na, i, j);
                     pids[2] = IDX(na, i, l);
                     pids[3] = IDX(na, j, k);
                     pids[4] = IDX(na, k, l);
                     pids[5] = IDX(na, j, l);
                  }
               } else {
                  if (q.v[1] < q.v[3]) {
                     // CACB -> jlik
                     qrp[0] = q.v[0];
                     qrp[1] = q.v[1];
                     qrp[3] = q.v[3];
                     q.v[0] = qrp[1];
                     q.v[1] = qrp[3];
                     q.v[3] = qrp[0];
                     pids[0] = IDX(na, j, l);
                     pids[1] = IDX(na, i, j);
                     pids[2] = IDX(na, j, k);
                     pids[3] = IDX(na, i, l);
                     pids[4] = IDX(na, k, l);
                     pids[5] = IDX(na, i, k);
                  } else {
                     // CBCA -> ljik
                     itmp = q.v[0];
                     q.v[0] = q.v[3];
                     q.v[3] = itmp;
                     pids[0] = IDX(na, j, l);
                     pids[1] = IDX(na, i, l);
                     pids[2] = IDX(na, k, l);
                     pids[3] = IDX(na, i, j);
                     pids[4] = IDX(na, j, k);
                     pids[5] = IDX(na, i, k);
                  }
               }
            }
         }
      }
   }
}

int GetSubType(Quadruplet& q) {
   int i, comp, count;
   comp = q.v[0];
   for (i = 1, count = 1; i < 4; i++) {
      if (q.v[i] == comp) count++;
   }
   if (count == 2) return 2;
   if (count == 4) return 0;
   if (q.v[1] != q.v[2]) return 3;
   return 1;
}

void GetPairTypes(Triplet& t3, vector<Pair>& t2, int *tps) {
   int i, j, k, l, tp;
   Pair p;
   for (i = 0, k = 0; i < 2; i++) {
      for (j = i + 1; j < 3; j++) {
         p.v[0] = t3.v[i];
         p.v[1] = t3.v[j];
         if (p.v[0] > p.v[1]) { l = p.v[0]; p.v[0] = p.v[1]; p.v[1] = l; }
         //auto t2pos = lower_bound(t2.begin(), t2.end(), p, PairCompare());
         auto t2pos = find(t2.begin(), t2.end(), p);
         tp = t2pos - t2.begin();
         tps[k++] = tp;
      }
   }
}

void GetPairTypes(Quadruplet& t4, vector<Pair>& t2, int *tps) {
   int i, j, k, l, tp;
   Pair p;
   for (i = 0, k = 0; i < 3; i++) {
      for (j = i + 1; j < 4; j++) {
         p.v[0] = t4.v[i];
         p.v[1] = t4.v[j];
         if (p.v[0] > p.v[1]) { l = p.v[0]; p.v[0] = p.v[1]; p.v[1] = l; }
         //auto t2pos = lower_bound(t2.begin(), t2.end(), p, PairCompare());
         auto t2pos = find(t2.begin(), t2.end(), p);
         tp = t2pos - t2.begin();
         tps[k++] = tp;
      }
   }
}

int GetSym(Quadruplet& q) {
   if (q.v[0] == q.v[1]) {
      if (q.v[2] == q.v[1]) {
         return 24;
      } else {
         return 4;
      }
   } else {
      if (q.v[2] == q.v[1]) {
         return 6;
      } else {
         return 2;
      }
   }
   return 0;
}

bool calculation_prepare(PType *params, Iterator<Mol>& mols, Descriptors& descs, int& nstr, int& n, double*& v, double*& Y, int nblocks[], char **blocks, double **gsp, int **p4sym) {
   double *dmin, *dmax, *rcoords;
   int *types2, *types3, *stypes3, *types4, *stypes4, *ngp, *ptypes3, *ptypes4;
   int i, j, k, l, m, o, na, nt1, nt2, nt3, tp, na3, na4, nt4, ii, jj, kk, nb_total;
   vector<string> tkns, tkns2;
   vector<int> t1;
   vector<Pair> t2;
   vector<Triplet> t3;
   vector<Quadruplet> t4;
   Pair p;
   Triplet t;
   Quadruplet q;
   vector<Pair>::iterator t2pos;
   vector<Triplet>::iterator t3pos;
   vector<Quadruplet>::iterator t4pos;
   bool res = true;
   int pids[6] = {};
   int prec = numeric_limits<double>::max_digits10;
   streamsize oldprec = cout.precision();
   string s = params->molFile;
   mols = ReadXYZFile(s.c_str());
   if (!mols) {
      cerr << "Error opening " << s << endl;
      return false;
   }
   nstr = mols.Size();
   cout << "Number of configurations = " << nstr << endl;
   na = mols[0].NAtoms();
   cout << "Number of atoms in configuration = " << na << endl;
   n = 3*na;
   j = nstr*n;
   v = new double[j];
   for (i = 0; i < nstr; i++) mols[i].GetCoordinates(v + i*n);
   mols.Resize(1);
   Y = new double[nstr];
   s = params->eFile;
   if (s.size() > 0) {
      if (!File2Array(s, Y)) {
         cerr << "Error reading " << s << endl;
         return false;
      }
   } else {
      cerr << "Input error" << endl;
      return false;
   }
   tkns = Tokenize(params->t1, "\r\n");
   nt1 = tkns.size();
   cout << "Number of atom types = " << nt1 << endl;
   cout << "Atom types:" << endl;
   for (i = 0; i < nt1; i++) {
      j = atoi(tkns[i].c_str());
      t1.push_back(j);
      cout << Atom::GetAtomicSymbol(j) << endl;
   }
   tkns = Tokenize(params->t2, "\r\n");
   nt2 = tkns.size();
   cout << "Number of atom type pairs = " << nt2 << endl;
   nblocks[0] = nt2;
   for (i = 0; i < nt2; i++) {
      tkns2 = Tokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      p.v[0] = j;
      j = atoi(tkns2[1].c_str());
      p.v[1] = j;
      t2.push_back(p);
   }
   dmin = new double[nt2];
   dmax = new double[nt2];
   types2 = new int[na*(na-1)/2];
   for (i = 0, m = 0; i < na-1; i++) {
      k = mols[0].GetAtomicNumber(i);
      for (j = i+1; j < na; j++) {
         p.v[0] = k;
         p.v[1] = mols[0].GetAtomicNumber(j);
         if (p.v[0] > p.v[1]) { l = p.v[0]; p.v[0] = p.v[1]; p.v[1] = l; }
         //t2pos = lower_bound(t2.begin(), t2.end(), p, PairCompare());
         t2pos = find(t2.begin(), t2.end(), p);
         tp = t2pos - t2.begin();
         types2[m++] = tp;
      }
   }
   tkns = Tokenize(params->t3, "\r\n");
   nt3 = tkns.size();
   cout << "Number of atom type triplets = " << nt3 << endl;
   nblocks[1] = nt3;
   for (i = 0; i < nt3; i++) {
      tkns2 = Tokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      t.v[0] = j;
      j = atoi(tkns2[1].c_str());
      t.v[1] = j;
      j = atoi(tkns2[2].c_str());
      t.v[2] = j;
      t3.push_back(t);
   }
   na3 = na*(na-1)*(na-2)/6;
   types3 = new int[4*na3];
   stypes3 = new int[nt3];
   ptypes3 = new int[nt3*3];
   for (i = 0; i < nt3; i++) {
      stypes3[i] = GetSubType(t3[i]);
      GetPairTypes(t3[i], t2, ptypes3+i*3);
   }
   for (i = 0, o = 0; i < na-2; i++) {
      l = mols[0].GetAtomicNumber(i);
      for (j = i+1; j < na-1; j++) {
         m = mols[0].GetAtomicNumber(j);
         for (k = j+1; k < na; k++) {
            t.v[0] = l;
            t.v[1] = m;
            t.v[2] = mols[0].GetAtomicNumber(k);
            analyze_triplet(t, na, i, j, k, pids);
            t3pos = find(t3.begin(), t3.end(), t);
            tp = t3pos - t3.begin();
            types3[o] = tp;
            types3[na3+3*o] = pids[0]; types3[na3+3*o+1] = pids[1]; types3[na3+3*o+2] = pids[2];
            o++;
         }
      }
   }
   tkns = Tokenize(params->t4, "\r\n");
   nt4 = tkns.size();
   cout << "Number of atom type quadruplets = " << nt4 << endl;
   nblocks[2] = nt4;
   for (i = 0; i < nt4; i++) {
      tkns2 = Tokenize(tkns[i], " ");
      j = atoi(tkns2[0].c_str());
      q.v[0] = j;
      j = atoi(tkns2[1].c_str());
      q.v[1] = j;
      j = atoi(tkns2[2].c_str());
      q.v[2] = j;
      j = atoi(tkns2[3].c_str());
      q.v[3] = j;
      t4.push_back(q);
   }
   na4 = na*(na-1)*(na-2)*(na-3)/24;
   types4 = new int[7*na4];
   stypes4 = new int[nt4];
   ptypes4 = new int[nt4*6];
   *p4sym = new int[nt4];
   for (i = 0; i < nt4; i++) {
      stypes4[i] = GetSubType(t4[i]);
      GetPairTypes(t4[i], t2, ptypes4+i*6);
      (*p4sym)[i] = GetSym(t4[i]);
   }
   for (i = 0, o = 0; i < na-3; i++) {
      ii = mols[0].GetAtomicNumber(i);
      for (j = i+1; j < na-2; j++) {
         jj = mols[0].GetAtomicNumber(j);
         for (k = j+1; k < na-1; k++) {
            kk = mols[0].GetAtomicNumber(k);
            for (l = k+1; l < na; l++) {
               q.v[0] = ii;
               q.v[1] = jj;
               q.v[2] = kk;
               q.v[3] = mols[0].GetAtomicNumber(l);
               analyze_quadruplet(q, na, i, j, k, l, pids);
               t4pos = find(t4.begin(), t4.end(), q);
               tp = t4pos - t4.begin();
               types4[o] = tp;
               types4[na4+6*o] = pids[0]; types4[na4+6*o+1] = pids[1]; types4[na4+6*o+2] = pids[2]; types4[na4+6*o+3] = pids[3]; types4[na4+6*o+4] = pids[4]; types4[na4+6*o+5] = pids[5];
               o++;
            }
         }
      }
   }
   ngp = new int[nt2];
   *gsp = nullptr;
   tkns = Tokenize(params->ngp, "\r\n");
   if (tkns.size() != nt2) {
      cerr << "Input error reading NUMBER_OF_GRID_POINTS" << endl;
      res = false;
      goto end;
   }
   for (i = 0; i < nt2; i++)
      ngp[i] = atoi(tkns[i].c_str());
   descs.Init(na, nt2, types2, nt3, stypes3, types3, ptypes3, nt4, stypes4, types4, ptypes4, ngp, params->nthreads);
   if (params->dist.size() == 0) {
      for (i = 0; i < nt2; i++) {
         dmin[i] = 1000.0;
         dmax[i] = 0.0;
      }
      for (i = 0; i < nstr; i++) {
         rcoords = v + i*n;
         descs.GetDistanceRanges(rcoords, dmin, dmax);
      }
   } else {
      tkns = Tokenize(params->dist, "\r\n");
      if (tkns.size() == nt2) {
         for (i = 0; i < nt2; i++) {
            tkns2 = Tokenize(tkns[i], " ");
            dmin[i] = atof(tkns2[0].c_str());
            dmax[i] = atof(tkns2[1].c_str());
         }
      } else {
         cerr << "Input error reading DISTANCES" << endl;
         res = false;
         goto end;
      }
   }
   *gsp = new double[nt2];
   tkns = Tokenize(params->gsp, "\r\n");
   if (tkns.size() != nt2) {
      cerr << "Input error reading GRID_START_POINTS" << endl;
      res = false;
      goto end;
   }
   for (i = 0; i < nt2; i++)
      (*gsp)[i] = atof(tkns[i].c_str());
   cout << "Atom type pairs, minimum and maximum distances, max-min and number of grid points:" << endl;
   for (i = 0; i < nt2; i++)
      cout << Atom::GetAtomicSymbol(t2[i].v[0]) << " " << Atom::GetAtomicSymbol(t2[i].v[1]) << " " << setprecision(prec) << dmin[i] << " " << dmax[i] << " " << setprecision(oldprec) << dmax[i] - dmin[i] << " " << ngp[i] << endl;
   cout << "Grid step size (inverse): " << params->gss << endl;
   nb_total = nt2 + nt3 + nt4;
   *blocks = new char[5*nb_total];
   for (i=0; i < nb_total; i++) {
      (*blocks)[5*i] = 'g';
      (*blocks)[5*i+4] = '\0';
   }
   for (i=0; i < nt2; i++) {
      (*blocks)[5*i+1] = '2';
      if (i < 10) {
         (*blocks)[5*i+2] = '0';
         (*blocks)[5*i+3] = '0' + i;
      } else {
         (*blocks)[5*i+2] = '1';
         (*blocks)[5*i+3] = '0' + (i - 10);
      }
   }
   for (i=nt2, j=0; i < nt2+nt3; i++, j++) {
      (*blocks)[5*i+1] = '3';
      if (j < 10) {
         (*blocks)[5*i+2] = '0';
         (*blocks)[5*i+3] = '0' + j;
      } else {
         (*blocks)[5*i+2] = '1';
         (*blocks)[5*i+3] = '0' + (j - 10);
      }
   }
   for (i=nt2+nt3, j=0; i < nb_total; i++, j++) {
      (*blocks)[5*i+1] = '4';
      if (j < 10) {
         (*blocks)[5*i+2] = '0';
         (*blocks)[5*i+3] = '0' + j;
      } else {
         (*blocks)[5*i+2] = '1';
         (*blocks)[5*i+3] = '0' + (j - 10);
      }
   }
   descs.SetGrid(dmin, dmax, *gsp, params->gss, *blocks, *p4sym);
   Descriptors::SetKernel2B(params->kernel_2b);
   Descriptors::SetKernel3B(params->kernel_3b);
   Descriptors::SetKernel4B(params->kernel_4b);
end:
   delete [] ngp;
   delete [] ptypes4;
   delete [] stypes4;
   delete [] types4;
   delete [] ptypes3;
   delete [] stypes3;
   delete [] types3;
   delete [] types2;
   delete [] dmax;
   delete [] dmin;
   return res;
}
