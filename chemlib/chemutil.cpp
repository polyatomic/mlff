#include "util.h"
#include "mathutil.h"

int GetSubType(Triplet& t) {
   if (t.v[0] == t.v[1]) {
      if (t.v[1] == t.v[2]) {
         return 0;
      }
   }
   return 1;
}

void analyze_triplet(Triplet& t, int na, int i, int j, int k, int pids[]) {
   int itmp;
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
      } else {
         // ABA -> jik
         itmp = t.v[0];
         t.v[0] = t.v[1];
         t.v[1] = itmp;
         pids[0] = IDX(na, i, j);
         pids[1] = IDX(na, j, k);
         pids[2] = IDX(na, i, k);
      }
   }
}

void analyze_quadruplet(Quadruplet& q, int na, int i, int j, int k, int l, int pids[]) {
   int itmp;
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
      } else {
         if (q.v[2] == q.v[3]) {
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
         } else {
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
   return 1;
}
