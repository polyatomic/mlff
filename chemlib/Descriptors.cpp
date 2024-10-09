#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "mathutil.h"
#include "util.h"

#include "Descriptors.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::string;
using std::sort;

Descriptors::Descriptors():
m_na(0),
m_nt2(0),
m_t2(0),
m_ngp2(0),
m_svp2(0),
m_ndesc(0),
m_na2(0),
m_dist(0),
m_t2b(0),
m_ngp2sum(0),
m_nt3(0),
m_na3(0),
m_t3(0),
m_st3(0),
m_t3b(0),
m_svp3(0),
m_svp3i(0),
m_ngp3(0),
m_nt4(0),
m_na4(0),
m_t4(0),
m_st4(0),
m_t4b(0),
m_ngp3sum(0),
m_svp4(0),
m_svp4i(0),
m_ngp4(0),
m_pt3(0),
m_pt4(0),
m_use_inverse_space(false) {
}

Descriptors::~Descriptors() {
   delete [] m_t2;
   delete [] m_ngp2;
   delete [] m_dist;
   delete [] m_t2b;
   delete [] m_t3;
   delete [] m_st3;
   delete [] m_t3b;
   delete [] m_t4;
   delete [] m_st4;
   delete [] m_t4b;
   delete [] m_pt3;
   delete [] m_pt4;
   ReleaseGrid();
}

bool Descriptors::Init(int na, int nt2, int *t2, int nt3, int *st3, int *t3, int *pt3, int nt4, int *st4, int *t4, int *pt4, int *ngp, int nthreads) {
   int i;
   ReleaseGrid();
   if (!t2) {
      return false;
   }
   m_na = na;
   m_na2 = na*(na-1)/2;
   m_nt2 = nt2;
   delete [] m_t2;
   m_t2 = new int[m_na2];
   for (i=0; i < m_na2; i++) m_t2[i] = t2[i];
   delete [] m_ngp2;
   m_ngp2 = new int[nt2];
   for (i=0; i < nt2; i++) m_ngp2[i] = ngp[i];
   delete [] m_dist;
   m_dist = new double[nthreads*m_na2];
   delete [] m_t2b;
   m_t2b = new int[nt2];
   m_nt3 = nt3;
   if (!nt3) return true;
   m_na3 = na*(na-1)*(na-2)/6;
   delete [] m_t3;
   m_t3 = new int[4*m_na3];
   for (i=0; i < 4*m_na3; i++) m_t3[i] = t3[i];
   delete [] m_st3;
   m_st3 = new int[nt3];
   for (i=0; i < nt3; i++) m_st3[i] = st3[i];
   delete [] m_pt3;
   m_pt3 = new int[3*nt3];
   for (i=0; i < 3*nt3; i++) m_pt3[i] = pt3[i];
   delete [] m_t3b;
   m_t3b = new int[nt3];
   m_nt4 = nt4;
   if (!nt4) return true;
   m_na4 = na*(na-1)*(na-2)*(na-3)/24;
   delete [] m_t4;
   m_t4 = new int[7*m_na4];
   for (i=0; i < 7*m_na4; i++) m_t4[i] = t4[i];
   delete [] m_st4;
   m_st4 = new int[nt4];
   for (i=0; i < nt4; i++) m_st4[i] = st4[i];
   delete [] m_pt4;
   m_pt4 = new int[6*nt4];
   for (i=0; i < 6*nt4; i++) m_pt4[i] = pt4[i];
   delete [] m_t4b;
   m_t4b = new int[nt4];
   return true;
}

void Descriptors::SetGrid(double *dmin, double *dmax, double *gmin, double step_size, const char *blocks, int *sym) {
   int i, j, dif;
   int *p;
   double b, stp;
   ReleaseGrid();
   m_svp2 = new double*[m_nt2];
   for (i=0; i < m_nt2; i++) m_svp2[i] = 0;
   stp = step_size;
   for (i=0; i < m_nt2; i++) {
      //b = 1.0/dmax[i];
      //e = 1.0/dmin[i];
      //stp = (e - b)/(m_ngp2[i] - 1);
      b = gmin[i];
      m_svp2[i] = new double[m_ngp2[i]];
      if (m_use_inverse_space) {
         for (j=0; j < m_ngp2[i]; j++)
            m_svp2[i][j] = b + j*stp;
      } else {
         for (j=0; j < m_ngp2[i]; j++)
            m_svp2[i][m_ngp2[i] - j - 1] = 1.0/(b + j*stp);
      }
   }
   for (i = 0, m_ndesc = 0; i < m_nt2; i++) {
      m_t2b[i] = m_ndesc;
      m_ndesc += m_ngp2[i];
   }
   m_ngp2sum = m_ndesc;
   if (m_nt3) {
      m_svp3 = new double*[m_nt3];
      for (i=0; i < m_nt3; i++) m_svp3[i] = 0;
      m_svp3i = new int*[m_nt3];
      for (i=0; i < m_nt3; i++) m_svp3i[i] = 0;
      m_ngp3 = new int[m_nt3];
      for (i = 0; i < m_nt3; i++) {
         p = m_pt3 + 3*i;
         if (!ReadGridPoints(p[0], p[1], p[2], i, (string(blocks+4*(m_nt2+i)) + ".txt").c_str())) AssignGridPoints(p[0], p[1], p[2], i);
      }
      for (i=0,j=0; i < m_nt3; i++) {
         m_t3b[i] = j;
         j += m_ngp3[i];
      }
      m_ndesc += j;
      m_ngp3sum = m_ndesc;
   }
   if (m_nt4) {
      m_svp4 = new double*[m_nt4];
      for (i=0; i < m_nt4; i++) m_svp4[i] = 0;
      m_svp4i = new int*[m_nt4];
      for (i=0; i < m_nt4; i++) m_svp4i[i] = 0;
      m_ngp4 = new int[m_nt4];
      for (i = 0; i < m_nt4; i++) {
         p = m_pt4 + 6*i;
         // TODO: Remove sym and induce it from m_st4
         if (!ReadGridPoints(p[0], p[1], p[2], p[3], p[4], p[5], i, (string(blocks+4*(m_nt2+m_nt3+i)) + ".txt").c_str())) AssignGridPoints(p[0], p[1], p[2], p[3], p[4], p[5], i, sym[i]);
      }
      for (i = 0, j = 0; i < m_nt4; i++) {
         m_t4b[i] = j;
         j += m_ngp4[i];
      }
      m_ndesc += j;
   }
   for (i = 0; i < m_nt2; i++) {
      ReadGridPoints(i, (string(blocks+4*i) + ".txt").c_str());
   }
   for (i = 0, j = 0; i < m_nt2; i++) {
      m_t2b[i] = j;
      j += m_ngp2[i];
   }
   dif = m_ngp2sum - j;
   m_ndesc -= dif;
   m_ngp2sum -= dif;
   m_ngp3sum -= dif;
}

int Descriptors::GetNDescriptors() {
   return m_ndesc;
}

int Descriptors::GetN2BDescriptors(int type) {
   if (type < 0) return m_ngp2sum;
   if (type == 0) return m_t2b[1];
   if (type == m_nt2-1) return m_ngp2sum - m_t2b[type];
   return m_t2b[type+1] - m_t2b[type];
}

int Descriptors::GetN3BDescriptors(int type) {
   if (!m_nt3) return 0;
   if (type < 0) return m_ngp3sum - m_ngp2sum;
   if (type == 0) return m_t3b[1];
   if (type == m_nt3-1) return m_ngp3sum - m_ngp2sum - m_t3b[type];
   return m_t3b[type+1] - m_t3b[type];
}

int Descriptors::GetN4BDescriptors(int type) {
   if (!m_nt4) return 0;
   if (type < 0) return m_ndesc - m_ngp3sum;
   if (type == 0) return m_t4b[1];
   if (type == m_nt4-1) return m_ndesc - m_ngp3sum - m_t4b[type];
   return m_t4b[type+1] - m_t4b[type];
}

void Descriptors::GetDistanceRanges(double *r, double *dmin, double *dmax) {
   int i, j, k, tp;
   double dist;
   double *rp;
   for (j = 0, i = 0; j < m_na-1; j++) {
      rp = r+3*j;
      for (k = j+1; k < m_na; k++) {
         tp = m_t2[i];
         dist = sqrt(sq_dist(rp, r+3*k, 3));
         if (dmin[tp] > dist) dmin[tp] = dist;
         if (dmax[tp] < dist) dmax[tp] = dist;
         i++;
      }
   }
}

void Descriptors::ReleaseGrid() {
   int i;
   if (!m_svp2) return;
   for (i=0; i < m_nt2; i++) delete [] m_svp2[i];
   delete [] m_svp2;
   m_svp2 = 0;
   if (!m_svp3) return;
   for (i=0; i < m_nt3; i++) {
      delete [] m_svp3[i];
      delete [] m_svp3i[i];
   }
   delete [] m_svp3;
   delete [] m_svp3i;
   delete [] m_ngp3;
   m_svp3 = 0;
   m_svp3i = 0;
   m_ngp3 = 0;
   if (!m_svp4) return;
   for (i = 0; i < m_nt4; i++) {
      delete [] m_svp4[i];
      delete [] m_svp4i[i];
   }
   delete [] m_svp4;
   delete [] m_svp4i;
   delete [] m_ngp4;
   m_svp4 = 0;
   m_svp4i = 0;
   m_ngp4 = 0;
}

void Descriptors::SetKernel2B(int kernel_type) {
   switch (kernel_type) {
      case RP_3_0:
         kernel2B = rp_3_0_kernel;
         cout << "Using 3,0 reciprocal power kernel for 2B terms" << endl;
         break;
      case RP_3_1:
         kernel2B = rp_3_1_kernel;
         cout << "Using 3,1 reciprocal power kernel for 2B terms" << endl;
         break;
      case RP_3_2:
         kernel2B = rp_3_2_kernel;
         cout << "Using 3,2 reciprocal power kernel for 2B terms" << endl;
         break;
      case RP_3_3:
         kernel2B = rp_3_3_kernel;
         cout << "Using 3,3 reciprocal power kernel for 2B terms" << endl;
         break;
      case RP_3_4:
         kernel2B = rp_3_4_kernel;
         cout << "Using 3,4 reciprocal power kernel for 2B terms" << endl;
         break;
      case RP_3_5:
         kernel2B = rp_3_5_kernel;
         cout << "Using 3,5 reciprocal power kernel for 2B terms" << endl;
         break;
      case RP_3_6:
         kernel2B = rp_3_6_kernel;
         cout << "Using 3,6 reciprocal power kernel for 2B terms" << endl;
         break;
   }
}

void Descriptors::SetKernel3B(int kernel_type) {
   switch (kernel_type) {
      case RP_3_0:
         kernel3B = rp_3_0_kernel;
         cout << "Using 3,0 reciprocal power kernel for 3B terms" << endl;
         break;
      case RP_3_1:
         kernel3B = rp_3_1_kernel;
         cout << "Using 3,1 reciprocal power kernel for 3B terms" << endl;
         break;
      case RP_3_2:
         kernel3B = rp_3_2_kernel;
         cout << "Using 3,2 reciprocal power kernel for 3B terms" << endl;
         break;
      case RP_3_3:
         kernel3B = rp_3_3_kernel;
         cout << "Using 3,3 reciprocal power kernel for 3B terms" << endl;
         break;
      case RP_3_4:
         kernel3B = rp_3_4_kernel;
         cout << "Using 3,4 reciprocal power kernel for 3B terms" << endl;
         break;
      case RP_3_5:
         kernel3B = rp_3_5_kernel;
         cout << "Using 3,5 reciprocal power kernel for 3B terms" << endl;
         break;
      case RP_3_6:
         kernel3B = rp_3_6_kernel;
         cout << "Using 3,6 reciprocal power kernel for 3B terms" << endl;
         break;
   }
}

void Descriptors::SetKernel4B(int kernel_type) {
   switch (kernel_type) {
      case RP_3_0:
         kernel4B = rp_3_0_kernel;
         cout << "Using 3,0 reciprocal power kernel for 4B terms" << endl;
         break;
      case RP_3_1:
         kernel4B = rp_3_1_kernel;
         cout << "Using 3,1 reciprocal power kernel for 4B terms" << endl;
         break;
      case RP_3_2:
         kernel4B = rp_3_2_kernel;
         cout << "Using 3,2 reciprocal power kernel for 4B terms" << endl;
         break;
      case RP_3_3:
         kernel4B = rp_3_3_kernel;
         cout << "Using 3,3 reciprocal power kernel for 4B terms" << endl;
         break;
      case RP_3_4:
         kernel4B = rp_3_4_kernel;
         cout << "Using 3,4 reciprocal power kernel for 4B terms" << endl;
         break;
      case RP_3_5:
         kernel4B = rp_3_5_kernel;
         cout << "Using 3,5 reciprocal power kernel for 4B terms" << endl;
         break;
      case RP_3_6:
         kernel4B = rp_3_6_kernel;
         cout << "Using 3,6 reciprocal power kernel for 4B terms" << endl;
         break;
   }
}

int Descriptors::getposition(int p[], int m[]) {
   int sum = p[0]*m[0] + p[1]*m[1] + p[2]*m[2] + p[3]*m[3] + p[4]*m[4] + p[5];
   return sum;
}

bool Descriptors::IsEmbeddable(double d1, double d2, double d3, double d4, double d5, double d6) {
   double U, V, W, u, v, w, X, x, Y, y, Z, z, a, b, c, d, tmp;
   U = d1;
   V = d2;
   W = d4;
   u = d6;
   v = d5;
   w = d3;
   X = (w - U + v)*(U + v + w);
   x = (U - v + w)*(v - w + U);
   Y = (u - V + w)*(V + w + u);
   y = (V - w + u)*(w - u + V);
   Z = (v - W + u)*(W + u + v);
   z = (W - u + v)*(u - v + W);
   tmp = x*Y*Z;
   if (tmp < 0.0) return false;
   a = sqrt(tmp);
   tmp = y*Z*X;
   if (tmp < 0.0) return false;
   b = sqrt(tmp);
   tmp = z*X*Y;
   if (tmp < 0.0) return false;
   c = sqrt(tmp);
   tmp = x*y*z;
   if (tmp < 0.0) return false;
   d = sqrt(tmp);
   tmp = (b + c + d - a)*(a - b + c + d)*(a + b - c + d)*(a + b + c - d);
   if (tmp < 0.0) return false;
   return true;
}

void Descriptors::reorder_idx(int idxsrc[], int idx[], int nperm, int perm[][6]) {
   int i, j;
   int *p;
   SixtupletCompare sc;
   Sixtuplet s_best, s;
   for (i=0; i < 6; i++) s_best.v[i] = idxsrc[i];
   for (i = 1; i < nperm; i++) {
      p = perm[i];
      for (j=0; j < 6; j++) s.v[j] = idxsrc[p[j]];
      if (sc(s, s_best)) {
         for (j=0; j < 6; j++) s_best.v[j] = s.v[j];
      }
   }
   for (i=0; i < 6; i++) idx[i] = s_best.v[i];
}

void Descriptors::Calculate(double *r, double *x, int wsi) {
   int i, j, k, l, tp, tpb, i1, i2, i3, i4;
   int *dloc;
   double *x1, *x2, *dist, *xp;
   double d1[6];
   double d2[6];
   double (*symmetrized_kernels3[3])(double*, double*, double(*)(double, double));
   symmetrized_kernels3[0] = symmetrized_kernel3_0;
   symmetrized_kernels3[1] = symmetrized_kernel3_1;
   symmetrized_kernels3[2] = symmetrized_kernel3_2;
   double (*symmetrized_kernels4[4])(double*, double*, double(*)(double, double));
   symmetrized_kernels4[0] = symmetrized_kernel4_0;
   symmetrized_kernels4[1] = symmetrized_kernel4_1;
   symmetrized_kernels4[2] = symmetrized_kernel4_2;
   symmetrized_kernels4[3] = symmetrized_kernel4_3;
   dist = m_dist + wsi*m_na2;
   for (i=0; i < m_ndesc; i++) x[i] = 0.0;
   for (i = 0, k = 0; i < m_na-1; i++) {
      x1 = r + 3*i;
      for (j = i + 1; j < m_na; j++) {
         x2 = r + 3*j;
         if (m_use_inverse_space)
            dist[k] = 1.0/sqrt(sq_dist(x1, x2, 3));
         else
            dist[k] = sqrt(sq_dist(x1, x2, 3));
         k++;
      }
   }
   for (i = 0; i < m_na2; i++) {
      tp = m_t2[i];
      tpb = m_t2b[tp];
      for (l = 0; l < m_ngp2[tp]; l++) {
         x[tpb+l] += kernel2B(m_svp2[tp][l], dist[i]);
      }
   }
   if (m_nt3) {
      xp = x + m_ngp2sum;
      dloc = m_t3 + m_na3;
      for (i1=0,k=0; i1 < m_na-2; i1++) {
         for (i2=i1+1; i2 < m_na-1; i2++) {
            for (i3=i2+1; i3 < m_na; i3++) {
               tp = m_t3[k];
               tpb = m_t3b[tp];
               d2[0] = dist[dloc[3*k]];
               d2[1] = dist[dloc[3*k+1]];
               d2[2] = dist[dloc[3*k+2]];
               for (l=0; l < m_ngp3[tp]; l++) {
                  d1[0] = m_svp3[tp][3*l];
                  d1[1] = m_svp3[tp][3*l+1];
                  d1[2] = m_svp3[tp][3*l+2];
                  xp[tpb+l] += symmetrized_kernels3[m_st3[tp]](d1, d2, kernel3B);
               }
               k++;
            }
         }
      }
   }
   if (m_nt4) {
      xp = x + m_ngp3sum;
      dloc = m_t4 + m_na4;
      for (i1=0,k=0; i1 < m_na-3; i1++) {
         for (i2=i1+1; i2 < m_na-2; i2++) {
            for (i3=i2+1; i3 < m_na-1; i3++) {
               for (i4=i3+1; i4 < m_na; i4++) {
                  tp = m_t4[k];
                  tpb = m_t4b[tp];
                  d2[0] = dist[dloc[6*k]];
                  d2[1] = dist[dloc[6*k+1]];
                  d2[2] = dist[dloc[6*k+2]];
                  d2[3] = dist[dloc[6*k+3]];
                  d2[4] = dist[dloc[6*k+4]];
                  d2[5] = dist[dloc[6*k+5]];
                  for (l=0; l < m_ngp4[tp]; l++) {
                     d1[0] = m_svp4[tp][6*l];
                     d1[1] = m_svp4[tp][6*l+1];
                     d1[2] = m_svp4[tp][6*l+2];
                     d1[3] = m_svp4[tp][6*l+3];
                     d1[4] = m_svp4[tp][6*l+4];
                     d1[5] = m_svp4[tp][6*l+5];
                     xp[tpb+l] += symmetrized_kernels4[m_st4[tp]](d1, d2, kernel4B);
                  }
                  k++;
               }
            }
         }
      }
   }
}

// Return closest position in 1D grid v to val
int Descriptors::GetPos(double *v, int n, double val) {
   auto lower = std::lower_bound(v, v+n, val);
   if (lower == v+n)
      return n-1;
   auto p = lower - v;
   if (p == 0)
      return 0;
   if ((*lower) - val < val - v[p-1])
      return p;
   return p - 1;
}

void Descriptors::GetIdx(double *d, int *idxs, int tp, int bo) {
   int j;
   if (bo == 2) {
      idxs[0] = GetPos(m_svp2[tp], m_ngp2[tp], d[0]);
   } else if (bo == 3) {
      for (j = 0; j < 3; j++) {
         idxs[j] = GetPos(m_svp2[m_pt3[3*tp+j]], m_ngp2[m_pt3[3*tp+j]], d[j]);
      }
   } else if (bo == 4) {
      for (j = 0; j < 6; j++) {
         idxs[j] = GetPos(m_svp2[m_pt4[6*tp+j]], m_ngp2[m_pt4[6*tp+j]], d[j]);
      }
   }
}

void Descriptors::FindClosestMultiplets(double *r, set<int> *pss, set<Triplet, TripletCompare> *tss, set<Sixtuplet, SixtupletCompare> *sss) {
   int i, j, k, tp, i1, i2, i3, i4, stp;
   int *dloc;
   int idx[6];
   int idxsrc[6];
   double *x1, *x2, *dist;
   double d2[6];
   Triplet t;
   Sixtuplet s;
   dist = m_dist;
   for (i=0,k=0; i < m_na-1; i++) {
      x1 = r + 3*i;
      for (j=i+1; j < m_na; j++) {
         tp = m_t2[k];
         x2 = r + 3*j;
         if (m_use_inverse_space)
            d2[0] = dist[k] = 1.0/sqrt(sq_dist(x1, x2, 3));
         else
            d2[0] = dist[k] = sqrt(sq_dist(x1, x2, 3));
         GetIdx(d2, idxsrc, tp, 2);
         idx[0] = idxsrc[0];
         pss[tp].insert(idx[0]);
         k++;
      }
   }
   if (!m_nt3) return;
   dloc = m_t3 + m_na3;
   for (i1=0,k=0; i1 < m_na-2; i1++) {
      for (i2=i1+1; i2 < m_na-1; i2++) {
         for (i3=i2+1; i3 < m_na; i3++) {
            tp = m_t3[k];
            stp = m_st3[tp];
            d2[0] = dist[dloc[3*k]];
            d2[1] = dist[dloc[3*k+1]];
            d2[2] = dist[dloc[3*k+2]];
            GetIdx(d2, idxsrc, tp, 3);
            idx[0] = idxsrc[0];
            idx[1] = idxsrc[1];
            idx[2] = idxsrc[2];
            if (stp == 0) sort(idx, idx+3);
            else if (stp == 1) sort(idx, idx+2);
            t.v[0] = idx[0];
            t.v[1] = idx[1];
            t.v[2] = idx[2];
            tss[tp].insert(t);
            k++;
         }
      }
   }
   if (!m_nt4) return;
   dloc = m_t4 + m_na4;
   for (i1=0,k=0; i1 < m_na-3; i1++) {
      for (i2=i1+1; i2 < m_na-2; i2++) {
         for (i3=i2+1; i3 < m_na-1; i3++) {
            for (i4=i3+1; i4 < m_na; i4++) {
               tp = m_t4[k];
               stp = m_st4[tp];
               d2[0] = dist[dloc[6*k]];
               d2[1] = dist[dloc[6*k+1]];
               d2[2] = dist[dloc[6*k+2]];
               d2[3] = dist[dloc[6*k+3]];
               d2[4] = dist[dloc[6*k+4]];
               d2[5] = dist[dloc[6*k+5]];
               GetIdx(d2, idxsrc, tp, 4);
               if (stp == 0) {
                  reorder_idx(idxsrc, idx, 24, m_perm0);
               } else if (stp == 1) {
                  reorder_idx(idxsrc, idx, 6, m_perm1);
               } else if (stp == 2) {
                  reorder_idx(idxsrc, idx, 4, m_perm2);
               } else {
                  reorder_idx(idxsrc, idx, 2, m_perm2);
               }
               s.v[0] = idx[0];
               s.v[1] = idx[1];
               s.v[2] = idx[2];
               s.v[3] = idx[3];
               s.v[4] = idx[4];
               s.v[5] = idx[5];
               sss[tp].insert(s);
               k++;
            }
         }
      }
   }
}

bool Descriptors::ReadGridPoints(int idst, const char *fn) {
   int i, ngp, gp;
   ifstream ifile;
   istringstream isstream;
   string line;
   ifile.open(fn);
   if (ifile) {
      getline(ifile, line, '\n');
      isstream.str(line);
      isstream >> ngp;
      isstream.clear();
      m_ngp2[idst] = ngp;
      for (i = 0; i < ngp; i++) {
         getline(ifile, line, '\n');
         isstream.str(line);
         isstream >> gp;
         m_svp2[idst][i] = m_svp2[idst][gp];
         isstream.clear();
      }
      cout << "Read " << ngp << " grid points for pairs of type " << idst << endl;
   } else {
      cerr << "Unable to open " << fn << endl;
      return false;
   }
   ifile.close();
   return true;
}

void Descriptors::AssignGridPoints(int isrc1, int isrc2, int isrc3, int idst) {
   int i1, i2, i3, count, inxt;
   double d1, d2, d3, d1o, d2o, d3o, hs;
   for (i1=0,count = 0; i1 < m_ngp2[isrc1]; i1++) {
      d1 = m_svp2[isrc1][i1];
      if (m_use_inverse_space) d1 = 1.0/d1;
      if (isrc1 == isrc2) inxt = i1; else inxt = 0;
      for (i2=inxt; i2 < m_ngp2[isrc2]; i2++) {
         d2 = m_svp2[isrc2][i2];
         if (m_use_inverse_space) d2 = 1.0/d2;
         if (isrc2 == isrc3) inxt = i2; else inxt = 0;
         for (i3=inxt; i3 < m_ngp2[isrc3]; i3++) {
            d3 = m_svp2[isrc3][i3];
            if (m_use_inverse_space) d3 = 1.0/d3;
            hs = 0.5*(d1 + d2 + d3);
            if ((hs - d1) >= 0.0 && (hs - d2) >= 0.0 && (hs - d3) >= 0.0) count++;
         }
      }
   }
   m_ngp3[idst] = count;
   m_svp3[idst] = new double[3*count];
   m_svp3i[idst] = new int[3*count];
   for (i1=0,count=0; i1 < m_ngp2[isrc1]; i1++) {
      d1o = m_svp2[isrc1][i1];
      if (m_use_inverse_space) d1 = 1.0/d1o;
      else d1 = d1o;
      if (isrc1 == isrc2) inxt = i1; else inxt = 0;
      for (i2=inxt; i2 < m_ngp2[isrc2]; i2++) {
         d2o = m_svp2[isrc2][i2];
         if (m_use_inverse_space) d2 = 1.0/d2o;
         else d2 = d2o;
         if (isrc2 == isrc3) inxt = i2; else inxt = 0;
         for (i3=inxt; i3 < m_ngp2[isrc3]; i3++) {
            d3o = m_svp2[isrc3][i3];
            if (m_use_inverse_space) d3 = 1.0/d3o;
            else d3 = d3o;
            hs = 0.5*(d1 + d2 + d3);
            if ((hs - d1) >= 0.0 && (hs - d2) >= 0.0 && (hs - d3) >= 0.0) {
               m_svp3[idst][3*count] = d1o; m_svp3[idst][3*count+1] = d2o; m_svp3[idst][3*count+2] = d3o;
               m_svp3i[idst][3*count] = i1; m_svp3i[idst][3*count+1] = i2; m_svp3i[idst][3*count+2] = i3;
               count++;
            }
         }
      }
   }
   cout << "Generated " << count << " grid points for triplets of type " << idst << endl;
}

void Descriptors::AssignGridPoints(int isrc1, int isrc2, int isrc3, int isrc4, int isrc5, int isrc6, int idst, int nperm) {
   int i, j, i1, i2, i3, i4, i5, i6, count, gptot, pos, countkeep;
   double d[6];
   int *p;
   int orig[6];
   int reordered[6];
   int (*perm)[6];
   bool *ignore;
   int m[5];
   if (nperm == 24) {
      perm = m_perm0;
   } else if (nperm == 4) {
      perm = m_perm2;
   } else if (nperm == 6) {
      perm = m_perm1;
   } else {
      perm = m_perm2;
   }
   gptot = m_ngp2[isrc1]*m_ngp2[isrc2]*m_ngp2[isrc3]*m_ngp2[isrc4]*m_ngp2[isrc5]*m_ngp2[isrc6];
   ignore = new bool[gptot];
   for (i = 0; i < gptot; i++) ignore[i] = false;
   m[4] = m_ngp2[isrc6];
   m[3] = m[4]*m_ngp2[isrc5];
   m[2] = m[3]*m_ngp2[isrc4];
   m[1] = m[2]*m_ngp2[isrc3];
   m[0] = m[1]*m_ngp2[isrc2];
   countkeep = 0;
   for (i1 = 0, count = 0; i1 < m_ngp2[isrc1]; i1++) {
      for (i2 = 0; i2 < m_ngp2[isrc2]; i2++) {
         for (i3 = 0; i3 < m_ngp2[isrc3]; i3++) {
            for (i4 = 0; i4 < m_ngp2[isrc4]; i4++) {
               for (i5 = 0; i5 < m_ngp2[isrc5]; i5++) {
                  for (i6 = 0; i6 < m_ngp2[isrc6]; i6++) {
                     if (!ignore[count]) {
                        orig[0] = i1; orig[1] = i2; orig[2] = i3; orig[3] = i4; orig[4] = i5; orig[5] = i6;
                        for (i=1; i < nperm; i++) {
                           p = perm[i];
                           for (j=0; j < 6; j++) reordered[j] = orig[p[j]];
                           pos = getposition(reordered, m);
                           if (pos > count) ignore[pos] = true;
                        }
                        if (m_use_inverse_space) {
                           d[0] = 1.0/m_svp2[isrc1][i1];
                           d[1] = 1.0/m_svp2[isrc2][i2];
                           d[2] = 1.0/m_svp2[isrc3][i3];
                           d[3] = 1.0/m_svp2[isrc4][i4];
                           d[4] = 1.0/m_svp2[isrc5][i5];
                           d[5] = 1.0/m_svp2[isrc6][i6];
                        } else {
                           d[0] = m_svp2[isrc1][i1];
                           d[1] = m_svp2[isrc2][i2];
                           d[2] = m_svp2[isrc3][i3];
                           d[3] = m_svp2[isrc4][i4];
                           d[4] = m_svp2[isrc5][i5];
                           d[5] = m_svp2[isrc6][i6];
                        }
                        if (IsEmbeddable(d[0], d[1], d[2], d[3], d[4], d[5])) {
                           countkeep++;
                        } else {
                           ignore[count] = true;
                        }
                     }
                     count++;
                  }
               }
            }
         }
      }
   }
   m_ngp4[idst] = countkeep;
   m_svp4[idst] = new double[6*countkeep];
   m_svp4i[idst] = new int[6*countkeep];
   countkeep = 0;
   for (i1 = 0, count = 0; i1 < m_ngp2[isrc1]; i1++) {
      for (i2 = 0; i2 < m_ngp2[isrc2]; i2++) {
         for (i3 = 0; i3 < m_ngp2[isrc3]; i3++) {
            for (i4 = 0; i4 < m_ngp2[isrc4]; i4++) {
               for (i5 = 0; i5 < m_ngp2[isrc5]; i5++) {
                  for (i6 = 0; i6 < m_ngp2[isrc6]; i6++) {
                     if (!ignore[count]) {
                        m_svp4[idst][6*countkeep] = m_svp2[isrc1][i1];
                        m_svp4[idst][6*countkeep+1] = m_svp2[isrc1][i2];
                        m_svp4[idst][6*countkeep+2] = m_svp2[isrc1][i3];
                        m_svp4[idst][6*countkeep+3] = m_svp2[isrc1][i4];
                        m_svp4[idst][6*countkeep+4] = m_svp2[isrc1][i5];
                        m_svp4[idst][6*countkeep+5] = m_svp2[isrc1][i6];
                        m_svp4i[idst][6*countkeep] = i1;
                        m_svp4i[idst][6*countkeep+1] = i2;
                        m_svp4i[idst][6*countkeep+2] = i3;
                        m_svp4i[idst][6*countkeep+3] = i4;
                        m_svp4i[idst][6*countkeep+4] = i5;
                        m_svp4i[idst][6*countkeep+5] = i6;
                        countkeep++;
                     }
                     count++;
                  }
               }
            }
         }
      }
   }
   cout << "Generated " << countkeep << " grid points for quadruplets of type " << idst << endl;
   delete [] ignore;
}

bool Descriptors::ReadGridPoints(int isrc1, int isrc2, int isrc3, int idst, const char *fn) {
   int i, ngp, gp1, gp2, gp3;
   ifstream ifile;
   istringstream isstream;
   string line;
   ifile.open(fn);
   if (ifile) {
      getline(ifile, line, '\n');
      isstream.str(line);
      isstream >> ngp;
      isstream.clear();
      m_ngp3[idst] = ngp;
      if (ngp) {
         m_svp3[idst] = new double[3*ngp];
         m_svp3i[idst] = new int[3*ngp];
      }
      for (i = 0; i < ngp; i++) {
         getline(ifile, line, '\n');
         isstream.str(line);
         isstream >> gp1 >> gp2 >> gp3;
         m_svp3[idst][3*i] = m_svp2[isrc1][gp1];
         m_svp3[idst][3*i+1] = m_svp2[isrc2][gp2];
         m_svp3[idst][3*i+2] = m_svp2[isrc3][gp3];
         m_svp3i[idst][3*i] = gp1;
         m_svp3i[idst][3*i+1] = gp2;
         m_svp3i[idst][3*i+2] = gp3;
         isstream.clear();
      }
      cout << "Read " << ngp << " grid points for triplets of type " << idst << endl;
   } else {
      cerr << "Unable to open " << fn << endl;
      return false;
   }
   ifile.close();
   return true;
}

bool Descriptors::ReadGridPoints(int isrc1, int isrc2, int isrc3, int isrc4, int isrc5, int isrc6, int idst, const char *fn) {
   int i, ngp, gp1, gp2, gp3, gp4, gp5, gp6;
   ifstream ifile;
   istringstream isstream;
   string line;
   ifile.open(fn);
   if (ifile) {
      getline(ifile, line, '\n');
      isstream.str(line);
      isstream >> ngp;
      isstream.clear();
      m_ngp4[idst] = ngp;
      if (ngp) {
         m_svp4[idst] = new double[6*ngp];
         m_svp4i[idst] = new int[6*ngp];
      }
      for (i = 0; i < ngp; i++) {
         getline(ifile, line, '\n');
         isstream.str(line);
         isstream >> gp1 >> gp2 >> gp3 >> gp4 >> gp5 >> gp6;
         m_svp4[idst][6*i] = m_svp2[isrc1][gp1];
         m_svp4[idst][6*i+1] = m_svp2[isrc2][gp2];
         m_svp4[idst][6*i+2] = m_svp2[isrc3][gp3];
         m_svp4[idst][6*i+3] = m_svp2[isrc4][gp4];
         m_svp4[idst][6*i+4] = m_svp2[isrc5][gp5];
         m_svp4[idst][6*i+5] = m_svp2[isrc6][gp6];
         m_svp4i[idst][6*i] = gp1;
         m_svp4i[idst][6*i+1] = gp2;
         m_svp4i[idst][6*i+2] = gp3;
         m_svp4i[idst][6*i+3] = gp4;
         m_svp4i[idst][6*i+4] = gp5;
         m_svp4i[idst][6*i+5] = gp6;
         isstream.clear();
      }
      cout << "Read " << ngp << " grid points for quadruplets of type " << idst << endl;
   } else {
      cerr << "Unable to open " << fn << endl;
      return false;
   }
   ifile.close();
   return true;
}

void Descriptors::Get3BGridPoints(set<Triplet, TripletCompare> *tss) {
   int i, j;
   Triplet t;
   for (i = 0; i < m_nt3; i++) {
      for (j = 0; j < m_ngp3[i]; j++) {
         t.v[0] = m_svp3i[i][3*j];
         t.v[1] = m_svp3i[i][3*j+1];
         t.v[2] = m_svp3i[i][3*j+2];
         tss[i].insert(t);
      }
   }
}

void Descriptors::Get4BGridPoints(set<Sixtuplet, SixtupletCompare> *sss) {
   int i, j;
   Sixtuplet s;
   for (i = 0; i < m_nt4; i++) {
      for (j = 0; j < m_ngp4[i]; j++) {
         s.v[0] = m_svp4i[i][6*j];
         s.v[1] = m_svp4i[i][6*j+1];
         s.v[2] = m_svp4i[i][6*j+2];
         s.v[3] = m_svp4i[i][6*j+3];
         s.v[4] = m_svp4i[i][6*j+4];
         s.v[5] = m_svp4i[i][6*j+5];
         sss[i].insert(s);
      }
   }
}

void Descriptors::SetUseInverseSpace(bool val) {
   m_use_inverse_space = val;
}

double (*Descriptors::kernel2B)(double, double) = 0;
double (*Descriptors::kernel3B)(double, double) = 0;
double (*Descriptors::kernel4B)(double, double) = 0;
int Descriptors::m_perm0[][6] =
{
   {0, 1, 2, 3, 4, 5},
   {0, 2, 1, 4, 3, 5},
   {1, 0, 2, 3, 5, 4},
   {1, 2, 0, 5, 3, 4},
   {2, 0, 1, 4, 5, 3},
   {2, 1, 0, 5, 4, 3},
   {0, 3, 4, 1, 2, 5},
   {0, 4, 3, 2, 1, 5},
   {3, 0, 4, 1, 5, 2},
   {3, 4, 0, 5, 1, 2},
   {4, 0, 3, 2, 5, 1},
   {4, 3, 0, 5, 2, 1},
   {1, 3, 5, 0, 2, 4},
   {1, 5, 3, 2, 0, 4},
   {3, 1, 5, 0, 4, 2},
   {3, 5, 1, 4, 0, 2},
   {5, 1, 3, 2, 4, 0},
   {5, 3, 1, 4, 2, 0},
   {2, 4, 5, 0, 1, 3},
   {2, 5, 4, 1, 0, 3},
   {4, 2, 5, 0, 3, 1},
   {4, 5, 2, 3, 0, 1},
   {5, 2, 4, 1, 3, 0},
   {5, 4, 2, 3, 1, 0}
};
int Descriptors::m_perm1[][6] =
{
   {0, 1, 2, 3, 4, 5},
   {0, 2, 1, 4, 3, 5},
   {1, 0, 2, 3, 5, 4},
   {1, 2, 0, 5, 3, 4},
   {2, 0, 1, 4, 5, 3},
   {2, 1, 0, 5, 4, 3}
};
int Descriptors::m_perm2[][6] =
{
   {0, 1, 2, 3, 4, 5},
   {0, 2, 1, 4, 3, 5},
   {0, 3, 4, 1, 2, 5},
   {0, 4, 3, 2, 1, 5}
};
