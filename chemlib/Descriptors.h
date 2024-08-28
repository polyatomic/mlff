#pragma once

#include <set>

using std::set;

#define RP_3_0 0
#define RP_3_1 1
#define RP_3_2 2
#define RP_3_3 3
#define RP_3_4 4
#define RP_3_5 5
#define RP_3_6 6

class Descriptors {
public:
   Descriptors();
   ~Descriptors();
   bool Init(int na, int nt2, int *t2, int nt3, int *st3, int *t3, int nt4, int *st4, int *t4, int *ngp, int nthreads = 1);
   void SetGrid(double *dmin, double *dmax, double *gmin, double step_size);
   int GetNDescriptors();
   int GetN2BDescriptors(int type = -1);
   int GetN3BDescriptors(int type = -1);
   int GetN4BDescriptors(int type = -1);
   void GetDistanceRanges(double *r, double *dmin, double *dmax);
   void Calculate(double *r, double *x, int wsi = 0);
   void FindClosestMultiplets(double *r, set<int> *pss, set<Triplet, TripletCompare> *tss, set<Sixtuplet, SixtupletCompare> *sss);
   bool ReadGridPoints(int idst, const char *fn);
   void AssignGridPoints(int isrc1, int isrc2, int isrc3, int idst);
   void AssignGridPoints(int isrc1, int isrc2, int isrc3, int isrc4, int isrc5, int isrc6, int idst);
   bool ReadGridPoints(int isrc1, int isrc2, int isrc3, int idst, const char *fn);
   bool ReadGridPoints(int isrc1, int isrc2, int isrc3, int isrc4, int isrc5, int isrc6, int idst, const char *fn);
   void Get3BGridPoints(set<Triplet, TripletCompare> *tss);
   void Get4BGridPoints(set<Sixtuplet, SixtupletCompare> *sss);
   void SetUseInverseSpace(bool val);

   static void SetKernel2B(int kernel_type);
   static void SetKernel3B(int kernel_type);
   static void SetKernel4B(int kernel_type);
   static int getposition(int p[], int m[]);
   static bool IsEmbeddable(double d1, double d2, double d3, double d4, double d5, double d6);

private:
   int m_na;         // Number of atoms
   int m_na2;        // Number of atom pairs
   int m_na3;        // Number of atom triplets
   int m_na4;        // Number of atom quadruplets
   int m_nt2;        // Number of atom type pairs
   int m_nt3;        // Number of atom type triplets
   int m_nt4;        // Number of atom type quadruplets
   int *m_t2;        // Vector of type ids for atom pairs
   int *m_t3;        // Vector of type ids for atom triplets and their pair ids
   int *m_t4;        // Vector of type ids for atom quadruplets and their pair ids
   int *m_st3;       // Subtype ids for each triplet type
   int *m_st4;       // Subtype ids for each quadruplet type
   int *m_ngp2;      // Number of grid points for each atom type pair
   int *m_ngp3;      // Number of grid points for each atom type triplet
   int *m_ngp4;      // Number of grid points for each atom type quadruplet
   double **m_svp2;  // Matrix of pair grid points
   double **m_svp3;  // Matrix of triplet grid points
   double **m_svp4;  // Matrix of quadruplet grid points
   int **m_svp3i;    // Triplet grid point indices
   int **m_svp4i;    // Quadruplet grid point indices
   int m_ndesc;      // Total number of descriptors
   double *m_dist;   // Workspace for interatomic distances
   int *m_t2b;       // Descriptor saving boundaries for different pair types
   int *m_t3b;       // Descriptor saving boundaries for different triplet types
   int *m_t4b;       // Descriptor saving boundaries for different quadruplet types
   int m_ngp2sum;    // Number of 2B descriptors
   int m_ngp3sum;    // Number of 2B and 3B descriptors
   bool m_use_inverse_space;

   void ReleaseGrid();
   int GetPos(double *v, int n, double val);
   void GetIdx(double *d, int *idxs, int tp, int bo);

   static double (*kernel2B)(double, double);
   static double (*kernel3B)(double, double);
   static double (*kernel4B)(double, double);
};
