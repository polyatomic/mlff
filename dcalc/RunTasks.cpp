#include <iostream>

#include "RunTasks.h"
#include "tasks.h"

extern PTypeExtended g_params;

using std::cout;
using std::endl;

void RunTasks() {
   int i, j, nstr, n, ndesc, nbtot, ndesc_block;
   double *v, *Y, *stdevs, *avgs, *gsp;
   int *p4sym;
   int nblocks[3] = {};
   char *blocks, *blockscp;
   char dummy[] = { 'd', '0', '0', '0' };
   Iterator<Mol> mols;
   Descriptors descs;
   Y = v = stdevs = avgs = gsp = 0;
   blocks = blockscp = 0;
   p4sym = 0;
   TFunctorDaDaI<Descriptors> dcalc;
   descs.SetUseInverseSpace(g_params.uis);
   if (calculation_prepare(&g_params, mols, descs, nstr, n, v, Y, nblocks, &blocks, &gsp, &p4sym)) {
      ndesc = descs.GetNDescriptors();
      stdevs = new double[ndesc];
      avgs = new double[ndesc];
      cout << "Number of descriptors: " << ndesc << endl;
      dcalc.init(&descs, &Descriptors::Calculate);
      nbtot = nblocks[0] + nblocks[1] + nblocks[2];
      blockscp = new char[5*nbtot];
      for (i=0; i < 5*nbtot; i++) blockscp[i] = '\0';
      auto select_block = [blocks, blockscp, dummy, &nbtot](int bl) {
         int i, j;
         for (i=0; i < nbtot; i++) {
            if (i == bl) {
               for (j=0; j < 4; j++)
                  blockscp[5*i + j] = blocks[5*i + j];
            } else {
               for (j=0; j < 4; j++)
                  blockscp[5*i + j] = dummy[j];;
            }
         }
      };
      // Calculate each descriptor block
      for (i = 0, j = 0; i < nblocks[0]; i++) {
         select_block(j);
         descs.SetGrid(nullptr, nullptr, gsp, g_params.gss, blockscp, p4sym);
         ndesc_block = descs.GetNDescriptors();
         cout << "Descriptors in block " << &blockscp[5*j] << ": " << ndesc_block << endl;
         calculate_descriptors(nstr, ndesc_block, n, v, &dcalc, stdevs, avgs, g_params.nthreads, g_params.maxrows, (string("mat") + &blockscp[5*j]).c_str());
         Matrix2File(stdevs, ndesc_block, 1, (string("stdevs") + &blockscp[5*j] + ".bin").c_str());
         Matrix2File(avgs, ndesc_block, 1, (string("avgs") + &blockscp[5*j] + ".bin").c_str());
         j++;
      }
      for (i = 0; i < nblocks[1]; i++) {
         select_block(j);
         descs.SetGrid(nullptr, nullptr, gsp, g_params.gss, blockscp, p4sym);
         ndesc_block = descs.GetNDescriptors();
         cout << "Descriptors in block " << &blockscp[5*j] << ": " << ndesc_block << endl;
         calculate_descriptors(nstr, ndesc_block, n, v, &dcalc, stdevs, avgs, g_params.nthreads, g_params.maxrows, (string("mat") + &blockscp[5*j]).c_str());
         Matrix2File(stdevs, ndesc_block, 1, (string("stdevs") + &blockscp[5*j] + ".bin").c_str());
         Matrix2File(avgs, ndesc_block, 1, (string("avgs") + &blockscp[5*j] + ".bin").c_str());
         j++;
      }
      for (i = 0; i < nblocks[2]; i++) {
         select_block(j);
         descs.SetGrid(nullptr, nullptr, gsp, g_params.gss, blockscp, p4sym);
         ndesc_block = descs.GetNDescriptors();
         cout << "Descriptors in block " << &blockscp[5*j] << ": " << ndesc_block << endl;
         calculate_descriptors(nstr, ndesc_block, n, v, &dcalc, stdevs, avgs, g_params.nthreads, g_params.maxrows, (string("mat") + &blockscp[5*j]).c_str());
         Matrix2File(stdevs, ndesc_block, 1, (string("stdevs") + &blockscp[5*j] + ".bin").c_str());
         Matrix2File(avgs, ndesc_block, 1, (string("avgs") + &blockscp[5*j] + ".bin").c_str());
         j++;
      }
   }
   delete [] blockscp;
   delete [] avgs;
   delete [] stdevs;
   delete [] blocks;
   delete [] p4sym;
   delete [] gsp;
   delete [] v;
   delete [] Y;
}
