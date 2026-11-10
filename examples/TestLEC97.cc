/**
 * In this example, we check that the results of the new LatMRG version
 * are consistent with the examples reported in the article [rLEC97c].
 * We found no discrepancy.  We can also compare the timings.
 */
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <algorithm>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/ReducerStatic.h"

// #include "latmrg/LCGComponent.h"
#include "latmrg/LCGCombined.h"
#include "latmrg/MRGCombined.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGLatticeLac.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/LCGLatticeLac.h"
#include "latmrg/MLCGLattice.h"
#include "latmrg/MLCGLatticeLac.h"

using namespace LatticeTester;
using namespace LatMRG;

/**
 * This function calculates the FOM of the dual lattice of an MRG with modulus 'm'
 * and vector of multipliers 'aa', for the successive coordinates in `lowDim` to `highDim` dimensions.
 */
template<typename Int, typename Real>
static void FOMSuccLattice(Int &m, IntVec &aa, int64_t lowDim, int64_t highDim) {
   int64_t k = aa.length() - 1;  // order
   NTL::Vec<int64_t> t;
   LatMRG::MRGLattice<Int, Real> mrg(m, aa, highDim);
   ReducerBB<Int, Real> red(mrg);
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MRGLattice: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);
   if (k == 1) {
      LatMRG::LCGLattice<Int, Real> lcg(m, aa[1], highDim, 1, highDim);
      red.setIntLattice(lcg);
      std::cout << "\nUsing LCGLattice: \n";
      std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
      fomdual.computeMeritSucc(lcg, lowDim, highDim);
   }
}

/**
 * Same as the previous function, but with lacunary indices defined by 'lac'.
 */
template<typename Int, typename Real>
static void FOMSuccLatticeLac(Int &m, IntVec &aa, IntVec lac, int64_t lowDim, int64_t highDim) {
   int64_t k = aa.length() - 1;  // order
   NTL::Vec<int64_t> t;
   LatMRG::MRGLatticeLac<Int, Real> mrg(m, aa, highDim);
   mrg.setLac(lac);
   ReducerBB<Int, Real> red(mrg);
   WeightsUniform weights(1.0);
   NormaBestLat normaDual(-log(m), k, highDim, L2NORM);
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red);
   fomdual.setVerbosity(3);
   std::cout << "\nUsing MRGLatticeLac: \n";
   std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
   fomdual.computeMeritSucc(mrg, lowDim, highDim);
   if (k == 1) {
      LatMRG::LCGLatticeLac<Int, Real> lcg(m, aa[1], highDim);
      lcg.setLac(lac);
      red.setIntLattice(lcg);
      std::cout << "\nUsing LCGLatticeLac: \n";
      std::cout << "lowDim = " << lowDim << ", highDim = " << highDim << "\n";
      fomdual.computeMeritSucc(lcg, lowDim, highDim);
   }
}

int main() {

   int64_t maxdim = 12;
   std::cout << "Types: NTL::ZZ, double \n";
   // NTL::Vec<int64_t> t; // The t-vector for the FoM.

   Int m;
   IntVec aa; // Vector a has size k+1 for an MRG of order k.

   // Example 1.
   std::cout << "\n=============================================================\n";
   std::cout << "An LCG with m = 4294967296, a = 1099087573. \n";
   maxdim = 30;
   m = 4294967296;
   m = m / 4;
   aa.SetLength(2);
   aa[1] = 1099087573;
   FOMSuccLattice<Int, Real>(m, aa, 2, maxdim);
   std::cout << "\nThe results and timings should be compared with Table I of [rLEC97c].\n";

   // Example 2.
   std::cout << "\n=============================================================\n";
   std::cout << "An LCG with m = 2147483647, a = 16807, lacunary with leap d = 131072, \n";
   std::cout << "Lacunary indices 1, 2, 3, d+1, d+2, d+3, 2*d+1, 2*d+2, ...\n";
   maxdim = 42;
   Int leap(131072);  // Value of d
   IntVec lac;
   lac.SetLength(maxdim);
   for (int64_t i = 0; i < maxdim; i++)
      lac[i] = 1 + (i % 3) + leap * (i / 3);
   m = 2147483647;
   aa.SetLength(2);
   aa[1] = 16807;
   FOMSuccLatticeLac<Int, Real>(m, aa, lac, 2, maxdim);
   std::cout << "\nThe results should be compared with Table II of [rLEC97c].\n";

   // Example 3.
   std::cout << "\n=============================================================\n";
   std::cout << "The combined LCG of L'Ecuyer (1988), lacunary with leap d = 2^{30}, \n";
   std::cout << "Lacunary indices 1, 2, 3, d+1, d+2, d+3, 2*d+1, 2*d+2, ..., as above.\n";
   // We first compute the equivalent LCG for the combination.
   LCGCombined<Int> comblcg88;
   comblcg88.addComponent(Int(2147483563), Int(40014));
   comblcg88.addComponent(Int(2147483399), Int(40692));
   comblcg88.computeCombination();
   m = comblcg88.getModulus();
   aa[1] = comblcg88.geta();
   std::cout << "Combined LCG has m = " << m << " and a = " << aa[1] << "\n";
   // Compute lacunary indices.
   leap = Int(1 << 30);  // Value of d = 2^{30}
   for (int64_t i = 0; i < maxdim; i++)
      lac[i] = 1 + (i % 3) + leap * (i / 3);
   FOMSuccLatticeLac<Int, Real>(m, aa, lac, 2, maxdim);
   std::cout << "\nThe results should be compared with Table III of [rLEC97c].\n";

   // Example 4.
   std::cout << "\n=============================================================\n";
   std::cout << "Combination of an MRG of order 2 with an LCG, for full lattice. \n";
   std::cout << "Equivalent to an MRG with k=2, m = 1059855887, a_1 = 919821343, a_2 = 650755204.\n";
   // We first compute the equivalent LCG for the combination.
   m = 1059855887;
   aa.SetLength(3);
   aa[1] = 919821343;
   aa[2] = 650755204;
   FOMSuccLattice<Int, Real>(m, aa, 3, 20);
   std::cout << "\nThe results should be compared with Table V of [rLEC97c].\n";

   // Example 5.
   std::cout << "\n=============================================================\n";
   std::cout << "An LCG with m = 2147483647 and a = 45991\n";
   m = 2147483647;
   aa.SetLength(2);
   aa[1] = 45991;
   FOMSuccLattice<Int, Real>(m, aa, 2, maxdim);
   std::cout << "\nThe results should be compared with Table VI of [rLEC97c].\n";

   // Example 6.
   std::cout << "\n=============================================================\n";
   std::cout << "An MRG with m = 9223372036854773561 and"
         << " aa = (1145902849652723, 0, -1184153554609676). \n";
   m = 9223372036854773561;
   aa.SetLength(4);
   aa[1] = 1145902849652723;
   aa[2] = 0;
   aa[3] = -1184153554609676;
   FOMSuccLattice<Int, Real>(m, aa, 4, maxdim);
   std::cout << "\nThe results should be compared with Table VII of [rLEC97c].\n";

   // Example 7.
   std::cout << "\n=============================================================\n";
   std::cout << "Combined MRGs with m_1 = 2^{63}-2247 and m_2 = 2^{63}-9609.\n";

   Int mrgmm[2] = { (Int(1) << 63) - Int(2247), (Int(1) << 63) - Int(9609) };
   Int mrgaa[6] = { Int(3866005879), Int(0), Int(-3472501966), Int(0), Int(48193584), Int(-3751984989) };
   MRGCombined<Int> combmrg63k3 (3, 2, mrgmm, mrgaa);
   combmrg63k3.computeCombination();
   m = combmrg63k3.getModulus();
   aa = combmrg63k3.getaa();
   std::cout << " m_1 = " << mrgmm[0] << ",  m_2 = " << mrgmm[1] << "\n";
   std::cout << "Combined LCG has m = " << m << " and \n aa = " << aa << "\n";

   FOMSuccLattice<Int, Real>(m, aa, 4, 12);
   std::cout << "\nThe results should be compared with Table VIII of [rLEC97c].\n";
}
