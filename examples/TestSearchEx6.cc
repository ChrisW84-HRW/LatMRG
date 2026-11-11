/**
 * In this example, we check that the results of the new LatMRG version
 * are consistent with Example 7.6 reported in the article [rLEC97c].
 * We found no relevant discrepancy which cannot be explained by
 * randomness but a significant reduction in runtime.
 */

#define TYPES_CODE ZD     // Int = ZZ, Real = double

#include <cstdint>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/SearchBase.h"
#include "latmrg/SearchMRG.h"

using namespace LatMRG;

int main()
{
    constexpr int k = 3;
    constexpr int64_t maxdim = 12;

    const NTL::ZZ m = NTL::to_ZZ("9223372036854773561");

    // Intended merit: M_8, using dimensions 4 through 8.
    NTL::Vec<int64_t> t;
    t.SetLength(3);
    t[0] = 8;
    t[1] = 0;
    t[2] = 0;

    ReducerBB<Int, Real> red(maxdim);
    red.maxNodesBB = 1000000;

    NormaBestLat norma(-log(m), k, maxdim);
    WeightsUniform weights(1.0);
    FigureOfMeritDualM<Int, Real> fom(t, weights, norma, &red, true);
    fom.setLowBound(0.0);

    IntVec lower, upper, boxSize;
    lower.SetLength(k + 1);
    upper.SetLength(k + 1);
    boxSize.SetLength(k + 1);

    lower[1] = NTL::to_ZZ("1125899906842624");
    upper[1] = NTL::to_ZZ("9223372036854765808");

    lower[2] = NTL::ZZ(0);
    upper[2] = NTL::ZZ(0);

    lower[3] = NTL::to_ZZ("-9223372036854765808");
    upper[3] = NTL::to_ZZ("-1125899906842624");

    boxSize[1] = NTL::ZZ(4);
    boxSize[2] = NTL::ZZ(1);
    boxSize[3] = NTL::ZZ(4);

    SearchMRG<Int, Real> searcher;
    searcher.setModulus(m);
    searcher.setOrder(k);
    searcher.setMaxDim(maxdim);
    searcher.setLowBoundaries(lower);
    searcher.setHighBoundaries(upper);
    searcher.setSubregionSize(boxSize);
    searcher.setNumRandomRegions(1000);

    searcher.setFOM(fom);
    searcher.setDualLattice(true);
    searcher.setCalculateFOM(true);
    searcher.setMaxPeriod(true);
    searcher.setRequireApproxFactoringCondition(true);
    searcher.setBestFOM(true);
    searcher.setNumBestGen(3);

    // There are at most 1000 * 4 * 1 * 4 completed vectors.
    // Keep this limit above that number so it does not end the run.
    searcher.setMaxNumGen(16001);
    searcher.setTimeLimit(7200.0);
    searcher.setProgress(true);
    searcher.setVerbose(1);
    searcher.setOutputToGenFile(false);

    const std::clock_t start = std::clock();
    searcher.performSearch(&SearchMRG<Int, Real>::searchLoopRandomApproxFactoringByAkBox);
    const std::clock_t end = std::clock();

    const double cpuSeconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    std::cout << std::setprecision(12)
              << "\nSearch CPU time: " << cpuSeconds << " seconds\n"
              << "Candidates examined: "
              << searcher.numCandidatesExamined << '\n'
              << "Maximal-period candidates: "
              << searcher.numMaxPeriodCandidates << '\n'
              << "Best merit: "
              << searcher.getBestMerit() << '\n';

    int rank = 0;
    for (const Candidate<Int>& candidate :
         searcher.bestCandidates.getList())
    {
        const IntVec& a = candidate.getAa();
        std::cout << ++rank
                  << ". a_1 = " << a[1]
                  << ", a_2 = " << a[2]
                  << ", a_3 = " << a[3]
                  << ", merit = " << candidate.getMerit()
                  << '\n';
    }

    return 0;
}