// TestSearchEx5.cc

#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <cmath>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/SearchBase.h"
#include "latmrg/SearchMRG.h"
#include "latmrg/SearchCombinedMRG.h"

using namespace LatMRG;

int main()
{
    constexpr int k = 1;
    constexpr int64_t maxdim = 10;

    const NTL::ZZ m = NTL::to_ZZ("2147483647");

    // Successive-coordinate FoM, starting at dimension 2
    // and ending at dimension 8.
    NTL::Vec<int64_t> t;
    t.SetLength(k);
    t[0] = 8;

    ReducerBB<Int, Real> red(maxdim);
    red.maxNodesBB = 1000000;

    NormaBestLat norma(-log(m), k, maxdim);
    WeightsUniform weights(1.0);
    FigureOfMeritDualM<Int, Real> fom(
        t, weights, norma, &red, true);
    fom.setLowBound(0.0);

    IntVec lower, upper;
    lower.SetLength(k + 1);
    upper.SetLength(k + 1);
    lower[1] = NTL::to_ZZ("40000");
    upper[1] = NTL::to_ZZ("1000000000");

    SearchMRG<Int, Real> searcher;
    searcher.setModulus(m);
    searcher.setOrder(k);
    searcher.setMaxDim(maxdim);
    searcher.setLowBoundaries(lower);
    searcher.setHighBoundaries(upper);

    searcher.setFOM(fom);
    searcher.setDualLattice(true);
    searcher.setCalculateFOM(true);
    searcher.setMaxPeriod(true);
    searcher.setRequireApproxFactoringCondition(true);
    searcher.setBestFOM(true);
    searcher.setNumBestGen(10);

    // With maxPeriod=true, this limits maximal-period candidates.
    // It is set above the number reported for this search.
    searcher.setMaxNumGen(13183);
    searcher.setTimeLimit(7200.0);
    searcher.setProgress(true);
    searcher.setVerbose(1);
    searcher.setOutputToGenFile(false);

    
    const std::clock_t start = std::clock();
    searcher.performSearch(
        &SearchMRG<Int, Real>::searchLoopExhaustive);
    const std::clock_t end = std::clock();
    const double cpuSeconds = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Search CPU time: " << cpuSeconds << " seconds\n";
    std::cout << "Best merit: " << searcher.getBestMerit() << "\n";
    

    std::cout << std::setprecision(12)
              << "\nCandidates passing coefficient condition: "
              << searcher.numCandidatesExamined << '\n'
              << "Maximal-period candidates: "
              << searcher.numMaxPeriodCandidates << '\n'
              << "Best merit: "
              << searcher.getBestMerit() << '\n';

    int rank = 0;
    for (const Candidate<Int>& candidate :
         searcher.bestCandidates.getList())
    {
        std::cout << ++rank
                  << ". a_1 = " << candidate.getAa()[1]
                  << ", merit = " << candidate.getMerit()
                  << '\n';
    }

    return 0;
}