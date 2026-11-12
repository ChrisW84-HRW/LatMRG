#define TYPES_CODE ZD     // Int = ZZ, Real = double

#include <cmath>
#include <cstdint>
#include <ctime>
#include <iomanip>
#include <iostream>

#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/SearchMRG.h"
#include "latmrg/SearchCombinedMRG.h"

using namespace LatMRG;

int main()
{
    constexpr long k = 3;
    constexpr int64_t maxdim = 12;

    const Int m1 = NTL::to_ZZ("9223372036854773561"); // 2^63 - 2247
    const Int m2 = NTL::to_ZZ("9223372036854766199"); // 2^63 - 9609

    IntVec low1, high1, box1;
    IntVec low2, high2, box2;

    low1.SetLength(k + 1);
    high1.SetLength(k + 1);
    box1.SetLength(k + 1);

    low2.SetLength(k + 1);
    high2.SetLength(k + 1);
    box2.SetLength(k + 1);

    // First component: 2^25 < a1 < 2^35, a2 = 0,
    //                  -2^35 < a3 < -2^25.
    low1[1]  = NTL::power(Int(2), 25) + 1;
    high1[1] = NTL::power(Int(2), 35) - 1;
    low1[2]  = high1[2] = Int(0);
    low1[3]  = -NTL::power(Int(2), 35) + 1;
    high1[3] = -NTL::power(Int(2), 25) - 1;

    box1[1] = Int(8);
    box1[2] = Int(1);
    box1[3] = Int(4);

    // Second component: a1 = 0, 2^25 < a2 < 2^32,
    //                   -2^32 < a3 < -2^25.
    low2[1]  = high2[1] = Int(0);
    low2[2]  = NTL::power(Int(2), 25) + 1;
    high2[2] = NTL::power(Int(2), 32) - 1;
    low2[3]  = -NTL::power(Int(2), 32) + 1;
    high2[3] = -NTL::power(Int(2), 25) - 1;

    box2[1] = Int(1);
    box2[2] = Int(8);
    box2[3] = Int(4);

    SearchMRG<Int, Real> first;
    SearchMRG<Int, Real> second;

    auto configureComponent = [&](SearchMRG<Int, Real>& search, const Int& modulus, const IntVec& low, const IntVec& high, const IntVec& box)
    {
        search.setModulus(modulus);
        search.setOrder(k);
        search.setLowBoundaries(low);
        search.setHighBoundaries(high);
        search.setSubregionSize(box);
        search.setNumRandomRegions(100);

        search.setRequireApproxFactoringCondition(true);
        search.setMaxPeriod(true);
        search.setCalculateFOM(false);

        // Each of 100 regions contains 8 * 1 * 4 = 32 vectors.
        // Leave room for every full-period vector found.
        search.setNumBestGen(3200);
        search.setMaxNumGen(3201);

        search.setTimeLimit(7200.0);
        search.setProgress(false);
        search.setVerbose(1);
    };

    configureComponent(first, m1, low1, high1, box1);
    configureComponent(second, m2, low2, high2, box2);

    first.performSearch(&SearchMRG<Int, Real>::searchLoopRandomApproxFactoringByAkBox);

    second.performSearch(&SearchMRG<Int, Real>::searchLoopRandomApproxFactoringByAkBox);

    const auto firstCount = first.bestCandidates.getList().size();
    const auto secondCount = second.bestCandidates.getList().size();

    std::cout << "Full-period candidates for component 1: " << firstCount << '\n'
              << "Full-period candidates for component 2: " << secondCount << '\n';

    SearchCombinedMRG<Int, Real> combined;
    combined.addComponentCandidates(first.bestCandidates);
    combined.addComponentCandidates(second.bestCandidates);

    // The combined performSearch() wrapper calculates the product modulus
    // before SearchBase constructs its projection lattice. Calculate it here
    // as well because the FOM normalizer needs the modulus now.
    combined.calculateModulus();
    combined.setMaxDim(maxdim);

    NTL::Vec<int64_t> t;
    t.SetLength(3);
    t[0] = 8;
    t[1] = 0;
    t[2] = 0;

    ReducerBB<Int, Real> red(maxdim);
    red.maxNodesBB = 1000000;

    NormaBestLat norma(-log(combined.getModulus()), k, maxdim);
    WeightsUniform weights(1.0);
    FigureOfMeritDualM<Int, Real> fom(t, weights, norma, &red, true);
    fom.setLowBound(0.0);

    combined.setFOM(fom);
    combined.setDualLattice(true);
    combined.setCalculateFOM(true);
    combined.setBestFOM(true);
    combined.setNumBestGen(2);

    // There can be at most 3200 * 3200 input combinations.
    // Set the generator limit above that bound so the odometer,
    // rather than this limit, ends a complete search.
    combined.setMaxNumGen(3200 * 3200 + 1);
    combined.setTimeLimit(14400.0);
    combined.setProgress(true);
    combined.setVerbose(1);

    const std::clock_t start = std::clock();

    combined.performSearch(&SearchCombinedMRG<Int, Real>::searchLoopDeterministic);

    const std::clock_t end = std::clock();

    std::cout << std::setprecision(12)
              << "\nCombined-search CPU time: "
              << static_cast<double>(end - start) / CLOCKS_PER_SEC
              << " seconds\n"
              << "Combinations examined: "
              << combined.numCandidatesExamined << '\n'
              << "Best merit: "
              << combined.getBestMerit() << '\n';

    int rank = 0;
    for (const Candidate<Int>& candidate : combined.bestCandidates.getList()) {
    std::cout << ++rank << ". Combined coefficients: " << candidate.getAa() << ", merit: " << candidate.getMerit() << '\n';
    const IntVec& a = candidate.getAa();
        for (long j = 1; j <= k; ++j) {
            Int a1 = a[j] % m1;
            Int a2 = a[j] % m2;
            if (a1 > m1 / 2) a1 -= m1;
            if (a2 > m2 / 2) a2 -= m2;
            std::cout << "   a_" << j
                  << ": component 1 = " << a1
                  << ", component 2 = " << a2 << '\n';
        }
    }

    return 0;
}