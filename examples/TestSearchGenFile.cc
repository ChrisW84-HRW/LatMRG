/**
 * In this example, we check that the results are the same
 * if we a) perform a search via a loop and store the best
 * candidates in a .gen file and b) we read in the data from
 * the .gen file and calculate the FoM using the stored data only.
 * We do this for MRGs as well as for combined MRGs.
 */

#define TYPES_CODE ZD     // Int = ZZ, Real = double

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>

#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/GenFileReader.h"
#include "latmrg/SearchMRG.h"
#include "latmrg/SearchCombinedMRG.h"

using namespace LatMRG;

namespace {

constexpr long k = 3;
constexpr int64_t maxdim = 16;

/* Keep each search's FOM and its dependencies alive for the entire search. */
struct FomSetup {ReducerBB<Int, Real> reducer;  NormaBestLat normalizer;  WeightsUniform weights;  FigureOfMeritDualM<Int, Real> fom;

    FomSetup(const Int& modulus, const NTL::Vec<int64_t>& t) : 
    reducer(maxdim),  normalizer(-log(modulus), k, maxdim), weights(1.0), fom(t, weights, normalizer, &reducer, true)
    {
        reducer.maxNodesBB = 1000000;
        fom.setLowBound(0.0);
    }
};

template<typename Search>
void configureSearch(Search& search, FomSetup& setup, int maxTests)
{
    search.setMaxDim(maxdim);
    search.setFOM(setup.fom);
    search.setDualLattice(true);
    search.setCalculateFOM(true);
    search.setMaxPeriod(false);
    search.setBestFOM(true);
    search.setNumBestGen(3);
    search.setMaxNumGen(maxTests);
    search.setTimeLimit(3600.0);
    search.setProgress(false);
    search.setVerbose(0);
}

void configureBounds(SearchMRG<Int, Real>& search, const Int& modulus)
{
    IntVec low, high;
    low.SetLength(k + 1);
    high.SetLength(k + 1);
    for (long j = 1; j <= k; ++j) {
        low[j] = 1;
        high[j] = modulus - 1;
    }
    search.setLowBoundaries(low);
    search.setHighBoundaries(high);
}

void report(const char* title, const CandidateList<Int>& candidates)
{
    std::cout << title << '\n';
    int rank = 0;
    for (const auto& candidate : candidates.getList()) {
        std::cout << "  " << ++rank
                  << ": merit = " << candidate.getMerit()
                  << ", coefficients = " << candidate.getAa() << '\n';
    }
}

} // namespace

int main()
{
    const Int m1 = NTL::to_ZZ(1048573);
    const Int m2 = NTL::to_ZZ(1048583);

    NTL::Vec<int64_t> t;
    t.SetLength(3);
    t[0] = 16;
    t[1] = 8;
    t[2] = 5;

    std::cout << "m1 = " << m1 << "\nm2 = " << m2 << '\n';

    // Stup 1: search each component and write the best three candidates to a file.
    FomSetup setup1(m1, t);
    SearchMRG<Int, Real> first;
    first.setModulus(m1);
    first.setOrder(k);
    configureBounds(first, m1);
    configureSearch(first, setup1, 8);
    first.setOutputToGenFile(true);
    first.setFilename("component1.gen");
    first.performSearch(&SearchMRG<Int, Real>::searchLoopRandom);
    report("Component 1: original search", first.bestCandidates);

    FomSetup setup2(m2, t);
    SearchMRG<Int, Real> second;
    second.setModulus(m2);
    second.setOrder(k);
    configureBounds(second, m2);
    configureSearch(second, setup2, 8);
    second.setOutputToGenFile(true);
    second.setFilename("component2.gen");
    second.performSearch(&SearchMRG<Int, Real>::searchLoopRandom);
    report("Component 2: original search", second.bestCandidates);

    // Reread the two one-component files and recompute their FOMs.
    GenFileReader<Int> reader1("component1.gen");
    FomSetup rereadSetup1(m1, t);
    SearchMRG<Int, Real> reread1;
    reread1.setGenFileReader(reader1);
    configureSearch(reread1, rereadSetup1, 3);
    reread1.performSearch(&SearchMRG<Int, Real>::searchLoopFromFile);
    report("Component 1: read from file", reread1.bestCandidates);

    GenFileReader<Int> reader2("component2.gen");
    FomSetup rereadSetup2(m2, t);
    SearchMRG<Int, Real> reread2;
    reread2.setGenFileReader(reader2);
    configureSearch(reread2, rereadSetup2, 3);
    reread2.performSearch(&SearchMRG<Int, Real>::searchLoopFromFile);
    report("Component 2: read from file", reread2.bestCandidates);

    // Step 2a: test all 3 * 3 combinations directly.
    SearchCombinedMRG<Int, Real> direct;
    direct.addComponentCandidates(first.bestCandidates);
    direct.addComponentCandidates(second.bestCandidates);
    direct.calculateModulus();

    FomSetup directSetup(direct.getModulus(), t);
    configureSearch(direct, directSetup, 9);
    direct.setOutputToGenFile(true);
    direct.setFilename("combined.gen");
    direct.performSearch(&SearchCombinedMRG<Int, Real>::searchLoopDeterministic);
    report("Combined: direct search", direct.bestCandidates);

    // Step 2b: read the retained combined generators from the .gen file.
    GenFileReader<Int> combinedReader("combined.gen");
    SearchCombinedMRG<Int, Real> fromFile;
    fromFile.setGenFileReader(combinedReader);
    fromFile.calculateModulus();

    FomSetup fileSetup(fromFile.getModulus(), t);
    configureSearch(fromFile, fileSetup, 3);
    fromFile.performSearch(&SearchCombinedMRG<Int, Real>::searchLoopFromFile);
    report("Combined: read from file", fromFile.bestCandidates);

    return 0;
}