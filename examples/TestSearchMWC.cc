#define TYPES_CODE ZD     // Int = ZZ, Real = double

#include <cmath>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/SearchMWC.h"

using namespace LatMRG;

namespace {

int64_t maxdim = 12;
int numMult = 1000000;

using SearchLoop = void (SearchMWC<Int, Real>::*)(IntLattice<Int, Real>*, int&);

void runExample(const std::string& title, long order, int64_t exponent, const IntVec& low, const IntVec& high, const NTL::Vec<int64_t>& t,
    SearchLoop loop, int numBestGen, const std::vector<int>& equalityGroups = {}, const std::vector<long>& numPowerTwo = {}, const std::vector<long>& maxPowerTwo = {})
{
    ReducerBB<Int, Real> red(maxdim);
    red.maxNodesBB = 1000000;

    WeightsUniform weights(1.0);

    // prepareCandidate() updates these bounds for each candidate modulus.
    NormaBestLat norma(-2, 1, maxdim);
    FigureOfMeritDualM<Int, Real> fom(t, weights, norma, &red, true);
    fom.setLowBound(0.0);

    SearchMWC<Int, Real> search;
    search.setOrder(order);
    search.setExponent(exponent);
    search.setLowBoundaries(low);
    search.setHighBoundaries(high);
    search.setEqualityGroups(equalityGroups);

    if (!numPowerTwo.empty()) {
        search.setNumPowerTwo(numPowerTwo);
        search.setMaxPowerTwo(maxPowerTwo);
    }

    search.setPeriodRequirement(PeriodRequirement::SafePrimeHalfPeriod);
    search.setMaxDim(maxdim);
    search.setFOM(fom);
    search.setDualLattice(true);
    search.setCalculateFOM(true);
    search.setMaxPeriod(true);
    search.setBestFOM(true);
    search.setNumBestGen(numBestGen);
    search.setMaxNumGen(numMult);
    search.setTimeLimit(7200.0);
    search.setProgress(true);
    search.setVerbose(0);

    std::cout << "\n #######################################"
              << "\n " << title
              << "\n #######################################\n";

    const std::clock_t start = std::clock();
    search.performSearch(loop);
    const std::clock_t end = std::clock();

    std::cout << "\nCPU time: "
              << static_cast<double>(end - start) / CLOCKS_PER_SEC
              << " seconds\n"
              << "Vectors examined: "
              << search.numCandidatesExamined << '\n'
              << "Safe-prime candidates: "
              << search.numMaxPeriodCandidates << '\n';

    int rank = 0;
    for (const Candidate<Int>& candidate :
         search.bestCandidates.getList()) {
        std::cout << ++rank
                  << ". FOM = " << candidate.getMerit()
                  << ", m = " << candidate.getModulus()
                  << ", aa = " << candidate.getAa() << '\n';
    }
}

} // namespace

int main()
{
    // Original order-1 example.
    NTL::Vec<int64_t> t1;
    t1.SetLength(5);
    t1[0] = 12;
    t1[1] = 12;
    t1[2] = 10;
    t1[3] = 10;
    t1[4] = 8;

    IntVec low1, high1;
    low1.SetLength(2);
    high1.SetLength(2);
    low1[0] = high1[0] = Int(-1);
    low1[1] = Int(0);
    high1[1] = NTL::power(Int(2), 62) - 1;

    runExample("Random MWCs of order 1 with a_0 = -1", 1, 64, low1, high1, t1, &SearchMWC<Int, Real>::searchLoopRandom, 1);

    // Original order-2 example: a_1 = 0; draw a_2.
    NTL::Vec<int64_t> t2 = t1;
    t2[1] = 0;

    IntVec low2, high2;
    low2.SetLength(3);
    high2.SetLength(3);
    low2[0] = high2[0] = Int(-1);
    low2[1] = high2[1] = Int(0);
    low2[2] = Int(0);
    high2[2] = NTL::power(Int(2), 62) - 1;

    runExample("Random MWCs of order 2 with a_0 = -1 and a_1 = 0", 2, 64, low2, high2, t2, &SearchMWC<Int, Real>::searchLoopRandom, 1);

    // a) Random search with a_1 = a_2 and a_0 = -1.
    IntVec lowA, highA;
    lowA.SetLength(3);
    highA.SetLength(3);
    lowA[0] = highA[0] = Int(-1);
    lowA[1] = lowA[2] = Int(0);
    highA[1] = highA[2] = NTL::power(Int(2), 62) - 1;

    runExample("Random MWCs of order 2 with a_0 = -1 and a_1 = a_2", 2, 16, lowA, highA, t2, &SearchMWC<Int, Real>::searchLoopRandom, 3, {1, 1});

    // b) Exhaustive search: 1000 choices of a_1 times
    //    1000 choices of a_2, with a_0 = -1.
    IntVec lowB, highB;
    lowB.SetLength(3);
    highB.SetLength(3);
    lowB[0] = highB[0] = Int(-1);
    lowB[1] = lowB[2] = Int(0);
    highB[1] = highB[2] = Int(999);

    runExample("Exhaustive MWCs of order 2: 1000 x 1000, a_0 = -1", 2, 16, lowB, highB, t2, &SearchMWC<Int, Real>::searchLoopExhaustive, 3);

    // c) Random signed sums of powers of two.
    IntVec lowC, highC;
    lowC.SetLength(3);
    highC.SetLength(3);
    lowC[0] = highC[0] = Int(-1);
    lowC[1] = lowC[2] = Int(-1023);
    highC[1] = highC[2] = Int(1023);

    runExample("Random power-of-two MWCs of order 2 with a_0 = -1", 2, 16, lowC, highC, t2, &SearchMWC<Int, Real>::searchLoopRandomPowerOfTwo, 3, {}, {1, 3, 3}, {0, 18, 18});

    return 0;
}