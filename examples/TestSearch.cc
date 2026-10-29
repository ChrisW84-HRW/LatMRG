// This file is part of LatMRG.
//
// Copyright (C) 2012-2026 The LatMRG authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0

#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <memory>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"

#include "latmrg/SearchBase.h"
#include "latmrg/SearchMRG.h"


using namespace LatMRG;

/**
 * Run one MRG search strategy and report its CPU time and best merit.
 *
 * The FoM lower bound is reset because the same FoM object is reused
 * for several independent searches.
 */
template<typename Int, typename Real>
void runSearch(const char* name, SearchMRG<Int, Real>& searcher,  FigureOfMeritM<Int, Real>& fom, bool (SearchMRG<Int, Real>::*generator)())
  {
    std::cout << "\n" << name << "\n";
    fom.setLowBound(0.0);
    const std::clock_t start = std::clock();
    searcher.performSearch(generator);
    const std::clock_t end = std::clock();
    const double cpuSeconds =  static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Search CPU time: " << cpuSeconds << " seconds\n";
    std::cout << "Best merit: " << searcher.getBestMerit() << "\n";
  }

/**
 * Compare four random-search strategies for MRGs of order 2.
 *
 * All methods use the same coefficient bounds, figure of merit,
 * and generator limit.
 */
int main()
{
    constexpr int k = 2;
    constexpr int64_t maxdim = 16;

    const NTL::ZZ m = NTL::to_ZZ("9223372036854773561");

    // Figure-of-merit parameters.
    NTL::Vec<int64_t> t;
    t.SetLength(k + 1);
    t[0] = 8;  // Successive coordinates, up to dimension 8.
    t[1] = 5;  // Pairs drawn from the first 5 coordinates.
    t[2] = 5;  // Triples drawn from the first 5 coordinates.

    ReducerBB<Int, Real> red(maxdim);
    NormaBestLat norma(-log(m), k, maxdim);
    WeightsUniform weights(1.0);
    FigureOfMeritDualM<Int, Real> fom(t, weights, norma, &red, true);

    // Inclusive bounds for coefficients a_1 and a_2.
    NTL::Vec<NTL::ZZ> lower;
    NTL::Vec<NTL::ZZ> upper;
    lower.SetLength(k + 1);
    upper.SetLength(k + 1);

    lower[1] = 1145902849652723;
    lower[2] = 5189485190151516;
    upper[1] = 1245902849752955;
    upper[2] = 5289485190251758;

    // Side lengths for box-based searches.
    NTL::Vec<NTL::ZZ> boxSize;
    boxSize.SetLength(k + 1);
    boxSize[1] = 10;
    boxSize[2] = 10;

    // Configure the search settings shared by all four strategies.
    SearchMRG<Int, Real> searcher;
    searcher.setModulus(m);
    searcher.setOrder(k);
    searcher.setLowBoundaries(lower);
    searcher.setHighBoundaries(upper);
    searcher.setMaxDim(maxdim);

    searcher.setFOM(fom);
    searcher.setDualLattice(true);
    searcher.setCalculateFOM(true);
    searcher.setMaxPeriod(true);
    searcher.setMaxNumGen(100000);

    searcher.setOutputToGenFile(false);
    searcher.setFilename("test.gen");
    searcher.setProgress(true);
    searcher.setVerbose(1);

    // Configure settings used by the box and fixed-a_k strategies.
    searcher.setSubregionSize(boxSize);
    searcher.setNumRandomRegions(10000);
    searcher.setNumRandomAk(100000);
    searcher.setNumRandomPointsPerAk(100);

    // Try the different standard search stategies
    runSearch("Pure random search", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandom);
    runSearch("Random search by a_k", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomByAk);
    runSearch("Random search by a_k and boxes", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomByAkBox);
    runSearch("Random search with boxes", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomBox);

    // Try search strategies with approximate factoring
    runSearch("Pure random search with approximate factoring", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomApproxFactoring);
    runSearch("Pure random search with approximate factoring + by Ak", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomApproxFactoringByAk);

    // Try search strategies with sums of powers of two
    searcher.setNumPowerTwo({2, 2});  // At most two terms for each coefficient
    searcher.setMaxPowerTwo({62, 62}); // Exponents from 0 through 62
    runSearch("Pure random search with two powers of two", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomPowersOfTwo);
    runSearch("Pure random search with two powers of two + by Ak", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandomPowersOfTwoByAk);
    
    // Impose that constraint that two coefficients are equal
    searcher.setEqualityGroups({1, 1});
    lower[2] = lower[1];
    upper[2] = upper[1];
    searcher.setLowBoundaries(lower);
    searcher.setHighBoundaries(upper);
    searcher.setEqualityGroups({1, 1});
    runSearch("Pure random search with first and second coefficient equal", searcher, fom, &SearchMRG<Int, Real>::nextGeneratorRandom);

    return 0;
}