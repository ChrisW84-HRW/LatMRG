// This file is part of LatMRG.
//
// Copyright (C) 2012-2026  The LatMRG authors, under the occasional supervision
// of Pierre L'Ecuyer at Universit? de Montr?al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <iostream>
#include <cstdint>
#include <memory>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"


#include "latmrg/SearchBase.h"
#include "latmrg/SearchMRG.h"
#include "latmrg/SearchMWC.h"
#include "latmrg/SearchCombinedMRG.h"
#include "latmrg/SearchCombinedMWC.h"

/**
 * This example implements a search for MRGs. It is far from being
 * the final version but should only be seen as a starting point.
 * The most important lines are the last lines of the code where the
 * actual search is performed.
 */

using namespace LatMRG;

// Everything which is needed for the configuration is set first
// Only the last two lines are really performing the seek.
int main() {
  
  cout << "The purpose is to test the search for good generators.\n";

  int k = 2;

  // Variables necessary for MRG calculations
  NTL::Vec<int64_t> t; // The t-vector for the FOM.
  t.SetLength(k+1);
  t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
  t[1] = 5;    // Then pairs and triples up to coord. 5.
  t[2] = 5;
  const NTL::ZZ m = NTL::to_ZZ("9223372036854773561");  // modulus  
  NTL::Vec<NTL::ZZ> b; // low boundaries for the multipliers
  b.SetLength(k+1);
  b[1] = 1145902849652723;  
  b[2] = 5189485190151516; 
  NTL::Vec<NTL::ZZ> c;
  c.SetLength(k+1); // upper boundaries for the multipliers
  c[1] = 1245902849752955;
  c[2] = 5289485190251758;
  
  NTL::Vec<NTL::ZZ> h;
  h.SetLength(k+1); // sizes for the box search
  h[1] = 10;
  h[2] = 10;
  const int64_t maxdim(16);
  ReducerBB<Int, Real> red(maxdim);
  NormaBestLat norma(-log(m), k, maxdim);
  WeightsUniform weights(1.0);

  
  
  // Define all the necessary details of the configuration  
  SearchMRG<Int, Real> searcher;
  searcher.setFOM(t, &norma, &red, &weights);  
  searcher.setDualLattice(true);
  searcher.setCalculateFOM(true);
  searcher.setMaxDim(maxdim); 
  searcher.setPermax(true);
  searcher.setOutputToGenFile(false);
  searcher.setFilename("test.gen");  
  searcher.setMaxGen(100000);  
  // Set paramters for box search
  searcher.setSubregionSize(h);
  searcher.setNumRandomRegions(10000);
  searcher.setNumRandomAk(100000);
  searcher.setNumRandomPointsPerAk(100);
  searcher.setCountBadAk(false);
  
  searcher.setModulus(m);
  searcher.setLowBoundaries(b);
  searcher.setHighBoundaries(c);
  searcher.setOrder(b.length() - 1);

  // Perform the actual search  
  
  std::cout << "\n \n Pure random search\n";
  std::clock_t startRandom = std::clock();
  searcher.performSearch(&SearchMRG<Int, Real>::nextGeneratorRandom);
  std::clock_t endRandom = std::clock();
  double cpuTimeRandom = double(endRandom - startRandom) / CLOCKS_PER_SEC;
  std::cout << "Search CPU time: " << cpuTimeRandom << " seconds\n";
  std::cout << "Best merit: " << searcher.getBestMerit() << "\n";
  
  
  std::cout << "\n \n Random search with Ak fixed\n";
  std::clock_t startAk = std::clock();
  searcher.performSearch(&SearchMRG<Int, Real>::nextGeneratorRandomAkFixed);
  std::clock_t endAk = std::clock();
  double cpuTimeAk = double(endAk - startAk) / CLOCKS_PER_SEC;
  std::cout << "Search CPU time: " << cpuTimeAk << " seconds\n";
  std::cout << "Best merit: " << searcher.getBestMerit() << "\n";

  
  std::cout << "\n \n Random search with Ak fixed + boxes \n";
  std::clock_t startBox = std::clock();
  searcher.performSearch(&SearchMRG<Int, Real>::nextGeneratorRandomAkFixedBox);
  std::clock_t endBox = std::clock();
  double cpuTimeBox = double(endBox - startBox) / CLOCKS_PER_SEC;
  std::cout << "Search CPU time: " << cpuTimeBox << " seconds\n";
  std::cout << "Best merit: " << searcher.getBestMerit() << "\n";

  
  std::cout << "\n \n Random search with boxes \n";
  std::clock_t startBox2 = std::clock();
  searcher.performSearch(&SearchMRG<Int, Real>::nextGeneratorRandomBox);
  std::clock_t endBox2 = std::clock();
  double cpuTimeBox2 = double(endBox2 - startBox2) / CLOCKS_PER_SEC;
  std::cout << "Search CPU time: " << cpuTimeBox2 << " seconds\n";
  std::cout << "Best merit: " << searcher.getBestMerit() << "\n";
  // searcher.printResults();
    
  
  
  
  return 0;
}
