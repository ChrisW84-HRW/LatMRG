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

#define LATMRG_SEARCH

#include <iostream>
#include <cstdint>
#include <memory>
#include <NTL/ZZ.h>

#include "latmrg/IntFactorization.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/ConfigSearchBase.h"
#include "latmrg/ConfigSearchMRG.h"
#include "latmrg/ConfigSearchMWC.h"
#include "latmrg/ConfigSearchCombinedMRG.h"
#include "latmrg/ConfigSearchCombinedMWC.h"
#include "latmrg/Search.h"
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
  
  cout << "The purpose is to test Seek.h.\n";

  // Variables necessary for MRG calculations
  NTL::Vec<int64_t> t; // The t-vector for the FOM.
  t.SetLength(3);
  t[0] = 8;    // We look at successive coordinates in up to t[0] dimensions.
  t[1] = 5;    // Then pairs and triples up to coord. 5.
  t[2] = 5;
  const NTL::ZZ m = NTL::to_ZZ("9223372036854773561");  // modulus  
  NTL::Vec<NTL::ZZ> b; // low boundaries for the multipliers
  b.SetLength(3);
  b[1] = 1145902849652723;  
  b[2] = 5189485190151516; 
  NTL::Vec<NTL::ZZ> c;
  c.SetLength(3); // upper boundaries for the multipliers
  c[1] = 1145902849652955;
  c[2] = 5189485190151758;
  
  NTL::Vec<NTL::ZZ> h;
  h.SetLength(3); // sizes for the box search
  h[1] = 2;
  h[2] = 3;
  const int64_t maxdim(16);
  ReducerBB<Int, Real> red(maxdim);
  NormaBestLat norma(log(m), 1, 16);
  WeightsUniform weights(1.0);
  
  // Define all the necessary details of the configuration  
  ConfigSearchMRG<Int, Real> conf;
  conf.setFOM(t, &norma, &red, &weights);  
  conf.setCalculateFOM(true);
  conf.setMaxDim(maxdim); 
  conf.setPermax(true);
  conf.setOutputToGenFile(false);
  conf.setFilename("test.gen");  
  conf.setMaxGen(100);  
  // Set paramters for box search
  conf.setSubRegionsSize(h);
  conf.setNumRandomRegions(10);
  
  conf.setModulus(m);
  conf.setLowBoundaries(b);
  conf.setHighBoundaries(c);
  conf.setOrder(b.length() - 1);
  // Perform the actual search
  SearchMRG<Int, Real> searcher(conf);
  searcher.performSearch(&SearchMRG<Int, Real>::nextGeneratorRandomBox);
  searcher.printResults();
    
  return 0;
}
