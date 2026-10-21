// This file is part of LatMRG.
//
// Copyright (C) 2012-2022  The LatMRG authors, under the occasional supervision
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

#ifndef LATMRG_SEARCHMRG_H
#define LATMRG_SEARCHMRG_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>

#include "latmrg/PrimesFinder.h"
#include "latmrg/ConfigSearchBase.h"
#include "latmrg/ConfigSearchMRG.h"
#include "latmrg/FigureOfMeritData.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsUniform.h"
#include "latmrg/Search.h"

#include "latmrg/LCGLattice.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"
#include "latmrg/MWCLattice.h"

namespace LatMRG {

  /**
  * \class SearchMRG
  *
  * This class implements the search procedure for MRG lattices. 
  * It derives from the abstract class Search and provides the 
  * lattice-specific implementation of the next generator methods 
  * required during the search. This class only defines how MRG 
  * candidate generators are created. 
  *
  * The method nextGenerator() performs an exhaustive enumeration of the
  * admissible MRG generators defined by the configuration. 
  *
  * The method nextGeneratorRandom() generates admissible MRG generators by
  * randomly selecting the recurrence coefficients within the ranges specified
  * in the configuration. 
  *
  * The constructor is mainly inherited from Search. It initializes the search using a
  * ConfigSearch object, which specifies the modulus, recurrence order,
  * admissible coefficient ranges, search limits, and the figure-of-merit
  * configuration.
  */

    template<typename Int, typename Real>
    class SearchMRG: public Search<MRGLattice<Int, Real>> {

        public:

        SearchMRG(ConfigSearchMRG<Int, Real>& config)        
        : Search<MRGLattice<Int, Real>>(config),
          mrgConf(config),
          mrg(config.getModulus(), config.getOrder())
        {
        }

        MRGLattice<Int, Real>* nextGenerator() override;

        MRGLattice<Int, Real>* nextGeneratorRandom() override;

        MRGLattice<Int, Real>* nextGeneratorRandomBox();

        bool checkMaxPeriod(MRGLattice<Int, Real>& lat) override;

        private:

        ConfigSearchMRG<Int, Real>& mrgConf;

        /**
        * Object for creating an MRG lattice
        */
        MRGComponent<Int> mrg;

        /*
        * Helper variables for box search -- maybe need to be initalized
        */
        long currentRandomRegion;
        Int currentSubIndex;
        Int currentSubregionTotal;
        IntVec alpha;
        bool hasLastMaxPeriodTest;
        NTL::Vec<NTL::ZZ> lastAA;
        bool lastMaxPeriodResult;

    };

//============================================================================
// Implementation


 /**
  * nextGenerator() performs an exhaustive enumeration of the admissible MRG 
  * generators defined by the configuration. The coefficients of the recurrence 
  * are obtained by decoding the current generator index into the corresponding 
  * multiplier values within the prescribed ranges. Each admissible set of 
  * coefficients is used to construct an MRG lattice that is subsequently analyzed 
  * by the search procedure defined in the base class.
  * If the search is restricted to generators with maximal period, then nextGenerator()
  * gives the next generator with this property (otherwise they are skipped).
  */
  template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGenerator()
  {
    const int k = mrgConf.getOrder();    
    Int range, val;
    Int tmp;
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);

    while (this->currentGen < this->conf.max_gen &&
           this->currentGen < mrgConf.getNoMultipliers())
    {
        // Set the vector of multipliers for the next generator
        tmp = NTL::to_ZZ(this->currentGen);
        for (int i = 1; i <= k; i++) {
            range = mrgConf.getHighBoundary(i) - mrgConf.getLowBoundary(i) + 1;
            val = tmp % range;
            tmp /= range;
            aa[i] = mrgConf.getLowBoundary(i) + val;
        }
        ++this->currentGen;
        
        // Preparation for being able to check primitivity -- only necessary once for MRGs but important to update for MWCs
        setModulusIntP<Int>(mrgConf.getModulus());
        if (this->conf.onlyMaxPeriod() && !mrg.maxPeriod(aa)) {continue;}

        return new MRGLattice<Int, Real>(mrgConf.getModulus(), aa, this->conf.maxdim);
    }
    return nullptr;
  }
  //===========================================================================

  /**
  * nextGeneratorRandom() generates admissible MRG generators by randomly selecting 
  * the recurrence coefficients within the ranges specified in the configuration. 
  * It provides the random search analagoue of nextGenerator() and can be used by 
  * performSearch() when a stochastic search. 
  * If the search is restricted to generators with maximal period, then 
  * nextGeneratorRandom() only yields generators having this property.
  */
  template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandom()
  {
    const int k = mrgConf.getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);

    while (this->currentGen < this->conf.max_gen)
    {
        for (int i = 1; i <= k; i++) aa[i] = randInt(mrgConf.getLowBoundary(i),  mrgConf.getHighBoundary(i));
        ++this->currentGen;
        // Preparation for being able to check primitivity -- not really necessary for MRGs but important for MWCs
        setModulusIntP<Int>(mrgConf.getModulus());
        if (this->conf.onlyMaxPeriod() && !mrg.maxPeriod(aa))  continue;

        return new MRGLattice<Int, Real>(mrgConf.getModulus(), aa, this->conf.maxdim);
    }
    return nullptr;
  }

  //===========================================================================

  /**
  * Description.
  */
  template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomBox()
  {
    const int k = mrgConf.getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);

    while (this->currentGen < this->conf.max_gen)
    {
        for (int i = 1; i <= k; i++) aa[i] = randInt(mrgConf.getLowBoundary(i),  mrgConf.getHighBoundary(i));
        ++this->currentGen;
        // Preparation for being able to check primitivity -- not really necessary for MRGs but important for MWCs
        setModulusIntP<Int>(mrgConf.getModulus());
        if (this->conf.onlyMaxPeriod() && !mrg.maxPeriod(aa))  continue;

        return new MRGLattice<Int, Real>(mrgConf.getModulus(), aa, this->conf.maxdim);
    }
    return nullptr;
  }

  /**
  * This method checks whether the generated lattice satisfies the
  * maximal-period requirement for an MRG.
  * Probably this function will not be used since it requires to pass a lattice object which has not been creater.
  * Delete it if it is not necessary.
  */
  template<typename Int, typename Real> bool SearchMRG<Int, Real>::checkMaxPeriod(MRGLattice<Int, Real>& lat)
  {
     return mrg.maxPeriod(lat.getaa());
  }

}// End namespace LatMRG
#endif
