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

#ifndef LATMRG_CONFIGSEARCHCOMBINEDMWC_H
#define LATMRG_CONFIGSEARCHCOMBINEDMWC_H

#include <vector>
#include <utility>
#include <cstdint>

#include "latmrg/ConfigSearchBase.h"
#include "latticetester/Util.h"

namespace LatMRG {

/**
 * Configuration for an MWC search.
 */
template<typename Int, typename Real>
struct ConfigSearchCombinedMWC : public ConfigSearch<Int, Real>
{
    /**
     * Configuration for one MWC component.
     */
    struct ComponentMWC
    {
        // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
        long powMod;   // Value of e such that the modulus is b = 2^e
        long order;    // The order k
        IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j
        IntVec highBoundaries; // Default values are 1 and b-1
        vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
        vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
        vector<int64_t> randomBits; // For random searches, the nonzero coefficients `a_j` are generated randomly with `e_j` random bits
        int64_t numaj; // number of positive coefficients 
        vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal

        long getOrder() const
        {
            return order;
        }

        Int getNoMultipliers() const
        {
            Int total(1);
            for (long i = 1; i <= order; i++) {
                total *= highBoundaries[i] - lowBoundaries[i] + 1;
            }
            return total;
        }

        Int getLowBoundary(int i) const
        {
            return lowBoundaries[i];
        }

        Int getHighBoundary(int i) const
        {
            return highBoundaries[i];
        }

        int64_t getRandomBits(int i) const
        {
            return randomBits[i];
        }

        long getPowMod() const
        {
            return powMod;
        }
    };

    std::vector<ComponentMWC> genComponents;

    void createComponents()
    {
        genComponents.resize(this->numComp);
    }

    const ComponentMWC & component(int i) const
    {
        return genComponents[i];
    }

    ComponentMWC & component(int i)
    {
        return genComponents[i];
    }

    Int getModulus() const override
    {
        Int mod(1);
        for (const auto& comp : genComponents) {
            mod *= comp.getModulus();
        }
        return mod;
    }
};

} // namespace LatMRG

#endif // LATMRG_CONFIGSEARCHCOMBINEDMWC_H