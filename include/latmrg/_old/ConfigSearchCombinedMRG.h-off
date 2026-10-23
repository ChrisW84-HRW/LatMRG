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

#ifndef LATMRG_CONFIGSEARCHCOMBINEDMRG_H
#define LATMRG_CONFIGSEARCHCOMBINEDMRG_H

#include <vector>
#include <utility>
#include <string>

#include "latmrg/ConfigSearchBase.h"
#include "latmrg/EnumTypes.h"
#include "latticetester/Util.h"

namespace LatMRG {

/**
 * Configuration for an MRG seek.
 */
template<typename Int, typename Real>
struct ConfigSearchCombinedMRG : public ConfigSearch<Int, Real>
{
    ConfigSearchCombinedMRG() : ConfigSearch<Int, Real>(MRG) {}

    /**
     * Configuration + methods for one MRG component.
     */
    struct ComponentMRG
    {
        Int modulus;   // The modulus m = b^e + r
        Int base;      // The base b
        Int exponent;  // The exponent e
        Int rest;      // The rest r
        long order;    // The order k
        IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j 
        IntVec highBoundaries; // Default values are 1 and m-1 
        bool appFac = false;   // true iff we enforce the approximate factoring condition (default = false)
        vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
        vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
        vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
        DecompType howFactorh; // Indicates how to factorize h=(m-1)/2
        IntVec factorsh;      // The factors of h, repeated according to multiplicity, when known
        DecompType howFactor;  // Indicates how to factorize r=(m^k-1)/(m-1)
        IntVec factorsr;      // The factors of r, repeated according to multiplicity, when known
        bool permaxPow2;   // True iff m is a power of 2, b=2, r=0, k=1, and we want max period for the LCG
                       // Default is false
        
        Int getModulus() const
        {
            return modulus;
        }

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
    };

    std::vector<ComponentMRG> genComponents;

    // This method resizes the vector of components to the right size
    void createComponents()
    {
        genComponents.resize(this->numComp);
    }

    // Gets the modulus of an MRG with multiple components
    Int getModulus() const override
    {
        Int mod(1);
        for (const auto& component : genComponents) {
            mod *= component.getModulus();
        }
        return mod;
    }

    // Gets a single component
    ComponentMRG & component(int i)
    {
        return genComponents[i];
    }

    // Gets the modulus of a single component
    Int getComponentModulus(int k) const override
    {
        return genComponents[k].getModulus();
    }

    // Gets the order of a single component
    long getComponentOrder(int k) const override
    {
        return genComponents[k].getOrder();
    }
};

} // namespace LatMRG

#endif // LATMRG_CONFIGSEARCHCOMBINEDMRG_H