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

#ifndef LATMRG_CONFIGSEARCHMWC_H
#define LATMRG_CONFIGSEARCHMWC_H

#include <vector>
#include <utility>
#include <cstdint>

#include "latmrg/ConfigSearchBase.h"
#include "latticetester/Util.h"

namespace LatMRG {

/**
 * Configuration for an MWC search.
 */
template<typename Int, typename Real> struct ConfigSearchMWC : public ConfigSearch<Int, Real>
{
    Int modulus;   // The modulus m of the equivalent LCG (to be computed)
    long powMod;   // Value of e such that the modulus is b = 2^e
    long order;    // The order k
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j
    IntVec highBoundaries; // Default values are 1 and b-1
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<int64_t> randomBits; // For random searches, the nonzero coefficients `a_j` are generated randomly with `e_j` random bits
    int64_t numaj; // number of positive coefficients 
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal

    // Get and set methods
    Int getModulus() const { return modulus; }
    void setModulus(const Int& value) { modulus = value;}
    long getPowMod() const { return powMod;}
    void setPowMod(long value) { powMod = value;}
    long getOrder() const { return order;}
    void setOrder(long value) { order = value; }
    const IntVec& getLowBoundaries() const {return lowBoundaries;}
    void setLowBoundaries(const IntVec& values) {lowBoundaries = values;}
    Int getLowBoundary(int i) const {return lowBoundaries[i];}
    void setLowBoundary(int i, const Int& value) {lowBoundaries[i] = value;}
    const IntVec& getHighBoundaries() const {return highBoundaries;}
    void setHighBoundaries(const IntVec& values) {highBoundaries = values;}
    Int getHighBoundary(int i) const {return highBoundaries[i];}
    void setHighBoundary(int i, const Int& value) {highBoundaries[i] = value;}
    const vector<long>& getNumPowerTwo() const {return numPowerTwo;}
    void setNumPowerTwo(const vector<long>& values) {numPowerTwo = values;}
    const vector<long>& getMaxPowerTwo() const {return maxPowerTwo;}
    void setMaxPowerTwo(const vector<long>& values) {maxPowerTwo = values;}
    const vector<int64_t>& getRandomBits() const {return randomBits;}
    void setRandomBits(const vector<int64_t>& values) {randomBits = values;}
    int64_t getRandomBits(int i) const {return randomBits[i];}
    void setRandomBits(int i, int64_t value) {randomBits[i] = value;}
    int64_t getNumaj() const {return numaj;}
    void setNumaj(int64_t value) {numaj = value;}
    const vector<vector<pair<long,long>>>& getEqualCoeffs() const {return equalCoeffs;}
    void setEqualCoeffs(const vector<vector<pair<long,long>>>& values) {equalCoeffs = values; }
    Int getNoMultipliers() const
    {
        Int total(1);
        for (long i = 1; i <= order; i++) {
            total *= highBoundaries[i] - lowBoundaries[i] + 1;
        }
        return total;
    }
    // There is only a single component
    Int getComponentModulus(int k = 0) const override {return getModulus();}
    long getComponentOrder(int k = 0) const override {return getOrder();}
};

} // namespace LatMRG

#endif // LATMRG_CONFIGSEARCHMWC_H