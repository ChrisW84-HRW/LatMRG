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

#ifndef LATMRG_CONFIGSEARCHMRG_H
#define LATMRG_CONFIGSEARCHMRG_H

#include <vector>
#include <utility>
#include <string>

#include "latmrg/ConfigSearchBase.h"
#include "latmrg/EnumTypes.h"
#include "latticetester/Util.h"

namespace LatMRG {

/**
 * Configuration for an MRG search.
 */
template<typename Int, typename Real>
struct ConfigSearchMRG : public ConfigSearch<Int, Real>
{
    ConfigSearchMRG() : ConfigSearch<Int, Real>(MRG) {}

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
    long numRandomRegions = 1;   // for search of random boxes (cooresponds to N_r in the guide)
    IntVec subregionSize;        // the allowed sizes of the boxes h_i, indexed from 1 to k
                       
    // Get and set methods
    Int getModulus() const {return modulus;}
    void setModulus(const Int& value) {modulus = value;}
    Int getBase() const {return base;}
    void setBase(const Int& value) {base = value;}
    Int getExponent() const {return exponent;}
    void setExponent(const Int& value) {exponent = value;}
    Int getRest() const {return rest;}
    void setRest(const Int& value) {rest = value;}
    long getOrder() const {return order;}
    void setOrder(long value) {order = value;}
    const IntVec& getLowBoundaries() const {return lowBoundaries;}
    void setLowBoundaries(const IntVec& values) {lowBoundaries = values;}
    Int getLowBoundary(int i) const {return lowBoundaries[i];}
    void setLowBoundary(int i, const Int& value) {lowBoundaries[i] = value;}
    const IntVec& getHighBoundaries() const {return highBoundaries;}
    void setHighBoundaries(const IntVec& values) {highBoundaries = values;}
    Int getHighBoundary(int i) const {return highBoundaries[i];}
    void setHighBoundary(int i, const Int& value) {highBoundaries[i] = value;}
    bool getAppFac() const {return appFac;}
    void setAppFac(bool value) {appFac = value;}
    const vector<long>& getNumPowerTwo() const {return numPowerTwo;}
    void setNumPowerTwo(const vector<long>& values) {numPowerTwo = values;}
    const vector<long>& getMaxPowerTwo() const {return maxPowerTwo;}
    void setMaxPowerTwo(const vector<long>& values) {maxPowerTwo = values;}
    const vector<vector<pair<long,long>>>& getEqualCoeffs() const {return equalCoeffs;}
    void setEqualCoeffs(const vector<vector<pair<long,long>>>& values) {equalCoeffs = values;}
    DecompType getHowFactorh() const {return howFactorh;}
    void setHowFactorh(DecompType value) {howFactorh = value;}
    const IntVec& getFactorsh() const {return factorsh;}
    void setFactorsh(const IntVec& values) {factorsh = values;}
    DecompType getHowFactor() const {return howFactor;}
    void setHowFactor(DecompType value) {howFactor = value;}
    const IntVec& getFactorsr() const {return factorsr;}
    void setFactorsr(const IntVec& values) {factorsr = values;}
    bool getPermaxPow2() const {return permaxPow2;}
    void setPermaxPow2(bool value) {permaxPow2 = value;}
    Int getNoMultipliers() const
    {
        Int total(1);
        for (long i = 1; i <= order; i++) {
            total *= highBoundaries[i] - lowBoundaries[i] + 1;
        }
        return total;
    }
    long getNumRandomRegions() {return numRandomRegions;}
    void setNumRandomRegions(long value) {numRandomRegions = value;}    
    Int getSubregionSize(int i) const {return subregionSize[i];}
    void setSubregionSize(IntVec value) {subregionSize = value;}
    // There is only a single component
    Int getComponentModulus(int k = 0) const override {return getModulus();}
    long getComponentOrder(int k = 0) const override {return getOrder();}
    
};

} // namespace LatMRG

#endif // LATMRG_CONFIGSEARCHMRG_H