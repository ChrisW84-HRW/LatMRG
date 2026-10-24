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

#ifndef LATMRG_SEARCHMRG_H
#define LATMRG_SEARCHMRG_H

#include <vector>
#include <utility>
#include <string>
#include <memory>

#include "latmrg/SearchBase.h"
#include "latmrg/EnumTypes.h"
#include "latticetester/Util.h"

namespace LatMRG {

/**
 * SearchMRG contain search procedures for MRG lattices.
 *
 * It derives from SearchBase with MRGLattice as the concrete lattice type and
 * stores the configuration data specific to MRG generators, such as the modulus,
 * recurrence order, coefficient bounds, and optional restrictions used by the
 * different search methods.
 *
 * The main feature of the class is that it  provides generator routines for different 
 * search methods such as exhaustive search and different types of random search.
 */
template<typename Int, typename Real>
struct SearchMRG : public SearchBase<Int, Real, MRGLattice<Int, Real>>
{
    SearchMRG() : SearchBase<Int, Real, MRGLattice<Int, Real>>(MRG) {}

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
    int numRandomPointsPerAk = 1;
    int numRandomAk = 1;

    bool countBadAk = false;
                       
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
    void setNumRandomPointsPerAk(int value) {numRandomPointsPerAk = value;}
    int getNumRandomPointsPerAk() const {return numRandomPointsPerAk;}
    void setNumRandomAk(int value) {numRandomAk = value;}
    int getNumRandomAk() const {return numRandomAk;}
    bool getCountBadAk() const {return countBadAk;};
    void setCountBadAk(bool value) {countBadAk = value;}
    // There is only a single component
    Int getComponentModulus(int k = 0) const override {return getModulus();}
    long getComponentOrder(int k = 0) const override {return getOrder();}
    
    MRGLattice<Int, Real>* nextGeneratorExhaustive() override;

    MRGLattice<Int, Real>* nextGeneratorRandom() override;

    MRGLattice<Int, Real>* nextGeneratorRandomBox();
    
    MRGLattice<Int, Real>* nextGeneratorRandomAkFixed();

    MRGLattice<Int, Real>* nextGeneratorRandomAkFixedBox();
    bool checkMaxPeriod(MRGLattice<Int, Real>& lat) override;

    void resetGeneratorState() override{
        currentRandomRegion = 0;
        currentSubIndex = 0;
        hasLastAk = false;
        lastMaxPeriodResult = false;
        mrg.reset();
        currentAkTrial = 0;
        currentAkSample = 0;
        hasCurrentAk = false;
        WorkingLattice = std::make_unique<MRGLattice<Int, Real>>(getModulus(), getLowBoundaries(), this->maxdim);
        aa.SetLength(getOrder() + 1);
        alpha.SetLength(getOrder() + 1);
    }

    private:

    /**
    * Object for creating an MRG lattice
    */
    std::unique_ptr<MRGComponent<Int>> mrg;

    /**
     * Helper function 
    */
    MRGComponent<Int>& getMRGComponent() {
        if (!mrg) {
            mrg.reset(new MRGComponent<Int>(getModulus(), getOrder()));
        }
        return *mrg;
    }

    /*
    * Helper Lattice
    */
    std::unique_ptr<MRGLattice<Int, Real>> WorkingLattice;

    /*
    * Vector for creating the multipliers
    */
    NTL::Vec<NTL::ZZ> aa;

    /*
    * Helper variables for box search
    */
    long currentRandomRegion = 0;
    Int currentSubIndex = Int(0);
    IntVec alpha;
    bool hasLastAk = false;
    NTL::ZZ lastAk;
    bool lastMaxPeriodResult = false;

    long currentAkTrial = 0;
    long currentAkSample = 0;
    bool hasCurrentAk = false;
    NTL::ZZ currentAk;
    
};

//============================================================================
// Implementation


 /**
  * nextGeneratorExhaustive() performs an exhaustive enumeration of the admissible MRG 
  * generators defined by the configuration. The coefficients of the recurrence 
  * are obtained by decoding the current generator index into the corresponding 
  * multiplier values within the prescribed ranges. Each admissible set of 
  * coefficients is used to construct an MRG lattice that is subsequently analyzed 
  * by the search procedure defined in the base class.
  * If the search is restricted to generators with maximal period, then nextGeneratorExhaustive()
  * gives the next generator with this property (otherwise they are skipped).
  */
  template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorExhaustive()
  {    
    Int range, val;
    Int tmp;
    while (this->currentGen < this->getMaxGen() &&
           this->currentGen < getNoMultipliers())
    {
        // Set the vector of multipliers for the next generator
        tmp = NTL::to_ZZ(this->currentGen);
        for (int i = 1; i <= getOrder(); i++) {
            range = getHighBoundary(i) - getLowBoundary(i) + 1;
            val = tmp % range;
            tmp /= range;
            aa[i] = getLowBoundary(i) + val;
        }
        ++this->currentGen;
        
        if (this->onlyMaxPeriod() && !getMRGComponent().maxPeriod(aa)) {
            if (!countBadAk) --this->currentGen;
            continue;
        }
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }
    return nullptr;
  }
  //===========================================================================

  /**
  * nextGeneratorRandom() generates admissible MRG generators by randomly selecting 
  * the recurrence coefficients within the ranges specified in the configuration. 
  * It provides the random search analagoue of nextGeneratorExhaustive() and can be used by 
  * performSearch() when a stochastic search. 
  * If the search is restricted to generators with maximal period, then 
  * nextGeneratorRandom() only yields generators having this property.
  */
  template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandom()
  {
    
    while (this->currentGen < this->getMaxGen())
    {
        for (int i = 1; i <= getOrder(); i++) aa[i] = randInt(getLowBoundary(i),  getHighBoundary(i));
        ++this->currentGen;
        if (this->onlyMaxPeriod() && !getMRGComponent().maxPeriod(aa)) {
            if (!countBadAk) --this->currentGen;
            continue;
        } 
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }
    return nullptr;
  }

  //===========================================================================

  /**
   * nextGeneratorRandomBox() performs a random-box search over the admissible
   * MRG coefficient region. It first selects random subregions, or boxes, within
   * the prescribed coefficient bounds. Then, each selected box is checked 
   * exhaustively. The side lengths of the boxes are specified by the subregion 
   * sizes h_i in the configuration.
   *
   * If the search is restricted to generators with maximal period, the maximal
   * period condition is checked before constructing the corresponding MRG
   * lattice. When possible, candidates sharing the same last coefficient a_k are
   * skipped together if this value cannot lead to maximal period, avoiding
    * unnecessary lattice constructions.
  */
  template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomBox()
  {
    Int subregionTotal(1);
    Int akStepSize(1);
    
    // Calculate subregion size
    for (int i = 1; i <= getOrder(); i++) {
       Int hi = getSubregionSize(i);
        subregionTotal *= hi;
        if (i < getOrder()) akStepSize *= hi;
    }

    while (this->currentGen < this->getMaxGen() &&
           currentRandomRegion < getNumRandomRegions())
    {
        // If the current box is exhausted, move to the next random box.
        if (currentSubIndex >= subregionTotal) {
            ++currentRandomRegion;
            currentSubIndex = 0;
            hasLastAk = false;
            continue;
        }

        // Draw a new random box.
        if (currentSubIndex == 0) {
            for (int i = 1; i <= getOrder(); i++) {
                alpha[i] = randInt(getLowBoundary(i), getHighBoundary(i) - getSubregionSize(i) + 1);
            }
            hasLastAk = false;
        }

        // Decode currentSubIndex into a vector inside the current box.
        Int tmp = currentSubIndex;
        Int offset;
        Int offsetK(0);
        for (int i = 1; i <= getOrder(); i++) {
            Int hi = getSubregionSize(i);
            offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
            // Save what is the current offset in the last variable
            if (i == getOrder()) offsetK = offset;
        }
        
        ++this->currentGen;
        // Check maximal period only when a_k has changed.
        if (this->onlyMaxPeriod()) {
            bool akChanged = !hasLastAk || aa[getOrder()] != lastAk;
            if (akChanged) {
                setModulusIntP<Int>(getModulus());
                lastMaxPeriodResult = getMRGComponent().maxPeriod(aa);
                lastAk = aa[getOrder()];
                hasLastAk = true;
            }
            // If this a_k fails, skip all remaining vectors with this a_k.
            if (!lastMaxPeriodResult) {
                // This ensures that we jump to the next value of a_k while 
                // but do not skip any other points in the box
                if (!countBadAk) --this->currentGen;
                currentSubIndex = (offsetK + 1) * akStepSize;
                continue;
            }
        }
        ++currentSubIndex;
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }
    
    return nullptr;
}

/**
 * nextGeneratorRandomAkFixed() performs a random search organized around the
 * last multiplier a_k. It first draws values of a_k at random and
 * retains only those that satisfy the maximal-period requirement, if this
 * requirement is chosen by the user.
 *
 * For each accepted value of a_k, the method generates a prescribed number of
 * random coefficient vectors by drawing the remaining coefficients
 * a_1, ..., a_{k-1} independently within their configured bounds. 
 */
template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomAkFixed()
{  
    while (this->currentGen < this->getMaxGen() &&
           currentAkTrial < getNumRandomAk())
    {
        /*
         * If we do not currently have a good a_k, or if we already generated
         * N_a vectors for the current a_k, draw a new a_k.
         */
        if (!hasCurrentAk || currentAkSample >= getNumRandomPointsPerAk()) {
            hasCurrentAk = false;
            currentAkSample = 0;

            
            for (int i = 1; i < getOrder(); i++) {
                aa[i] = randInt(getLowBoundary(i), getHighBoundary(i));
            }
        
            while (currentAkTrial < getNumRandomAk() &&
                   this->currentGen < this->getMaxGen())
            {
                currentAk = randInt(getLowBoundary(getOrder()), getHighBoundary(getOrder()));
                aa[getOrder()] = currentAk;
                ++currentAkTrial;

                /*
                 * Check whether this a_k has maximal period.
                 */
                if (this->onlyMaxPeriod()) {
                    if (!getMRGComponent().maxPeriod(aa)) {                        
                        if (countBadAk) ++this->currentGen;
                        continue;
                    }
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return nullptr;
        }
        ++currentAkSample;
        ++this->currentGen;
        
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }
    return nullptr;
}

/**
 * nextGeneratorRandomAkFixedBox() performs a hybrid random-box search organized
 * around the las mulitpiliers a_k. It first draws values of a_k at
 * random and retains only those that satisfy the maximal-period requirement, if
 * wanted by the user.
 *
 * For each accepted value of a_k, the method draws a random box for the
 * remaining coefficients a_1, ..., a_{k-1}. The box is defined by the configured
 * subregion sizes h_i for i = 1, ..., k-1, while a_k remains fixed, corresponding
 * to h_k = 1. The points in this box are then checked exhaustively.
 */
template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomAkFixedBox()
{
    /*
     * h_k is assumed to be 1. Thus the box is only over the coordinates
     * a_1, ..., a_{k-1}.
     */
    Int boxSize(1);
    for (int i = 1; i < getOrder(); i++) {
        Int hi = getSubregionSize(i);
        boxSize *= hi;
    }

    while (this->currentGen < this->getMaxGen() &&
           (hasCurrentAk || currentAkTrial < getNumRandomAk()))
    {
        /*
         * If there is no current good a_k, or if the current box has been
         * exhausted, draw a new a_k and a new random box for the other
         * coordinates.
         */
        if (!hasCurrentAk || currentSubIndex >= boxSize) {
            hasCurrentAk = false;
            currentSubIndex = 0;

            while (currentAkTrial < getNumRandomAk() &&
                   this->currentGen < this->getMaxGen())
            {
                /*
                 * Draw the lower corner of the box for coordinates
                 * a_1, ..., a_{k-1}.
                 */
                for (int i = 1; i < getOrder(); i++) {
                    alpha[i] = randInt(getLowBoundary(i), getHighBoundary(i) - getSubregionSize(i) + 1);
                    aa[i] = alpha[i];
                }

                /*
                 * Draw a_k independently.
                 */
                currentAk = randInt(getLowBoundary(getOrder()), getHighBoundary(getOrder()));
                aa[getOrder()] = currentAk;

                ++currentAkTrial;

                /*
                 * Check whether this a_k is good. The check is performed
                 * only once, using the first point of the box.
                 */
                if (this->onlyMaxPeriod()) {
                    if (!getMRGComponent().maxPeriod(aa)) {                        
                        if (countBadAk) ++this->currentGen;
                        continue;
                    }
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return nullptr;
        }

        /*
         * Enumerate the current box exhaustively for a_1, ..., a_{k-1},
         * while keeping a_k fixed.
         */
        Int tmp = currentSubIndex;

        for (int i = 1; i < getOrder(); i++) {
            Int hi = getSubregionSize(i);
            Int offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
        }
        aa[getOrder()] = currentAk;

        ++currentSubIndex;
        ++this->currentGen;
        
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
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
    return getMRGComponent().maxPeriod(lat.getaa());
  }

} // namespace LatMRG

#endif // LATMRG_SEARCHMRG_H