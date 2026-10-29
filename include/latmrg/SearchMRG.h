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

    // Definition of search range    
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j 
    IntVec highBoundaries; // Default values are 1 and m-1 

    // Settings for random next generator algorithms
    long numRandomRegions = 1;   // For search of random boxes (corresponds to N_r in the guide)
    IntVec subregionSize;        // The allowed sizes of the boxes h_i, indexed from 1 to k
    int numRandomPointsPerAk = 1; // Number of random values for a_1, ..., a_{k-1} per a_k which satisfies the maximal period condition; default = 1
    int numRandomAk = 1000; // Maximal number of random values fixed a_k
    bool countBadAk = false; // Count a_k's which do not have maximal period in the total number of generators checked

    // Specific constraints on the generators
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    std::vector<int> equalityGroups; // Entry i-1 is the group of a_i.
                       
    // Get and set methods
    const IntVec& getLowBoundaries() const {return lowBoundaries;}
    void setLowBoundaries(const IntVec& values) {lowBoundaries = values;}
    Int getLowBoundary(int i) const {return lowBoundaries[i];}
    void setLowBoundary(int i, const Int& value) {lowBoundaries[i] = value;}
    const IntVec& getHighBoundaries() const {return highBoundaries;}
    void setHighBoundaries(const IntVec& values) {highBoundaries = values;}
    Int getHighBoundary(int i) const {return highBoundaries[i];}
    void setHighBoundary(int i, const Int& value) {highBoundaries[i] = value;}

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

    const vector<long>& getNumPowerTwo() const {return numPowerTwo;}
    void setNumPowerTwo(const vector<long>& values) {numPowerTwo = values;}
    const vector<long>& getMaxPowerTwo() const {return maxPowerTwo;}
    void setMaxPowerTwo(const vector<long>& values) {maxPowerTwo = values;}
    const std::vector<int>& getEqualityGroups() const {return equalityGroups;}
    void setEqualityGroups(const std::vector<int>& groups) {equalityGroups = groups;}

    Int getNumMultipliers() const
    {
        Int total(1);
        for (long i = 1; i <= this->order; i++) {
            total *= highBoundaries[i] - lowBoundaries[i] + 1;
        }
        return total;
    }

    // These functions are not called at all but need to be there,
    // because we need them for combined MRGs and thus virtual methods
    // have been defined in the SearchBase class.
    virtual Int getComponentModulus(int k) const { return this->modulus;}
    virtual long getComponentOrder(int k) const { return this->order;}

    /**
    * Vector of multipliers
    */
    NTL::Vec<NTL::ZZ> aa;

    /**
     * These methods generate the next MRG candidate according to the selected search strategy
     * Each method returns the next admissible lattice, or nullptr when its
     * search limits have been reached.
     */
    MRGLattice<Int, Real>* nextGeneratorExhaustive() override;
    MRGLattice<Int, Real>* nextGeneratorRandom() override;
    MRGLattice<Int, Real>* nextGeneratorRandomBox();    
    MRGLattice<Int, Real>* nextGeneratorRandomAkFixed();
    MRGLattice<Int, Real>* nextGeneratorRandomAkFixedBox();
    MRGLattice<Int, Real>* nextGeneratorRandomApproxFactoring();
    MRGLattice<Int, Real>* nextGeneratorRandomApproxFactoringAkFixed();
    MRGLattice<Int, Real>* nextGeneratorRandomPowersOfTwo();
    MRGLattice<Int, Real>* nextGeneratorRandomPowersOfTwoAkFixed();

    
    /*
    * These methods draw random coefficient with additional constaints.
    * drawRandomCoefficients() respects the set equal coefficients.
    * drawPowerOfTwoCoefficient draws a single coefficient which is the sum of powers of two
    */
    void drawRandomCoefficients(int fixedIndex = 0);
    Int drawPowerOfTwoCoefficient(int i);
    
    /**
     * This method checks whether the generated lattice satisfies the
     * maximal-period requirement for an MRG.
     */
    bool checkMaxPeriod(MRGLattice<Int, Real>& lat) override;

    /*
    * Helper function which calculates the groups of equal coefficients 
    * according to the vector equalityGroups
    */
    void prepareEqualityGroups();

    /**
     * Resets generator-specific state before a new search, so repeated calls
     * to performSearch() start with a fresh enumeration or random search.
     */
    void resetGeneratorState() override{
        currentRandomRegion = 0;
        currentSubIndex = 0;
        hasLastAk = false;
        lastMaxPeriodResult = false;
        mrg.reset();
        currentAkTrial = 0;
        currentAkSample = 0;
        hasCurrentAk = false;
        WorkingLattice = std::make_unique<MRGLattice<Int, Real>>(this->modulus, lowBoundaries, this->maxdim);
        aa.SetLength(this->order + 1);
        alpha.SetLength(this->order + 1);
        if (lowBoundaries.length() == this->order + 1 && highBoundaries.length() == this->order + 1)  prepareEqualityGroups();
        exhaustiveIndex = Int(0);
    }

    /**
     * Returns the vector of multipliers
     */
    const IntVec& getCurrentCoefficients() const override {return aa;}
    
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
            mrg.reset(new MRGComponent<Int>(this->modulus, this->order));
        }
        return *mrg;
    }

    /*
    * This lattice is re-used as output of nextGenerator type functions
    */
    std::unique_ptr<MRGLattice<Int, Real>> WorkingLattice;

    /*
    * Helper variables for box searches
    */
    long currentRandomRegion = 0;
    Int currentSubIndex = Int(0);
    IntVec alpha;
    bool hasLastAk = false;
    Int lastAk;
    bool lastMaxPeriodResult = false;

    long currentAkTrial = 0;
    long currentAkSample = 0;
    bool hasCurrentAk = false;
    Int currentAk;    

    Int exhaustiveIndex = Int(0);
    
    /*
    * Helper variables for constraints on coefficients
    */
    int numEqualityGroups = 0;
    std::vector<Int> groupLow;
    std::vector<Int> groupHigh;
    std::vector<Int> groupValue;
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
    while (this->currentGen < this->maxNumGen &&
           exhaustiveIndex < getNumMultipliers())
    {
        // Set the vector of multipliers for the next generator
        tmp = exhaustiveIndex;
        ++exhaustiveIndex; // Always advance, even if the period test fails.
        
        for (int i = 1; i <= this->order; i++) {
            range = highBoundaries[i] - lowBoundaries[i] + 1;
            val = tmp % range;
            tmp /= range;
            aa[i] = lowBoundaries[i] + val;
        }
        ++this->currentGen;
        
        if (this->maxPeriod && !getMRGComponent().maxPeriod(aa)) {
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
    while (this->currentGen < this->maxNumGen)
    {
        drawRandomCoefficients();
        ++this->currentGen;
        if (this->maxPeriod && !getMRGComponent().maxPeriod(aa)) {
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
    for (int i = 1; i <= this->order; i++) {
       Int hi = subregionSize[i];
        subregionTotal *= hi;
        if (i < this->order) akStepSize *= hi;
    }

    while (this->currentGen < this->maxNumGen &&
           currentRandomRegion < numRandomRegions)
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
            for (int i = 1; i <= this->order; i++) {
                alpha[i] = randInt(lowBoundaries[i], highBoundaries[i] - subregionSize[i] + 1);
            }
            hasLastAk = false;
        }

        // Decode currentSubIndex into a vector inside the current box.
        Int tmp = currentSubIndex;
        Int offset;
        Int offsetK(0);
        for (int i = 1; i <= this->order; i++) {
            Int hi = subregionSize[i];
            offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
            // Save what is the current offset in the last variable
            if (i == this->order) offsetK = offset;
        }
        
        ++this->currentGen;
        // Check maximal period only when a_k has changed.
        if (this->maxPeriod) {
            bool akChanged = !hasLastAk || aa[this->order] != lastAk;
            if (akChanged) {
                setModulusIntP<Int>(this->modulus);
                lastMaxPeriodResult = getMRGComponent().maxPeriod(aa);
                lastAk = aa[this->order];
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
//===========================================================================

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
    while (this->currentGen < this->maxNumGen &&
           currentAkTrial < numRandomAk)
    {
        /*
         * If we do not currently have a good a_k, or if we already generated
         * N_a vectors for the current a_k, draw a new a_k.
         */
        if (!hasCurrentAk || currentAkSample >= numRandomPointsPerAk) {
            hasCurrentAk = false;
            currentAkSample = 0;
            
            while (currentAkTrial < numRandomAk &&
                   this->currentGen < this->maxNumGen)
            {
                currentAk = randInt(lowBoundaries[this->order], highBoundaries[this->order]);
                aa[this->order] = currentAk;
                ++currentAkTrial;

                /*
                 * Check whether this a_k has maximal period.
                 */
                if (this->maxPeriod) {
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
        // Keeps the class of coefficients which equal a_k fixed and draws only the other coefficients
        drawRandomCoefficients(this->order);
        ++currentAkSample;
        ++this->currentGen;
        
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }
    return nullptr;
}
//===========================================================================

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
    for (int i = 1; i < this->order; i++) {
        Int hi = subregionSize[i];
        boxSize *= hi;
    }

    while (this->currentGen < this->maxNumGen &&
           (hasCurrentAk || currentAkTrial < numRandomAk))
    {
        /*
         * If there is no current good a_k, or if the current box has been
         * exhausted, draw a new a_k and a new random box for the other
         * coordinates.
         */
        if (!hasCurrentAk || currentSubIndex >= boxSize) {
            hasCurrentAk = false;
            currentSubIndex = 0;

            while (currentAkTrial < numRandomAk &&
                   this->currentGen < this->maxNumGen)
            {
                /*
                 * Draw the lower corner of the box for coordinates
                 * a_1, ..., a_{k-1}.
                 */
                for (int i = 1; i < this->order; i++) {
                    alpha[i] = randInt(lowBoundaries[i], highBoundaries[i] - subregionSize[i] + 1);
                    aa[i] = alpha[i];
                }

                /*
                 * Draw a_k independently.
                 */
                currentAk = randInt(lowBoundaries[this->order], highBoundaries[this->order]);
                aa[this->order] = currentAk;

                ++currentAkTrial;

                /*
                 * Check whether this a_k is good. The check is performed
                 * only once, using the first point of the box.
                 */
                if (this->maxPeriod) {
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

        for (int i = 1; i < this->order; i++) {
            Int hi = subregionSize[i];
            Int offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
        }
        aa[this->order] = currentAk;

        ++currentSubIndex;
        ++this->currentGen;
        
        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }

    return nullptr;
}
//===========================================================================

/**
 * Generates MRG candidates whose coefficients satisfy a sufficient condition
 * for approximate factoring.
 *
 * For each coefficient a_i, the method draws an integer j uniformly from
 * 0 < j < sqrt(m). It then randomly chooses |a_i| to be either j or
 * floor(m / j), and independently chooses the sign. The usual coefficient
 * boundaries and equality constraints are not used.
 */
template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomApproxFactoring()
{
    Int j, maxJ;
    NTL::SqrRoot(maxJ, this->modulus);
    if (maxJ * maxJ == this->modulus) --maxJ;

    while (this->currentGen < this->maxNumGen) {
        for (int i = 1; i <= this->order; ++i) {
            j = randInt(NTL::ZZ(1), maxJ);
            aa[i] = randInt(0,1) == 0  ? j : this->modulus / j;
            if (randInt(0,1) == 1) aa[i] = -aa[i];
        }
        ++this->currentGen;
        if (this->maxPeriod && !getMRGComponent().maxPeriod(aa)) {
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
 * Like nextGeneratorRandomApproxFactoring but keeps an
 * accepted a_k fixed for several random points.
 */
template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomApproxFactoringAkFixed()
{
    Int j, maxJ;

    NTL::SqrRoot(maxJ, this->modulus);
    if (maxJ * maxJ == this->modulus) --maxJ;

    while (this->currentGen < this->maxNumGen &&
           (hasCurrentAk || currentAkTrial < numRandomAk))
    {
        if (!hasCurrentAk ||
            currentAkSample >= numRandomPointsPerAk) {
            hasCurrentAk = false;
            currentAkSample = 0;

            while (currentAkTrial < numRandomAk &&
                   this->currentGen < this->maxNumGen)
            {
                j = randInt(Int(1), maxJ);
                currentAk = randInt(0, 1) == 0 ? j : this->modulus / j;
                if (randInt(0, 1) == 1) currentAk = -currentAk;
                aa[this->order] = currentAk;
                ++currentAkTrial;

                if (this->maxPeriod &&
                    !getMRGComponent().maxPeriod(aa)) {
                    if (countBadAk)
                        ++this->currentGen;
                    continue;
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return nullptr;
        }
        // Draw fresh coefficients while keeping the max period a_k fixed.
        for (int i = 1; i < this->order; ++i) {
            j = randInt(Int(1), maxJ);
            aa[i] = randInt(0, 1) == 0 ? j: this->modulus / j;
            if (randInt(0, 1) == 1) aa[i] = -aa[i];
        }
        aa[this->order] = currentAk;
        ++currentAkSample;
        ++this->currentGen;

        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }
    return nullptr;
}
//===========================================================================

/**
 * Generates random MRG coefficients as signed sums of powers of two.
 *
 * For each a_i, chooses between 1 and n_i distinct exponents from
 * 0, ..., p_i, then independently chooses the sign of each term.
 * Equality constraints and coefficient boundaries are not used.
 */
template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomPowersOfTwo()
{
    Int power, value;
    long ni, pi, choice, maxTerms, terms;

    while (this->currentGen < this->maxNumGen) {
        for (int i = 1; i <= this->order; ++i) {
            // These std::vectors are indexed from 0: entry i-1 describes a_i.
            ni = numPowerTwo[i - 1];
            pi = maxPowerTwo[i - 1];
            maxTerms = std::min(ni, pi + 1);
            terms = randInt(1,maxTerms);

            // Draw distinct exponents without replacement.
            std::vector<long> exponents(pi + 1);
            for (long e = 0; e <= pi; ++e)
                exponents[e] = e;

            value = 0;
            for (long j = 0; j < terms; ++j) {
                choice = randInt(j, pi);
                std::swap(exponents[j], exponents[choice]);
                power = NTL::power(NTL::ZZ(2), exponents[j]);
                if (randInt(0,1) == 0) value += power;
                else value -= power;
            }
            aa[i] = value;
        }

        ++this->currentGen;

        if (this->maxPeriod &&
            !getMRGComponent().maxPeriod(aa)) {
            if (!countBadAk)
                --this->currentGen;
            continue;
        }

        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }

    return nullptr;
}
//===========================================================================

/**
 * Like nextGeneratorRandomPowerofTwo but for each accepted a_k 
 * it draws several independent vectors for the remaining
 * coefficients while keeping a_k fixed.
 */
template<typename Int, typename Real> MRGLattice<Int, Real>* SearchMRG<Int, Real>::nextGeneratorRandomPowersOfTwoAkFixed()
{
    while (this->currentGen < this->maxNumGen &&
           (hasCurrentAk || currentAkTrial < numRandomAk))
    {
        // Find a new acceptable a_k when needed.
        if (!hasCurrentAk ||
            currentAkSample >= numRandomPointsPerAk) {
            hasCurrentAk = false;
            currentAkSample = 0;

            while (currentAkTrial < numRandomAk &&
                   this->currentGen < this->maxNumGen)
            {
                currentAk = drawPowerOfTwoCoefficient(this->order);
                aa[this->order] = currentAk;
                ++currentAkTrial;

                if (this->maxPeriod &&
                    !getMRGComponent().maxPeriod(aa)) {
                    if (countBadAk)
                        ++this->currentGen;
                    continue;
                }
                hasCurrentAk = true;
                break;
            }
            if (!hasCurrentAk)
                return nullptr;
        }
        // Draw a fresh vector for the accepted a_k.
        for (int i = 1; i < this->order; ++i)
            aa[i] = drawPowerOfTwoCoefficient(i);

        aa[this->order] = currentAk;

        ++currentAkSample;
        ++this->currentGen;

        WorkingLattice->setaa(aa);
        return WorkingLattice.get();
    }

    return nullptr;
}
//===========================================================================


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
//===========================================================================

/**
 * Prepares the equality constraints for a search.
 *
 * Groups coefficients with the same equality label and determines the common
 * admissible range for each group by intersecting bounds. An empty 
 * equalityGroups vector treats all coefficients as independent.
 */
template<typename Int, typename Real> void SearchMRG<Int, Real>::prepareEqualityGroups()
{
    const int k = this->order;

    if (!equalityGroups.empty() && equalityGroups.size() != static_cast<std::size_t>(k))
        throw std::invalid_argument("One equality-group label is required per coefficient");

    numEqualityGroups = 0;
    groupLow.resize(k + 1);
    groupHigh.resize(k + 1);
    groupValue.resize(k + 1);

    std::vector<bool> seen(k + 1, false);

    for (int i = 1; i <= k; ++i) {
        const int group = equalityGroups.empty() ? i : equalityGroups[i - 1];

        if (group < 1 || group > k)
            throw std::invalid_argument("Invalid equality-group label");

        if (group > numEqualityGroups)
            numEqualityGroups = group;

        if (!seen[group]) {
            groupLow[group] = lowBoundaries[i];
            groupHigh[group] = highBoundaries[i];
            seen[group] = true;
        } else {
            if (lowBoundaries[i] > groupLow[group])
                groupLow[group] = lowBoundaries[i];
            if (highBoundaries[i] < groupHigh[group])
                groupHigh[group] = highBoundaries[i];
        }
    }

    for (int group = 1; group <= numEqualityGroups; ++group) {
        if (!seen[group] || groupLow[group] > groupHigh[group])
            throw std::invalid_argument(
                "Invalid or incompatible equality group");
    }
}
//===========================================================================

/**
 * Draws one random value for each prepared equality group and copies it into
 * every corresponding position of aa. This produces a coefficient vector
 * satisfying the equality constraints.
 *
 * prepareEqualityGroups() must be called before this method.
 */
template<typename Int, typename Real> void SearchMRG<Int, Real>::drawRandomCoefficients(int fixedIndex)
{
    int fixedGroup = 0; // 0 means no coefficient is fixed.

    if (fixedIndex != 0) {fixedGroup = equalityGroups.empty() ? fixedIndex : equalityGroups[fixedIndex - 1];}

    for (int group = 1; group <= numEqualityGroups; ++group) {
        if (group == fixedGroup) groupValue[group] = aa[fixedIndex];
        else groupValue[group] = randInt(groupLow[group], groupHigh[group]);
    }

    for (int i = 1; i <= this->order; ++i) {
        const int group = equalityGroups.empty() ? i : equalityGroups[i - 1];
        aa[i] = groupValue[group];
    }
}
//===========================================================================

/**
 * Draws coefficient a_i as a signed sum of distinct powers of two.
 *
 * Randomly selects between 1 and numPowerTwo[i-1] terms, with distinct
 * exponents from 0 through maxPowerTwo[i-1]. The sign of each term is chosen
 * independently.
 */
template<typename Int, typename Real> Int SearchMRG<Int, Real>::drawPowerOfTwoCoefficient(int i)
{
    const long ni = numPowerTwo[i - 1];
    const long pi = maxPowerTwo[i - 1];
    const long terms = randInt(1, std::min(ni, pi + 1));

    std::vector<long> exponents(pi + 1);
    for (long e = 0; e <= pi; ++e)
        exponents[e] = e;

    Int value(0);

    for (long j = 0; j < terms; ++j) {
        const long choice = randInt(j, pi);
        std::swap(exponents[j], exponents[choice]);

        const Int power = NTL::power(Int(2), exponents[j]);
        if (randInt(0, 1) == 0)
            value += power;
        else
            value -= power;
    }

    return value;
}
//===========================================================================



} // namespace LatMRG

#endif // LATMRG_SEARCHMRG_H