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
    long numRandomRegions = 1;   // For search of random boxes (cooresponds to N_r in the guide)
    IntVec subregionSize;        // The allowed sizes of the boxes h_i, indexed from 1 to k
    int numRandomPointsPerAk = 1; // Number of random points per a_k which satisfies the maximal period conditon; default = 1
    int numRandomAk = 1000; // Maximal number of random values for a_1, ..., a_{k-1} drawn for fixed a_k
    bool countBadAk = false; // Count a_k's which do not have maximal period in the total number of generators checked

    // Specific constraints on the generators
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal
                       
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
    const vector<vector<pair<long,long>>>& getEqualCoeffs() const {return equalCoeffs;}
    void setEqualCoeffs(const vector<vector<pair<long,long>>>& values) {equalCoeffs = values;}

    Int getNoMultipliers() const
    {
        Int total(1);
        for (long i = 1; i <= this->getOrder(); i++) {
            total *= highBoundaries[i] - lowBoundaries[i] + 1;
        }
        return total;
    }

    // These functions are not called at all but need to be there,
    // because we need them for combined MRGs and thus virtual methods
    // have been defined in the SearchBase class.
    virtual Int getComponentModulus(int k) const { return this->getModulus();}
    virtual long getComponentOrder(int k) const { return this->getOrder();}

    /*
    * Vector for creating the multipliers
    */
    NTL::Vec<NTL::ZZ> aa;
    
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
        WorkingLattice = std::make_unique<MRGLattice<Int, Real>>(this->getModulus(), getLowBoundaries(), this->maxdim);
        aa.SetLength(this->getOrder() + 1);
        alpha.SetLength(this->getOrder() + 1);
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
            mrg.reset(new MRGComponent<Int>(this->getModulus(), this->getOrder()));
        }
        return *mrg;
    }

    /*
    * This lattice is re-used as output of nextGenerator-type functions
    */
    std::unique_ptr<MRGLattice<Int, Real>> WorkingLattice;

    /*
    * Helper variables for box searches
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
    while (this->currentGen < this->getMaxNumGen() &&
           this->currentGen < getNoMultipliers())
    {
        // Set the vector of multipliers for the next generator
        tmp = NTL::to_ZZ(this->currentGen);
        for (int i = 1; i <= this->getOrder(); i++) {
            range = getHighBoundary(i) - getLowBoundary(i) + 1;
            val = tmp % range;
            tmp /= range;
            aa[i] = getLowBoundary(i) + val;
        }
        ++this->currentGen;
        
        if (this->getMaxPeriod() && !getMRGComponent().maxPeriod(aa)) {
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
    
    while (this->currentGen < this->getMaxNumGen())
    {
        for (int i = 1; i <= this->getOrder(); i++) aa[i] = randInt(getLowBoundary(i),  getHighBoundary(i));
        ++this->currentGen;
        if (this->getMaxPeriod() && !getMRGComponent().maxPeriod(aa)) {
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
    for (int i = 1; i <= this->getOrder(); i++) {
       Int hi = getSubregionSize(i);
        subregionTotal *= hi;
        if (i < this->getOrder()) akStepSize *= hi;
    }

    while (this->currentGen < this->getMaxNumGen() &&
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
            for (int i = 1; i <= this->getOrder(); i++) {
                alpha[i] = randInt(getLowBoundary(i), getHighBoundary(i) - getSubregionSize(i) + 1);
            }
            hasLastAk = false;
        }

        // Decode currentSubIndex into a vector inside the current box.
        Int tmp = currentSubIndex;
        Int offset;
        Int offsetK(0);
        for (int i = 1; i <= this->getOrder(); i++) {
            Int hi = getSubregionSize(i);
            offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
            // Save what is the current offset in the last variable
            if (i == this->getOrder()) offsetK = offset;
        }
        
        ++this->currentGen;
        // Check maximal period only when a_k has changed.
        if (this->getMaxPeriod()) {
            bool akChanged = !hasLastAk || aa[this->getOrder()] != lastAk;
            if (akChanged) {
                setModulusIntP<Int>(this->getModulus());
                lastMaxPeriodResult = getMRGComponent().maxPeriod(aa);
                lastAk = aa[this->getOrder()];
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
    while (this->currentGen < this->getMaxNumGen() &&
           currentAkTrial < getNumRandomAk())
    {
        /*
         * If we do not currently have a good a_k, or if we already generated
         * N_a vectors for the current a_k, draw a new a_k.
         */
        if (!hasCurrentAk || currentAkSample >= getNumRandomPointsPerAk()) {
            hasCurrentAk = false;
            currentAkSample = 0;

            
            for (int i = 1; i < this->getOrder(); i++) {
                aa[i] = randInt(getLowBoundary(i), getHighBoundary(i));
            }
        
            while (currentAkTrial < getNumRandomAk() &&
                   this->currentGen < this->getMaxNumGen())
            {
                currentAk = randInt(getLowBoundary(this->getOrder()), getHighBoundary(this->getOrder()));
                aa[this->getOrder()] = currentAk;
                ++currentAkTrial;

                /*
                 * Check whether this a_k has maximal period.
                 */
                if (this->getMaxPeriod()) {
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
    for (int i = 1; i < this->getOrder(); i++) {
        Int hi = getSubregionSize(i);
        boxSize *= hi;
    }

    while (this->currentGen < this->getMaxNumGen() &&
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
                   this->currentGen < this->getMaxNumGen())
            {
                /*
                 * Draw the lower corner of the box for coordinates
                 * a_1, ..., a_{k-1}.
                 */
                for (int i = 1; i < this->getOrder(); i++) {
                    alpha[i] = randInt(getLowBoundary(i), getHighBoundary(i) - getSubregionSize(i) + 1);
                    aa[i] = alpha[i];
                }

                /*
                 * Draw a_k independently.
                 */
                currentAk = randInt(getLowBoundary(this->getOrder()), getHighBoundary(this->getOrder()));
                aa[this->getOrder()] = currentAk;

                ++currentAkTrial;

                /*
                 * Check whether this a_k is good. The check is performed
                 * only once, using the first point of the box.
                 */
                if (this->getMaxPeriod()) {
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

        for (int i = 1; i < this->getOrder(); i++) {
            Int hi = getSubregionSize(i);
            Int offset = tmp % hi;
            tmp /= hi;
            aa[i] = alpha[i] + offset;
        }
        aa[this->getOrder()] = currentAk;

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