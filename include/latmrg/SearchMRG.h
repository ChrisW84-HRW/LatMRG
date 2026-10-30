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
 * different search functions.
 *
 * The main feature of the class is that it  provides generator routines for different 
 * search functions such as exhaustive search and different types of random search.
 */
template<typename Int, typename Real>
struct SearchMRG : public SearchBase<Int, Real, MRGLattice<Int, Real>>
{
    SearchMRG() : SearchBase<Int, Real, MRGLattice<Int, Real>>(MRG) {}

    // Definition of search range    
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j 
    IntVec highBoundaries; // Default values are 1 and m-1 

    // Settings for random search loop algorithms
    long numRandomRegions = 1;   // For search of random boxes (corresponds to N_r in the guide)
    IntVec subregionSize;        // The allowed sizes of the boxes h_i, indexed from 1 to k
    int numRandomPointsPerAk = 1; // Number of random values for a_1, ..., a_{k-1} per a_k which satisfies the maximal period condition; default = 1
    int numRandomAk = 1000; // Maximal number of random values fixed a_k

    // Specific constraints on generators
    vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
    vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
    std::vector<int> equalityGroups; // Entry i-1 is the group of a_i.
                       
    // Get and set functions
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
    // because we need them for combined MRGs and thus virtual functions
    // have been defined in the SearchBase class.
    virtual Int getComponentModulus(int k) const { return this->modulus;}
    virtual long getComponentOrder(int k) const { return this->order;}


    /**
    * Search strategies. Each function generates and tests candidate coefficient
    * vectors until its search region or a configured limit is reached.
    *
    * Exhaustive: enumerates the full coefficient range.
    * Random: draws independent vectors within the coefficient bounds.
    * RandomBox: exhaustively enumerates randomly selected boxes.
    * RandomByAk: draws several random vectors for each accepted a_k.
    * RandomByAkBox: enumerates a random box for each accepted a_k.
    * RandomApproxFactoring: draws vectors in approximate-factoring form.
    * RandomApproxFactoringByAk: keeps each accepted a_k fixed for several such vectors.
    * RandomPowerOfTwo: draws vectors as signed sums of powers of two.
    * RandomPowerOfTwoByAk: keeps each accepted a_k fixed for several such vectors.
    */
    void searchLoopExhaustive(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandom(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomBox(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomByAk(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomByAkBox(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomApproxFactoring(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomApproxFactoringByAk(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomPowerOfTwo(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomPowerOfTwoByAk(IntLattice<Int, Real>* proj, int& old);

    /**
    * Exhaustively tests the coefficient vectors in the specified box.
    *
    * 'boxLow' and 'boxHigh' give the lower and upper bounds for each coefficient
    * a_i, indexed from 1 through order. 'proj' is the projection lattice used
    * for FOM calculation; it may be null when FOM calculation is disabled.
    * 'old' tracks the last progress-bar position and is updated when progress
    * is printed.
    *
    * If 'akAlreadyAccepted' is true, the box must have a fixed a_k
    * (boxLow[order] == boxHigh[order]) that has already passed the
    * maximal-period test. Otherwise, the test is performed once per distinct
    * a_k in the box. The search stops at the box boundary or a search limit.
    */
    void searchLoopBox(const IntVec& boxLow, const IntVec& boxHigh, IntLattice<Int, Real>* proj, int& old, bool akAlreadyAccepted = false);
    
    /*
    * These functions draw random coefficient with additional constaints.
    * drawRandomCoefficients() respects the set equal coefficients.
    * drawPowerOfTwoCoefficient draws a single coefficient which is the sum of powers of two
    */
    void drawRandomCoefficients(int fixedIndex = 0);
    Int drawPowerOfTwoCoefficient(int i);
    
    /**
     * This function checks whether the generated lattice satisfies the
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
    void resetSearchState() override {
        mrg.reset();
        WorkingLattice.reset();
        aa.SetLength(this->order + 1);
        if (lowBoundaries.length() == this->order + 1 && highBoundaries.length() == this->order + 1)  prepareEqualityGroups();
    }

    /**
     * Returns the vector of multipliers
     */
    const IntVec& getCurrentCoefficients() const override {return aa;}

    
    /**
    * Vector of multipliers
    */
    NTL::Vec<NTL::ZZ> aa;

    
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
    * This lattice is re-used for storing the current lattice if needed.
    * The function gets the lattice.
    */
    std::unique_ptr<MRGLattice<Int, Real>> WorkingLattice;
    MRGLattice<Int, Real>& getWorkingLattice() override;

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
 * Returns the working lattice for the current candidate.
 *
 * Constructs the lattice on first use and updates it with the current
 * coefficients on subsequent calls. This function is called only when
 * computing a figure of merit, so period-only searches do not construct
 * a lattice.
 */
template<typename Int, typename Real> MRGLattice<Int, Real>& SearchMRG<Int, Real>::getWorkingLattice()
{
    if (!WorkingLattice) {
        WorkingLattice = std::make_unique<MRGLattice<Int, Real>>(
            this->modulus, aa, this->maxdim
        );
    } else {
        WorkingLattice->setaa(aa);
    }
    return *WorkingLattice;
}


//===========================================================================

/**
* Exhaustively tests the coefficient vectors in the specified box.
* It assumes valid coefficient bounds.
*/
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopBox(const IntVec& boxLow, const IntVec& boxHigh, IntLattice<Int, Real>* proj, int& old, bool akAlreadyAccepted)
{
    Int boxTotal(1);
    Int akStepSize(1);
    Int tmp, sideLength, offset;
    Int offsetK(0);
    for (int i = 1; i <= this->order; ++i) {
        const Int sideLength = boxHigh[i] - boxLow[i] + 1;
        boxTotal *= sideLength;
        if (i < this->order)
            akStepSize *= sideLength;
    }
    Int index(0);
    bool hasTestedAk = akAlreadyAccepted && this->maxPeriod;
    Int lastTestedAk;
    if (hasTestedAk) lastTestedAk = boxLow[this->order];
    bool lastPeriodResult = hasTestedAk;
    while (index < boxTotal && !this->reachedNumGeneratorLimit() &&  !this->timer.timeOver(this->timeLimit))
    {
        tmp = index;
        // Decode the index as a mixed-radix number. The first coefficient
        // changes fastest, so all vectors with one a_k form a contiguous block.
        for (int i = 1; i <= this->order; ++i) {
            sideLength = boxHigh[i] - boxLow[i] + 1;
            offset = tmp % sideLength;
            tmp /= sideLength;
            aa[i] = boxLow[i] + offset;
            if (i == this->order)  offsetK = offset;
        }
        ++index;
        ++this->numCandidatesExamined;
        if (this->maxPeriod) {
            if (!hasTestedAk || aa[this->order] != lastTestedAk) {
                lastPeriodResult = getMRGComponent().maxPeriod(aa);
                lastTestedAk = aa[this->order];
                hasTestedAk = true;
            }
            if (!lastPeriodResult) {
                // Move to the first vector with the next a_k.
                // Skipped vectors are not counted as examined.
                index = (offsetK + 1) * akStepSize;
                if (this->progress) old = this->printProgress(old);
                continue;
            }
            ++this->numMaxPeriodCandidates;
        }
        this->testCandidate(aa, proj);
        if (this->progress) old = this->printProgress(old);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}


//===========================================================================

/**
* sarchLoopRandom() generates admissible MRG generators by randomly selecting 
* the recurrence coefficients within the ranges specified in the configuration.
* If the search is restricted to generators with maximal period, then 
* searchLoopRandom() only yields generators having this property.
*/
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandom(IntLattice<Int, Real>* proj, int& old)
{
    // Keep drawing until the time limit or the maximal number of generators 
    // is reached.
    while (!this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
    {
        drawRandomCoefficients();
        ++this->numCandidatesExamined;
        if (this->progress) old = this->printProgress(old);
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(aa)) continue;
            ++this->numMaxPeriodCandidates;
        }
        this->testCandidate(aa, proj);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}

//===========================================================================


/**
* searchLoopExhaustive() performs an exhaustive enumeration of the admissible MRG 
* generators defined by the configuration. Passes the full coefficient bounds to 
* searchLoopBox(), which enumerates the coefficient vectors. If maximal period is 
* required, vectors with an a_k that fails the period test are skipped. Accepted vectors 
* are passed to testCandidate(); a lattice is constructed only when the FOM
* is calculated.
*/
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopExhaustive(IntLattice<Int, Real>* proj, int& old)
{
    searchLoopBox(lowBoundaries, highBoundaries, proj, old);
}

//===========================================================================

/**
 * saarchLoopRandomByAk() performs a random search organized around the
 * last multiplier a_k. It first draws values of a_k at random and
 * retains only those that satisfy the maximal-period requirement, if this
 * requirement is chosen by the user.
 *
 * For each accepted value of a_k, the function generates a prescribed number of
 * random coefficient vectors by drawing the remaining coefficients
 * a_1, ..., a_{k-1} independently within their configured bounds.
 */
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandomByAk(IntLattice<Int, Real>* proj, int& old)
{
    long akTrials = 0;
    long pointsForCurrentAk = 0;
    bool hasAcceptedAk = false;
    Int acceptedAk;
    while (!this->reachedNumGeneratorLimit() && (hasAcceptedAk || akTrials < numRandomAk) && !this->timer.timeOver(this->timeLimit))
    {
        if (!hasAcceptedAk ||
            pointsForCurrentAk >= numRandomPointsPerAk) {
            hasAcceptedAk = false;
            pointsForCurrentAk = 0;            
            // Keep drawing until the time limit or the maximal number of generators 
            // or the maximal number of trials for a_k is reached.
            while (akTrials < numRandomAk && !this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
            {
                acceptedAk = randInt(lowBoundaries[this->order], highBoundaries[this->order]);
                aa[this->order] = acceptedAk;
                ++akTrials;
                if (this->maxPeriod && !getMRGComponent().maxPeriod(aa)) {
                    ++this->numCandidatesExamined;
                    if (this->progress) old = this->printProgress(old);
                    continue;
                }
                hasAcceptedAk = true;
                break;
            }
            if (!hasAcceptedAk) break;
        }
        // Draw a fresh vector for each point, keeping the accepted a_k.
        aa[this->order] = acceptedAk;
        drawRandomCoefficients(this->order);
        ++pointsForCurrentAk;
        ++this->numCandidatesExamined;
        if (this->maxPeriod) ++this->numMaxPeriodCandidates;
        this->testCandidate(aa, proj);
        if (this->progress) old = this->printProgress(old);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}

//===========================================================================


/**
 * searchLoopRandomByAkBox() performs a hybrid random-box search organized
 * around the las mulitpiliers a_k. It first draws values of a_k at
 * random and retains only those that satisfy the maximal-period requirement, if
 * wanted by the user.
 *
 * For each accepted value of a_k, the function draws a random box for the
 * remaining coefficients a_1, ..., a_{k-1}. The box is defined by the configured
 * subregion sizes h_i for i = 1, ..., k-1, while a_k remains fixed, corresponding
 * to h_k = 1. The points in this box are then checked exhaustively.
 */
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandomByAkBox(IntLattice<Int, Real>* proj, int& old)
{
    IntVec boxLow, boxHigh;
    Int ak;
    bool acceptedAk;
    boxLow.SetLength(this->order + 1);
    boxHigh.SetLength(this->order + 1);
    long akTrials = 0;
    // Draw up to numRandomAk values of a_k, stopping earlier at a search limit.
    // Each accepted a_k has its entire box enumerated in this iteration.
    while (!this->reachedNumGeneratorLimit() && akTrials < numRandomAk && !this->timer.timeOver(this->timeLimit))
    {
        acceptedAk = false;
        // Draw boxes and values of a_k until one passes the period test,
        // or until a search limit is reached.
        while (akTrials < numRandomAk &&  !this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
        {
            for (int i = 1; i < this->order; ++i) {
                boxLow[i] = randInt(lowBoundaries[i], highBoundaries[i] - subregionSize[i] + 1);
                boxHigh[i] = boxLow[i] + subregionSize[i] - 1;
                // Preserve the original behavior: the period test receives
                // the first coefficient vector in the drawn box.
                aa[i] = boxLow[i];
            }
            ak = randInt(lowBoundaries[this->order], highBoundaries[this->order]);
            aa[this->order] = ak;
            ++akTrials;
            if (this->maxPeriod && !getMRGComponent().maxPeriod(aa)) {
                ++this->numCandidatesExamined;
                if (this->progress) old = this->printProgress(old);
                continue;
            }
            boxLow[this->order] = ak;
            boxHigh[this->order] = ak;
            acceptedAk = true;
            break;
        }
        if (!acceptedAk) break;
        // The entire accepted box is processed here, including when this
        // was the final permitted a_k draw.
        searchLoopBox(boxLow, boxHigh, proj, old, true);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}
//===========================================================================

/**
 * searchLoopRandomBox() searches randomly selected boxes within the 
 * coefficient bounds. For each box, the function draws a lower corner for every coefficient,
 * including a_k. Its upper corner is determined by the configured
 * subregion sizes. searchLoopBox() then enumerates the box exhaustively,
 * testing the maximal-period condition once per distinct a_k in that box.
 * 
 * The search stops when all requested boxes have been visited or a
 * generator, time, or period-only list limit is reached.
 */
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandomBox(IntLattice<Int, Real>* proj, int& old)
{
    IntVec boxLow, boxHigh;
    boxLow.SetLength(this->order + 1);
    boxHigh.SetLength(this->order + 1);
    for (long region = 0;  region < numRandomRegions && !this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit);
         ++region)
    {
        for (int i = 1; i <= this->order; ++i) {
            boxLow[i] = randInt(lowBoundaries[i], highBoundaries[i] - subregionSize[i] + 1);
            boxHigh[i] = boxLow[i] + subregionSize[i] - 1;
        }
        searchLoopBox(boxLow, boxHigh, proj, old);
        if (!this->calculateFOM && this->bestCandidates.full())
            break;
    }
}

//===========================================================================

/**
 * Searches for MRGs whose coefficients satisfy a sufficient condition
 * for approximate factoring.
 * 
 * For each coefficient a_i, the function draws an integer j uniformly from
 * 0 < j < sqrt(m). It then randomly chooses |a_i| to be either j or
 * floor(m / j), and independently chooses the sign. The usual coefficient
 * boundaries and equality constraints are not used.
 */
template<typename Int, typename Real> 
void SearchMRG<Int, Real>::searchLoopRandomApproxFactoring(IntLattice<Int, Real>* proj, int& old)
{
    Int maxJ, j, draw;
    long bits;
    NTL::SqrRoot(maxJ, this->modulus);
    if (maxJ * maxJ == this->modulus) --maxJ;
    while (!this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
    {
        for (int i = 1; i <= this->order; ++i) {
            // A single uniform draw encodes both j and two choice bits.
            // Dividing by 4 gives j in [1, maxJ]; the remainder selects
            // either j or floor(modulus / j) and determines its sign.
            draw = randInt(Int(0), 4 * maxJ - 1);
            j = draw / 4 + 1;
            bits = NTL::conv<long>(draw % 4);
            aa[i] = (bits & 1) ? this->modulus / j : j;
            if (bits & 2) aa[i] = -aa[i];
        }
        ++this->numCandidatesExamined;
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(aa)) {
                if (this->progress) old = this->printProgress(old);
                continue;
            }
            ++this->numMaxPeriodCandidates;
        }
        this->testCandidate(aa, proj);
        if (this->progress) old = this->printProgress(old);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}

//===========================================================================

/**
 * Like searchLoopRandomApproxFactoring() but keeps an
 * accepted a_k fixed for several random points.
 */
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandomApproxFactoringByAk(IntLattice<Int, Real>* proj, int& old)
{
    Int maxJ, draw, j;
    long bits;
    NTL::SqrRoot(maxJ, this->modulus);
    if (maxJ * maxJ == this->modulus) --maxJ;
    long akTrials = 0;
    long pointsForCurrentAk = 0;
    bool hasAcceptedAk = false;
    Int acceptedAk;
    while (!this->reachedNumGeneratorLimit() && (hasAcceptedAk || akTrials < numRandomAk) && !this->timer.timeOver(this->timeLimit))
    {
        if (!hasAcceptedAk ||
            pointsForCurrentAk >= numRandomPointsPerAk)
        {
            hasAcceptedAk = false;
            pointsForCurrentAk = 0;
            // Draw a_k until one passes the period test or a limit is reached.
            while (akTrials < numRandomAk && !this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
            {
                // A single uniform draw encodes both j and two choice bits.
                // Dividing by 4 gives j in [1, maxJ]; the remainder selects
                // either j or floor(modulus / j) and determines its sign.
                draw = randInt(Int(0), 4 * maxJ - 1);
                j = draw / 4 + 1;
                bits = NTL::conv<long>(draw % 4);
                acceptedAk = (bits & 1) ? this->modulus / j : j;
                if (bits & 2) acceptedAk = -acceptedAk;
                aa[this->order] = acceptedAk;
                ++akTrials;
                if (this->maxPeriod && !getMRGComponent().maxPeriod(aa))
                {
                    ++this->numCandidatesExamined;
                    if (this->progress) old = this->printProgress(old);
                    continue;
                }
                hasAcceptedAk = true;
                break;
            }
            if (!hasAcceptedAk) break;
        }
        // Draw the remaining coefficients while keeping a_k fixed.
        for (int i = 1; i < this->order; ++i) {
            draw = randInt(Int(0), 4 * maxJ - 1);
            j = draw / 4 + 1;
            bits = NTL::conv<long>(draw % 4);
            aa[i] = (bits & 1) ? this->modulus / j : j;
            if (bits & 2) aa[i] = -aa[i];
        }
        aa[this->order] = acceptedAk;
        ++pointsForCurrentAk;
        ++this->numCandidatesExamined;
        if (this->maxPeriod) ++this->numMaxPeriodCandidates;
        this->testCandidate(aa, proj);
        if (this->progress) old = this->printProgress(old);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}

//===========================================================================

/**
 * Searches by drawing independent coefficient vectors as signed sums
 * of distinct powers of two.
 * 
 * For each a_i, chooses between 1 and n_i distinct exponents from
 * 0, ..., p_i, then independently chooses the sign of each term.
 * Equality constraints and coefficient boundaries are not used.
 */
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandomPowerOfTwo(IntLattice<Int, Real>* proj, int& old)
{
    while (!this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
    {
        for (int i = 1; i <= this->order; ++i)
            aa[i] = drawPowerOfTwoCoefficient(i);
        ++this->numCandidatesExamined;
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(aa)) {
                if (this->progress) old = this->printProgress(old);
                continue;
            }
            ++this->numMaxPeriodCandidates;
        }
        this->testCandidate(aa, proj);
        if (this->progress) old = this->printProgress(old);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}


//===========================================================================

/**
 * Like searchLoopRandomPowerofTwo but for each accepted a_k 
 * it draws several independent vectors for the remaining
 * coefficients while keeping a_k fixed.
 */
template<typename Int, typename Real>
void SearchMRG<Int, Real>::searchLoopRandomPowerOfTwoByAk(IntLattice<Int, Real>* proj, int& old)
{
    long akTrials = 0;
    long pointsForCurrentAk = 0;
    bool hasAcceptedAk = false;
    Int acceptedAk;
    while (!this->reachedNumGeneratorLimit() && (hasAcceptedAk || akTrials < numRandomAk) && !this->timer.timeOver(this->timeLimit))
    {
        if (!hasAcceptedAk ||
            pointsForCurrentAk >= numRandomPointsPerAk)
        {
            hasAcceptedAk = false;
            pointsForCurrentAk = 0;
            // Draw a_k until one passes the period test or a limit is reached.
            while (akTrials < numRandomAk && !this->reachedNumGeneratorLimit() && !this->timer.timeOver(this->timeLimit))
            {
                acceptedAk = drawPowerOfTwoCoefficient(this->order);
                aa[this->order] = acceptedAk;
                ++akTrials;
                if (this->maxPeriod && !getMRGComponent().maxPeriod(aa))
                {
                    ++this->numCandidatesExamined;
                    if (this->progress) old = this->printProgress(old);
                    continue;
                }
                hasAcceptedAk = true;
                break;
            }
            if (!hasAcceptedAk) break;
        }
        // Draw fresh coefficients while keeping the accepted a_k fixed.
        for (int i = 1; i < this->order; ++i)
            aa[i] = drawPowerOfTwoCoefficient(i);
        aa[this->order] = acceptedAk;
        ++pointsForCurrentAk;
        ++this->numCandidatesExamined;
        if (this->maxPeriod) ++this->numMaxPeriodCandidates;
        this->testCandidate(aa, proj);
        if (this->progress) old = this->printProgress(old);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}

/**
* This function checks whether the generated lattice satisfies the
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
 * prepareEqualityGroups() must be called before this function.
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