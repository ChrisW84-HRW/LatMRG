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

#ifndef LATMRG_SEARCHMWC_H
#define LATMRG_SEARCHMWC_H

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

#include "latmrg/SearchBase.h"
#include "latmrg/LCGLattice.h"
#include "latmrg/MWCComponent.h"

namespace LatMRG {

template<typename Int, typename Real>
struct SearchMWC : public SearchBase<Int, Real, LCGLattice<Int, Real>>
{
    SearchMWC() : SearchBase<Int, Real, LCGLattice<Int, Real>>(MWC) {this->order = 0;}

    // MWC base b = 2^e. The associated LCG modulus is candidate-dependent
    // and is stored separately in this->modulus.
    int64_t exponent = 0;
    Int base = Int(0);

     // Definition of search range    
    IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j
    IntVec highBoundaries;

    // Optional coefficient constraints.
    std::vector<int> equalityGroups;
    std::vector<long> numPowerTwo;
    std::vector<long> maxPowerTwo;

    PeriodRequirement periodRequirement = PeriodRequirement::SafePrimeHalfPeriod;
    int64_t primalityTrials = 100;

    // Get and set functions
    int64_t getExponent() const {return exponent;}
    void setExponent(int64_t value) {base = NTL::power(Int(2), value);}
    const Int& getBase() const {return base;}
    const IntVec& getLowBoundaries() const {return lowBoundaries;}
    void setLowBoundaries(const IntVec& values) {lowBoundaries = values;}
    const IntVec& getHighBoundaries() const {return highBoundaries;}
    void setHighBoundaries(const IntVec& values) {highBoundaries = values;}
    const std::vector<int>& getEqualityGroups() const {return equalityGroups;}
    void setEqualityGroups(const std::vector<int>& values) { equalityGroups = values;}
    const std::vector<long>& getNumPowerTwo() const {return numPowerTwo;}
    void setNumPowerTwo(const std::vector<long>& values) {numPowerTwo = values;}
    const std::vector<long>& getMaxPowerTwo() const {return maxPowerTwo;}
    void setMaxPowerTwo(const std::vector<long>& values) {maxPowerTwo = values;}
    PeriodRequirement getPeriodRequirement() const {return periodRequirement;}
    void setPeriodRequirement(PeriodRequirement value) {periodRequirement = value;}
    int64_t getPrimalityTrials() const { return primalityTrials; }
    void setPrimalityTrials(int64_t value) { primalityTrials = value;}
    Int getComponentModulus(int) const override {return this->modulus;}
    long getComponentOrder(int) const override {return this->order;}
    const IntVec& getCurrentCoefficients() const override {return this->aa;}

    /**
    * Search strategies. Each function generates and tests candidate coefficient
    * vectors until its search region or a configured limit is reached.
    *
    * Exhaustive: enumerates the full coefficient range.
    * Random: draws independent vectors within the coefficient bounds.
    * RandomPowerOfTwo: draws vectors as signed sums of powers of two.
    * FromFile reads a .gen file and tests all generators stored therein
    */
    void searchLoopExhaustive(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandom(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandomPowerOfTwo(IntLattice<Int, Real>* proj, int& old);
    void searchLoopFromFile(IntLattice<Int, Real>* proj, int& old);

    /**
     * This function checks whether the generated lattice satisfies the
     * maximal-period requirement for an MRG.
     */
    bool checkMaxPeriod(LCGLattice<Int, Real>& lat) override;
    
    /*
    * Helper function which calculates the groups of equal coefficients 
    * according to the vector equalityGroups
    */
    void prepareEqualityGroups();

private:

    /*
    * This lattice is re-used for storing the current lattice if needed.
    * The function gets the lattice.
    */
    std::unique_ptr<LCGLattice<Int, Real>> workingLattice;
    LCGLattice<Int, Real>& getWorkingLattice() override;

   /*
    * These functions draw random coefficient with additional constaints.
    * drawRandomCoefficients() respects the set equal coefficients.
    * drawPowerOfTwoCoefficient draws a single coefficient which is the sum of powers of two
    */
    void drawRandomCoefficients();
    Int drawPowerOfTwoCoefficient(int i);

    /**
    * Computes the associated LCG modulus for the current MWC coefficients
    * and prepares the lattice, projection, and FOM normalizer for evaluation.
    * 
    * Returns false for inadmissible generator vectors.
    */
    bool prepareCandidate(IntLattice<Int, Real>* proj);
    
    /**
    * Returns whether the current MWC candidate satisfies the configured
    * period requirement. Called only when maximal period is requested.
    */
    bool currentCandidateHasRequiredPeriod();

    int numEqualityGroups = 0;
    IntVec groupLow;
    IntVec groupHigh;
    IntVec groupValue;

};


//===========================================================================


/**
* searchLoopExhaustive performs an exhaustive enumeration of the admissible MWC generators 
* generators defined by the configured bounds.
*/
template<typename Int, typename Real>
void SearchMWC<Int, Real>::searchLoopExhaustive(IntLattice<Int, Real>* proj, int& old)
{
    this->aa.SetLength(this->order + 1);
    for (long j = 0; j <= this->order; ++j) this->aa[j] = lowBoundaries[j];
    unsigned iterationsUntilTimeCheck = 0;
    bool finished = false;
    while (!finished && this->numCandidatesExamined < this->maxNumGen) {
        if (this->timeExpired(iterationsUntilTimeCheck)) break;
        ++this->numCandidatesExamined;
        if (prepareCandidate(proj)) {
            if (!this->maxPeriod || currentCandidateHasRequiredPeriod()) {
                if (this->maxPeriod) ++this->numMaxPeriodCandidates;
                this->testCandidate(proj);
                if (!this->calculateFOM && this->bestCandidates.full()) break;
            }
        }
        if (this->progress)     old = this->printProgress(old);
        // Advance to the next vector; a_0 changes fastest.
        for (long j = 0; j <= this->order; ++j) {
            if (this->aa[j] < highBoundaries[j]) {
                ++this->aa[j];
                break;
            }
            this->aa[j] = lowBoundaries[j];
            if (j == this->order)
                finished = true;
        }
    }
}

//===========================================================================

/**
* sarchLoopRandom generates admissible MRG generators by randomly selecting 
* the recurrence coefficients within the ranges specified in the configuration.
* If the search is restricted to generators with maximal period, then 
* searchLoopRandom only yields generators having this property.
*/
template<typename Int, typename Real>
void SearchMWC<Int, Real>::searchLoopRandom(IntLattice<Int, Real>* proj, int& old)
{
    this->aa.SetLength(this->order + 1);
    prepareEqualityGroups();
    unsigned iterationsUntilTimeCheck = 0;
    while (this->numCandidatesExamined < this->maxNumGen) {
        if (this->timeExpired(iterationsUntilTimeCheck)) break;
        drawRandomCoefficients();
        ++this->numCandidatesExamined;
        if (prepareCandidate(proj)) {
            if (!this->maxPeriod || currentCandidateHasRequiredPeriod()) {
                if (this->maxPeriod)  ++this->numMaxPeriodCandidates;
                this->testCandidate(proj);
                if (!this->calculateFOM && this->bestCandidates.full()) break;
            }
        }
        if (this->progress) old = this->printProgress(old);
    }
}

/**
 * Searches by drawing independent coefficient vectors as signed sums
 * of distinct powers of two.
 * 
 * For each a_i, chooses between 1 and n_i distinct exponents from
 * 0, ..., p_i, then independently chooses the sign of each term.
 * Equality constraints and coefficient boundaries are not used.
 */
template<typename Int, typename Real>
void SearchMWC<Int, Real>::searchLoopRandomPowerOfTwo(IntLattice<Int, Real>* proj, int& old)
{
    this->aa.SetLength(this->order + 1);
    unsigned iterationsUntilTimeCheck = 0;
    while (this->numCandidatesExamined < this->maxNumGen) {
        if (this->timeExpired(iterationsUntilTimeCheck)) break;
        for (int i = 0; i <= this->order; ++i) this->aa[i] = drawPowerOfTwoCoefficient(i);
        ++this->numCandidatesExamined;
        if (prepareCandidate(proj)) {
            if (!this->maxPeriod || currentCandidateHasRequiredPeriod()) {
                if (this->maxPeriod) ++this->numMaxPeriodCandidates;
                this->testCandidate(proj);
                if (!this->calculateFOM && this->bestCandidates.full()) break;
            }
        }        
        if (this->progress) old = this->printProgress(old);
    }
}

//===========================================================================

/**
 * Computes the associated LCG modulus for the current MWC coefficients.
 * Returns false if the candidate is inadmissible. Otherwise, updates the
 * modulus-dependent state used for period and FOM calculations.
 */
template<typename Int, typename Real>
bool SearchMWC<Int, Real>::prepareCandidate(IntLattice<Int, Real>* proj)
{
    // computeLCGModulusMWC requires gcd(a_0, b) == 1.
    if (NTL::GCD(this->aa[0], base) != Int(1)) return false;
    const Int m = computeLCGModulusMWC(base, this->aa);
    if (m <= Int(1)) return false;
    this->modulus = m;
    setModulusIntP<Int>(this->modulus);
    if (this->calculateFOM) {
        proj->setModulus(this->modulus);
        if (this->dualLattice) this->fom->getNormalizer()->computeBounds(-log(this->modulus), 1);
        else this->fom->getNormalizer()->computeBounds(log(this->modulus), 1);
    }
    return true;
}

//===========================================================================

/**
 * Returns the working lattice for the current candidate.
 *
 * Constructs the lattice on first use and updates it with the current
 * coefficients on subsequent calls. This function is called only when
 * computing a figure of merit, so period-only searches do not construct
 * a lattice.
 */
template<typename Int, typename Real>
LCGLattice<Int, Real>& SearchMWC<Int, Real>::getWorkingLattice()
{
    if (!workingLattice) {
        workingLattice = std::make_unique<LCGLattice<Int, Real>>(this->maxdim, this->norm);
    }
    workingLattice->setModulus(this->modulus);
    workingLattice->seta(base);
    return *workingLattice;
}

//===========================================================================

/**
 * Returns whether the current candidate meets the configured MWC period
 * requirement: the safe-prime condition or the full-period test.
 * Assumes prepareCandidate() has already computed the associated modulus.
 */
template<typename Int, typename Real>
bool SearchMWC<Int, Real>::currentCandidateHasRequiredPeriod()
{
    // prepareCandidate() has already computed this->modulus.
    if (periodRequirement == PeriodRequirement::SafePrimeHalfPeriod)
        return mIsSafePrime(this->modulus, primalityTrials);

    return maxPeriodMWC(this->modulus, base, primalityTrials);
}

//===========================================================================

/**
 * Draws one random value for each prepared equality group and copies it into
 * every corresponding position of aa. This produces a coefficient vector
 * satisfying the equality constraints.
 *
 * prepareEqualityGroups() must be called before this function.
 */
template<typename Int, typename Real>
void SearchMWC<Int, Real>::drawRandomCoefficients()
{
    this->aa[0] = randInt(lowBoundaries[0], highBoundaries[0]);
    for (int group = 1; group <= numEqualityGroups; ++group)
        groupValue[group] = randInt(groupLow[group], groupHigh[group]);
    for (int i = 1; i <= this->order; ++i) {
        const int group = equalityGroups.empty() ? i : equalityGroups[i - 1];
        this->aa[i] = groupValue[group];
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
template<typename Int, typename Real> Int SearchMWC<Int, Real>::drawPowerOfTwoCoefficient(int i)
{
    const long ni = numPowerTwo[i];
    const long pi = maxPowerTwo[i];
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

/**
 * Prepares the equality constraints for a search.
 *
 * Groups coefficients with the same equality label and determines the common
 * admissible range for each group by intersecting bounds. An empty 
 * equalityGroups vector treats all coefficients as independent.
 */
template<typename Int, typename Real> void SearchMWC<Int, Real>::prepareEqualityGroups()
{
    const int k = this->order;
    if (!equalityGroups.empty() && equalityGroups.size() != static_cast<std::size_t>(k))
        throw std::invalid_argument("One equality-group label is required per coefficient");
    numEqualityGroups = 0;
    groupLow.SetLength(k + 1);
    groupHigh.SetLength(k + 1);
    groupValue.SetLength(k + 1);
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
* This function checks whether the generated lattice satisfies the
* maximal-period requirement for an MWC.
* Probably this function will not be used since it requires to pass a lattice object which has not been creater.
* Delete it if it is not necessary.
*/
template<typename Int, typename Real>
bool SearchMWC<Int, Real>::checkMaxPeriod(LCGLattice<Int, Real>&)
{
    // Uses the current MWC candidate, prepared by prepareCandidate().
    return currentCandidateHasRequiredPeriod();
}

//===========================================================================


} // namespace LatMRG
#endif // LATMRG_SEARCHMWC_H