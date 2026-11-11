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

#ifndef LATMRG_SEARCHCOMBINEDMRG_H
#define LATMRG_SEARCHCOMBINEDMRG_H

#include <stdexcept>
#include <vector>
#include <NTL/ZZ.h>

#include "latmrg/SearchBase.h"

namespace LatMRG {

template<typename Int, typename Real>
struct SearchCombinedMRG : public SearchBase<Int, Real, MRGLattice<Int, Real>>
{

    SearchCombinedMRG() : SearchBase<Int, Real, MRGLattice<Int, Real>>(COMBO) {}

    // We need to make sure the combined modulus is calculated before performSearch is executed.
    template<typename SearchType> int performSearch(void (SearchType::*searchLoop)(IntLattice<Int, Real>*, int&))
    {
        calculateModulus();
        WorkingLattice.reset();
        return SearchBase<Int, Real, MRGLattice<Int, Real>>::template performSearch<SearchType>(searchLoop);
    }

    // These functions the the modulus, the order and the coefficients of a specific component
    Int getComponentModulus(int i) const override  {return components.at(i).getModulus(); }
    long getComponentOrder(int i) const override {return components.at(i).getOrder();}
    const IntVec& getCurrentCoefficients() const override {return this->aa;}

    /**
     * Search strategies. Each function generates and tests candidate coefficient
     * vectors until its search region or a configured limit is reached.
     * searchLoopDeterministic calculates all combinations of the stored list of candidates.
     * SearchLoopRandom draws random combinations of the stored lists of candidates.
     */
    void searchLoopDeterministic(IntLattice<Int, Real>* proj, int& old);
    void searchLoopRandom(IntLattice<Int, Real>* proj, int& old);

    // Required by SearchBase; combined searching is not implemented yet.
    bool checkMaxPeriod(MRGLattice<Int, Real>&) override;

    /**
     * These functions add a candidate list per components, delete the
     * candidate lists for all components or for a specific component
     */
    // 
    void addComponentCandidates(const CandidateList<Int>& candidates);
    void clearComponents();
    void removeComponent(int64_t index);

    /*
    * Calculates the modulus of the combined MRG.
    */
    void calculateModulus();

    /**
     * Computes the combined coefficients for one choice of component candidates.
    *
    * selectedCandidates[c] is the index of the candidate chosen
    * from the components. The result is stored in aa.
    */
    void computeCombination(const IntVec& selectedCandidates);
    
    /*
    * This lattice is re-used for storing the current lattice if needed.
    * The function gets the lattice.
    */
    std::unique_ptr<MRGLattice<Int, Real>> WorkingLattice;

    // Input components; each has its own, separate candidate list.
    std::vector<CandidateList<Int>> components;

    
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

    // Number of selected candidates per component in random search
    int64_t numSelectedPerComponent = 1000;    
    int64_t getNumSelectedPerComponent() const {return numSelectedPerComponent;}
    void setnNmSelectedPerComponent(int64_t & values) {numSelectedPerComponent = values;}

private:
    // Returns the current combined lattice
    MRGLattice<Int, Real>& getWorkingLattice() override
    {
        if (!WorkingLattice) {
            WorkingLattice = std::make_unique<MRGLattice<Int, Real>>(this->modulus, this->aa, this->maxdim);
        } else {
            WorkingLattice->setaa(this->aa);
            WorkingLattice->setModulus(this->modulus);
        }
        return *WorkingLattice;
    }
};

//============================================================================
// Implementation

template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::addComponentCandidates(const CandidateList<Int>& candidates)
{
    const auto& list = candidates.getList();
    if (list.empty()) throw std::invalid_argument("Cannot add an empty candidate list");

    const Candidate<Int>& first = list.front();
    const Int componentModulus = first.getModulus();
    const long componentOrder = first.getAa().length() - 1;

    if (!components.empty() && components.front().getOrder() != componentOrder)
        throw std::invalid_argument("All components must have the same order");

    for (const CandidateList<Int>& component : components) {
        if (NTL::GCD(componentModulus, component.getModulus()) != Int(1))
            throw std::invalid_argument("Component moduli must be pairwise coprime");
    }

    components.push_back(candidates);
    this->numComp = static_cast<long>(components.size());
    this->order = componentOrder;
    this->modulus = Int(0); // Recalculate before starting a search.
}

//============================================================================

template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::clearComponents()
{
    components.clear();
    this->numComp = 0;
    this->order = 0;
    this->modulus = Int(0);
}

//============================================================================

template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::removeComponent(int64_t index)
{
    if (index >= components.size())
        throw std::out_of_range("Component index is out of range");

    components.erase(components.begin() + index);
    this->numComp = static_cast<long>(components.size());
    this->order = components.empty() ? 0 : components.front().getOrder();
    this->modulus = Int(0); // Recalculate before starting a search.
}

//============================================================================

/*
 * Calculates the modulus of the combined MRG.
*/
template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::calculateModulus()
{
    Int product(1);
    for (std::size_t c = 0; c < components.size(); ++c) {
        const CandidateList<Int>& component = components[c];
        if (component.getOrder() != this->order || component.getModulus() <= Int(1))
            throw std::invalid_argument("Invalid component order or modulus");
        product *= component.getModulus();
    }
    this->modulus = product;
    setModulusIntP<Int>(this->modulus);
}

//============================================================================

/**
 * Computes the coefficients of one combined MRG candidate.
 * 'selectedCandidates' contains the number of the selected
 * candidate in the list for each component.
*/
template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::computeCombination(const IntVec& selectedCandidates)
{
    const Int otherModuli, inverse, weight;
    if (selectedCandidates.length() != static_cast<long>(components.size()))
        throw std::invalid_argument("One index is required per component");

    if (components.empty() || this->order <= 0)
        throw std::invalid_argument("Combined components and order must be configured");

    Int product(1);
    for (long c = 0; c < selectedCandidates.length(); ++c) {
        const CandidateList<Int>& component = components[c];
        if (component.getOrder() != this->order || component.getModulus() <= Int(1))
        throw std::invalid_argument("Invalid component order or modulus");
         product *= component.getModulus();
    }
    this->modulus = product;

    this->aa.SetLength(this->order + 1);
    for (long j = 0; j <= this->order; ++j) this->aa[j] = 0;

    // Adapted from computeCombination() from MRGCombined: 
    // For each component, get the candidate indexed by selectedCandidates[c].
    // Compute the weight for that component's modulus, then add the selected
    // contribution to each combined coefficient.
    for (std::size_t c = 0; c < components.size(); ++c) {
        // Use references here
        const CandidateList<Int>& component = components[c];
        const std::list<Candidate<Int>>& list = component.getList();
        const Int& indexValue = selectedCandidates[c];
        auto it = list.begin();
        for (Int n(0); n < indexValue && it != list.end(); ++n) ++it;
        if (it == list.end()) throw std::out_of_range("Candidate index is out of range");
        const Candidate<Int>& candidate = *it;
        otherModuli = this->modulus / component.getModulus();
        inverse = NTL::InvMod(otherModuli % component.getModulus(), component.getModulus());
        weight = otherModuli * inverse;
        for (long j = 1; j <= this->order; ++j) {
            this->aa[j] = (this->aa[j] + candidate.getAa()[j] * weight) % this->modulus;
            if (this->aa[j] < 0) this->aa[j] += this->modulus;
        }
    }
}

//============================================================================

/**
* This function checks whether the generated lattice satisfies the
* maximal-period requirement for an MRG.
* Probably this function will not be used since it requires to pass a lattice object which has not been creater.
* Delete it if it is not necessary.
*/
template<typename Int, typename Real> bool SearchCombinedMRG<Int, Real>::checkMaxPeriod(MRGLattice<Int, Real>& lat)
{
  return getMRGComponent().maxPeriod(lat.getaa());
}

//============================================================================

/**
 * Tests combinations of the candidate lists supplied for the components.
 *
 * Enumerates the Cartesian product of the lists, selecting one candidate
 * from each component. For every selection, computes the combined
 * coefficient vector, evaluates its figure of merit, and retains it if
 * it ranks among the best combined generators.
 *
 * Assumes the supplied lists contain full-period component candidates;
 * it does not search for components or recheck their periods.
 */
template<typename Int, typename Real>
void SearchCombinedMRG<Int, Real>::searchLoopDeterministic(IntLattice<Int, Real>* proj, int& old)
{
    // positions is an odometer: one position per component list.
    std::vector<std::size_t> positions(components.size(), 0);
    IntVec selectedCandidates;
    selectedCandidates.SetLength(static_cast<long>(components.size()));
    unsigned iterationsUntilTimeCheck = 0;

    while (!this->reachedNumGeneratorLimit()) {
        // Check before the first combination, then once every 1,000.
        if (iterationsUntilTimeCheck == 0) {
            if (this->timer.timeOver(this->timeLimit))
                break;
            iterationsUntilTimeCheck = 1000;
        }
        --iterationsUntilTimeCheck;
        for (std::size_t c = 0; c < components.size(); ++c)
            selectedCandidates[static_cast<long>(c)] = NTL::conv<Int>(static_cast<long>(positions[c]));
        computeCombination(selectedCandidates);
        ++this->numCandidatesExamined;
        // The input lists are assumed to contain only full-period component
        // candidates, as supplied by SearchMRG.
        bool accepted = true;
        if (this->maxPeriod) {
            accepted = getMRGComponent().maxPeriod(this->aa);
            if (accepted) ++this->numMaxPeriodCandidates;
        }
        if (this->progress) old = this->printProgress(old);
        if (accepted) {
            this->testCandidate(proj);
            if (!this->calculateFOM && this->bestCandidates.full()) break;
        }
        // Advance the odometer. The last component varies fastest.
        std::size_t c = components.size();
        while (c > 0) {
            --c;
            if (++positions[c] < components[c].getList().size()) break;
            positions[c] = 0;
        }
        // Every combination has been visited.
        if (c == 0 && positions[0] == 0)  break;
    }
}
//============================================================================

/**
 * Randomly tests combinations of the supplied component candidate lists.
 *
 * On each iteration, draws one candidate from each component list,
 * computes the combined coefficient vector, evaluates its figure of merit,
 * and retains it if it ranks among the best combined generators. Draws are 
 * made with replacement, so a combination may be tested more than once. 
 *
 * Assumes the supplied lists already contain full-period components;
 * it does not generate or rank component candidates.
 */
template<typename Int, typename Real>
void SearchCombinedMRG<Int, Real>::searchLoopRandom(
    IntLattice<Int, Real>* proj, int& old)
{
    IntVec selectedCandidates;
    selectedCandidates.SetLength(components.size());
    unsigned iterationsUntilTimeCheck = 0;

    std::size_t combinationLimit = 1;
    const auto generatorLimit = static_cast<std::size_t>(this->maxNumGen);

    for (std::size_t c = 0; c < components.size(); ++c) {
        if (combinationLimit > generatorLimit / numSelectedPerComponent) {
            combinationLimit = generatorLimit;
            break;
        }
        combinationLimit *= numSelectedPerComponent;
    }

    while (!this->reachedNumGeneratorLimit() && this->numCandidatesExamined < combinationLimit) {
        // Check before the first draw, then every 1,000 draws.
        if (iterationsUntilTimeCheck == 0) {
            if (this->timer.timeOver(this->timeLimit))
                break;
            iterationsUntilTimeCheck = 1000;
        }
        --iterationsUntilTimeCheck;
        // The loop randomly selects one candidate from the list of each compoent
        // and stores its index in selectedCandidates
        for (std::size_t c = 0; c < components.size(); ++c) {
            const long lastIndex = components[c].getList().size() - 1;
            selectedCandidates[c] = randInt(0, lastIndex);
        }
        computeCombination(selectedCandidates);
        ++this->numCandidatesExamined;
        // Assumes the supplied lists contain full-period components and
        // that no additional combined-period check is required.
        if (this->maxPeriod) {
            if (!getMRGComponent().maxPeriod(this->aa)) {
                if (this->progress) old = this->printProgress(old);
                continue;
            }
            ++this->numMaxPeriodCandidates;
        }
        if (this->progress) old = this->printProgress(old);
        testCandidate(proj);
        if (!this->calculateFOM && this->bestCandidates.full()) break;
    }
}//============================================================================

} // namespace LatMRG

#endif // LATMRG_SEARCHCOMBINEDMRG_H