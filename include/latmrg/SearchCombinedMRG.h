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

    struct ComponentMRG
    {
        Int modulus;
        long order;
        CandidateList<Int> bestCandidates;
    };

    // Input components; each has its own, separate candidate list.
    std::vector<ComponentMRG> components;

    // Coefficients of a combined candidate, for the future search loop.
    IntVec aa;

    // These functions the the modulus, the order and the coefficients of a specific component
    Int getComponentModulus(int i) const override  {return components.at(i).modulus; }
    long getComponentOrder(int i) const override {return components.at(i).order;}
    const IntVec& getCurrentCoefficients() const override {return aa;}

    // Required by SearchBase; combined searching is not implemented yet.
    bool checkMaxPeriod(MRGLattice<Int, Real>&) override
    {
        throw std::logic_error("Combined MRG search is not implemented");
    }

    /**
     * These functions add a candidate list per components, delete the
     * candidate lists for all components or for a specific component
     */
    // 
    void addComponentCandidates(const CandidateList<Int>& candidates);
    void clearComponents();
    void removeComponent(int64_t index);

    /*
    * Calculates the modulus of the combined MRG
    */
    void calculateModulus();

    /**
     * Computes the combined coefficients for one choice of component candidates.
    *
    * selectedCandidates[c] is the index of the candidate chosen
    * from components[c].bestCandidates. The result is stored in aa.
    */
    void computeCombination(const IntVec& selectedCandidates);

private:
    // Required by SearchBase; combined searching is not implemented yet.
    MRGLattice<Int, Real>& getWorkingLattice() override
    {
        throw std::logic_error("Combined MRG search is not implemented");
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

    if (!components.empty() && components.front().order != componentOrder)
        throw std::invalid_argument("All components must have the same order");

    for (const ComponentMRG& component : components) {
        if (NTL::GCD(componentModulus, component.modulus) != Int(1))
            throw std::invalid_argument("Component moduli must be pairwise coprime");
    }

    components.push_back({componentModulus, componentOrder, candidates});
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
    this->order = components.empty() ? 0 : components.front().order;
    this->modulus = Int(0); // Recalculate before starting a search.
}

//============================================================================

template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::calculateModulus()
{
    Int product(1);
    for (std::size_t c = 0; c < components.size(); ++c) {
        const ComponentMRG& component = components[c];
        if (component.order != this->order || component.modulus <= Int(1))
            throw std::invalid_argument("Invalid component order or modulus");
        product *= component.modulus;
    }
    this->modulus = product;
}

//============================================================================

template<typename Int, typename Real> void SearchCombinedMRG<Int, Real>::computeCombination(const IntVec& selectedCandidates)
{
    const Int otherModuli, inverse, weight;
    if (selectedCandidates.length() != static_cast<long>(components.size()))
        throw std::invalid_argument("One index is required per component");

    if (components.empty() || this->order <= 0)
        throw std::invalid_argument("Combined components and order must be configured");

    Int product(1);
    for (long c = 0; c < selectedCandidates.length(); ++c) {
        const ComponentMRG& component =  components[c];
        if (component.order != this->order || component.modulus <= Int(1))
        throw std::invalid_argument("Invalid component order or modulus");
         product *= component.modulus;
    }
    this->modulus = product;

    aa.SetLength(this->order + 1);
    for (long j = 0; j <= this->order; ++j) aa[j] = 0;

    // Adapted from computeCombination() from MRGCombined: 
    // For each component, get the candidate indexed by selectedCandidates[c].
    // Compute the weight for that component's modulus, then add the selected
    // contribution to each combined coefficient.
    for (std::size_t c = 0; c < components.size(); ++c) {
        // Use references here
        const ComponentMRG& component = components[c];
        const std::list<Candidate<Int>>& list = component.bestCandidates.getList();
        const Int& indexValue = selectedCandidates[c];
        auto it = list.begin();
        for (Int n(0); n < indexValue && it != list.end(); ++n) ++it;
        if (it == list.end()) throw std::out_of_range("Candidate index is out of range");
        const Candidate<Int>& candidate = *it;
        otherModuli = this->modulus / component.modulus;
        inverse = NTL::InvMod(otherModuli % component.modulus, component.modulus);
        weight = otherModuli * inverse;
        for (long j = 1; j <= this->order; ++j) {
            aa[j] = (aa[j] + candidate.getAa()[j] * weight) % this->modulus;
            if (aa[j] < 0) aa[j] += this->modulus;
        }
    }
}
} // namespace LatMRG

#endif // LATMRG_SEARCHCOMBINEDMRG_H