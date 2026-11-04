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

    Int getComponentModulus(int i) const override  {return components.at(i).modulus; }

    long getComponentOrder(int i) const override {return components.at(i).order;}

    const IntVec& getCurrentCoefficients() const override {return aa;}

    // Required by SearchBase; combined searching is not implemented yet.
    bool checkMaxPeriod(MRGLattice<Int, Real>&) override
    {
        throw std::logic_error("Combined MRG search is not implemented");
    }

private:
    // Required by SearchBase; combined searching is not implemented yet.
    MRGLattice<Int, Real>& getWorkingLattice() override
    {
        throw std::logic_error("Combined MRG search is not implemented");
    }
};

} // namespace LatMRG

#endif // LATMRG_SEARCHCOMBINEDMRG_H