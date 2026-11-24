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

#ifndef LATMRG_SEARCHGENFILE_H
#define LATMRG_SEARCHGENFILE_H

#include <string>
#include <utility>

#include "latmrg/IntFactorization.h"
#include "latmrg/GenFileReader.h"
#include "latmrg/SearchMRG.h"
#include "latmrg/SearchCombinedMRG.h"

namespace LatMRG {

/**
 * Coordinates a search over generators stored in a .gen file.
 *
 * GenFileReader reads the header during construction. 'run' selects
 * SearchMRG for one component or SearchCombinedMRG for multiple
 * components and then calls the corresponding search function. The same
 * reader object is passed on so generator records can be read without
 * reopening the file.
 */
template<typename Int, typename Real> class SearchGenFile
{
public:
    explicit SearchGenFile(const std::string& filename) : reader(filename) {}

    template<typename SimpleSearch, typename CombinedSearch>
    int run(SimpleSearch&& searchSimple, CombinedSearch&& searchCombined)
    {
        if (reader.getNumComponents() == 1) {
            SearchMRG<Int, Real> search;  
            return std::forward<SimpleSearch>(searchSimple)(search, reader);
        }
        SearchCombinedMRG<Int, Real> search;
        return std::forward<CombinedSearch>(searchCombined)(search, reader);
    }

private:
    GenFileReader<Int> reader;
};

} // namespace LatMRG

#endif // LATMRG_SEARCHGENFILE_H