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

#ifndef LATMRG_GENFILEREADER_H
#define LATMRG_GENFILEREADER_H

#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "latticetester/FlexTypes.h"

namespace LatMRG {

/**
 * Reads generators from the old LatMRG .gen format.
 *
 * Header:
 *   number of components
 *   modulus of component 1
 *   order of component 1
 *   ...
 *   number of generators
 *
 * Each generator then contains the coefficients of component 1, followed
 * by those of component 2, and so on. Coefficients are read as Int;
 * blank lines between them are ignored.
 */
template<typename Int>
class GenFileReader
{
public:
    explicit GenFileReader(const std::string& filename): file(filename)
    {
        if (!file) throw std::runtime_error("Failed to open .gen file: " + filename);
        long componentCount = 0;
        file >> componentCount;
        components.reserve(static_cast<std::size_t>(componentCount));
        for (long c = 0; c < componentCount; ++c) {
            ComponentInfo component;
            if (!(file >> component.modulus) || component.modulus <= Int(1)) throw std::runtime_error("Invalid component modulus in .gen file");
            // If there is explanatory text to the modulus
            // line, the reader needs only the initial integer.
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            if (!(file >> component.order) || component.order <= 0) throw std::runtime_error("Invalid component order in .gen file");
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            components.push_back(std::move(component));
        }
        long generatorCount;
        if (!(file >> generatorCount) || generatorCount < 0) throw std::runtime_error("Invalid generator count in .gen file");
        numGenerators = static_cast<std::size_t>(generatorCount);
    }

    std::size_t getNumComponents() const
    {
        return components.size();
    }

    const Int& getComponentModulus(std::size_t c) const
    {
        return components.at(c).modulus;
    }

    long getComponentOrder(std::size_t c) const
    {
        return components.at(c).order;
    }

    std::size_t getNumGenerators() const
    {
        return numGenerators;
    }

    std::size_t getNumGeneratorsRead() const
    {
        return generatorsRead;
    }

    /**
     * Reads the next complete generator into componentCoefficients.
     *
     * Returns false after the number of generators declared in the header
     * has been read. Throws if the file ends in the middle of a generator.
     * Each vector is indexed from 1 through its component's order; entry
     * zero is initialized to zero.
     */
    bool readNext(std::vector<IntVec>& componentCoefficients)
    {
        if (generatorsRead == numGenerators)
            return false;
        std::vector<IntVec> result;
        result.reserve(components.size());
        for (const ComponentInfo& component : components) {
            IntVec coefficients;
            coefficients.SetLength(component.order + 1);
            coefficients[0] = Int(0);
            for (long i = 1; i <= component.order; ++i) {
                if (!(file >> coefficients[i]))  throw std::runtime_error("Incomplete generator in .gen file");
            }
            result.push_back(std::move(coefficients));
        }
        componentCoefficients = std::move(result);
        ++generatorsRead;
        return true;
    }

private:
    struct ComponentInfo
    {
        Int modulus;
        long order = 0;
    };
    std::ifstream file;
    std::vector<ComponentInfo> components;
    std::size_t numGenerators = 0;
    std::size_t generatorsRead = 0;
};

} // namespace LatMRG

#endif // LATMRG_GENFILEREADER_H