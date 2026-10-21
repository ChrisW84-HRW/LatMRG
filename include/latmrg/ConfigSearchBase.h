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

#ifndef LATMRG_CONFIGSEARCHBASE_H
#define LATMRG_CONFIGSEARCHBASE_H

#include <ctime>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cstdint>
#include <memory>

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"

#include "latmrg/EnumTypes.h"
#include "latmrg/MWCLattice.h"
#include "latmrg/MRGLattice.h"
#include "latmrg/MRGComponent.h"

using LatticeTester::IntLatticeExt;
using LatticeTester::Normalizer;
using namespace std;
using namespace LatMRG;
using namespace LatticeTester::Random;

namespace LatMRG {

/**
 * ConfigMerit stores the information on the figure of merit
 * and how it is computed.
 */
template<typename Int, typename Real>
struct ConfigMerit
{
    bool calculateFOM = false;

    NTL::Vec<int64_t> t;

    LatticeTester::Weights* weights = nullptr;
    LatticeTester::Normalizer* norma = nullptr;
    LatticeTester::ReducerBB<Int, Real>* red = nullptr;

    bool includeFirst = true;
    bool dualLattice = false;

    LatticeTester::NormType norm = LatticeTester::L2NORM;

    double minMerit = 0.0;
    double maxMerit = 1.0;

    #ifdef LATMRG_SEARCH
        bool best = true;
        int num_gen = 0;
        double currentMerit = double(0);
    #endif
    
    int no_bestGen = 5; // Number of best generators to be kept 
};


/**
 * Common configuration details for all seearch operations.
 *
 * This struct contains only settings that are independent of
 * the generator type.
 */
template<typename Int, typename Real>
struct ConfigSearch
{
    /**
     * We need to make explicit the type of the Configuration to have it available for necessary case distinctions.
     */
    explicit ConfigSearch(GenType type) : genType(type) {}

    GenType genType;   // Type of generator, currently MRG or MWC.

    long numComp = 1; // Number of components.  If > 1, we have a combined generator.
    int64_t maxdim = 0; // Maximal dimension of the lattice

    bool permax = true; // When true, we retain only max-period generators. Default is true.

    bool outputToGenFile = false; // When true, the output is written to a .gen file instead of the terminal
    std::string filename; // Name of the .gen file

    int max_gen = 10; // Maximal number of generators to be tested

    LatticeType latticeType; // The type of lattice that is constructed (see EnumTypes).
    SearchMethod searchMethod;  // The method of search for the search (exhaustive or random).

    double timeLimit = 600.0; // CPU time limit in seconds.

    #ifdef LATMRG_SEARCH
        bool progress = true; // Prints the program progress between generators
    #endif

    ConfigMerit<Int, Real> configFOM;

    // Get and set methods
    GenType getGenType() const {return genType;}
    void setGenType(GenType value) {genType = value;}
    long getNumComp() const {return numComp;}
    void setNumComp(long value) {numComp = value;}
    int64_t getMaxDim() const {return maxdim;}
    void setMaxDim(int64_t value) {maxdim = value;}
    bool getPermax() const {return permax;}
    void setPermax(bool value) {permax = value;}
    bool getOutputToGenFile() const {return outputToGenFile;}
    void setOutputToGenFile(bool value) {outputToGenFile = value;}
    const std::string& getFilename() const {return filename;}
    void setFilename(const std::string& value) {filename = value;}
    int getMaxGen() const {return max_gen;}
    void setMaxGen(int value) {max_gen = value;}
    LatticeType getLatticeType() const {return latticeType;}
    void setLatticeType(LatticeType value) {latticeType = value;}
    SearchMethod getSearchMethod() const {return searchMethod;}
    void setSearchMethod(SearchMethod value) {searchMethod = value;}
    double getTimeLimit() const {return timeLimit;}
    void setTimeLimit(double value) {timeLimit = value;}
    #ifdef LATMRG_SEARCH
        bool getProgress() const {return progress;}
        void setProgress(bool value) {progress = value;}
    #endif
    bool onlyMaxPeriod() const { return permax;}

    // Get and set functions for the figure of merit
    void setFOM(const NTL::Vec<int64_t>& tValue, LatticeTester::Normalizer* normaValue, 
        LatticeTester::ReducerBB<Int, Real>* redValue, LatticeTester::Weights* weightsValue)
    {
        configFOM.t = tValue;
        configFOM.norma = normaValue;
        configFOM.red = redValue;
        configFOM.weights = weightsValue;
    }
    void setCalculateFOM(bool value) {configFOM.calculateFOM = value;}
    bool getCalculateFOM() const {return configFOM.calculateFOM;}
    const NTL::Vec<int64_t>& getT() const {return configFOM.t;}
    void setT(const NTL::Vec<int64_t>& value) {configFOM.t = value;}
    LatticeTester::Weights* getWeights() const {return configFOM.weights;}
    void setWeights(LatticeTester::Weights* value) {configFOM.weights = value;}
    LatticeTester::Normalizer* getNorma() const {return configFOM.norma;}
    void setNorma(LatticeTester::Normalizer* value) {configFOM.norma = value;}
    LatticeTester::ReducerBB<Int, Real>* getRed() const {return configFOM.red;}
    void setRed(LatticeTester::ReducerBB<Int, Real>* value) {configFOM.red = value;}
    bool getIncludeFirst() const {return configFOM.includeFirst;}
    void setIncludeFirst(bool value) {configFOM.includeFirst = value;}
    bool getDualLattice() const {return configFOM.dualLattice;}
    void setDualLattice(bool value) {configFOM.dualLattice = value;}
    LatticeTester::NormType getNorm() const {return configFOM.norm;}
    void setNorm(LatticeTester::NormType value) {configFOM.norm = value;}
    double getMinMerit() const {return configFOM.minMerit;}
    void setMinMerit(double value) {configFOM.minMerit = value;}
    double getMaxMerit() const {return configFOM.maxMerit;}
    void setMaxMerit(double value) {configFOM.maxMerit = value;}
    #ifdef LATMRG_SEARCH
        bool getBest() const {return configFOM.best;}
        void setBest(bool value) {configFOM.best = value;}
        int getNumGen() const {return configFOM.num_gen;}
        void setNumGen(int value) {configFOM.num_gen = value;}
        double getCurrentMerit() const {return configFOM.currentMerit;}
        void setCurrentMerit(double value) {configFOM.currentMerit = value;}
    #endif
    int getNoBestGen() const {return configFOM.no_bestGen;}
    void setNoBestGen(int value) {configFOM.no_bestGen = value;}

    /**
     * These are different virtual methods to be specified by the subclasses.
     * ~ ConfigSearch is the virtual destructor
     * createComponents() initializes the correct number of components as defined by numComp
     * getModulus() gets the modulus
     * getComponentModulus() gets the modulus of a single component
     * getComponentOrder() gets the order of a single component
     */
    virtual ~ConfigSearch() = default;
    virtual Int getModulus() const = 0;
    virtual Int getComponentModulus(int k) const = 0;
    virtual long getComponentOrder(int k) const = 0;
};

} // namespace LatMRG

#endif // LATMRG_CONFIGSEARCH_H