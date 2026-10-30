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

#ifndef LATMRG_SEARCHBASE_H
#define LATMRG_SEARCHBASE_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <functional>
#include <ctime>
#include <cstdlib>
#include <list>
#include <vector>
#include <string>
#include <cstdint>
#include <memory>

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsUniform.h"

#include "latmrg/Candidate.h"
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
 * SearchBase contains the common configuration data and search logic used by
 * all (RNG type specific) search classes.
 *
 * Generator-specific subclasses provide functions that generate candidate
 * coefficient vectors. SearchBase evaluates their figures of merit when
 * requested, constructing a working lattice only for that purpose.
 */
template<typename Int, typename Real, typename Lat>
struct SearchBase
{
    /**
     * We need to make explicit the type of the Configuration to have it available for necessary case distinctions.
     */
    explicit SearchBase(GenType type) : genType(type), bestCandidates(numBestGen, bestFOM) {}

    // Lattice specific parameters
    GenType genType;   // Type of generator, currently MRG, MWC or COMBO.
    long numComp = 1; // Number of components.  If > 1, we have a combined generator.    
    bool dualLattice = false; // Work with primal/dual lattice
    NormType norm = L2NORM; // Norm Type
    int64_t maxdim = 0; // Maximal dimension of the lattice
    Int modulus;   // The modulus of the lattice
    long order;    // The order k of the generator

    // Settings for the search
    bool bestFOM = true; // Want to keep the best/worst FoMs, default is true
    bool maxPeriod = true; // When true, we retain only max-period generators. default is true.
    int numBestGen = 5; // Number of best generators to be kept, default is 5
    bool calculateFOM = true; // Decides whether we want to calculate the FOM
    int maxNumGen = 1000; // Maximal number of generators to be tested, default = 1000

    // Technical setting
    double timeLimit = 600.0; // CPU time limit in seconds. default = 600
    bool print_time = true; // Decides if running time is printed in printResults
    bool progress = false; // Prints the program progress between generators, default = false    
    bool outputToGenFile = false; // When true, the relevant output is written to a .gen file, default false
    std::string filename = "outputSearch.gen"; // Name of the .gen file, default = "outputSearch.gen"
    int verbose = 1; // Level of verbosity of the printed output, default = 1
    
    /**
     * Pointer to a figure of merit object
     */
    FigureOfMeritM<Int, Real>* fom = nullptr;
    

    // Get and set functions    
    GenType getGenType() const {return genType;}
    void setGenType(GenType value) {genType = value;}
    long getNumComp() const {return numComp;}
    void setNumComp(long value) {numComp = value;}
    bool getDualLattice() const {return dualLattice;}
    void setDualLattice(bool value) {dualLattice = value;}
    NormType getNorm() const { return norm;}
    void setNorm(NormType value) { norm = value;}
    int64_t getMaxDim() const {return maxdim;}
    void setMaxDim(int64_t value) {maxdim = value;}    
    Int getModulus() const {return modulus;}
    void setModulus(Int value) {modulus = value;}
    long getOrder() const {return order;}
    void setOrder(long value) {order = value;}

    bool getBestFOM() const {return bestFOM;}
    void setBestFOM(bool value) {bestFOM = value;}
    bool getMaxPeriod() const {return maxPeriod;}
    void setMaxPeriod(bool value) {maxPeriod = value;}
    int getNumBestGen() const {return numBestGen;}
    void setNumBestGen(int value) {numBestGen = value;}
    bool getCalculateFOM() const {return calculateFOM;}
    void setCalculateFOM(bool value) {calculateFOM = value;}
    int getMaxNumGen() const {return maxNumGen;}
    void setMaxNumGen(int value) {maxNumGen = value;}

    double getTimeLimit() const {return timeLimit;}
    void setTimeLimit(double value) {timeLimit = value;}
    bool getPrintTime() const { return print_time;}
    void setPrintTime(bool value) { print_time = value;}
    bool getProgress() const {return progress;}
    void setProgress(bool value) {progress = value;}
    bool getOutputToGenFile() const {return outputToGenFile;}
    void setOutputToGenFile(bool value) {outputToGenFile = value;}
    const std::string& getFilename() const {return filename;}
    void setFilename(const std::string& value) {filename = value;}
    int getVerbose() const {return verbose;}
    void setVerbose(int value) {verbose = value;}

    void setFOM(FigureOfMeritM<Int, Real> & value) {fom = &value;}

    /**
     * These are different virtual functions to be specified by the subclasses.
     * ~ SearhBase is the virtual destructor
     * getComponentModulus() gets the modulus of a single component
     * getComponentOrder() gets the order of a single component
     * getCurrentCoefficients() allows to make an independent copy of the vector of multipliers.
     */
    virtual ~SearchBase() = default;
    virtual Int getComponentModulus(int k) const = 0;
    virtual long getComponentOrder(int k) const = 0;
    virtual const IntVec& getCurrentCoefficients() const = 0;
    
   /**
   * Resets generator type specific search state before a new search is started.
   *
   * The base implementation does nothing. Subclasses override this function
   * to reset state variables used by their generator routines, such as
   * enumeration counters.
   */
   virtual void resetSearchState() {}

  /**
  * Performs a search using a search-loop function from a derived class.
  *
  * The search loop generates and counts candidates, checks period and stopping
  * conditions, and calls testCandidate() for candidates to be retained.
  */
  template<typename SearchType> int performSearch(void (SearchType::*searchLoop)(IntLattice<Int, Real>*, int&));

  /**
  * Evaluates and retains a candidate coefficient vector.
  *
  * If FOM calculation is enabled, constructs or updates the working lattice,
  * computes its merit, and retains it if it ranks among the best candidates.
  * Otherwise, retains the coefficients without constructing a lattice.
  */ 
   void testCandidate(const IntVec& coefficients, IntLattice<Int, Real>* proj);
   
   /**
   * This function checks whether the generated lattice satisfies the
   * maximal-period requirement (if applicable). The implementation
   * defers per lattice type and is provided in the subclasses.
   */
   virtual bool checkMaxPeriod(Lat& lat) = 0;

   /*
   * Checks if we the upper bound on the number of generators to be checked has been
   * reached. Depends on whether we restrict to maximal period generators.
   */
   bool reachedNumGeneratorLimit() const {return (maxPeriod ? numMaxPeriodCandidates : numCandidatesExamined) >= maxNumGen;}
    
    /**
    * This function reports the progress of the current search. 
    */
    int printProgress(int old);

    /**
    * This function prints the results of the random search. It
    * should only be called after PerformSearch() has been executed. 
    */
    void printResults();

    /**
     * Writes the output to a .gen file.
     */
    void writeGenFile() const;

    /**
    * Returns the merit of the best retained lattice.
    * Returns NaN if no lattice has been retained.
    */
    double getBestMerit();

    // Technical variables

    /*
    * This lattice is re-used for FoM calculation.
    */
    virtual Lat& getWorkingLattice() = 0;

    /**
    * Counter for the number of RNGs examined.
    */
    int numCandidatesExamined = 0;   
    
    /**
    * Counter for the number of RNGs examined which have maximal period.
    */
    int numMaxPeriodCandidates = 0;
    
    /**
     * Object to store the information of the lattices with the
     * best figures of merit.
     */    
    CandidateList<Int> bestCandidates;

    /**
    * Program time.
    */
    Chrono timer; 
   
};
    
//============================================================================
// Implementation

  
  /**
  * This function is supposed to report the progress of the current  
  * search. 
  */
  template<typename Int, typename Real, typename Lat> int SearchBase<Int, Real, Lat>::printProgress(int old) {  
    // Progress is measured with respect to both stopping criteria. Since the
    // search terminates as soon as either the time limit or the generator limit
    // is reached, we display the larger of the two relative progresses.
    double time_progress = double(timer.val(Chrono::SEC)) / timeLimit;
    double gen_progress = double(maxPeriod ? numMaxPeriodCandidates : numCandidatesExamined) / maxNumGen;
    int per_80 = static_cast<int>(80 * std::max(time_progress, gen_progress));
    if (per_80 > 80) per_80 = 80;
    if (per_80 < 0) per_80 = 0;  
    // We do not print for no reason as this slows the program a lot.
    if (per_80 <= old) return old;
    std::cout << "Program progress: [";
    for (int i = 0; i < per_80; i++) std::cout << "#";
    for (int i = per_80; i < 80; i++) std::cout << " ";
    std::cout << "] ";std::cout << std::setw(3) << int(per_80/80.0*100) << " %";
    if (per_80 == 80)
        std::cout << '\n';
    else
        std::cout << '\r';
    std::cout << std::flush;
    return per_80;
    
  }

  //===========================================================================
  
  /** 
   * This function prints the results of the search.
   * The function is not yet complete but rather a prototype
   * to which additional information can be added later.
  */
  template<typename Int, typename Real, typename Lat> void SearchBase<Int, Real, Lat>::printResults() { 
    std::cout << "\nSearch: A search program for Random Number Generators\n";
    std::cout << ((numComp>1)?"Combined generators":"Simple generator")
      << " configuration" << ((numComp>1)?"s":"") << "\n\n";
    if (genType == MRG) {
        std::cout << "Modulo:         m = " << modulus;
        std::cout << "\n";
        std::cout << "Order:          k = " << order << "\n";
    } else if (genType == MWC) {
    } else if (genType == MMRG) {
    } else if (genType == COMBO) {
      for (int k = 0; k < numComp; k++) {
       if (k > 0) std::cout << "\n";
       if (numComp >1) std::cout << "Component " << k+1 << ":\n";
       std::cout << "Modulo:         m = " << getComponentModulus(k);
       std::cout << "\n";
       std::cout << "Order:          k = " << getComponentOrder(k) << "\n";
      //std::cout << "Generator type: " << toStringGenType(getGenType()) << "\n";
      }
    }
    std::cout << (maxPeriod?"Full":"Any") << " period length\n";
    std::cout << "\nTest:\n" << (bestFOM?"Best":"Worst") << " generators \n";
    if (calculateFOM) std::cout << "Figure of Merit M with vector t=" << fom->getTVector() << "\n";
    std::cout << "for the " << (dualLattice?"dual":"primal") << " lattice \n";
    std::cout << "Number of generators kept: " << numBestGen << "\n";
    std::cout << "Number of generators tested: " << maxNumGen << "\n\n"; 
    if (print_time) {
      std::cout << "Allowed running time: " << timeLimit << "s.\n";
      std::cout << "Actual CPU time: " << timer.toString() << "\n\n";
    }    
    std::cout << "Retained generators (from best to worst):\n";
    for (const auto& candidate : bestCandidates.getList()) {
      std::cout << "Multipliers: " << candidate.getAa() << "\n";
      std::cout << "Merit:      " << candidate.getMerit() << "\n\n";
    }
  }

  //===========================================================================

/**
 * Writes the retained candidates to a .gen file.
 * The format is provisional and needs to be checked against the old .gen format.
 */
template<typename Int, typename Real, typename Lat>
void SearchBase<Int, Real, Lat>::writeGenFile() const
{
    std::ofstream genFile(filename, std::ios::out | std::ios::trunc);
    if (!genFile) {
        std::cerr << "Failed to open " << filename << '\n';
        return;
    }
    for (const auto& candidate : bestCandidates.getList()) {
        genFile << toStringGenType(genType) << '\n';
        genFile << candidate.getModulus() << '\n';
        genFile << candidate.getAa() << '\n';
    }
}


  //===========================================================================

/**
 * If a figure of merit is requested, the function obtains the working lattice, computes
 * the merit using proj, and retains the candidate if it belongs among the
 * best found so far. Otherwise, stores the candidate with merit zero while
 * space remains in the list.
 *
 * The search loop is responsible for period checks, candidate counters, and
 * stopping conditions; this function does not perform those tasks.
 * The coefficients must match those of the current working lattice.
 */
template<typename Int, typename Real, typename Lat> 
void SearchBase<Int, Real, Lat>::testCandidate(const IntVec& coefficients, IntLattice<Int, Real>* proj)
{
    Candidate<Int> candidate{modulus, coefficients, 0.0};
    if (!calculateFOM) {
        bestCandidates.addWithoutMerit(candidate);
        return;
    }
    assert(fom != nullptr && proj != nullptr);
    Lat& lattice = getWorkingLattice();
    candidate.merit = fom->computeMerit(lattice, *proj);
    // For now, this ranking condition supports best-FOM searches.
    assert(bestFOM);
    const bool retained = !bestCandidates.full() || candidate.getMerit() > bestCandidates.getCutoffMerit();
    if (retained) {
        bestCandidates.add(candidate);
        if (verbose > 2) {
            std::cout << "New element in the list of best candidates:\n"
                      << "Multipliers: " << candidate.getAa() << '\n'
                      << "Merit:      " << candidate.getMerit() << "\n\n";
        }
    }
    if (bestCandidates.full()) fom->setLowBound(bestCandidates.getCutoffMerit());
}

  //===========================================================================

template<typename Int, typename Real, typename Lat> template<typename SearchType> 
int SearchBase<Int, Real, Lat>::performSearch( void (SearchType::*searchLoop)(IntLattice<Int, Real>*, int&))
{
    assert(numComp > 0);
    assert(calculateFOM || maxPeriod);
    numCandidatesExamined = 0;
    numMaxPeriodCandidates = 0;
    resetSearchState();
    bestCandidates = CandidateList<Int>(numBestGen, bestFOM);
    std::unique_ptr<IntLattice<Int, Real>> proj;
    if (calculateFOM) {
        assert(fom != nullptr);
        proj = std::make_unique<IntLattice<Int, Real>>(modulus, fom->getTVector().length(), norm);
    }
    // Initialize the modulus in NTL modular integer arithmetic.
    // used for checking maximal period condition.
    setModulusIntP<Int>(modulus);
    timer.init();
    int old = 0;
    if (progress)
        old = printProgress(-1);
    (static_cast<SearchType*>(this)->*searchLoop)(proj.get(), old);
    if (progress) printProgress(old);
    if (verbose > 1) printResults();
    if (outputToGenFile) writeGenFile();
    return 0;
}
  //===========================================================================

  /**
  * This function returns the merit of the best retained lattice.
  * Returns NaN if no lattice has been retained.
  */
 template<typename Int, typename Real, typename Lat> double SearchBase<Int, Real, Lat>::getBestMerit()
  {
    const auto& list = bestCandidates.getList();
    if (list.empty())
        return std::numeric_limits<double>::quiet_NaN();
    return list.front().merit;
  }

} // namespace LatMRG

#endif // LATMRG_SEARCHBASE_H