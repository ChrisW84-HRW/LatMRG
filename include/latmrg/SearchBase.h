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

#include "latmrg/FigureOfMeritData.h"
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

    bool best = true;
    int num_gen = 0;
    double currentMerit = double(0);
    
    
    int no_bestGen = 5; // Number of best generators to be kept 
};


/**
 * SearchBase contains the common configuration data and search logic used by
 * all (RNG type specific) search classes.
 *
 * Generator-specific subclasses provide the methods that construct candidate
 * lattices, while SearchBase implements the shared functionality such as the
 * main facility to perform searches, storage of the best retained lattices, 
 * progress reporting, and result printing.
 */
template<typename Int, typename Real, typename Lat>
struct SearchBase
{
    /**
     * We need to make explicit the type of the Configuration to have it available for necessary case distinctions.
     */
    explicit SearchBase(GenType type) : genType(type), bestLattices(configFOM.no_bestGen, configFOM.best) {}

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

    bool progress = true; // Prints the program progress between generators
    
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
    bool getProgress() const {return progress;}
    void setProgress(bool value) {progress = value;}
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
    bool getBest() const {return configFOM.best;}
    void setBest(bool value) {configFOM.best = value;}
    int getNumGen() const {return configFOM.num_gen;}
    void setNumGen(int value) {configFOM.num_gen = value;}
    double getCurrentMerit() const {return configFOM.currentMerit;}
    void setCurrentMerit(double value) {configFOM.currentMerit = value;}
    int getNoBestGen() const {return configFOM.no_bestGen;}
    void setNoBestGen(int value) {configFOM.no_bestGen = value;}

    /**
     * These are different virtual methods to be specified by the subclasses.
     * ~ SearhBase is the virtual destructor
     * createComponents() initializes the correct number of components as defined by numComp
     * getModulus() gets the modulus
     * getComponentModulus() gets the modulus of a single component
     * getComponentOrder() gets the order of a single component
     */
    virtual ~SearchBase() = default;
    virtual Int getModulus() const = 0;
    virtual Int getComponentModulus(int k) const = 0;
    virtual long getComponentOrder(int k) const = 0;

  /**
   * This abstract method looks for the next generator based on the chosen configuration.
   * It uses the variable currentGen to calculate the next generator in the list.
   * nextGeneratorExhaustive() is performing an exhaustive search. This method needs to
   * be implemented separately for the different types of lattices.
   */
   virtual Lat* nextGeneratorExhaustive() = 0;

   /**
   * As nextGeneratorExhaustive() but it chooses a random next generator which is within
   * the range defined by the configuration.
   */
   virtual Lat* nextGeneratorRandom() = 0;  
   
   /**
   * Resets generator type specific search state before a new search is started.
   *
   * The base implementation does nothing. Subclasses override this method
   * to reset state variables used by their generator routines, such as
   * enumeration counters.
   */
   virtual void resetGeneratorState() {}

   /**
   * This method performs the search based on the current configuration.
   * The parameter *generator defines the method of the search, e.g.,
   * a random or an exhaustive search.
   * The extra template parameter Generator allows to pass a pointer to a member function 
   * of the derived class, for example SearchMRG::nextGeneratorRandomBox, without hard-coding 
   * the derived type in the base class. Thus the base algorithm remains generic and reusable, 
   * while the generator-specific enumeration logic stays in the subclasses.
   */
   template<typename Generator> int performSearch(Lat* (Generator::*generator)());
   
   /**
   * This method checks whether the generated lattice satisfies the
   * maximal-period requirement (if applicable). The implementation
   * defers per lattice type and is provided in the subclasses.
   */
   virtual bool checkMaxPeriod(Lat& lat) = 0;
    
    /**
    * This method reports the progress of the current search. 
    */
    int printProgress(int old);

    /**
    * This method prints the results of the random search. It
    * should only be called after PerformSearch() has been executed. 
    */
    void printResults();

    /**
    * Returns the merit of the best retained lattice.
    * Returns NaN if no lattice has been retained.
    */
    double getBestMerit();

    /**
    * Counter for the current RNG.
    */
    int currentGen = 0;         
    
    /**
     * Object to store the information of the lattices with the
     * best figures of merit.
     */    
    MeritList<Lat> bestLattices;

    
    private:
        
    /**
    * The current lattice which is analyzed.
    */
   Lat* lat = nullptr;

    /**
    * Figure of merit objects for the primal and dual lattice.
    */        
    // FigureOfMeritM<Int, Real> fomPrimal;
    // FigureOfMeritDualM<Int, Real> fomDual;

    /**
    * Program time.
    */
    Chrono timer; 

    /**
     * Decides if runnning time is printed in printResults().
     */
    bool print_time = true;

    /**
     * Output streams for terminal and file
     */
    std::ostream* out = &std::cout;
    std::ofstream file;
   

};
    
//============================================================================
// Implementation

  
  /**
  * This method is supposed to report the progress of the current  
  * search. 
  */
  template<typename Int, typename Real, typename Lat> int SearchBase<Int, Real, Lat>::printProgress(int old) {  
    // Progress is measured with respect to both stopping criteria. Since the
    // search terminates as soon as either the time limit or the generator limit
    // is reached, we display the larger of the two relative progresses.
    double time_progress = double(timer.val(Chrono::SEC)) / getTimeLimit();
    double gen_progress  = double(currentGen) / getMaxGen();

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
   * This method prints the results of the search.
   * The function is not yet complete but rather a prototype
   * to which additional information can be added later.
  */
  template<typename Int, typename Real, typename Lat> void SearchBase<Int, Real, Lat>::printResults() { 

    if (getOutputToGenFile()) {
        file.open(getFilename(), std::ios::out | std::ios::trunc);
        if (!file) {
            std::cerr << "Failed to open " << getFilename() << "\n";
        }
        out = &file;
    }
    
    *out << "\nSearch: A search program for Random Number Generators\n";
    *out << ((getNumComp()>1)?"Combined generators":"Simple generator")
      << " configuration" << ((getNumComp()>1)?"s":"") << "\n\n";

      
    
    for (int k = 0; k < getNumComp(); k++) {
      if (k > 0) *out << "\n";
      if (getNumComp() >1) *out << "Component " << k+1 << ":\n";
      *out << "Generator type: " << toStringGenType(getGenType()) << "\n";
      if (getGenType() == MRG) {
        *out << "Modulo:         m = " << getComponentModulus(k); //" = " << conf.genComponents[k]->getB() << "^"
          //<< conf.fact[k]->getE();
        // if (conf.fact[k]->getR() > 0) *out << "+" << conf.genComponents[k]->getR();
        // if (conf.fact[k]->getR() < 0) *out << conf.genComponents[k]->getR();

        *out << "\n";
        *out << "Order:          k = " << getComponentOrder(k) << "\n";
      } else if (getGenType() == MWC) {
      } else if (getGenType() == MMRG) {
      }
      *out << (onlyMaxPeriod()?"Full":"Any") << " period length\n";
      }
    
    *out << "\nTest:\n" << (configFOM.best?"Best":"Worst") << " generators \n";
    *out << "Figure of Merit M with vector t=" << configFOM.t << "\n";
    *out << "for the " << (configFOM.dualLattice?"dual":"primal") << " lattice \n";
    *out << "Number of generators kept: " << configFOM.no_bestGen << "\n";
    *out << "Number of generators tested: " << getMaxGen() << "\n\n"; 
    if (print_time) {
      *out << "Allowed running time: " << getTimeLimit() << "s.\n";
      *out << "Actual CPU time: " << timer.toString() << "\n\n";
    }
    

    *out << "Retained generators (from best to worst):\n";
    for (auto it = bestLattices.getList().begin(); it!= bestLattices.getList().end(); it++) {
      *out << (*it).getBasis() << "\n";
      *out << (*it).toStringMerit() << "\n";
    }
    
  
  }

  //===========================================================================

  /**
  * This method performs the search based on the current configuration.
  * The parameter *generator defines the method of the search, e.g.,
  * a random or an exhaustive search.
  */
  template<typename Int, typename Real, typename Lat> template<typename Generator>
  int SearchBase<Int, Real, Lat>::performSearch(Lat* (Generator::*generator)())  {  
    
    assert(getNumComp() > 0);
    assert(getCalculateFOM() || getPermax());

    currentGen = 0;
    resetGeneratorState();
    bestLattices = MeritList<Lat>(configFOM.no_bestGen, configFOM.best);

    int old = 0;
    // Launching the tests
    if (progress) {
      old = printProgress(-1);
    }
    timer.init();
    
    IntLattice<Int, Real> proj(getModulus(), configFOM.t.length(), configFOM.norm);    
    FigureOfMeritData<Lat> fomData;
    
    FigureOfMeritM<Int, Real> fomPrimal(configFOM.t, *configFOM.weights, *configFOM.norma, configFOM.red, configFOM.includeFirst);
    FigureOfMeritDualM<Int, Real> fomDual(configFOM.t, *configFOM.weights, *configFOM.norma, configFOM.red, configFOM.includeFirst);
    // Variable which depends on whether the FoM is calculated for the primal / dual lattice.
    // Avoids duplicate code.
    
    auto& fom = configFOM.dualLattice ? fomDual : fomPrimal;  
    
    setModulusIntP<Int>(getModulus());   

    // Loop through all lattices
    do {      
      lat = (static_cast<Generator*>(this)->*generator)();
      if (!lat) continue;   

      fomData.setLattice(lat);
      
      /** If the FOM is calculated, then the lattice is only added if its FOM is among 
       *  the best ones. If we do not calculate the FOM, we just add the lattices with 
       * maximal period to the list.
      */
      if (getCalculateFOM()) { 
        fomData.setMerit(fom.computeMerit(*lat, proj));  
        if (configFOM.dualLattice) lat->buildDualBasis(maxdim);
        else lat->buildBasis(maxdim);      
        fomData.setBasis(configFOM.dualLattice ? lat->getDualBasis() : lat->getBasis());
        fomData.setMeritProj(fom.getMinMeritProj());
        fomData.setMeritSqlen(fom.getMinMeritSqlen());
        bestLattices.add(fomData);
        if (configFOM.best)
          fom.setLowBound(bestLattices.getSmallestMerit());
        configFOM.currentMerit = bestLattices.getMerit();
      }
      else {
        // We need to buiild the basis for the FOM because we currently want to store it.
        // We might want to change this later.
        if (configFOM.dualLattice) lat->buildDualBasis(getMaxDim());
        else lat->buildBasis(getMaxDim());
        fomData.setBasis(configFOM.dualLattice ? lat->getDualBasis() : lat->getBasis());
        bestLattices.addWithoutMerit(fomData);
          if (bestLattices.full())
          break;
      }
      configFOM.num_gen++;
      // The output is currently only for test purposes
      // *out << fom.computeMerit(*lat, proj) << "\n";
      if (progress) old = printProgress(old);
    } while (!timer.timeOver(getTimeLimit()) && lat);
    old = printProgress(old);
    //std::cout << "\n \n \n !!! END !!! \n \n \n";
    //std::cout << "CurrentGen: " << currentGen;
    //std::cout << "\n \n \n !!! END !!! \n \n \n";
    return 0;
  }

  /**
  * This method returns the merit of the best retained lattice.
  * Returns NaN if no lattice has been retained.
  */
  template<typename Int, typename Real, typename Lat> double SearchBase<Int, Real, Lat>::getBestMerit()
  {
    auto& list = bestLattices.getList();
    if (list.empty())
        return std::numeric_limits<double>::quiet_NaN();
    return list.front().getMerit();
  }


} // namespace LatMRG

#endif // LATMRG_SEARCHBASE_H