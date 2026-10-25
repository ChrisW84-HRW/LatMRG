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
    explicit SearchBase(GenType type) : genType(type), bestLattices(numBestGen, bestFOM) {}

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
    bool progress = true; // Prints the program progress between generators, default = true    
    bool outputToGenFile = false; // When true, the output is written to a .gen file instead of the terminal, default false
    std::string filename = "outputSearch.gen"; // Name of the .gen file, default = "outputSearch.gen"
    
    /**
     * Pointer to a figure of merit object
     */
    FigureOfMeritM<Int, Real>* fom = nullptr;
    

    // Get and set methods    
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

    void setFOM(FigureOfMeritM<Int, Real> & value) {fom = &value;}

    /**
     * These are different virtual methods to be specified by the subclasses.
     * ~ SearhBase is the virtual destructor
     * getComponentModulus() gets the modulus of a single component
     * getComponentOrder() gets the order of a single component
     */
    virtual ~SearchBase() = default;
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

    // Technical variables

    /**
    * Counter for the current RNG.
    */
    int currentGen = 0;         
    
    /**
     * Object to store the information of the lattices with the
     * best figures of merit.
     */    
    MeritList<Lat> bestLattices;

    /**
    * The current lattice which is analyzed.
    */
    Lat* lat = nullptr;

    /**
    * Program time.
    */
    Chrono timer; 

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
    double gen_progress  = double(currentGen) / getMaxNumGen();

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
    
    if (getGenType() == MRG) {
        *out << "Modulo:         m = " << getModulus();
        *out << "\n";
        *out << "Order:          k = " << getOrder() << "\n";
    } else if (getGenType() == MWC) {
    } else if (getGenType() == MMRG) {
    } else if (getGenType() == COMBO) {
      for (int k = 0; k < getNumComp(); k++) {
       if (k > 0) *out << "\n";
       if (getNumComp() >1) *out << "Component " << k+1 << ":\n";
       *out << "Modulo:         m = " << getComponentModulus(k);
       *out << "\n";
       *out << "Order:          k = " << getComponentOrder(k) << "\n";
      //*out << "Generator type: " << toStringGenType(getGenType()) << "\n";
      }
    }
    *out << (getMaxPeriod()?"Full":"Any") << " period length\n";
    *out << "\nTest:\n" << (bestFOM?"Best":"Worst") << " generators \n";
    *out << "Figure of Merit M with vector t=" << fom.getTVector() << "\n";
    *out << "for the " << (dualLattice?"dual":"primal") << " lattice \n";
    *out << "Number of generators kept: " << numBestGen << "\n";
    *out << "Number of generators tested: " << getMaxNumGen() << "\n\n"; 
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
    assert(getCalculateFOM() || getMaxPeriod());

    currentGen = 0;
    resetGeneratorState();
    bestLattices = MeritList<Lat>(numBestGen, bestFOM);

    int old = 0;
    // Launching the tests
    if (progress) {
      old = printProgress(-1);
    }
    timer.init();
    
    IntLattice<Int, Real> proj(getModulus(), fom->getTVector().length(), norm);    
    FigureOfMeritData<Lat> fomData;
    
    // Variable which depends on whether the FoM is calculated for the primal / dual lattice.
    // Avoids duplicate code.
        
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
        fomData.setMerit(fom->computeMerit(*lat, proj));  
        if (dualLattice) lat->buildDualBasis(maxdim);
        else lat->buildBasis(maxdim);      
        fomData.setBasis(dualLattice ? lat->getDualBasis() : lat->getBasis());
        fomData.setMeritProj(fom->getMinMeritProj());
        fomData.setMeritSqlen(fom->getMinMeritSqlen());
        bestLattices.add(fomData);
        if (bestFOM)
          fom->setLowBound(bestLattices.getSmallestMerit());
      }
      else {
        // We need to buiild the basis for the FOM because we currently want to store it.
        // We might want to change this later.
        if (dualLattice) lat->buildDualBasis(getMaxDim());
        else lat->buildBasis(getMaxDim());
        fomData.setBasis(dualLattice ? lat->getDualBasis() : lat->getBasis());
        bestLattices.addWithoutMerit(fomData);
          if (bestLattices.full())
          break;
      }
      // The output is currently only for test purposes
      // *out << fom.computeMerit(*lat, proj) << "\n";
      if (progress) old = printProgress(old);
    } while (!timer.timeOver(getTimeLimit()) && lat);
    old = printProgress(old);
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