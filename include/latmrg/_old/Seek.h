// This file is part of LatMRG.
//
// Copyright (C) 2012-2022  The LatMRG authors, under the occasional supervision
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

#ifndef LATMRG_SEEK_H
#define LATMRG_SEEK_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <functional>

#include "latmrg/PrimesFinder.h"
#include "latmrg/ConfigSearchBase.h"
#include "latmrg/FigureOfMeritData.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsUniform.h"

#include "latmrg/LCGLattice.h"
#include "latmrg/LCGComponent.h"
#include "latmrg/MWCComponent.h"
#include "latmrg/MWCLattice.h"

namespace LatMRG {

 /**
 * \class Seek
 *
 * This abstract class performs an exhaustive or random search for RNGs.
 *
 * A search consists of generating candidate lattices one by one,
 * optionally verifying that the corresponding generator has maximal
 * period, and computing a figure of merit for each admissible lattice.
 * Depending on the configuration, the figure of merit is evaluated
 * either for the primal or the dual lattice. The best generators
 * encountered during the search are retained throughout the execution.
 *
 * The constructor initializes the search from a ConfigSearch object.
 * In particular, it stores the search configuration, constructs the
 * required figure of merit objects, and initializes the auxiliary
 * objects needed during the search.
 *
 * The method performSeek() is the main facility of the class. It
 * repeatedly calls a next generator function supplied by the user, evaluates
 * the figure of merit of the resulting lattices, and keeps the best
 * generators found during the search. The generator function is passed
 * as a pointer to a member function, see below, allowing the search procedure to
 * use either exhaustive enumeration or random generation while keeping
 * the search algorithm independent of the lattice type.
 *
 * The methods nextGenerator() and nextGeneratorRandom() define how
 * candidate generators are produced. The former performs an exhaustive
 * enumeration of the admissible generator space defined by the
 * configuration, whereas the latter generates candidates randomly within
 * the prescribed ranges. These methods are pure virtual and must be
 * implemented in the subclasses corresponding to the different lattice
 * types (e.g., MRG or MWC lattices).
 *
 * Thus, the subclasses are responsible only for the construction of
 * admissible candidate lattices, while the common search logic, including
 * merit computation, period checking, ranking, and progress reporting,
 * is implemented in this class.
 */

  template<typename Lat> class Seek {
  
         
    public:

    
   /**
   * Initialization of objects is based on the ConfigSearch parameter.
   */
   Seek(ConfigSearch<Int, Real>& config) : conf(config)
                                        , bestLattices(conf.configFOM.no_bestGen, conf.configFOM.best) 
                                        , fomPrimal(config.configFOM.t, *config.configFOM.weights, *config.configFOM.norma, config.configFOM.red, config.configFOM.includeFirst)
                                        , fomDual(config.configFOM.t, *config.configFOM.weights, *config.configFOM.norma, config.configFOM.red, config.configFOM.includeFirst)
                                        {            
                                          fomPrimal.setLowBound(config.configFOM.minMerit);  
                                          fomDual.setLowBound(config.configFOM.minMerit);     
                                          // Set if the output is written to a file (or the terminal) 
                                          if (config.outputToGenFile) {
                                             file.open(config.filename, std::ios::out | std::ios::trunc);
                                            if (!file) {
                                              std::cerr << "Failed to open " << config.filename << "\n";
                                            }
                                           out = &file;
                                          }                            
                                        }
    
   /**
    * Destructor
    */
   virtual ~Seek() = default;
                                          
   /**
   * This abstract method looks for the next generator based on the chosen configuration.
   * It uses the variable currentGen to calculate the next generator in the list.
   * nextGenerator() is performing an exhaustive search. This method needs to
   * be implemented separately for the different types of lattices.
   */
   virtual Lat* nextGenerator() = 0;

   /**
   * As nextGenerator() but it chooses a random next generator which is within
   * the range defined by the configuration.
   */
   virtual Lat* nextGeneratorRandom() = 0;   

   /**
   * This method performs the seek based on the current configuration.
   * The parameter *generator defines the method of the seek, e.g.,
   * a random or an exhaustive search.
   */
   template<typename Generator> int performSeek(Lat* (Generator::*generator)());
   
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
    * should only be called after PerformSeek() has been executed. 
    */
    void printResults();

    /**
    * The configuration for the seek.
    */
    ConfigSearch<Int, Real> & conf;
    
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
    std::unique_ptr<Lat> lat;

    /**
    * Figure of merit objects for the primal and dual lattice.
    */        
    FigureOfMeritM<Int, Real> fomPrimal;
    FigureOfMeritDualM<Int, Real> fomDual;

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
  template<typename Lat> int Seek<Lat>::printProgress(int old) {  
    // Progress is measured with respect to both stopping criteria. Since the
    // search terminates as soon as either the time limit or the generator limit
    // is reached, we display the larger of the two relative progresses.
    double time_progress = double(timer.val(Chrono::SEC)) / conf.timeLimit;
    double gen_progress  = double(currentGen) / conf.max_gen;

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
  template<typename Lat> void Seek<Lat>::printResults() {
    
    *out << "\nSeek: A search program for Random Number Generators\n";
    *out << ((conf.numComp>1)?"Combined generators":"Simple generator")
      << " configuration" << ((conf.numComp>1)?"s":"") << "\n\n";

      
    
    for (int k = 0; k < conf.numComp; k++) {
      if (k > 0) *out << "\n";
      if (conf.numComp >1) *out << "Component " << k+1 << ":\n";
      *out << "Generator type: " << toStringGenType(conf.genType) << "\n";
      if (conf.genType == MRG) {
        *out << "Modulo:         m = " << conf.getComponentModulus(k); //" = " << conf.genComponents[k]->getB() << "^"
          //<< conf.fact[k]->getE();
        // if (conf.fact[k]->getR() > 0) *out << "+" << conf.genComponents[k]->getR();
        // if (conf.fact[k]->getR() < 0) *out << conf.genComponents[k]->getR();

        *out << "\n";
        *out << "Order:          k = " << conf.getComponentOrder(k) << "\n";
      } else if (conf.genType == MWC) {
      } else if (conf.genType == MMRG) {
      }
      *out << (conf.onlyMaxPeriod()?"Full":"Any") << " period length\n";
      }

    
    *out << "\nTest:\n" << (conf.configFOM.best?"Best":"Worst") << " generators \n";
    *out << "Figure of Merit M with vector t=" << conf.configFOM.t << "\n";
    *out << "for the " << (conf.configFOM.dualLattice?"dual":"primal") << " lattice \n";
    *out << "Number of generators kept: " << conf.configFOM.no_bestGen << "\n";
    *out << "Number of generators tested: " << conf.max_gen << "\n\n"; 
    if (print_time) {
      *out << "Allowed running time: " << conf.timeLimit << "s.\n";
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
  * This method performs the seek based on the current configuration.
  * The parameter *generator defines the method of the seek, e.g.,
  * a random or an exhaustive search.
  */
  template<typename Lat> template<typename Generator> int Seek<Lat>::performSeek(Lat* (Generator::*generator)())  {  
    assert(conf.numComp > 0);
    int old = 0;
    // Launching the tests
    if (conf.progress) {
      old = printProgress(-1);
    }
    timer.init();
    
    IntLattice<Int, Real> proj(conf.getModulus(), conf.configFOM.t.length(), conf.configFOM.norm);
    
    FigureOfMeritData<Lat> fomData;

    // Preparation for being able to check primitivity
    setModulusIntP<Int>(conf.getModulus());
    
    // Loop through all lattices
    do {      
      lat.reset((static_cast<Generator*>(this)->*generator)());
      if (!lat) continue;   
      // If the period must be maximal, test if the specific component has max period.
      // Otherwise continue.
      if (conf.onlyMaxPeriod())
      { 
        if (!checkMaxPeriod(*lat))
          continue;
      }
      
      // Variable which depends on whether the FoM is calculated for the primal / dual lattice.
      // Avoids duplicate code.
      auto& fom = conf.configFOM.dualLattice ? fomDual : fomPrimal;       
      fomData.setMerit(fom.computeMerit(*lat, proj));
      fomData.setLattice(lat.get());
      fomData.setBasis(conf.configFOM.dualLattice ? lat->getDualBasis() : lat->getBasis());
      fomData.setMeritProj(fom.getMinMeritProj());
      fomData.setMeritSqlen(fom.getMinMeritSqlen());

      bestLattices.add(fomData); 
      // Update the lower bound for the FoM to equal the smallest stored FoM
      if (conf.configFOM.best)
        fom.setLowBound(bestLattices.getSmallestMerit());

      conf.configFOM.num_gen++;
      conf.configFOM.currentMerit = bestLattices.getMerit(); 
      // The output is currently only for test purposes
      // *out << fom.computeMerit(*lat, proj) << "\n";
      if (conf.progress) old = printProgress(old);
    } while (!timer.timeOver(conf.timeLimit) && lat);
    old = printProgress(old);
    return 0;
  }

}// End namespace LatMRG
#endif
