  //===========================================================================

  /**
  * This function performs the search based on the current configuration.
  * The parameter *nextGenerator defines the search strategy, e.g.,
  * a random or an exhaustive search.
  * 
  * Each call updates the current coefficient vector and returns true when
  * a candidate is available, or false when the search is exhausted.
  * The template parameter SearchType allows SearchBase to call that member
  * function without knowing the concrete derived search class.
  */
  template<typename Int, typename Real, typename Lat> template<typename SearchType>
  int SearchBase<Int, Real, Lat>::performSearchNextGenerator(bool (SearchType::*nextGenerator)())  {      
    assert(numComp > 0);
    assert(calculateFOM || maxPeriod);
    numCandidatesExamined = 0;
    numMaxPeriodCandidates = 0;
    resetSearchState();
    bestCandidates = CandidateList<Int>(numBestGen, bestFOM);
    int old = 0;
    // Launching the tests
    if (progress) {
      old = printProgress(-1);
    }
    timer.init();    
    std::unique_ptr<IntLattice<Int, Real>> proj;
    if (calculateFOM) {
      assert(fom != nullptr);
      proj = std::make_unique<IntLattice<Int, Real>>(modulus, fom->getTVector().length(), norm);
    }     
    // Set the modulus used by the modular-arithmetic in NTL, which are
    // called when checking whether a candidate has maximal period.
    setModulusIntP<Int>(modulus);   
    // Loop through all lattices
    do {      
      // Generate the next coefficient vector. False means the search is exhausted..
      if (!(static_cast<SearchType*>(this)->*nextGenerator)())
        break;
      Candidate<Int> candidate{getModulus(), getCurrentCoefficients(), 0.0}; // Copies aa; does not retain a reference to it.
      /** If the FOM is calculated, then the lattice is only added if its FOM is among 
      *  the best ones. If we do not calculate the FOM, we just add the lattices with 
      * maximal period to the list.
      */
      if (calculateFOM) {
        Lat& lat = getWorkingLattice();
        candidate.merit = fom->computeMerit(lat, *proj);
        const bool retained = !bestCandidates.full() || candidate.getMerit() > bestCandidates.getCutoffMerit();
        if (retained) {
          bestCandidates.add(candidate);
          if (verbose > 2) {
            std::cout << "New element in the list of best candidates:\n";
            std::cout << "Multipliers: " << candidate.getAa() << "\n";
            std::cout << "Merit:      " << candidate.getMerit() << "\n\n";
          }
        }
        // For a best-FOM search, the cutoff is the worst merit retained,
        // and should be used only once the list is full.
        if (bestCandidates.full())
          fom->setLowBound(bestCandidates.getCutoffMerit());          
      } else {
        bestCandidates.addWithoutMerit(candidate);
        if (bestCandidates.full())
          break;
      }
      if (progress) old = printProgress(old);
    } while (!timer.timeOver(getTimeLimit()) && !this->reachedNumGeneratorLimit());
    if (progress) old = printProgress(old);
    if (verbose > 1) printResults();
    if (outputToGenFile) writeGenFile();
    return 0;
  }
