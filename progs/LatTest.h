#ifndef LATMRG_LATTEST_H
#define LATMRG_LATTEST_H

extern std::ostream* out;
extern bool print_time;

template<typename Int, typename Real> struct LatTest {
  typedef NTL::vector<Int> IntVec;

  ConfigLat<Int, Real> conf;  // This is defined in ExecCommon.h
  // MeritList  is defined in FigureOfMerit.h.

  Chrono timer;
  std::vector<bool> full_period;

  /**
   * Prints the results of the program execution.
   */
  template<typename Lat>
    void printResults(MeritList<Lat>& bestLattice) {
      *out << "LatTest: A program to test the lattice structure of a linear RNG.\n";
      *out << delim;
      *out << ((conf.num_comp>1)?"Combined generators":"Simple generator")
        << " configuration" << ((conf.num_comp>1)?"s":"") << "\n\n";
      for (int k = 0; k < conf.num_comp; k++) {
        if (k > 0) *out << "\n";
        if (conf.num_comp >1) *out << "Component " << k+1 << ":\n";
        *out << "Generator type: " << toStringGen(conf.fact[k]->get_type()) << "\n";
        *out << "Modulo:         m = " ;
        if (conf.fact[k]->get_type() == MWC) *out << conf.fact[k]->m_MWCb;
        else *out << conf.fact[k]->getM();
        *out << " = " << conf.fact[k]->getB() << "^" << conf.fact[k]->getE();
        if (conf.fact[k]->getR() > 0) *out << "+" << conf.fact[k]->getR();
        if (conf.fact[k]->getR() < 0) *out << conf.fact[k]->getR();
        *out << "\n";
        *out << "Order:          k = " << conf.fact[k]->getK() << "\n";
        *out << (conf.period[0]?"Check":"Don't check") << " full period length\n";
      }
      *out << "\nTest:\n";
      if (conf.criterion == LatticeTester::SPECTRAL) {
        *out << "Spectral Test";
        if (conf.normaType != LatticeTester::NONE) *out << "\nNormalizer used: "
          << toStringNorma(conf.normaType);
      } else if (conf.criterion == LatticeTester::BEYER) *out << "Beyer quotient";
      else if (conf.criterion == LatticeTester::LENGTH) *out << "Shortest vector length";
      *out << "\n\n";
      *out << "Dimensions and projections:\n";
      *out << conf.proj->toString();
      if (print_time) {
        *out << delim;
        *out << "Allowed running time: " << conf.timeLimit << "s.\n";
        *out << "Actual CPU time: " << timer.toString() << "\n";
      }
      for (auto it = bestLattice.getList().begin(); it!= bestLattice.getList().end(); it++) {
        *out << delim;
        *out << (*it).getLattice() << "\n";
        if (conf.num_comp > 1) {
          bool print = false;
          for (int i = 0; i<conf.num_comp; i++) {
            if (conf.period[i]){
              *out << "Component " << i+1
              << ((full_period[i])?" has":" does not have") << " full period.\n";
              print = true;
            }
          }
          if (print) *out << "\n";
        } else {
          if (conf.period[0]) *out << "Full period: " << (full_period[0]?"yes":"no") << "\n\n";
        }
        if (conf.detail == 0) {
          *out << (*it).toStringMerit();
        } else if (conf.detail == 1) {
          *out << (*it).toStringDim();
        } else if (conf.detail == 2) {
          *out << (*it).toStringProjections();
        }
      }
    }

  //==========================================================================

  int TestLat () {
    if (!conf.gen_set) {
      std::cerr << "No generator set for in lattest tag. Aborting.\n";
      return 1;
    }
    if (!(conf.test_set)) {
      std::cerr << "No test set for in lattest tag. Aborting.\n";
      return 1;
    }
    if (!conf.proj_set) {
      std::cerr << "No projections set for in lattest tag. Aborting.\n";
      return 1;
    }
    // Initializing values
    timer.init();

    // Testing the generator(s)
    // Generators are initiated with dim = conf.proj->numProj() because this is
    // the best way to make sure the normalizer construction will work while
    // being efficient in the execution. This is so that we will precompute all
    // bounds we will use more than once, but do not instantiate the object with
    // a dimension so big it throws an error.
    if (conf.num_comp > 1) {
      full_period.resize(conf.num_comp);   // Number of components.
      // Checking full period of components that require it
      for (int i = 0; i < conf.num_comp; i++) {
        IntVec temp(conf.fact[i]->getK()+1);
        temp[0] = Int(0);
        for (int j = 1; j < conf.fact[i]->getK()+1; j++) temp[j] = conf.coeff[i][j-1];
        if (conf.period[i]) full_period[i] = conf.fact[i]->maxPeriod(temp);
        // This is not set because it is not needed in Seek
        conf.fact[i]->setA(conf.coeff[i]);
      }
      // Combined generators case
      MeritList<ComboLattice<Int, Real>> bestLattice(conf.max_gen, true);
      MRGLattice<Int, Real>* mrg = getLatCombo<Int, Real>(conf.fact, conf.proj->numProj());
      ComboLattice<Int, Real> combolat(conf.fact, *mrg);
      bestLattice.add(test_lat(combolat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MRG || conf.fact[0]->get_type() == LCG) {
      // Single MRG component.
      full_period.resize(1);  // This becomes a single boolean variable.
      MeritList<MRGLattice<Int, Real>> bestLattice(conf.max_gen, true);
      IntVec temp(conf.fact[0]->getK()+1);
      temp[0] = Int(0);
      for (int i = 1; i < conf.fact[0]->getK()+1; i++) temp[i] = conf.coeff[0][i-1];
      if (conf.fact[0]->get_type() == MRG && conf.period[0]) {
        full_period[0] = conf.fact[0]->maxPeriod(temp);
      } else if (conf.period[0]) {
        full_period[0] = conf.fact[0]->maxPeriod(conf.coeff[0][0]);
      }
      MRGLattice<Int, Real> mrglat(conf.fact[0]->getM(), temp, conf.proj->numProj(), conf.fact[0]->getK(), FULL, conf.norm);
      bestLattice.add(test_lat(mrglat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MWC) {
      full_period.resize(1);
      MeritList<MWCLattice<Int, Real>> bestLattice(conf.max_gen, true);
      MWCLattice<Int, Real> mwclat(conf.fact[0]->m_MWCb, conf.fact[0]->getM(), conf.proj->numProj());
      if (conf.period[0]) full_period[0] = conf.fact[0]->maxPeriod(mwclat.getCoef());
      bestLattice.add(test_lat(mwclat, conf));
      printResults(bestLattice);
    } else if (conf.fact[0]->get_type() == MMRG) {
      full_period.resize(1);
      if (conf.period[0]) full_period[0] = conf.fact[0]->maxPeriod(conf.fact[0]->getMatrix());
      MeritList<MMRGLattice<Int, Real>> bestLattice(conf.max_gen, true);
      MMRGLattice<Int, Real> mmrglat(conf.fact[0]->getM(), conf.fact[0]->getMatrix(), conf.proj->numProj(), conf.fact[0]->getK(), conf.norm);
      bestLattice.add(test_lat(mmrglat, conf));
      printResults(bestLattice);
    }
    delete conf.proj;
    for (int i = 0; i < conf.num_comp; i++) delete conf.fact[i];
    return 0;
  }

}; // end struct LatTest
#endif
