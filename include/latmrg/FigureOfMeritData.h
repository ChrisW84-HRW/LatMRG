#ifndef LATMRG_FIGUREOFMERIT_H
#define LATMRG_FIGUREOFMERIT_H

#include <list>
#include <sstream>
#include <iomanip>

#include "latticetester/Coordinates.h"

namespace LatMRG {

  /**
   * This class implements the basic structures for computing and storing a general
   * figure of merit defined as a minimum or maximum over a set of projections.
   * It is intended as a return type for functions that actually perform the computations.
   * 
   * In the class, the figure of merit, the worst projection and the shortest vector
   * length for this projection are stored.
   */
template<typename Lat>
class FigureOfMeritData {


   private:
      /*
      * Pointer to the lattice that is being tested.
      */
      const Lat* m_lattice = nullptr;

      /*
      * A copy of the basis of the lattice.
      */
      IntMat m_basis;

      /*
      * Used to store the merit of the lattice
      **/
      Real m_merit = 0;

      /*
      * Used to store the minimal squared vector length
      */
      double m_minMeritSqlen = 0.0;

      /* 
       * Projection for which the worst FoM is achieved
      */
      Coordinates m_minMeritProj;

      /*
      * Index of the best or worst value in `m_merits`.
      */
      long m_best_worst = 0;

      std::string m_latticeString;

      // Some variables are missing.  For example, we need beta_1, we need a normalizer, etc.   *********

   public:

      /**
       * get and set functions the relevant variables.
       * */
      Real getMerit() const { return m_merit; }

      void setMerit(Real m) { m_merit = m; }

      IntMat getBasis() const { return m_basis; }

      void setBasis(IntMat basis) { m_basis = basis; }

      double getMeritSqlen() const { return m_minMeritSqlen;}

      void setMeritSqlen(double len) {m_minMeritSqlen = len;} 

      Coordinates getMeritProj() const { return m_minMeritProj;}

      void setMeritProj(const Coordinates coord) {m_minMeritProj = coord;}

      /**
       * Returns a string containing the merit of the generator,  the projection
       * for which this is obtained and the shortest vector for this projection.
       * */
      std::string toStringMerit() const{
        std::ostringstream stream;
        stream << "Merit: " << getMerit() << "\n" << "Worst Projection: "
          << getMeritProj() << "\n" << "Shortest Vector for this projection: "
          << getMeritSqlen() << "\n";
        return stream.str();
      }

      /// Returns the string associated with this test.
      std::string getLattice() const { return m_latticeString;}

      void setLattice(Lat* lattice) { m_lattice = lattice;}    
      
      void storeLatticeAsString() { m_latticeString = m_lattice->toString(); }

    };


  /**
   * This class stores a certain number of `FigureOfMeritData` objects.
   * New elements are only added if the list is not full or if they have a better
   * merit than what is on the list.
   * */
template<typename Lat>
class MeritList{
      public:
      /**
       * Creates a `MeritList` that will hold at most `maxLength` `FigureOfMeritData`.
       * */
      MeritList(unsigned long maxLength, bool best){
        m_max = maxLength;
        m_merit = best?Real(0):Real(1);
        m_best = best;
      }

      /**
       * Adds `test` in `m_tests` if there is space left or if `test` as a
       * higher (lower) merit.
       * */
      void add(FigureOfMeritData<Lat> test) {
        if ((m_tests.size() >= m_max) && !(m_best ^ (test.getMerit() > m_tests.back().getMerit()))) {
          m_tests.pop_back();
        } else if (m_tests.size() >= m_max) return;
        posInsert(test);
      }

      unsigned long getSize(){
        return m_tests.size();
      }

      std::list<FigureOfMeritData<Lat>>& getList() {return m_tests;}

      Real getMerit(){
        return m_merit;
      }
      
      
      /**
      * Get the smallest merit in the list. 
      */
      Real getSmallestMerit() const {
         if (m_tests.empty()) {
          throw std::runtime_error("MeritList is empty");
         }
        return m_best ? m_tests.front().getMerit()
                      : m_tests.back().getMerit();
      }

    bool full() const {return m_tests.size() >= m_max;}

    void addWithoutMerit(FigureOfMeritData<Lat> test) {
      if (full())
        return;
      m_tests.push_back(test);
    }


      private:

      /**
       * The tests stored in this object.
       * */
      std::list<FigureOfMeritData<Lat>> m_tests;

      /**
       * The maximal number of tests this object will store.
       * */
      unsigned long m_max;

      /**
       * The minimum merit stored in this object.
       * */
      Real m_merit;

      /**
       * Indicates wether or not this object keeps best or worst lattices.
       * */
      bool m_best;

      /**
       * Finds the right position in the vector and inserts the new test.
       * */
      void posInsert(FigureOfMeritData<Lat>& test){
        if (m_tests.size() == 0) {
          m_tests.push_back(test);
          return;
        }
        for (auto rit = m_tests.crbegin(); rit != m_tests.crend(); rit++) {
          if (m_best) {
            if (test.getMerit() < rit->getMerit()) {
              test.storeLatticeAsString();
              m_tests.insert(rit.base(), test);
              if (rit == m_tests.crbegin() && m_tests.size() == m_max) m_merit = test.getMerit();
              return;
            }
          } else {
            if (test.getMerit() > rit->getMerit()) {
              test.storeLatticeAsString();
              m_tests.insert(rit.base(), test);
              if (rit == m_tests.crbegin() && m_tests.size() == m_max) m_merit = test.getMerit();
              return;
            }
          }
        }
        test.storeLatticeAsString();
        m_tests.insert(m_tests.cbegin(), test);
      }
    };

} // end namespace LatMRG
#endif
