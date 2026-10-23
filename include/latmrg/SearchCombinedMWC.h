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

#ifndef LATMRG_CONFIGSEARCHCOMBINEDMWC_H
#define LATMRG_CONFIGSEARCHCOMBINEDMWC_H

#include <vector>
#include <utility>
#include <cstdint>

#include "latmrg/SearchBase.h"
#include "latticetester/Util.h"

namespace LatMRG {

/**
 * Configuration for a combined MWC search.
 * This is only a DUMMY CLASS.
 */
template<typename Int, typename Real>
struct SearchCombinedMWC : public SearchBase<Int, Real, MWCLattice<Int, Real>>
{
    SearchCombinedMWC() : SearchBase<Int, Real, MWCLattice<Int, Real>>(MWC) {}
    /**
     * Configuration for one MWC component.
     */
    struct ComponentMWC
    {
        // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
        long powMod;   // Value of e such that the modulus is b = 2^e
        long order;    // The order k
        IntVec lowBoundaries;  // The low and high boundaries for the coefficients a_j
        IntVec highBoundaries; // Default values are 1 and b-1
        vector<long> numPowerTwo;  // Values of n_i (max number of powers of 2)
        vector<long> maxPowerTwo;  // Values of p_i (max power of 2)
        vector<int64_t> randomBits; // For random searches, the nonzero coefficients `a_j` are generated randomly with `e_j` random bits
        int64_t numaj; // number of positive coefficients 
        vector<vector<pair<long,long>>> equalCoeffs;  // List of pairs of coefficients that must be equal

        long getOrder() const
        {
            return order;
        }

        Int getNoMultipliers() const
        {
            Int total(1);
            for (long i = 1; i <= order; i++) {
                total *= highBoundaries[i] - lowBoundaries[i] + 1;
            }
            return total;
        }

        Int getLowBoundary(int i) const
        {
            return lowBoundaries[i];
        }

        Int getHighBoundary(int i) const
        {
            return highBoundaries[i];
        }

        int64_t getRandomBits(int i) const
        {
            return randomBits[i];
        }

        long getPowMod() const
        {
            return powMod;
        }
    };

    std::vector<ComponentMWC> genComponents;

    void createComponents()
    {
        genComponents.resize(this->numComp);
    }

    const ComponentMWC & component(int i) const
    {
        return genComponents[i];
    }

    ComponentMWC & component(int i)
    {
        return genComponents[i];
    }

    Int getModulus() const override
    {
        Int mod(1);
        for (const auto& comp : genComponents) {
            mod *= comp.getModulus();
        }
        return mod;
    }

    
    MWCLattice<Int, Real>* nextGenerator() override;

    MWCLattice<Int, Real>* nextGeneratorRandom() override;
        
    bool checkMaxPeriod(MWCLattice<Int, Real>& lat) override;

    private: 

    /**
    * Helper class used to find admissible MWC moduli.
    *
    * Searches around 2^e - 1 +/- c and returns a modulus m such that:
    *
    *   - m is prime
    *   - (m-1)/2 is prime
    *   - 2^((m-1)/2) != 1 (mod m)
    *
    * These conditions provide the safe-prime properties required
    * for MWC maximal-period behaviour.
    */
    class Modulus {
        public:

        /**
        * Search will begin at `2^e±c` and will increase of decrease depending on
        * `increase`.
        */
        Modulus(long e, Int c, bool increase) {
            m = (Int(1)<<e) - 1;
            this->increase = increase;
            if (increase) m+=c;
            else m-=c;
        }

        /**
        * This function finds the next value for `m`. This can return 1 if
        * `increase` is `false` and `m` gets to 1.
        */
        Int next() {
            while (m > 1) {
                nextM();
                LatMRG::PrimeType status = LatMRG::IntFactor<Int>::isPrime (m, KTRIALS);
                if (status == LatMRG::PRIME || status == LatMRG::PROB_PRIME) {
                  if (1 == m % 4) continue;
                  Int m1 = (m - 1)/2;
                  status = LatMRG::IntFactor<Int>::isPrime (m1, KTRIALS);
                  if (status != LatMRG::PRIME && status != LatMRG::PROB_PRIME) continue;
                  // For MWC we check 2^{(m-1)/2} \neq 1 mod m.
                  NTL::ZZ_p::init(NTL::ZZ(m));
                  NTL::ZZ_p a = NTL::power(NTL::ZZ_p(2), m1);
                  if (NTL::IsOne(a)) {
                    continue;
                  }
                  return m;
                }
              }
            return Int(-1);
        }

        private:

        static const long KTRIALS = 200;  //  The user must be allowed to select this constant!

        /**
        * Last valid modulus found (or 2^e).
        */
        Int m;

        /**
        * As passed to constructor.
        */
        bool increase;

        /**
        * Increment/decrement m
        */
        void nextM () {
            if (increase) m += 2;
            else m -= 2;
            if (0 == m % 5){
                if (increase) m += 2;
                else m -= 2;
            }
        }
    }; 

};

//============================================================================
// Implementation

  /**
  * nextGenerator() performs an exhaustive enumeration of the admissible MWC 
  * generators defined by the configuration. The multiplier coefficients are 
  * obtained by decoding the current generator index into the corresponding 
  * coefficient values within the prescribed ranges. The resulting coefficients 
  * are then used to construct the associated MWC lattice.
  */
  template<typename Int, typename Real> MWCLattice<Int, Real>* SearchCombinedMWC<Int, Real>::nextGenerator()   {
    auto* comp = this->conf.genComponents[0];
    
    Int range, val;
    Int tmp = NTL::to_ZZ(this->currentGen);
    
    int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);

    
    // decode currentGen into multi-index
    for (int i = 1; i <= k; i++) {
      range = comp->getHighBoundary(i) - comp->getLowBoundary(i) + 1;
      val = tmp % range;
      tmp /= range;
      aa[i] = comp->getLowBoundary(i) + val;
    }

    Int b = NTL::power(Int(2), comp->getPowMod());

    if (this->currentGen < this->conf.max_gen && this->currentGen < comp->getNoMultipliers())
    {
      ++this->currentGen;
      return new MWCLattice<Int, Real>(b, aa, this->conf.maxdim, this->conf.maxdim, this->conf.maxdim);
    }

    return nullptr;
  }

  //===========================================================================

  /**
  * nextGeneratorRandom() generates admissible MWC generators by randomly 
  * selecting the multiplier coefficients according to the bit-size restrictions
  * and configuration parameters. It provides the random search analagoue of 
  * nextGenerator() and can be used by performSearch() when a stochastic search.
  */
  template<typename Int, typename Real> MWCLattice<Int, Real>* SearchCombinedMWC<Int, Real>::nextGeneratorRandom() {
    auto* comp = this->conf.genComponents[0];
    int k = comp->getOrder();
    NTL::Vec<NTL::ZZ> aa;
    aa.SetLength(k + 1);
    // Initialized
    aa[0] = -1;
    for (int64_t j = 1; j < k; j++)
      aa[j] = 0;
    Int b = NTL::power(Int(2), comp->getPowMod());

    // Get random MWC generator
    if (comp->getRandomBits(0) > 0) {
         aa[0] = conv<Int>(LatticeTester::RandBits(comp->getRandomBits(0)));
         if ((aa[0] % 2) == 0) aa[0] += 1;  // a_0 must be odd.
      }
    aa[k] = conv<Int>(LatticeTester::RandBits(comp->getRandomBits(k)));  // `ek` random bits for a_k.
    for (int64_t j = 1; j < min(asMWC(comp)->numaj, k); j++)
        aa[k-j] = conv<Int>(LatticeTester::RandBits(comp->getRandomBits(j)));  // `ej` random bits for a_{k-j}.
    for (int64_t j = 1; j < -asMWC(comp)->numaj; j++)  // To impose equal coefficients, for testing.
        aa[k-j] = aa[k];
    
    if (this->currentGen < this->conf.max_gen)
    {
      ++this->currentGen;
      return new MWCLattice<Int, Real>(b, aa, this->conf.maxdim, this->conf.maxdim, this->conf.maxdim);
    }

    return nullptr;
  }

  //===========================================================================
  
  /**
  * This method checks whether the generated lattice satisfies the
  * maximal-period requirement for an MWC lattice.
  */
  template<typename Int, typename Real> bool SearchCombinedMWC<Int, Real>::checkMaxPeriod(MWCLattice<Int, Real>& lat) {
    return mIsSafePrime(lat.getModulus(), 100) && maxPeriodHalfMWC(lat.getModulus(), this->conf.genComponents[0]->getPowMod());
 }

} // namespace LatMRG

#endif // LATMRG_CONFIGSEARCHCOMBINEDMWC_H