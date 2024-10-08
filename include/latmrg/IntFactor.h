// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
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

#ifndef LATMRG_INTFACTOR_H
#define LATMRG_INTFACTOR_H

#include <string>
#include <iomanip>
#include <cstdint>
#include <sstream>
#include <NTL/ZZ.h>

#include "latticetester/EnumTypes.h"
#include "latmrg/EnumTypes.h"

namespace LatMRG {

static constexpr uint64_t NB_PRIMES = 6543;
const std::array<uint64_t, NB_PRIMES> PRIMES_ARRAY = {
#include "../data/primes.dat"
      };

/**
 * An object of this class represents a factor in the decomposition of a positive integer.
 * It is usually a prime factor, but not always.
 * The class also contains very basic functions to test whether an integer is prime,
 * probably prime, or composite.
 */
template<typename Int> class IntFactor {

public:

   /**
    * Constructs a factor \f$x\f$ of multiplicity `mult` with given `PrimeType status`.
    */
   IntFactor(const Int &x, int64_t mult = 1, PrimeType status = UNKNOWN) :
         m_factor(x), m_multiplicity(mult), m_status(status) {
   }

   /**
    * Returns the numerical value of this factor.
    */
   Int getFactor() const {
      return m_factor;
   }

   /**
    * Sets the value of this factor to \f$x\f$.
    */
   void setFactor(const Int &x) {
      m_factor = x;
   }

   /**
    * Returns the multiplicity of this factor.
    */
   int64_t getMultiplicity() const {
      return m_multiplicity;
   }

   /**
    * Sets the multiplicity of this factor to \f$m\f$.
    */
   void setMultiplicity(int64_t m) {
      m_multiplicity = m;
   }

   /**
    * Returns the `PrimeType` of this factor.
    */
   PrimeType getStatus() const {
      return m_status;
   }

   /**
    * Sets the `PrimeType` of this factor to \f$s\f$.
    */
   void setStatus(PrimeType s) {
      m_status = s;
   }

   /**
    * Tests whether \f$y\f$ is prime or not. First tests whether \f$y\f$ is
    * divisible by all small primes \f$p\f$ (\f$p < 2^{16}\f$) that are
    * kept in file `prime.dat`. If no factor is found, the Miller-Rabin probability
    * test from NTL is applied with \f$k\f$ trials.
    */
   static PrimeType isPrime(const Int &y, std::int64_t k);

   /**
    * Same as `isPrime (y, k)` with the current object in place of `y`.
    */
   PrimeType isPrime(std::int64_t k);

   /**
    * Applies the Miller-Rabin probability test to \f$y\f$ with \f$k\f$ trials.
    */
   static PrimeType isProbPrime(const Int &y, std::int64_t k);

   /**
    * Transforms the status `status` to a string and returns it.
    */
   static std::string toString(PrimeType status);

   /**
    * Returns a string that represents this object.
    */
   std::string toString() const;

   //===========================================================================

private:

   /**
    * This factor as an `Int`.
    */
   Int m_factor;

   /**
    * The multiplicity of this factor.
    */
   int64_t m_multiplicity;

   /**
    * The status of this factor.
    */
   PrimeType m_status;

};
// class IntFactor

//===========================================================================

template<typename Int>
std::string IntFactor<Int>::toString() const {
   char c;
   switch (m_status) {
   case PRIME:
      c = 'P';
      break;
   case PROB_PRIME:
      c = 'Q';
      break;
   case COMPOSITE:
      c = 'C';
      break;
   default:
      c = 'U';
      break;
   }
   std::ostringstream sortie;
   sortie << m_factor << std::setw(10) << m_multiplicity << std::setw(10) << c;
   return sortie.str();
}

//===========================================================================

template<typename Int>
inline std::string IntFactor<Int>::toString(PrimeType status) {
   return toStringPrimeType(status);
}

//===========================================================================

template<typename Int>
PrimeType IntFactor<Int>::isPrime(const Int &y, std::int64_t k) {
   // NbPrem has to be instantiated if we use NTL types
   Int NbPrem(2);
   // NbPrem = 2;
   NTL::ZZ LIM = NTL::conv<NTL::ZZ>(4295098369);  // A bit more than 2^{32}.
   // LIM = NTL::conv<NTL::ZZ> ("4295098369");
   Int ys = NTL::SqrRoot(y);
   uint64_t i = 1;
   while (i < NB_PRIMES && (NbPrem <= ys)) {
      if (y % NbPrem == 0) return COMPOSITE;
      NbPrem = PRIMES_ARRAY[i];
      i++;
   }
   if (y <= LIM) return PRIME;
   // This y is divisible by none of the prime numbers smaller than 2^16.
   return isProbPrime(y, k);
}

//===========================================================================

template<typename Int>
PrimeType IntFactor<Int>::isPrime(std::int64_t k) {
   return isPrime(m_factor, k);
}

//===========================================================================

template<typename Int>
PrimeType IntFactor<Int>::isProbPrime(const Int &y, std::int64_t k) {
   PrimeType status;
   std::int64_t res = NTL::ProbPrime(y, k);
   switch (res) {
   case 0:
      status = COMPOSITE;
      break;
   case 1:
      status = PROB_PRIME;
      break;
   default:
      status = UNKNOWN;
   }
   return status;
}

template class IntFactor<NTL::ZZ> ;
template class IntFactor<std::int64_t> ;

} // namespace LatMRG

#endif
