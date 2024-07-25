#ifndef LATMRG_PRIMITIVEINT_H
#define LATMRG_PRIMITIVEINT_H

#include "latmrg/IntFactorization.h"

namespace LatMRG {

/**
 * This file deals with primitive roots and primitive elements modulo an
 * integer. Suppose that \f$a\f$, \f$e\f$, and \f$p\f$ are integers, \f$p\f$ is a
 * prime number, and \f$a\f$ is relatively prime with \f$m=p^e\f$.
 * The smallest positive integer \f$\lambda(m,a)\f$ for which
 * \f$a^{\lambda(m,a)}= 1 \ \bmod\ m \f$ is called the order of \f$a\f$ modulo
 * \f$m\f$. Any \f$a\f$ which has the maximum possible order for a given
 * \f$m\f$ is called a primitive element modulo \f$m\f$.
 * This file provides a functions to test for this property.
 */

/**
 * Returns `true` iff `a` is a primitive element modulo \f$p^e\f$.
 * The prime factor decomposition of \f$p-1\f$ must be given in `f`
 * and the list of inverse factors in there must be up to date.
 */
template<typename Int>
static bool isPrimitiveElement(const Int &a, const IntFactorization<Int> &f,
      const Int &p, long e = 1);


//============================================================================
// Implementation

template<typename Int>
bool isPrimitiveElement(const Int &a, const IntFactorization<Int> &f,
      const Int &p, long e) {
   if (0 == p) throw std::range_error("PrimitiveInt::isPrimitiveElement: p = 0");
   if (0 == a) return false;
   Int t1, t2;
   Int m;
   m = NTL::power(p, e);
   t1 = a;
   if (t1 < 0) t1 += m;
   const std::vector<Int> invList = f.getInvFactorList();
   assert (!(invList.empty ()));
   for (auto it = invList.begin(); it != invList.end(); it++) {
      if (*it == (m - 1)) continue;
      t2 = NTL::PowerMod(t1, *it, m);
      if (t2 == 1) return false;
   }
   return true;
}

}  // namespace LatMRG
#endif