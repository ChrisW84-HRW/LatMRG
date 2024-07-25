#ifndef LATMRG_PRIMITIVEPOLY_H
#define LATMRG_PRIMITIVEPOLY_H

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>

#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "NTL/ZZ_pXFactoring.h"
#include "NTL/ZZ_pE.h"
#include "NTL/lzz_p.h"
#include "NTL/lzz_pE.h"
#include "NTL/lzz_pX.h"
#include "NTL/lzz_pXFactoring.h"

#include "latticetester/Util.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/IntFactor.h"
#include "latmrg/PrimitiveInt.h"

namespace LatMRG {

template<typename Int>
class ModInt {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};

template<>
class ModInt<std::int64_t> {
public:
   typedef NTL::zz_p IntP;
   typedef NTL::zz_pE PolE;
   typedef NTL::vec_zz_p IntVecP;
   typedef NTL::mat_zz_p IntMatP;
   typedef NTL::zz_pX PolX;
};

template<>
class ModInt<NTL::ZZ> {
public:
   typedef NTL::ZZ_p IntP;
   typedef NTL::ZZ_pE PolE;
   typedef NTL::vec_ZZ_p IntVecP;
   typedef NTL::mat_ZZ_p IntMatP;
   typedef NTL::ZZ_pX PolX;
};

/**
 * This class deals with polynomials \f$P(x)\f$ in \f$\mathbb Z_m[X]\f$
 * defined as
 * \anchor REF__PrimitivePoly_eq_poly1
 * \f[
 *   P(x) = c_0 + c_1x^1 + c_2 x^2 + \cdots + c_n x^n \tag{eq.poly1}
 * \f]
 * with degree \f$n\f$ and integer coefficients \f$c_i\f$ in
 * \f$\mathbb Z_m\f$. The arithmetic operations on these polynomials are
 * done modulo \f$m\f$ and modulo a polynomial \f$f(x)\f$ of degree \f$k\f$.
 * Thus all polynomials will be reduced modulo \f$f(x)\f$. In LatMRG, the
 * modulus polynomial \f$f(x)\f$ is usually written in the form
 * \anchor REF__PrimitivePoly_eq_poly2
 * \f[
 *   f(x) = x^k - a_1x^{k-1} - \cdots- a_{k-1} x - a_k, \tag{eq.poly2}
 * \f]
 * and is associated with the recurrence
 * \anchor REF__PrimitivePoly_eq_rec2
 * \f[
 *   x_n = (a_1 x_{n-1} + a_2 x_{n-2} + \cdots+ a_k x_{n-k}) \bmod m. \tag{eq.rec2}
 * \f]
 * The two functions `setM` and `setF` *must* be called to initialize the
 * modulus \f$m\f$ and the modulus polynomial \f$f(x)\f$ before doing any
 * arithmetic operations on `PrimitivePoly` objects.
 *
 * We recall that the polynomial \f$f(x)\f$ in {@link REF__PrimitivePoly_eq_poly2 (eq.poly2)}
 * is a primitive polynomial modulo \f$m\f$ if and only if  the three
 * following conditions are satisfied:
 * \anchor REF__PrimitivePoly_isprimi
 * <dl> <dt>None</dt>
 * <dd>
 * \f$[(-1)^{k+1} a_k]^{(m-1)/q} \bmod m \neq1\f$ for each prime
 * factor \f$q\f$ of \f$m - 1\f$;
 * </dd>
 * <dt>None</dt>
 * <dd>
 * \f$x^r \bmod(f(x),m) =  (-1)^{k+1} a_k \bmod m\f$;
 * </dd>
 * <dt>None</dt>
 * <dd>
 * \f$x^{r/q} \bmod(f(x), m) \f$ has positive degree for each prime
 * factor \f$q\f$ of \f$r\f$, with \f$1<q< r\f$;
 * </dd>
 * </dl> where \f$r = (m^k - 1)/(m - 1)\f$.
 * Condition 1 is the same as saying that \f$(-1)^{k+1} a_k\f$ is a
 * primitive root of \f$m\f$. Condition 3 is automatically satisfied
 * when \f$r\f$ is prime.
 *
 * Type `Int` is used to represent the polynomial coefficients.
 * The possible associated types `IntVec` are `int64_t*`         // *****  Not NTL::vector<ZZ> ???
 *  and <tt>vec_ZZ</tt>. Type `PolE` for the polynomials may be chosen
 * as <tt>zz_pE</tt> when \f$m\f$ is small enough, or may be set to
 * <tt>ZZ_pE</tt> which is implemented with the big integer type <tt>NTL::ZZ_p</tt>.
 *
 */

/**
 * Returns `true` if the polynomial \f$f(x)\f$ with coefficients in C is a primitive polynomial
 * modulo \f$m\f$. The factorizations of
 * \f$m-1\f$ and \f$r\f$ must be in `fm` and `fr` respectively.
 */
template<typename Int>
static bool isPrimitive(const IntVec &C, const Int &m, const IntFactorization<Int> &fm,
      const IntFactorization<Int> &fr);

// I think that if we choose to make the following public because they willbe used elsewhere,
//  then the functions should have more significant names!                          *************


/**
 * This method returns `true` iff the primitivity conditions 2 and 3 given in Knuth (to be added in the guide)
 * are satisfied by the modulus \f$f(x)\f$. It does not check condition 1, assuming it to be `true`.
 */
template<typename Int>
static bool isPrimitive(const IntVec &C, const Int &m, const IntFactorization<Int> &fr);

/**
 * Returns the modulus polynomial \f$f(x)\f$ for given coefficients 'C' and modulus 'm'.
 * The output polynomail is stored in the variable 'f'.
 */
static void getF(const IntVec &C, const Int &m, ModInt<Int>::PolX f);

/**
 * Returns the polynomial to \f$x^j \mod f(x) (\bmod m)\f$ where \f$f(x)\f$ is given according to the
 * coeffcients 'C' and modulus 'm'. The output polynomail is stored in the variable 'fj'.
 */
static void powerMod(const Int &j, const IntVec &C, const Int &m, ModInt<Int>::PolX fj);

//===========================================================================
// Implementation

template<typename Int>
static bool isPrimitive(const IntVec &C, const Int &m, const IntFactorization<Int> &fm,
      const IntFactorization<Int> &fr) {
   Int a0;
   static int64_t k;
   typename ModInt<Int>::PolX f;
   getF(C, m, f);
   k = C.length() - 1;
   // rep is the NTL::ZZ equivalent of the NTL::ZZ_p element.
   a0 = -rep(ConstTerm(f));
   if ((k & 2) == 0) a0 = -a0;
   if (!isPrimitiveElement(a0, fm, m)) return false;
   return isPrimitive(C, m, fr);
}

template<typename Int>
static bool isPrimitive(const IntVec &C, const Int &m, const IntFactorization<Int> &fr) {
   typename ModInt<Int>::PolX f;
   getF(C, m, f);
   static int64_t k;
   k = C.length() - 1;

   if (1 == k) return true;
   // Is f irreducible ?
   typename ModInt<Int>::PolX Q;
   Q = f;
   //  if (!isIrreducible())      // slow
   //  if (0 == DetIrredTest(Q))   // medium slow
   if (0 == IterIrredTest(Q))   // fastest
   return false;

   // ---- Test Condition 2
   Int r0;
   r0 = fr.getNumber();
   Q = PowerXMod(r0, f);
   if (0 != deg(Q)) return false;
   Int T1;
   T1 = rep(ConstTerm(f));
   if ((k & 1) == 1) T1 = -T1;
   if (T1 < 0) T1 += m;
   if (rep(ConstTerm(Q)) != T1) return false;

   // ---- Test Condition 3
   if (fr.getStatus() == LatticeTester::PRIME) return true;
   std::vector<Int> invFactorList = fr.getInvFactorList();
   assert(!invFactorList.empty());
   typename std::vector<Int>::const_iterator it = invFactorList.begin();
   while (it != invFactorList.end()) {
      Q = PowerXMod(*it, f);
      if (0 == deg(Q)) return false;
      ++it;
   }
   return true;
}

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&   ACHTUNG!  &&&&&&&&&

// Here two new objects f and cP are created inside this function and are
// never destroyed. So more and more objects are created and never recycled.
// This could lead to a memory leak!
// If ther are destroyed automatically on exit, then the f that is returned
// longer exists.  So there is a problem both ways.

// I think at least the f should rather be passed and returned as a parameter of the function.
//                                                                           ************
static void getF(const IntVec &C, const Int &m, ModInt<Int>::PolX f) {
   static int64_t k;
   typename ModInt<Int>::IntVecP cP;
   conv(cP, C);
   k = cP.length() - 1;
   for (int64_t i = 0; i <= k; i++)
      SetCoeff(f, i, cP[i]);
   f.normalize();
   ModInt<Int>::PolE::init(f);
}

// Same:  The result should rather be returned in one of the arguments.   ************
static void powerMod(const Int &j, const IntVec &C, const Int &m, ModInt<Int>::PolX fj) {
   typename ModInt<Int>::PolX f;
   getF(C, m, f);
   power(fj, f, j);
}

}
#endif