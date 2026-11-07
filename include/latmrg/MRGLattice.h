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

#ifndef LATMRG_MRGLATTICE_H
#define LATMRG_MRGLATTICE_H

#include <cassert>
#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/FlexTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/BasisConstruction.h"
#include "latmrg/FlexModInt.h"

namespace LatMRG {

/**
 * This subclass of `IntLatticeExt` handles lattices associated with multiple recursive generators
 * (MRGs) with modulus \f$m\f$, order \f$k\f$, and vector of multipliers \f$\mathbf{a} = (a_1, . . . , a_k)\f$.
 * In this implementation, instead of computing and storing the full basis and/or its m-dual,
 * we store only a vector \f$\mathbf{y} = (y_0,y_1,\dots,y_t)\f$ that contains the successive values
 * produced by the MRG recurrence when the \f$k\f$ initial values are \f$(0,\dots,0,1)\f$.
 * The other entries of this vector are computed using the recurrence.
 * Its coordinates are sufficient to build a basis or its m-dual
 * for the lattice or for a projection, as explained in the guide.
 */

template<typename Int, typename Real>
class MRGLattice: public IntLatticeExt<Int, Real> {

public:

   /**
    * This constructor takes as input the modulus `m`, the order `k` or the vector of multipliers `aa`,
    * the maximal dimension allowed for the full lattice, a projection, and a coordinate index,
    * and the norm used to measure the vector lengths.
    * When the vector `aa` is given, the order `k` is deduced its length via `k = aa.length()-1`,
    * and the MRG coefficients are set to \f$a_j=\f$`aa[j]` for \f$j=1,...,k\f$.
    * This vector `aa` can be set or changed separately (several times) by `setaa`,
    * but the order `k` cannot be changed.  The modulus `m` can be changed via `setModulus` from `IntLattice.
    * The parameters `maxDim`,  `maxDimProj`, and`maxCoord` give the maximal dimension of a basis (primal or m-dual),
    * the maximal dimension of a projection (for `buildProjection`),
    * and the maximal coordinate index that we may consider for the lattice and for any projection.
    * These values are used by the constructor to allocate space for matrix bases and for the vector \f$\mathbf{y}\f$.
    * This vector has length `maxCoord`\f$+k-1\f$ and it is filled by the function `setaa`
    * each time we set of change the vector `aa`.
    * The constructor does not build a basis.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord,
         NormType norm = L2NORM);

   /**
    * This version takes `maxCoord = maxDimProj = maxDim'.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Same as the first constructor, except that the vector `aa` is not given, only the order `k` is given.
    * The vector `aa` must be set afterwards via `setaa`.
    */
   MRGLattice(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord,
         NormType norm = L2NORM);

   MRGLattice(const Int &m, int64_t k, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Here the modulus `m` is not specified. It must be set later via `setModulus`.
    */
   MRGLattice(int64_t k, int64_t maxDim, int64_t maxDimProj, int64_t maxCoord,
         NormType norm = L2NORM);

   MRGLattice(int64_t k, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Destructor.
    */
   ~MRGLattice();

   /**
    * Sets the vector of multipliers to `aa`.  This vector must have length \f$k+1\f$
    * where  \f$k\f$ is the order `k` of this MRG.
    * One will have \f$a_j=\f$`aa[j]` for \f$j=1,...,k\f$.
    * The function also recomputes the vector \f$\mathbf{y} = (y_0,\dots,y_{t+k-2})\f$
    * where \f$t=\f$`maxCoord`\f. This vector has length $t+k-1\f$.
    * If `maxCoord` is given, `m_maxCoord` is changed to this new value,
    * otherwise its current value is taken.
    * This function is useful when the same `MRGLattice` object is used for several
    * different vectors `aa` of multipliers.
    */
   virtual void setaa(const IntVec &aa, int64_t maxCoord);

   virtual void setaa(const IntVec &aa);

   /**
    * Returns the vector of multipliers `aa`.
    */
   virtual IntVec getaa();

   /**
    * Compute the entries of \f$(y_0, y_1, ..., y_{t+k-2})\f$, where \f$\t=\f$`maxCoord` is the largest
    * coordinate index that we want to use for the lattice or for a projection.
    * If `maxCoord` is given, `m_maxCoord` is updated to this new value, otherwise the old
    *  `m_maxCoord` is used and kept.
    */
   virtual void computey(int64_t maxCoord);

   virtual void computey();

   /**
    * Builds a primal basis in `dim` dimensions. This `dim` must be at least the order `k` and
    * must not exceed `m_maxDim` or `m_maxCoord`.
    * This initial primal basis will be upper triangular.
    */
   virtual void buildBasis(int64_t dim);

   /**
    * Builds a lower-triangular m-dual basis directly in `dim` dimensions.
    * This `dim` must be at least the order `k`and  must not exceed `m_maxDim` or `m_maxCoord`.
    */
   virtual void buildDualBasis(int64_t dim);

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `m_maxDim` or `m_maxCoord`.
    */
   virtual void incDimBasis();

   /**
    * Increases the current dimension of the m-dual basis by 1.
    * The new increased dimension must not exceed `m_maxDim` or `m_maxCoord`.
    * This function uses the simplified method for MRG lattices given in the lattice tester guide.
    * It increases by 1 the number of columns of `m_bV0` that are computed.
    * This matrix is initialized upon the first call of `incDimDualBasis`.
    */
   virtual void incDimDualBasis();

   /**
    * Builds a basis for the projection of this `MRGLattice` onto the coordinates
    * in `coordSet` and puts it as the `m_basis` of `projLattice`.
    * The construction uses the vector \f$\mathbf{y}\f$, as explained in the guide.
    * The largest coordinate of the projection must not exceed `m_maxCoord` and
    * the number of coordinates must not exceed `m_maxDimProj`.
    * If the first `k` coordinates are in the projection, the construction is direct,
    * otherwise we may have to build a basis from a set of up to \f$k+s\f$ generating vectors,
    * where \f$s\f$ is the dimension of the projection.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection and puts it
    * as the `m_dualbasis` in `projLattice`.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet)
         override;

   /**
    * Returns the vector 'm_y'.
    */
   // IntVec gety();
protected:

   /**
    * This function initializes the first `dim` columns of `m_bV0`.
    * One must have `m_order <= dim <= maxDim`.
    * It is called by `buildDualBasis`.
    */
   void buildV0(int64_t dim);

   // Order of this MRG.
   int m_order;

   /**
    * The coefficients \f$a_1, ..., a_k\f$ of the MRG recurrence, a_j stored in `m_aa[j]`.
    */
   IntVec m_aa;

   /**
    * Largest index of a coordinate that can be considered.
    */
   int64_t m_maxCoord;

   /**
    * A vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used in the matrix \f$V^{(p)}\f$,
    * with \f$t=\f$`m_maxCoord`.
    */
   IntVec m_y;

   /**
    * The current dimension of m_bV0
    */
   int64_t m_dimbV0 = 0;

   /**
    * Used to construct an m-dual basis and when increasing its dimension.
    */
   IntMat m_bV0;

   /**
    * Maximal number of coordinates allowed in a projection.
    */
   int64_t m_maxDimProj;

   /**
    * This auxiliary matrix is used to store the set of generating vectors of a projection
    * before reducing them into a triangular basis.
    * Its size must be at least `(d + k) x d`, where `d = maxDimProj` is the largest dimension of a projection.
    * It is used only in `buildProjection`.
    */
   IntMat m_genTemp;

};

//============================================================================
// IMPLEMENTTION

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(int64_t k, int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      IntLatticeExt<Int, Real>(maxDim, norm) {
   this->m_order = k;
   m_aa.SetLength(k + 1);
   m_maxDimProj = maxDimProj;
   m_maxCoord = maxCoord;
   m_y.SetLength(maxCoord + k - 1);
   m_genTemp.SetDims(maxDimProj + k, maxDimProj);
   m_bV0.SetDims(k, maxDim);
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(int64_t k, int64_t maxDim, NormType norm) :
      MRGLattice<Int, Real>(k, maxDim, maxDim, maxDim, norm) {
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, int64_t k, int64_t maxDim, int64_t maxDimProj,
      int64_t maxCoord, NormType norm) :
      MRGLattice<Int, Real>(k, maxDim,  maxDimProj, maxCoord, norm) {
   this->setModulus(m);
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, int64_t k, int64_t maxDim, NormType norm) :
      MRGLattice<Int, Real>(m, k, maxDim, maxDim, maxDim, norm) {
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim,
      int64_t maxDimProj, int64_t maxCoord, NormType norm) :
      MRGLattice<Int, Real>(m, aa.length() - 1, maxDim, maxDimProj, maxCoord, norm) {
   setaa(aa);
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm) :
      MRGLattice<Int, Real>(m, aa.length() - 1, maxDim, maxDim, maxDim, norm) {
   setaa(aa);
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::~MRGLattice() {
   m_aa.kill();
   m_y.kill();
   m_genTemp.kill();
   m_bV0.kill();
}

//============================================================================

/*
 template<typename Int, typename Real>
 MRGLattice<Int, Real>& MRGLattice<Int, Real>::operator=(const MRGLattice<Int, Real> &lat) {
 if (this == &lat) return *this;
 // this->copy(lat); // copy constructor currently not available
 m_aa = lat.m_aa;
 m_order = lat.m_order;
 m_maxdim_m_y = lat.m_maxdim_m_y;
 return *this;
 }

 //============================================================================

 template<typename Int, typename Real>
 MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
 IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(), lat.getNormType()) {
 // this->copy(lat); // copy constructor currently not available
 m_aa = lat.m_aa;
 m_order = lat.m_order;
 m_maxdim_m_y = lat.m_maxdim_m_y;
 }
 */

//============================================================================
template<typename Int, typename Real>
void MRGLattice<Int, Real>::setaa(const IntVec &aa, int64_t maxCoord) {
   m_maxCoord = maxCoord;
   m_y.SetLength(maxCoord + m_order - 1);
   setaa(aa);
}

template<typename Int, typename Real>
void MRGLattice<Int, Real>::setaa(const IntVec &aa) {
   assert(aa.length() == m_order + 1);  // `aa` must have the right length.
   m_aa = aa;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   computey(m_maxCoord);
   m_dimbV0 = 0;    // This one is also invalid.
}

//============================================================================
template<typename Int, typename Real>
IntVec MRGLattice<Int, Real>::getaa() {
   return m_aa;
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension is set to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   int64_t k = m_order;
   assert((d >= k) && (d <= this->m_maxDim) && (d <= m_maxCoord));
   this->setDim(d);
   int64_t i, j;  // Row i and column j.
   for (i = 0; i < k; i++) {
      for (j = 0; j < d; j++)
         this->m_basis[i][j] = m_y[j - i + k - 1];
   }
   for (i = k; i < d; i++) {
      for (j = 0; j < d; j++)
         this->m_basis[i][j] = (i == j) * this->m_modulo;
   }
   // std::cout << "buildBasis, basis = \n" << this->m_basis << "\n";
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis(int64_t d) {
   int64_t k = this->m_order;
   assert((d >= k) && (d <= this->m_maxDim) && (d <= m_maxCoord));
   this->setDimDual(d);
   int64_t i, j;
   if (k == 1) {
      // Fill first column.
      this->m_dualbasis[0][0] = this->m_modulo;
      for (i = 1; i < d; i++)
         this->m_dualbasis[i][0] = -m_y[i] % this->m_modulo;
      // Then the other columns.
      for (i = 0; i < d; i++)
         for (j = 1; j < d; j++)
            this->m_dualbasis[i][j] = (i == j);
   } else {
      buildBasis(d);  // This basis is upper-triangular.
      mDualUpperTriangular(this->m_dualbasis, this->m_basis, this->m_modulo, d);
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
   assert((this->m_dim < min(this->m_maxDim, m_maxCoord)));
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   int64_t k = m_order;
   int64_t i, j, jj;
   // We compute the entries of the new column by using the recurrence.
   for (i = 0; i < d - 1; i++) {
      this->m_basis[i][d - 1] = 0;
      if (d > k) {
         for (jj = 1; jj <= k; jj++)
            this->m_basis[i][d - 1] += m_aa[jj] * this->m_basis[i][d - 1 - jj] % this->m_modulo;
      }
   }
   // Add `m e_k` as new row to the primal basis.
   for (j = 0; j < d - 1; j++)
      this->m_basis[d - 1][j] = 0;
   this->m_basis[d - 1][d - 1] = this->m_modulo;
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis() {
   assert((this->m_dimdual < min(this->m_maxDim, m_maxCoord)));
   int64_t d = 1 + this->getDimDual();  // New current dimension. 
   this->setDimDual(d);
   int64_t i;
   // If this is the first time we call `incDimDualBasis` for this basis,
   // we need to compute the initial part of m_bV0, otherwise we just add one column.
   if (m_dimbV0 < d - 1) buildV0(d);
   else if (m_dimbV0 == d - 1) {
      // Compute a new column of m_bV0.
      for (i = 0; i < this->m_order; i++) {
         m_bV0[i][d - 1] = 0;
         for (int jj = 1; jj <= m_order; jj++)
            m_bV0[i][d - 1] += m_aa[jj] * m_bV0[i][d - 1 - jj] % this->m_modulo;
      }
      m_dimbV0++;
   }
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d - 1; i++) {
      this->m_dualbasis[i][d - 1] = 0;
   }
   // Use `m_bV0` to add a new row to the m-dual basis.
   for (i = 0; i < this->m_order; i++) {
      this->m_dualbasis[d - 1][i] = -m_bV0[i][d - 1];
   }
   for (i = this->m_order; i < d - 1; i++)
      this->m_dualbasis[d - 1][i] = 0;
   this->m_dualbasis[d - 1][d - 1] = 1;
}

//============================================================================

// The number of generating vectors will be this->m_dim.
// The dimension of the projection will be equal to the cardinality of `coordSet`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   int64_t d = coordSet.size();
   assert(d <= m_maxDimProj);
   assert(*coordSet.end() <= uint64_t(m_maxCoord));  // Vector y is large enough.
   projLattice.setDim(d);
   IntMat* pbasis = &projLattice.getBasis(); // Pointer to basis of projection.
   // This `pbasis` will refer to either the basis of projLattice or to `m_genTemp`,
   // depending on the value of `projCase1`. This is to avoid copying.
   int64_t i, j;
   int64_t k = this->m_order;
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `coordSet`.
   
   // Check if we are in case 1.
   // This assumes that the coordinates of each projection are always in increasing order!  ***   
   // Note that we do not have direct access to coordinate `k` to check if it equals `k`.
   // This is why we need to loop over the coordinates.
   if (d < (unsigned) m_order) projCase1 = false;
   else {
      j = 1;
      for (auto it = coordSet.begin(); it != coordSet.end(), j <= k; it++, j++) {
         // if (j <= k) {
         if (*it != unsigned(j)) {
            projCase1 = false;
            break;
         }
      }
   }
   int64_t iadd = 0;
   if (!projCase1) {
      // In this case, we make `pbasis` point to `m_genTemp`, in which we will put
      // a set of generating vectors that will be reduced to a triangular basis.
      pbasis = &m_genTemp;  // Pointer will point a set of gen vectors.
      iadd = k;            // We will have k more rows in pbasis.
   }
   // Recall that if (projCase1), pbasis points to the basis of `projLattice`.
   // We first copy the selected coordinates for the first k rows.
   j = 0;
   for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++)
      for (i = 0; i < k; i++)
         (*pbasis)[i][j] = m_y[*it - i + k - 2];
   // Then the other rows.
   for (i = k; i < d + iadd; i++)
      for (j = 0; j < d; j++)
         (*pbasis)[i][j] = this->m_modulo * (i == j+iadd);
   // If not case1, we must reduce the set of gen vectors.
   if (!projCase1)
      upperTriangularBasis(projLattice.getBasis(), m_genTemp, this->m_modulo, this->m_dim, d);
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   int64_t d = coordSet.size();     // The dimension of this projection.
   assert(d <= m_maxDimProj);
   assert(*coordSet.end() <= uint64_t(m_maxCoord));  // Vector y is large enough.
   int64_t i, j;
   projLattice.setDimDual(d);
   IntMat & pbasis = projLattice.getBasis();  // Reference to basis of projection.
   IntMat & pdualBasis = projLattice.getDualBasis();   // And its m-dual.
   // In the special case where i_1 = k = 1, we get the m-dual basis directly.
   // so there is not need to build the primal basis.
   if (this->m_order == 1 && *coordSet.begin() == 1) {
      j = 0;
      for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++) {
         pdualBasis[j][0] = -m_y[*it - 1];
      }
      for (j = 0; j < d; j++)
         pdualBasis[0][j] = (i == j) * this->m_modulo;
      for (i = 1; i < d; i++) {
         for (j = 1; j < d; j++)
            pdualBasis[i][j] = (i == j);
      }
   } else {
      // Otherwise, we build a basis for the primal basis of the projection,
      // and we compute its m-dual.
      this->buildProjection(projLattice, coordSet);
      mDualUpperTriangular(pdualBasis, pbasis, this->m_modulo, d);
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::computey(int64_t maxCoord) {
   m_maxCoord = maxCoord;
   m_y.SetLength(maxCoord + m_order - 1);
   computey();
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::computey() {
   int64_t k = m_order;
   int64_t j, jj;
   for (j = 0; j < k - 1; j++)
      m_y[j] = 0;
   m_y[k - 1] = 1;
   for (j = k; j < m_maxCoord + k - 1; j++) {
    m_y[j] = 0;
    for (jj = 1; jj <= k; jj++)
        m_y[j] += m_aa[jj] * m_y[j - jj] % this->m_modulo;
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildV0(int64_t d) {
   int64_t k = m_order;
   assert((d >= k) && (d <= this->m_maxDim) && (d <= m_maxCoord));
   // int64_t dk = min(d, k);
   int64_t i, j, jj;
   for (i = 0; i < k; i++) {
      for (j = 0; j < k; j++)
         m_bV0[j][i] = (i == j);
      for (j = k; j < d; j++) {
         m_bV0[i][j] = 0;
         for (jj = 1; jj <= k; jj++)
            m_bV0[i][j] += m_aa[jj] * m_bV0[i][j - jj] % this->m_modulo;
      }
   }
   m_dimbV0 = d;
}

//============================================================================

template class MRGLattice<std::int64_t, double> ;
template class MRGLattice<NTL::ZZ, double> ;
template class MRGLattice<NTL::ZZ, xdouble> ;
template class MRGLattice<NTL::ZZ, quad_float> ;
template class MRGLattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatMRG

#endif
