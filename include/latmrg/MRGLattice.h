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
#include "latticetester/IntLatticeExt.h"
#include "latticetester/BasisConstruction.h"
#include "latmrg/FlexModInt.h"

namespace LatMRG {

/**
 * This subclass of `IntLatticeExt` defines an MRG lattice and is similar to Rank1Lattice,
 * but constructs lattices associated with multiple recursive generators (MRGs) with 
 * modulus m, order k, and vector of multipliers a = (a_1, . . . , a_k).
 *
 * The functions `buildBasis` and `incDimBasis` always build and update the primal basis.
 * To build the m-dual, use `buildDualBasis`.
 */

template<typename Int, typename Real>
class MRGLattice: public IntLatticeExt<Int, Real> {

private:
   typedef NTL::vector<Int> IntVec;
   typedef NTL::matrix<Int> IntMat;
   typedef NTL::vector<Real> RealVec;

public:

   /**
    * This constructor takes as input the modulus `m`, the vector of multipliers `aa`,
    * and the norm used to measure the vector lengths.
    * The vector `aa` must have length `k+1`, with \f$a_j\f$ in `aa[j]`.
    * The maximal dimension `maxDim` will be the maximal dimension of the basis.
    * This constructor does not build the basis, so we can build it for a smaller
    * number of dimensions or only for selected projections.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /*
    * Copy constructor.
    */
   MRGLattice(const MRGLattice<Int, Real> &Lat);

   /**
    * Assigns `Lat` to this object.
    */
   MRGLattice& operator=(const MRGLattice<Int, Real> &Lat);

   /**
    * Destructor.
    */
   ~MRGLattice();

   /**
    * Sets the vector of multipliers. The order `k` of the lattice is set equal to
    * the length of this vector, minus 1.  `aa[j]` must contain \f$a_j\f$ for j=1,...,k.
    */
   void setaa(const IntVec &aa);
   
   /**
    * Builds the vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used in the matrix V^{(p)}.
    */
   void buildyPol(int64_t dim);
   
   /**
    * Return the current dimesnion of the vector y.
    */
   int64_t getDimy() { return m_y.length();};
   
   /**
    * Builds a basis in `dim` dimensions. This `dim` must not exceed `this->maxDim()`.
    * This initial primal basis will be upper triangular - CW: based on V^{(0)}.
    */
   void buildBasis(int64_t dim);
   
        /**
    * Builds a basis in `dim` dimensions. This `dim` must not exceed `this->maxDim()`.
    * This initial primal basis will be upper triangular - CW: based on V^{(p)}.
    */
   void buildBasisPol(int64_t dim);

   /**
    * Builds both the primal and an m-dual lower triangular basis directly
    * in `dim` dimensions.  This `dim` must not exceed `maxDim`.
    */
   void buildDualBasis(int64_t dim);

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `maxDim`.
    */
   void incDimBasis();

   /**
    * Increases the current dimension of both the primal and m-dual basis by 1.
    * The new increased dimension must not exceed `maxDim`.
    * This function uses the simplified method for MRG lattices given in the lattice tester guide.
    * It requires the original primal basis to update the m-dual.
    */
   void incDimDualBasis();

   /**
    * This method overrides its namesake in `IntLattice`. A basis for the projection of this
    * `MRGLattice` over the coordinates in `proj` is returned in `projLattice`.
    * The implementation used here exploits the rank-k lattice structure and it
    * is faster than the general one. See Section 5.5 of the LatMRG guide.
    */
   void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Overrides the same function from `IntLattice`.
    * When all the first k coordinates {1,...,k} belong to the projection, both the primal
    * and m-dual constructions are direct, just by selecting the rows and columns
    * whose indices are in `proj`. Otherwise an upper triangular basis is constructed
    * for the basis and a lower-triangular basis for the m-dual.
    */
   void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Returns the first `dim` components of the generating vector \f$\ba\f$ as a string,
    * where `dim` is the current lattice dimension.
    */
   std::string toStringCoef() const;

   // Order of this MRG.
   int m_order;

   /**
    * The coefficients \f$a_1, ..., a_k\f$ of the MRG recurrence, a_j stored in `m_aCoeff[j]`.
    */
   IntVec m_aCoeff;
   
   /**
    * Boolean variable which decides whether the polynomial basis V^{(p)} is used instead of V^{(0)}
    */   
   bool use_pol_basis = true;

protected:

   /**
    * The following protected functions take the basis as a parameter for more flexibility.
    * They are used inside buildBasis, buildBasisDual, incDimBasis, etc., with either `m_basis` or `m_y`.
    */
   void buildBasis0(IntMat &basis, int64_t d);
   
   void buildBasis0Pol(IntMat &basis, int64_t d);

   void incDimBasis0(IntMat &basis, int64_t d);

   bool buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis, const Coordinates &proj);
   
   /**
    * Takes the polynomial `pcol` and returns in `col` the corresponding column in the
    * matrix of generating vectors.
    */
   void polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol);

   /**
    * A vector to store \f$y_0, y_1, ..., y_{t+k-2}\f$ used in the matrix V^{(p)}.
    */
   IntVec m_y;

   /**
    * This auxiliary matrix is used to store the generating vectors of a projections
    * before reducing them into a triangular basis.
    */
   IntMat m_genTemp;
   
   typename FlexModInt<Int>::PolX m_Pz;   // Maybe not needed.  
   typename NTL::zz_pX m_Pz_test;


};

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_maxDim = maxDim;
   setaa(aa);
   this->m_dim = 0;
   m_genTemp.resize(maxDim, maxDim);  
   
   // Build the vector y
   if (use_pol_basis) {
      // Polynomial case
      FlexModInt<Int>::mod_init(m);
      typename FlexModInt<Int>::PolE polDegOne;
      typename FlexModInt<Int>::PolE polPower;  
   
      // Set the characteristic polynomial of the recurrence   
      for (int64_t j = 1; j < this->m_aCoeff.length(); j++) {
        SetCoeff(m_Pz, this->m_aCoeff.length() - j - 1, FlexModInt<Int>::to_Int_p(-this->m_aCoeff[j]));
      }   
      SetCoeff(m_Pz, this->m_aCoeff.length() - 1, 1);    
      FlexModInt<Int>::PolE::init(m_Pz);
   
      buildyPol(maxDim + m_order - 1);
   }
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::~MRGLattice() {
   m_aCoeff.kill();
}

//============================================================================
template<typename Int, typename Real>
MRGLattice<Int, Real>& MRGLattice<Int, Real>::operator=(const MRGLattice<Int, Real> &lat) {
   if (this == &lat) return *this;
   this->copy(lat);
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_y = lat.m_y;
   return *this;
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(), lat.getNormType()) {
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_y = lat.m_y;
   // Should also copy the basis and all other variables???  ******
}

//============================================================================

/**
 * Sets the generating vector to `aa`.
 */
template<typename Int, typename Real>
void MRGLattice<Int, Real>::setaa(const IntVec &aa) {
   m_aCoeff = aa;
   m_order = aa.length() - 1;
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
}

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildyPol(int64_t dim) {
   int64_t k = m_order;
   int64_t j, i;
   int64_t n;
   typename FlexModInt<Int>::PolE polDegOne;
   typename FlexModInt<Int>::PolE polPower;  
   IntVec col;
   
   // Set the characteristic polynomial of the recurrence   
   for (j = 1; j < this->m_aCoeff.length(); j++) {
     SetCoeff(m_Pz, this->m_aCoeff.length() - j - 1, FlexModInt<Int>::to_Int_p(-this->m_aCoeff[j]));
   }   
   SetCoeff(m_Pz, this->m_aCoeff.length() - 1, 1);    
   FlexModInt<Int>::PolE::init(m_Pz);
   
   
   // Auxilliary variables to calculate powers p^n(z)
   std::string str = "[0 1]";
   std::istringstream in (str);
   in >> polDegOne;   
  
   m_y.resize(dim);
   // Temporary
   for (j = 0; j < dim; j++)
     m_y[j] = 0;
     
   for (j = 0; j < min(dim, k-1); j++)
      m_y[j] = 0;
   m_y[k-1] = 1;
   n = ceil(dim/k);
   for (j = 1; j < n+1; j++) {
      // Calculate powers p^\mu-1 
      // And fill the first k rows
      power(polPower, polDegOne, j*k);
      polyToColumn(col, polPower);
      //std::cout << j << ": " << col << "\n";
      for (i = 0; i < k; i++) {
         if (j*k+i < dim) {
            m_y[j*k+i] = col[k-i-1];
         }
      }
   }
   
}
//============================================================================
/*
// Builds an upper-triangular basis directly in `d` dimensions, as explained in Section 4.1 of
// the guide of LatMRG, puts it in `basis`.  Must have d <= m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis00(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t dk = min(d, m_order);
   int64_t i, j, k;
   // Put the identity matrix in the upper left corner
   for (i = 0; i < dk; i++) {
      for (j = 0; j < dk; j++)
         basis[i][j] = (i == j);  // Avoid "if" statements.
   }
   if (d > m_order) {
      // Put m times the identity matrix to the lower right part and the zero matrix to the lower left part.
      for (i = m_order; i < d; i++)
         for (j = 0; j < d; j++)
            basis[i][j] = this->m_modulo * (i == j);
      // Fill the rest of the first m_order rows
      for (i = 0; i < m_order; i++) {
         for (j = m_order; j < d; j++) {
            basis[i][j] = 0;
            // Calculate the components of v_{i,j}. The first component of the coefficient is m_aCoeff[0] here
            for (k = 1; k <= m_order; k++)
               basis[i][j] += m_aCoeff[k] * basis[i][j - k] % this->m_modulo;
         }
      }
   }
}
*/

//============================================================================

// Builds an upper-triangular basis directly in `d` dimensions, as explained in Section 4.1 of
// the guide of LatMRG, using the matrix V^{(p)} that contains y_k,...,y_{k+t-2}.
// Puts this matrix in `basis`.  Must have d <= m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t k = m_order;
   int64_t dk = min(d, k);
   int64_t i, j, jj;
   for (j = 0; j < k-1; j++)
      m_y[j] = 0;
   m_y[k-1] = 1;
   // Put the lower-triangular part of the identity matrix in the upper left corner
   for (i = 0; i < dk; i++) {
      for (j = 0; j <= i; j++)
         basis[i][j] = (i == j);  // Avoid "if" statements.
   }
   // Put m times the identity matrix to the lower right part and the zero matrix to the lower left part.
   for (i = k; i < d; i++)  // If d <= m_order, this does nothing.
      for (j = 0; j < d; j++)
         basis[i][j] = this->m_modulo * (i == j);
   // Fill the rest of the first m_order rows
   for (i = 0; i < dk; i++) {
     for (j = dk; j < d; j++) {
        basis[i][j] = 0;
        for (jj = 1; jj <= m_order; jj++)
           basis[i][j] += m_aCoeff[jj] * basis[i][j - jj] % this->m_modulo;  
     }     
   }
   
}

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0Pol(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t k = this->m_order;
   int64_t dk = min(d, k);
   int64_t i, j;
   
   for (j = 0; j < d; j++) {
      for (i = 0; i < dk; i++)
         basis[i][j] = m_y[j-i+k-1];
   }
   // Fill the rest of the rows
   for (i = dk; i < d; i++) { 
      for (j = 0; j < d; j++)
         basis[i][j] = (i == j) * this->m_modulo;
   }
   


}



//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension isset to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   this->setDim(d);
   this->buildBasis0(this->m_basis, d);
   this->setNegativeNorm();
}

template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasisPol(int64_t d) {
   this->setDim(d);
   this->buildBasis0Pol(this->m_basis, d);
   this->setNegativeNorm();
}

//============================================================================

// Builds the m-dual basis in a direct way in d dimensions.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis(int64_t d) {
   if (d > this->m_dim) {  // Compute the initial primal basis up to d dimensions.
      this->m_dim = d;
   }
   this->setDimDual(d);
   int64_t dk = min(d, m_order);
   int64_t i, j;
   for (i = 0; i < dk; i++)
      for (j = 0; j < d; j++)
         this->m_dualbasis[i][j] = (i == j) * this->m_modulo;
   for (i = dk; i < d; i++) {
      for (j = 0; j < dk; j++) {
         this->m_dualbasis[i][j] = -m_y[i-j+m_order-1];
         }
      for (j = dk; j < d; j++)
         this->m_dualbasis[i][j] = (i == j);
   }
   this->setDualNegativeNorm();
}

//============================================================================

/*
// Increases the dimension of given basis from d-1 to d dimensions.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis00(IntMat &basis, int64_t d) {
   // int64_t d = 1 + this->getDim();  // New current dimension.
   assert(d <= this->m_maxDim);
   int64_t i, j, jj;
   // Add new row and new column of the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   if (d > m_order) basis[d - 1][d - 1] = this->m_modulo;
   else basis[d - 1][d - 1] = 1;
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= m_order) {
         for (jj = 1; jj <= m_order; jj++)
            basis[i][d - 1] += m_aCoeff[k] * basis[i][d - 1 - jj] % this->m_modulo;
      }
   }
}
*/

//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   // int64_t d = 1 + this->getDim();  // New current dimension.
   assert(d <= this->m_maxDim);
   int64_t k = m_order;
   int64_t i, j, jj;
   // Add new row to the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   if (d > k) basis[d - 1][d - 1] = this->m_modulo;
   else basis[d - 1][d - 1] = 1;
   /**
   j = d + k - 2;
   // Add new column.
   basis[0][d - 1] = m_y[j];
   for (i = 1; i < min(d, k); i++)
      basis[i][d - 1] = basis[i - 1][d - 2];
   for (i = k; i < d - 1; i++)
      basis[i][d - 1] = 0;
   */
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= m_order) {
         for (jj = 1; jj <= m_order; jj++)
            basis[i][d - 1] += m_aCoeff[k] * basis[i][d - 1 - jj] % this->m_modulo;
      }
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   this->incDimBasis0(this->m_basis, d);
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis() {
   int64_t d = 1 + this->getDimDual();  // New current dimension.
   this->setDimDual(d);
   while (this->m_dim < d) {  // Increase dimension if needed.
      this->m_dim++;
   }
   int64_t i;
   // Add one extra 0 coordinate to each vector of the m-dual basis.
   for (i = 0; i < d - 1; i++) {
      this->m_dualbasis[i][d - 1] = 0;
      this->m_dualbasis[d - 1][i] = 0;
   }
   // Add the new vector the dual basis
   if (d <= m_order) this->m_dualbasis[d - 1][d - 1] = this->m_modulo;
   else {
      this->m_dualbasis[d - 1][d - 1] = 1;
      for (i = 0; i < m_order; i++)
         this->m_dualbasis[d - 1][i] = -m_y[(d-1)-i+m_order-1];
   }
   //  this->setDualNegativeNorm();   // just for testing ...  not needed.
}

//============================================================================

// We use the columns of basis to construct generating vectors and a pbasis for the projection.
// This function returns the value of `projCase`, which is `true` iff the first m_order coordinates
// are all in the projection.
template<typename Int, typename Real>
bool MRGLattice<Int, Real>::buildProjection0(IntMat &basis, int64_t dimbasis, IntMat &pbasis,
      const Coordinates &proj) {
   int64_t d = proj.size();
   int64_t i, j;
   // projCase1 = false;
   // Check if we are in case 1.
   // This assumes that the coordinates of each projection are always in increasing order!  ***
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.
   if (d < (unsigned) m_order) projCase1 = false;
   else {
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         if (j < m_order) {
            if (*it != unsigned(j + 1)) projCase1 = false;
         } else break;
      }
   }
   if (projCase1) {
      // We first compute the first m_order rows of the projection basis.
      for (i = 0; i < m_order; i++) {
         j = 0;
         for (auto it = proj.begin(); it != proj.end(); it++, j++)
            pbasis[i][j] = basis[i][*it - 1];
      }
      // Then the other rows.
      for (i = m_order; i < d; i++)
         for (j = 0; j < d; j++)
            pbasis[i][j] = this->m_modulo * (i == j);
   } else {
      // In this case we need to use the more general algorithm.
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         // Set column j of all generating vectors, for (j+1)-th coordinate of proj.
         for (i = 0; i < dimbasis; i++)
            m_genTemp[i][j] = basis[i][*it - 1];
      }
      // std::cout << " Generating vectors: \n" << m_genTemp << "\n";
      upperTriangularBasis(m_genTemp, pbasis, this->m_modulo, dimbasis, d);
   }
   return projCase1;
}

//============================================================================

// The number of generating vectors will be this->m_dim.
// The dimension of the projection will be equal to the cardinality of `proj`
// (see the last statement).
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // int64_t d = proj.size();
   assert(*proj.end() <= uint64_t(this->m_dim));
   // IntMat &pbasis = projLattice.getBasis(); // We want to construct this basis of the projection.
   projLattice.setDim(proj.size());
   this->buildProjection0(this->m_basis, this->m_dim, projLattice.getBasis(), proj);
}

//============================================================================
// The number of generating vectors here will be m_dim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t d = proj.size();     // The dimension of this projection.
   // If dimension is not large enough, we increase it.
   while (*proj.end() > uint64_t(this->m_dim)) {
      this->m_dim++;
   }
   IntMat &pbasis = projLattice.getBasis();  // Reference to basis of projection.
   IntMat &pdualBasis = projLattice.getDualBasis();   // And its m-dual.
   projLattice.setDim(d);
   projLattice.setDimDual(d);

   // We first build a basis for the primal basis of the projection.
   bool projCase1 = this->buildProjection0(this->m_basis, this->m_dim, pbasis, proj);
   // After this function call, the class variable `projCase1` is set.

   // Then we compute its m-dual.
   int64_t i, j, k;
   if (projCase1) { // Get the m-dual basis directly. ** Is it really faster and worthwhile?   ****
      // Get the first k columns directly from m_y.
      for (j = 0; j < m_order; j++) {
         pdualBasis[j][j] = this->m_modulo;
         i = m_order;
         auto it = proj.begin();
         for (k = 0; k < m_order - 1; k++)
            it++;
         for (it++; it != proj.end(); ++it, ++i) {
            pdualBasis[i][j] = --m_y[(*it-1)-i+m_order-1];
         }
      }
      // The other columns are trivial.
      for (i = 0; i < d; i++)
         for (j = m_order; j < d; j++)
            pdualBasis[i][j] = (i == j);
   } else {
      // We invert pbasis to get the m-dual.
      mDualUpperTriangular(pbasis, pdualBasis, this->m_modulo, d);
   }
}

// This applies phi inverse as described in the guide, and reverses the coordinates.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::polyToColumn(IntVec &col, typename FlexModInt<Int>::PolE &pcol) {
   int i, j, k;
   k = this->m_order;
   IntVec c;
   c.SetLength(k);
   col.SetLength(k);
   for (j = 1; j < k+1; j++) {
      c[j-1] = 0;
      NTL::conv(c[j-1], coeff(rep(pcol), k-j));
      for (i = 1; i < j; i ++) {
         c[j-1] += this->m_aCoeff[i]*c[j - 1 - i];
      }
      c[j-1] = c[j-1] % this->m_modulo;
   }
   for (j = 0; j < k; j++) {
      col[j] = c[k-j-1];
   }
}


//============================================================================
template<typename Int, typename Real>
std::string MRGLattice<Int, Real>::toStringCoef() const {
   return toString(m_aCoeff, 1, this->getDim() + 1);
}

//============================================================================

template class MRGLattice<std::int64_t, double> ;
template class MRGLattice<NTL::ZZ, double> ;
template class MRGLattice<NTL::ZZ, xdouble> ;
template class MRGLattice<NTL::ZZ, quad_float> ;
template class MRGLattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatMRG

#endif
