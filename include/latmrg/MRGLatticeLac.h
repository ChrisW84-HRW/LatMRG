#ifndef LATMRG_MRGLATTICELAC_H
#define LATMRG_MRGLATTICELAC_H

#include <string>
#include "latticetester/EnumTypes.h"
#include "latticetester/Lacunary.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/Types.h"
// #include "latticetester/Const.h"
// #include "latticetester/Lacunary.h"
#include "latticetester/MRGLattice.h"
// #include "latmrg/MRGComponent.h"

namespace LatMRG {

/**
 * This subclass of `MRGLattice` constructs and handles lattice bases built from MRGs as in `MRGLattice`,
 * but with arbitrary lacunary indices that are regularly spaced by packets of the same size.
 *
 * Perhaps the new functions offered here could be integrated into MRGLattice,
 * to reduce the number of classes.  ???            ************
 *
 */
class MRGLatticeLac: public LatticeTester::MRGLattice {
public:

   // Parent:
   // MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Constructor with modulus of congruence \f$m\f$, order of the recurrence
    * \f$k\f$, multipliers \f$A\f$, maximal dimension `maxDim`, and lattice type
    * `lattype`. Vector and matrix indices vary from 1 to `maxDim`. The length of
    * the basis vectors is computed with `norm`. The bases are built using the
    * *lacunary indices* `lac`.
    * The basis is built for the lacunary
    * indices `lac`.  `aa` has to be a vector of k+1 components
    * with `a[i]`=\f$a_i\f$ for compatibility with other classes.
    */
   MRGLatticeLac(const Int &m, const IntVec &aa, int64_t maxDim, IntVec &lac,
         NormType norm = L2NORM);

   /**
    * Copy constructor. The maximal dimension of the new basis is set to
    * <tt>Lat</tt>’s current dimension.
    */
   // MRGLatticeLac(const MRGLatticeLac &Lat);

   /**
    * Assigns `Lat` to this object. The maximal dimension of this basis is
    * set to <tt>Lat</tt>’s current dimension.
    */
   // MRGLatticeLac& operator=(const MRGLatticeLac &Lat);

   /**
    * Destructor.
    */
   virtual ~MRGLatticeLac();

   /**
    * Sets the lacunary indices for this lattice to `lat`.
    */
   void setLac(const IntVec &lat);

   /**
    * Returns the \f$j\f$-th lacunary index.
    */
   Int& getLac(int j);

protected:

   void initStates();

   /**
    * Increases the dimension of the basis by 1.
    */
   void incDimBasis(int);

   /**
    * The lacunary indices.
    */
   IntVec m_lac;

   //===========================================================================

   /** Max order for lacunary case in this class; otherwise, it takes too much memory.
    For order > ORDERMAX, use subclass MRGLatticeLac instead.
    This means that we can have short lacunary indices, supported here,
    and also long lacunary indices (e.g., for multiple streams), supported in MRGLatticeLac.  ******* ???
    */
#define ORDERMAX 100

};

//===========================================================================
// Implementation:



//============================================================================

// Builds an upper-triangular basis directly in `d` dimensions with the lacunary indices,
// as explained in Section 4.1 of the guide of LatMRG, and puts it in `basis`.
// Must have d <= m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0(IntMat &basis, int64_t d) {
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



//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   // int64_t d = 1 + this->getDim();  // New current dimension.
   assert(d <= this->m_maxDim);
   int64_t i, j, k;
   // Add new row and new column of the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   if (d - 1 >= m_order) basis[d - 1][d - 1] = this->m_modulo;
   else basis[d - 1][d - 1] = 1;
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= m_order) {
         for (k = 1; k <= m_order; k++)
            basis[i][d - 1] += m_aCoeff[k] * basis[i][d - 1 - k] % this->m_modulo;
      }
   }
}


//============================================================================

// We must be able to do this with basis equal to either m_basis or m_basis0.
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

// The old implementation.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildLaBasis(int64_t d) {

   // NOT USED, see: MRGLatticeLac::buildBasis

   if (this->m_order > ORDERMAX)
      LatticeTester::MyExit(1, "MRGLattice::buildLaBasis:   k > ORDERMAX");

   initStates();
   int64_t IMax = m_lac.getSize();

   IntVec b;
   b.SetLength(this->m_order + 1);
   LatticeTester::Invert(m_aCoeff, b, this->m_order);

   // b is the characteristic polynomial of the MRG.
   PolyPE < Int > ::setModulus(this->m_modulo);
   PolyPE < Int > ::setPoly(b);
   PolyPE<Int> pol;
   int64_t ord = 0;

   // Construction d'un systeme generateur modulo m.
   for (int64_t k = 0; k < IMax; k++) {
      // pour chaque indice lacunaire
      NTL::conv(m_e, m_lac[k]);

      // x^m_e Mod f(x) Mod m
      pol.powerMod(m_e);
      pol.toVector(m_xi);

      ord = 0;
      for (int64_t i = 1; i <= this->m_order; i++) {
         if (m_ip[i]) {
            ++ord;
            m_t5 = 0;
            for (int64_t j = 1; j <= this->m_order; j++)
               m_t5 += m_sta[i][j] * m_xi[j - 1];
            this->m_wSI[ord][k] = m_t5;
         }
      }
   }

   //  From here we can use BasisConstruction.  *********

   /* On veut s'assurer que la base m_v soit triangulaire (pour satisfaire
    * les conditions de l'article \cite{rLEC94e} [sec. 3, conditions sur
    * V_i >= i]) et de plein rang (on remplace les lignes = 0 par lignes
    * avec m sur la diagonale).
    * */
   LatticeTester::Triangularization<IntMat>(this->m_wSI, this->m_vSI, ord, IMax, this->m_modulo);
   LatticeTester::CalcDual < IntMat > (this->m_vSI, this->m_wSI, IMax, this->m_modulo);

   // Construire la base de dimension 1
   this->m_basis[0][0] = this->m_vSI[0][0];
   this->m_dualbasis[0][0] = this->m_wSI[0][0];
   this->setDim(1);

   this->setNegativeNorm();
   this->setDualNegativeNorm();

   //  This approach could be slow!  *****
   for (int64_t i = 2; i <= d; i++)
      incDimLaBasis(IMax);

   // for debugging
   // trace("ESPION_1", 1);
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimLaBasis(int64_t IMax) {

   LatticeTester::IntLatticeExt<Int, Real>::incDim();
   const int64_t dim = this->getDim(); // new dimension (dim++)

   /*
    if (dim >= IMax) {
    MyExit (0,
    "Dimension of the basis is too big:\nDim > Number of lacunary indices.");
    }
    */

   IntVec tempLineBasis(dim);
   IntVec tempColBasis(dim);

   for (int64_t i = 0; i < dim - 1; i++) {

      // tempLineBasis <- m_basis[i]
      for (int64_t k = 0; k < dim - 1; k++)
         tempLineBasis[k] = this->m_basis[i][k];

      for (int64_t i1 = 0; i1 < dim - 1; i1++) {

         Int tempScalDual;
         LatticeTester::ProdScal<Int>(tempLineBasis, this->m_wSI[i1], dim, tempScalDual);
         LatticeTester::Quotient(tempScalDual, this->m_modulo, tempScalDual);
         this->m_t1 = tempScalDual * this->m_vSI[i1][dim - 1];
         tempColBasis[i] += this->m_t1;
      }
      LatticeTester::Modulo(tempColBasis[i], this->m_modulo, tempColBasis[i]);
      this->m_basis[i][dim - 1] = tempColBasis[i];
   }

   for (int64_t j = 0; j < dim - 1; j++)
      this->m_basis[dim - 1][j] = 0;
   this->m_basis[dim - 1][dim - 1] = this->m_vSI[dim - 1][dim - 1];

   for (int64_t i = 0; i < dim - 1; i++)
      this->m_dualbasis[i][dim - 1] = 0;

   for (int64_t j = 0; j < dim - 1; j++) {

      Int tempScalDualBis;

      for (int64_t i = 0; i < dim - 1; i++) {
         this->m_t1 = this->m_dualbasis[i][j];
         this->m_t1 *= tempColBasis[i];
         tempScalDualBis += this->m_t1;
      }
      if (tempScalDualBis != 0) tempScalDualBis = -tempScalDualBis;

      LatticeTester::Quotient(tempScalDualBis, this->m_vSI[dim - 1][dim - 1], tempScalDualBis);
      this->m_dualbasis[dim - 1][j] = tempScalDualBis;
   }

   LatticeTester::Quotient(this->m_modulo, this->m_vSI[dim - 1][dim - 1], this->m_t1);
   this->m_dualbasis[dim - 1][dim - 1] = this->m_t1;

   this->setNegativeNorm();
   this->setDualNegativeNorm();

}

}
#endif

