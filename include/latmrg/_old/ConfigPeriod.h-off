#ifndef LATMRG_CONFIGPERIOD_H
#define LATMRG_CONFIGPERIOD_H

#include <ctime>
#include <cstdlib>
#include <list>

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/Random.h"
#include "latticetester/Reducer.h"

#include "latmrg/MWCLattice.h"
#include "latmrg/MMRGLattice.h"
#include "latmrg/ComboLattice.h"
#include "latmrg/Chrono.h"
#include "latmrg/ParamReaderExt.h"
#include "latmrg/Projections.h"
#include "latmrg/FigureOfMerit.h"
#include "latmrg/EnumTypes.h"
#include "latmrg/Projections.h"
#include "latmrg/Test.h"

using LatticeTester::IntLatticeExt;
using LatticeTester::Normalizer;
using namespace std;
using namespace LatMRG;
using namespace LatticeTester::Random;

typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

/**
 * This `struct` stores all the information required for the program that verifies
 * if an MRG has maximal period or not.
 * The information is organized in the same way as in a XML file for the `period` mode.
 */
template<typename Int, typename Real> struct ConfigPeriod {

    GenType genType;          // Type of generator, currently MRG or MWC.
    ConfigPeriodComponent<Int>* genComponent; // Configuration for the MRG or MWC component.

    OutputType outType;       // Type of output.
    long verboseLevel = 1;    // Level of verbosity in the output.  Integer from 0 to 3.
    string resFile;           // If present, the results will be printed in file `resFile.res`.
    bool showTimes = true;    // If true, prints CPU times for various operations.
};

/**
 * This `struct` stores the information required for a max period test for a MRG component.
 */
template<typename Int, typename Real> struct ConfigPeriodMRG : ConfigPeriodComponent {
    Int modulus;   // The modulus m = b^e + r
    Int base;      // The base b
    Int exponent;  // The exponent e
    Int rest;      // The rest r
    long order;    // The order k
    Intvect coefficients;  // The coefficients a_j
    DecompType howFactorh; // Indicates how to factorize h=(m-1)/2
    IntVect factorsh;      // The factors of h, repeated according to multiplicity, when known
    DecompType howFactor;  // Indicates how to factorize r=(m^k-1)/(m-1)
    IntVect factorsr;      // The factors of r, repeated according to multiplicity, when known
    bool permaxPow2;   // True iff m is a power of 2, b=2, r=0, k=1, and we want max period for the LCG
                       // Default is false
};

/**
 * This stores the information required for a max period test for a MWC component.
 */
template<typename Int, typename Real> struct ConfigPeriodMWC : ConfigPeriodComponent {
    // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
    long powMod;   // Value of e such that the modulus is b = 2^e
    long order;    // The order k
    Intvect coefficients;  // The coefficients a_j
};


/**
 * A generic type of struct for storing the information required for a max period test.
 */
template<typename Int, typename Real> abstract struct ConfigPeriodComponent;





