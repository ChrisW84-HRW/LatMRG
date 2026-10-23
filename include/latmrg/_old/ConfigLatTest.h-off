#ifndef LATMRG_CONFIGLATTEST_H
#define LATMRG_CONFIGLATTEST_H

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
 * This `struct` stores all the information required for a lattice test operation.
 * The information is organized in the same way as in a XML file for a lattest.
 * The variables names correspond approximately to the tag names.
 */
template<typename Int, typename Real> struct ConfigLatTest {

    GenType genType;   // Type of generator, currently MRG or MWC.
    long numComp;      // Number of components.  If > 1, we have a combined generator.

    // List of configurations for the `numComp` components. They can be MRG or MWC.
    // The size of this list must be equal to numComp.
    // Each entry is a configuration for one component.
    vector<ConfigLatTestComponent<Int>*> genComponents;

    bool leapBlocks;  // True iff we use lacunary indices with blocks at large leaps.
    Int leapSize;     // The leap size d when leapBlocks = true.
    Int blockSize;    // The block size s when leapBlocks = true.

    LatticeType latticeType; // The type of lattice that is constructed (see EnumTypes).
    ConfigMerit configFOM;  // The configuration for the FOM and its computation.

    OutputType outType;       // Type of output.
    long verboseLevel = 1;    // Level of verbosity in the output.  Integer from 0 to 3.
    string resFile;           // If present, the results will be printed in file `resFile.res`.
    string genFile;           // If present, the results will be printed .gen format in file `genFile.gen`.
    bool showTimes = true;    // If true, prints CPU times for various operations.
};


/**
 * Abstract structure to store the information required for a component in a lattice test.
 */
template<typename Int, typename Real> abstract struct ConfigLatTestComponent;


/**
 * This `struct` stores the information required for lattest for a MRG component.
 */
template<typename Int, typename Real> struct ConfigLatTestMRG : ConfigLatTestComponent {
    Int modulus;   // The modulus m = b^e + r
    Int base;      // The base b
    Int exponent;  // The exponent e
    Int rest;      // The rest r
    long order;    // The order k
    Intvect coefficients;  // The coefficients a_j
};

/**
 * This stores the information required for lattest for a MWC component.
 */
template<typename Int, typename Real> struct ConfigLatTestMWC : ConfigLatTestComponent {
    // Int modulus;   // The modulus m of the equivalent LCG (to be computed)
    long powMod;   // Value of e such that the modulus is b = 2^e
    long order;    // The order k
    Intvect coefficients;  // The coefficients a_j
};


/**
 * This stores the information required for a seek from a .gen file.
 */
template<typename Int, typename Real> struct ConfigLatTestGenFile : ConfigLatTestComponent {
    string genFile;
};



