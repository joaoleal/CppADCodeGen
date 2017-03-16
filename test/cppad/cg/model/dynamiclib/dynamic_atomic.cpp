/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
#include "CppADCGDynamicAtomicTest.hpp"

namespace CppAD {
namespace cg {

#define MODEL1

class CppADCGDynamicAtomicModel1Test : public CppADCGDynamicAtomicTest {
protected:
    static const size_t n;
    static const size_t m;
public:

    inline CppADCGDynamicAtomicModel1Test(bool verbose = false, bool printValues = false) :
        CppADCGDynamicAtomicTest("dynamicAtomic", verbose, printValues) {
        this->verbose_ = false;
    }

    virtual std::vector<ADCGD> model(const std::vector<ADCGD>& x) {
        std::vector<ADCGD> y(m);
#ifdef MODEL1
        y[0] = cos(x[0]);
        y[1] = x[1] * x[2] + sin(x[0]);
        y[2] = x[2] * x[2] + sin(x[1]);
        y[3] = x[0] / x[2] + x[1] * x[2] + 5.0;
#else
        y[0] = 1.0 / x[0];
#endif
        return y;
    }

};

#ifdef MODEL1
const size_t CppADCGDynamicAtomicModel1Test::n = 3;
const size_t CppADCGDynamicAtomicModel1Test::m = 4;
#else
const size_t CppADCGDynamicAtomicModel1Test::n = 1;
const size_t CppADCGDynamicAtomicModel1Test::m = 1;
#endif

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;
using namespace std;

TEST_F(CppADCGDynamicAtomicModel1Test, DynamicForRev) {
    using namespace std;
    using CppAD::vector;

    vector<Base> x(n);
#ifdef MODEL1
    for (size_t j = 0; j < n; j++)
        x[j] = j + 2;
#else
    x[0] = 0.5;
#endif

    this->testADFunAtomicLib(x);
    
    this->testAtomicSparsities(x);
}
