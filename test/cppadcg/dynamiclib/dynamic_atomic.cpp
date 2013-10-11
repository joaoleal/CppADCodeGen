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

        virtual std::vector<ADCGD> model(const std::vector<ADCGD>& u) {
            std::vector<ADCGD> Z(m);
#ifdef MODEL1
            Z[0] = cos(u[0]);
            Z[1] = u[1] * u[2] + sin(u[0]);
            Z[2] = u[2] * u[2] + sin(u[1]);
            Z[3] = u[0] / u[2] + u[1] * u[2] + 5.0;
#else
            Z[0] = 1.0 / u[0];
#endif
            return Z;
        }

    };

#ifdef MODEL1
    const size_t CppADCGDynamicAtomicModel1Test::n = 3;
    const size_t CppADCGDynamicAtomicModel1Test::m = 4;
#else
    const size_t CppADCGDynamicAtomicModel1Test::n = 1;
    const size_t CppADCGDynamicAtomicModel1Test::m = 1;
#endif

}

using namespace CppAD;
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
}
