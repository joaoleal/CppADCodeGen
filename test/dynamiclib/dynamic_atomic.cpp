/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
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

TEST_F(CppADCGDynamicAtomicModel1Test, multMatrixMatrixSparsityTrans) {
    std::vector<std::set<size_t> > aT(4); // a: 3 x 4
    aT[0].insert(0);
    aT[1].insert(1);
    aT[2].insert(1);
    aT[3].insert(0);
    aT[3].insert(2);
    std::vector<std::set<size_t> > b(4); // b: 4 x 2
    b[0].insert(0);
    b[2].insert(0);
    b[2].insert(1);
    b[3].insert(1);
    CppAD::vector<std::set<size_t> > rT(2); // r: 3 x 2 

    multMatrixMatrixSparsityTrans(aT, b, rT, 4, 2, 3);

    CppAD::vector<std::set<size_t> > rTExpected(2);
    rTExpected[0].insert(0);
    rTExpected[0].insert(1);
    rTExpected[1].insert(0);
    rTExpected[1].insert(1);
    rTExpected[1].insert(2);

    compareVectorSetValues(rT, rTExpected);
}

TEST_F(CppADCGDynamicAtomicModel1Test, multMatrixMatrixSparsityTrans2) {
    size_t m = 2;
    size_t n = 7;
    size_t q = 33;
    std::vector<std::set<size_t> > sT(m);
    sT[0].insert(29);
    sT[1].insert(30);

    std::vector<std::set<size_t> > jac(m);
    jac[0].insert(0);
    jac[0].insert(2);
    jac[1].insert(1);
    jac[1].insert(3);

    CppAD::vector<std::set<size_t> > rT(n);

    multMatrixMatrixSparsityTrans(sT, jac, rT, m, n, q);

    CppAD::vector<std::set<size_t> > rTExpected(n);
    rTExpected[0].insert(29);
    rTExpected[1].insert(30);
    rTExpected[2].insert(29);
    rTExpected[3].insert(30);

    compareVectorSetValues(rT, rTExpected);
}