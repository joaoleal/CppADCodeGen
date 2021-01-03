/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2020 Joao Leal
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

class CppADCGDynamicAtomicCustomSparsityTest : public CppADCGDynamicAtomicTest {
protected:
    static const size_t n = 3;
    static const size_t m = 4;
    CppAD::vector<Base> x;
    CppAD::vector<Base> par;
    const CppAD::vector<Base> xNorm;
    const CppAD::vector<Base> eqNorm;
    std::vector<std::set<size_t> > jacInner, hessInner;
    std::vector<std::set<size_t> > jacOuter, hessOuter;
public:

    explicit CppADCGDynamicAtomicCustomSparsityTest(bool verbose = false,
                                                    bool printValues = false) :
        CppADCGDynamicAtomicTest("model_atomic", verbose, printValues),
        x(n),
        jacInner(m), hessInner(n),
        jacOuter(m - 1), hessOuter(n) {
        this->verbose_ = false;

        using namespace std;
        using CppAD::vector;

        for (size_t j = 0; j < n; j++)
            x[j] = j + 2;

        /*
         * Elements for the custom Jacobian/Hessian tests
         */
        jacInner[0] = {0};
        jacInner[1] = {0, 1, 2};
        jacInner[2] = {1, 2};
        jacInner[3] = {0, 1, 2};

        // lower left side (with 1 exception)
        hessInner[0] = {0, 2};
        hessInner[1] = {1};
        hessInner[2] = {1, 2}; // flipped

        jacOuter[0] = jacInner[0];
        jacOuter[1] = jacInner[1];
        jacOuter[2] = jacInner[2];
        jacOuter[2].insert(jacInner[3].begin(), jacInner[3].end());

        hessOuter = hessInner; // only lower left side
        hessOuter[2].erase(1);
        hessOuter[1].insert(2);
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& x,
                             const std::vector<ADCGD>& par) override {
        std::vector<ADCGD> y(m);

        y[0] = cos(x[0]);
        y[1] = x[1] * x[2] + sin(x[0]);
        y[2] = x[2] * x[2] + sin(x[1]);
        y[3] = x[0] / x[2] + x[1] * x[2] + 5.0;

        return y;
    }

};

TEST_F(CppADCGDynamicAtomicCustomSparsityTest, AtomicLibModelBridgeCustomRev2) {
    this->testAtomicLibModelBridgeCustom(x, par,
                                         xNorm,
                                         eqNorm,
                                         jacInner, hessInner,
                                         jacOuter, hessOuter,
                                         true,
                                         1e-14, 1e-13);
}

TEST_F(CppADCGDynamicAtomicCustomSparsityTest, AtomicLibModelBridgeCustomDirect) {
    this->testAtomicLibModelBridgeCustom(x, par,
                                         xNorm,
                                         eqNorm,
                                         jacInner, hessInner,
                                         jacOuter, hessOuter,
                                         false,
                                         1e-14, 1e-13);
}

} // END cg namespace
} // END CppAD namespace
