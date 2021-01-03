/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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

/**
 * Nonlinear inner model with a single equation and an outer model which
 * just returns the result of the inner model.
 */
class SingleVarAtomicGenericModelTest : public CppADCGDynamicAtomicTest {
protected:
    static const size_t n = 1;
    static const size_t m = 1;
public:

    inline explicit SingleVarAtomicGenericModelTest(bool verbose = false,
                                                    bool printValues = false) :
            CppADCGDynamicAtomicTest("dynamicAtomic1", verbose, printValues) {
        this->verbose_ = false;
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& x,
                             const std::vector<ADCGD>& par) override {
        std::vector<ADCGD> y(m);

        y[0] = 1.0 / x[0];

        return y;
    }

    std::vector<ADCGD> modelOuter(const std::vector<ADCGD>& y,
                                  const std::vector<ADCGD>& par) override {
        std::vector<ADCGD> z(1);

        z[0] = y[0];

        return z;
    }

};

/**
 * Nonlinear inner model and an outer model which performs linear
 * operations with the result of the inner model.
 */
class CppADCGDynamicAtomicModel2Test : public CppADCGDynamicAtomicTest {
protected:
    static const size_t n = 3;
    static const size_t m = 4;
public:

    inline explicit CppADCGDynamicAtomicModel2Test(bool verbose = false,
                                                   bool printValues = false) :
        CppADCGDynamicAtomicTest("dynamicAtomic2", verbose, printValues) {
        this->verbose_ = false;
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

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;
using namespace std;

/**
 * @test Nonlinear inner model and an outer model which
 *       just returns the result of the inner model
 */
TEST_F(SingleVarAtomicGenericModelTest, DynamicForRevNonlinearLinear) {
    using namespace std;
    using CppAD::vector;

    std::vector<Base> x(n);
    x[0] = 0.5;
    std::vector<Base> par;

    // simple wrap (no outer model used)
    this->testADFunAtomicLibSimple(x, par); // one compiled model used as an atomic function by an ADFun

    this->testAtomicSparsities(x, par);

    // use outer model
    this->testADFunAtomicLib(x, par); // 1 compiled inner model used by CppAD

    this->testAtomicLibAtomicLib(x, par); // 2 models in 2 dynamic libraries
}

/**
 * @test Nonlinear inner model and an outer model which performs linear
 *       operations with the result of the inner model
 */
TEST_F(CppADCGDynamicAtomicModel2Test, DynamicForRevNonlinearLinear) {
    using namespace std;
    using CppAD::vector;

    std::vector<Base> x(n);
    for (size_t j = 0; j < n; j++)
        x[j] = j + 2;

    std::vector<Base> par;

    // simple wrap (no outer model used)
    this->testADFunAtomicLibSimple(x, par); // one compiled model used as an atomic function by an ADFun

    this->testAtomicSparsities(x, par);

    // use outer model
    this->testADFunAtomicLib(x, par); // 1 compiled inner model used by CppAD

    this->testAtomicLibAtomicLib(x, par); // 2 models in 2 dynamic libraries
}
