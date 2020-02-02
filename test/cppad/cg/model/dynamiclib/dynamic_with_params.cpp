/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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
#include "CppADCGDynamicTest.hpp"

namespace CppAD {
namespace cg {

class CppADCGDynamicTest1 : public CppADCGDynamicTest {
public:

    inline CppADCGDynamicTest1(bool verbose = false,
                               bool printValues = false) :
        CppADCGDynamicTest("dynamic_with_params", verbose, printValues) {

        // independent variables
        _xTape = {1, 1, 1};
        _parTape = {2, 3};
        _xRun = {1, 2, 1};
        _parRun = {3, 4};

        this->_forwardOne = false;
        this->_reverseOne = false;
        this->_reverseTwo = false;
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ax,
                             const std::vector<ADCGD>& ap) override {
        std::vector<ADCGD> ay(2);

        ay[0] = cos(ax[0]) * ap[0];
        ay[1] = ax[1] * ax[2] + sin(ax[0]) + ap[1];

        return ay;
    }

};

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;
using namespace std;

TEST_F(CppADCGDynamicTest1, ForwardZero1Assign) {
    _maxAssignPerFunc = 1;
    this->testForwardZero();
}

TEST_F(CppADCGDynamicTest1, ForwardZero) {
    _maxAssignPerFunc = 1000;
    this->testForwardZero();
}

TEST_F(CppADCGDynamicTest1, DenseJacobian) {
    this->testDenseJacobian();
}

TEST_F(CppADCGDynamicTest1, DenseHessian) {
    this->testDenseHessian();
}

TEST_F(CppADCGDynamicTest1, Jacobian) {
    this->testJacobian();
}

TEST_F(CppADCGDynamicTest1, Hessian) {
    this->testHessian();
}
