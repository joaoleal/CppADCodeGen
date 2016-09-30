/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2016 Ciengis
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

class CppADCGThreadPoolTest : public CppADCGDynamicTest {
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;
protected:
    std::vector<ADCG> u;
    std::vector<double> x;
public:

    inline CppADCGThreadPoolTest(bool verbose = true) :
        CppADCGDynamicTest("pool", verbose, false),
        u(9),
        x(u.size()) {
        this->_multithread = MultiThreadingType::PTHREADS;

        // independent variables
        for (auto& ui : u)
            ui = 1;

        for (auto& xi : x)
            xi = 1.5;
    }

    virtual std::vector<ADCGD> model(const std::vector<ADCGD>& x) {
        std::vector<ADCGD> y(6);

        for (size_t i = 0; i < 3; ++i) {
            size_t i0 = i * 2;
            size_t j0 = i * 3;

            y[i0] = cos(x[j0]);
            y[i0 + 1] = x[j0 + 1] * x[j0 + 2] + sin(x[j0]);
        }

        return y;
    }

};

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;
using namespace std;

TEST_F(CppADCGThreadPoolTest, DisabledFullVars) {
    this->_multithreadDisabled = true;

    this->_reverseOne = true;
    this->_reverseTwo = true;
    this->_denseJacobian = false;
    this->_denseHessian = false;

    this->testDynamicFull(u, x, 1000);
}

TEST_F(CppADCGThreadPoolTest, SingleJobFullVars) {
    this->_multithreadDisabled = false;
    this->_multithreadScheduler = ThreadPoolScheduleStrategy::SINGLE_JOB;

    this->_reverseOne = true;
    this->_reverseTwo = true;
    this->_denseJacobian = false;
    this->_denseHessian = false;

    this->testDynamicFull(u, x, 1000);
}

TEST_F(CppADCGThreadPoolTest, MultiJobFullVars) {
    this->_multithreadDisabled = false;
    this->_multithreadScheduler = ThreadPoolScheduleStrategy::MULTI_JOB;

    this->_reverseOne = true;
    this->_reverseTwo = true;
    this->_denseJacobian = false;
    this->_denseHessian = false;

    this->testDynamicFull(u, x, 1000);
}

TEST_F(CppADCGThreadPoolTest, StaticFullVars) {
    this->_multithreadDisabled = false;
    this->_multithreadScheduler = ThreadPoolScheduleStrategy::STATIC;

    this->_reverseOne = true;
    this->_reverseTwo = true;
    this->_denseJacobian = false;
    this->_denseHessian = false;

    this->testDynamicFull(u, x, 1000);
}

TEST_F(CppADCGThreadPoolTest, DynamicCustomElements) {

    std::vector<size_t> jacRow(3), jacCol(3); // all elements except 1
    jacRow[0] = 0;
    jacCol[0] = 0;
    jacRow[1] = 1;
    jacCol[1] = 0;
    jacRow[2] = 1;
    jacCol[2] = 2;

    std::vector<size_t> hessRow(2), hessCol(2); // all elements except 1
    hessRow[0] = 0;
    hessCol[0] = 0;
    hessRow[1] = 2;
    hessCol[1] = 1;

    this->_reverseOne = true;
    this->_reverseTwo = true;
    this->_denseJacobian = false;
    this->_denseHessian = false;

    this->testDynamicCustomElements(u, x, jacRow, jacCol, hessRow, hessCol);
}