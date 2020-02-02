/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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
#include "../CppADCGDynamicTest.hpp"
#include "cstr.hpp"

namespace CppAD {
namespace cg {

class CstrDynamicTest : public CppADCGDynamicTest {
public:

    inline explicit CstrDynamicTest(bool verbose = false,
                                    bool printValues = false) :
        CppADCGDynamicTest("cstr", verbose, printValues) {
        // independent variable vector
        _xTape = std::vector<double>(28, 1.0);
        _xNorm.resize(28);
        _xNorm[0] = 0.3; // h
        _xNorm[1] = 7.82e3; // Ca
        _xNorm[2] = 304.65; // Tr
        _xNorm[3] = 301.15; // Tj

        _xNorm[4] = 2.3333e-04; // u1
        _xNorm[5] = 6.6667e-05; // u2

        _xNorm[6] = 6.2e14; //
        _xNorm[7] = 10080; //
        _xNorm[8] = 2e3; //
        _xNorm[9] = 10e3; //
        _xNorm[10] = 1e-11; //
        _xNorm[11] = 6.6667e-05; //
        _xNorm[12] = 294.15; //
        _xNorm[13] = 294.15; //
        _xNorm[14] = 1000; //
        _xNorm[15] = 4184; //Cp
        _xNorm[16] = -33488; //deltaH
        _xNorm[17] = 299.15; // Tj0
        _xNorm[18] = 302.65; //   Tj2
        _xNorm[19] = 7e5; // cwallj
        _xNorm[20] = 1203; // csteam
        _xNorm[21] = 3.22; //dsteam
        _xNorm[22] = 950.0; //Ug
        _xNorm[23] = 0.48649427192323; //vc6in
        _xNorm[24] = 1000; //rhoj
        _xNorm[25] = 4184; //Cpj
        _xNorm[26] = 0.014; //Vj
        _xNorm[27] = 1e-7; //cwallr

        _xRun = _xTape;

        _eqNorm.resize(4);
        _eqNorm[0] = 0.3;
        _eqNorm[1] = 7.82e3;
        _eqNorm[2] = 304.65;
        _eqNorm[3] = 301.15;
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ind,
                             const std::vector<ADCGD>& par) override {
        return CstrFunc<CG<double> >(ind);
    }

};

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;

TEST_F(CstrDynamicTest, ForwardZero) {
    this->testForwardZero();
}

TEST_F(CstrDynamicTest, DenseJacobian) {
    this->testDenseJacobian();
}

TEST_F(CstrDynamicTest, DenseHessian) {
    this->testDenseHessian();
}

TEST_F(CstrDynamicTest, Jacobian) {
    this->testJacobian();
}

TEST_F(CstrDynamicTest, Hessian) {
    this->testHessian();
}

namespace CppAD {
namespace cg {

class CstrDynamicWithParamsTest : public CppADCGDynamicTest {
public:

    inline CstrDynamicWithParamsTest(bool verbose = false,
                                     bool printValues = false) :
            CppADCGDynamicTest("cstr_params", verbose, printValues) {
        // independent variable vector
        _xTape = std::vector<double>(6, 1.0);
        _xNorm.resize(6);
        _xNorm[0] = 0.3; // h
        _xNorm[1] = 7.82e3; // Ca
        _xNorm[2] = 304.65; // Tr
        _xNorm[3] = 301.15; // Tj

        _xNorm[4] = 2.3333e-04; // u1
        _xNorm[5] = 6.6667e-05; // u2

        _xRun = _xTape;

        _parTape.resize(22);
        _parTape[0] = 6.2e14; //
        _parTape[1] = 10080; //
        _parTape[2] = 2e3; //
        _parTape[3] = 10e3; //
        _parTape[4] = 1e-11; //
        _parTape[5] = 6.6667e-05; //
        _parTape[6] = 294.15; //
        _parTape[7] = 294.15; //
        _parTape[8] = 1000; //
        _parTape[9] = 4184; //Cp
        _parTape[10] = -33488; //deltaH
        _parTape[11] = 299.15; // Tj0
        _parTape[12] = 302.65; //   Tj2
        _parTape[13] = 7e5; // cwallj
        _parTape[14] = 1203; // csteam
        _parTape[15] = 3.22; //dsteam
        _parTape[16] = 950.0; //Ug
        _parTape[17] = 0.48649427192323; //vc6in
        _parTape[18] = 1000; //rhoj
        _parTape[19] = 4184; //Cpj
        _parTape[20] = 0.014; //Vj
        _parTape[21] = 1e-7; //cwallr

        _parRun = _parTape;

        _eqNorm.resize(4);
        _eqNorm[0] = 0.3;
        _eqNorm[1] = 7.82e3;
        _eqNorm[2] = 304.65;
        _eqNorm[3] = 301.15;
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ind,
                             const std::vector<ADCGD>& par) override {
        assert(!par.empty());
        return CstrFunc<CG<double> >(ind, par);
    }

};

} // END cg namespace
} // END CppAD namespace


TEST_F(CstrDynamicWithParamsTest, ForwardZero) {
    this->testForwardZero();
}

TEST_F(CstrDynamicWithParamsTest, DenseJacobian) {
    this->testDenseJacobian();
}

TEST_F(CstrDynamicWithParamsTest, DenseHessian) {
    this->testDenseHessian();
}

TEST_F(CstrDynamicWithParamsTest, Jacobian) {
    this->testJacobian();
}

TEST_F(CstrDynamicWithParamsTest, Hessian) {
    this->testHessian();
}
