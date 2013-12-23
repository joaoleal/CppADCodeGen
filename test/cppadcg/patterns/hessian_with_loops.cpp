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
#include "CppADCGPatternTest.hpp"

typedef double Base;
typedef CppAD::CG<Base> CGD;
typedef CppAD::AD<CGD> ADCGD;

namespace CppAD {

    class CppADCGHessLoopTest : public CppADCGPatternTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    public:

        inline CppADCGHessLoopTest(bool verbose = false, bool printValues = false) :
            CppADCGPatternTest(verbose, printValues) {
            this->testJacobian_ = false;
        }
    };
}

using namespace CppAD;

/**
 * @test Linear model (no Hessian)
 */
std::vector<ADCGD> modelLinear(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = 3 * x[i * n];
        y[i * m + 1] = x[i * n + 1] + x[i * n];
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelLinear) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelLinear);
    testLibCreation("modelLinear", m, n, 6);
}

/**
 * @test Bi-linear model of indexed variables
 */
std::vector<ADCGD> modelBiLinearIndexed(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = 3 * x[i * n];
        y[i * m + 1] = x[i * n + 1] * x[i * n];
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelBiLinearIndexed) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelBiLinearIndexed);
    testLibCreation("modelBiLinearIndexed", m, n, 6);
}

/**
 * @test Bi-linear model of indexed-nonIndexed variables
 */
std::vector<ADCGD> modelBiLinearIndexedNonIndexed(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = 3 * x[i * n];
        y[i * m + 1] = x[i * n + 1] * x[0];
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelBiLinearIndexedNonIndexed) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelBiLinearIndexedNonIndexed);
    testLibCreation("modelBiLinearIndexedNonIndexed", m, n, 6);
}

/**
 * @test Model with a bilinear term with a temporary (which depends linearly
 *       on a non-indexed) and an indexed variable
 */
std::vector<ADCGD> modelIndexedTemporary(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = 5 * x[0];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = 3 * x[i * n];
        y[i * m + 1] = x[i * n + 1] * tmp;
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelIndexedTemporary) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelIndexedTemporary);
    testLibCreation("modelIndexedTemporary", m, n, 6);
}

/**
 * @test Model with two temporary variables (which depend linearly on non-indexed)
 */
std::vector<ADCGD> modelTemporaryTemporary(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp1 = 5 * x[0];
    ADCGD tmp2 = 4 * x[2];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = 3 * tmp1 * x[i * n] * tmp2;
        y[i * m + 1] = x[i * n + 1];
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryTemporary) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelTemporaryTemporary);
    testLibCreation("modelTemporaryTemporary", m, n, 6);
}

/**
 * @test Model with a non-indexed and a temporary variable (which depends
 *       linearly on a non-indexed)
 */
std::vector<ADCGD> modelTemporaryNonIndexed(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = 4 * x[2];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = x[0] * x[i * n] * tmp;
        y[i * m + 1] = x[i * n + 1];
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryNonIndexed) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelTemporaryNonIndexed);
    testLibCreation("modelTemporaryNonIndexed", m, n, 6);
}

/**
 * @test Model with a temporary variable (which depends non-linearly on a
 *       non-indexed)
 */
std::vector<ADCGD> modelTemporary(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = 4 * x[2] * x[0];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = tmp;
        y[i * m + 1] = x[i * n + 1];
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporary) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelTemporary);
    testLibCreation("modelTemporary", m, n, 6);
}

/**
 * @test Model with bilinear term with a temporary variable (which depends 
 *       non-linearly on a non-indexed) and an indexed variable
 */
std::vector<ADCGD> modelTemporaryIndexed(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 1;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = cos(x[2]);
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = x[i * n] * tmp;
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryIndexed) {
    size_t m = 1;
    size_t n = 2;

    setModel(modelTemporaryIndexed);
    testLibCreation("modelTemporaryIndexed", m, n, 6);
}

/**
 * @test Model with 2 bilinear term with a temporary variable (which depends 
 *       non-linearly on a non-indexed) and an indexed variable
 */
std::vector<ADCGD> modelTemporaryIndexed2(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp1 = cos(x[2]);
    ADCGD tmp2 = cos(x[0]);
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = x[i * n] * tmp1;
        y[i * m + 1] = x[i * n + 1] * tmp2;
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryIndexed2) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelTemporaryIndexed2);
    testLibCreation("modelTemporaryIndexed2", m, n, 6);
}

/**
 * @test Model with a temporary variable (which depends non-linearly on a 
 *       non-indexed) and a non-linear term with an indexed variable
 */
std::vector<ADCGD> modelTemporaryIndexed3(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 1;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = cos(x[0]);
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = x[i * n + 1] / x[i * n] * tmp;
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryIndexed3) {
    size_t m = 1;
    size_t n = 2;

    setModel(modelTemporaryIndexed3);
    testLibCreation("modelTemporaryIndexed3", m, n, 6);
}

/**
 * @test Model with two temporary variables (which depends non-linearly on the
 *       same non-indexed variable) and an indexed variable
 */
std::vector<ADCGD> modelTemporaryIndexed4(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 1;
    size_t n = 1;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp1 = 2 * x[0];
    ADCGD tmp2 = 3 * x[0];
    ADCGD tmp3 = 4 * x[0];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = tmp1 * x[i * n] + x[i * n] * tmp2 + x[i * n] * tmp3;
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryIndexed4) {
    size_t m = 1;
    size_t n = 1;

    setModel(modelTemporaryIndexed4);
    testLibCreation("modelTemporaryIndexed4", m, n, 6);
}

/**
 * @test Model with 2 equations each with a bilinear term with a temporary 
 *       variable (which depends non-linearly on a non-indexed) and an indexed
 *       variable. Tests the grouping of expressions in the same if-else branch
 */
std::vector<ADCGD> modelTemporaryIndexed5(const std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp1 = cos(x[0]);
    ADCGD tmp2 = sin(x[1]);
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = x[i * n] * tmp1;
        y[i * m + 1] = x[i * n + 1] * tmp2;
    }

    return y;
}

TEST_F(CppADCGHessLoopTest, modelTemporaryIndexed5) {
    size_t m = 2;
    size_t n = 2;

    setModel(modelTemporaryIndexed5);
    testLibCreation("modelTemporaryIndexed5", m, n, 6);
}
