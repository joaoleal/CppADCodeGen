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
            this->jacobian_ = IGNORE;
        }
    };
}

using namespace CppAD;

/**
 * @test Linear model (no Hessian)
 */
std::vector<ADCGD> modelLinear(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelLinear", modelLinear, m, n, 6);
}

/**
 * @test Bi-linear model of indexed variables
 */
std::vector<ADCGD> modelBiLinearIndexed(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelBiLinearIndexed", modelBiLinearIndexed, m, n, 6);
}

/**
 * @test Bi-linear model of indexed-nonIndexed variables
 */
std::vector<ADCGD> modelBiLinearIndexedNonIndexed(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelBiLinearIndexedNonIndexed", modelBiLinearIndexedNonIndexed, m, n, 6);
}

/**
 * @test Model with a bilinear term with a temporary (which depends linearly
 *       on a non-indexed) and an indexed variable
 */
std::vector<ADCGD> modelIndexedTemporary(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelIndexedTemporary", modelIndexedTemporary, m, n, 6);
}

/**
 * @test Model with two temporary variables (which depend linearly on non-indexed)
 */
std::vector<ADCGD> modelTemporaryTemporary(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelTemporaryTemporary", modelTemporaryTemporary, m, n, 6);
}

/**
 * @test Model with a non-indexed and a temporary variable (which depends
 *       linearly on a non-indexed)
 */
std::vector<ADCGD> modelTemporaryNonIndexed(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelTemporaryNonIndexed", modelTemporaryNonIndexed, m, n, 6);
}

/**
 * @test Model with a temporary variable (which depends non-linearly on a
 *       non-indexed)
 */
std::vector<ADCGD> modelTemporary(std::vector<ADCGD>& x, size_t repeat) {
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

    testLibCreation("modelTemporary", modelTemporary, m, n, 6);
}