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

using namespace CppAD;

/**
 * @test All variables are indexed, no temporaries
 */
std::vector<ADCGD> model0(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * n]);
        y[i * m + 1] = x[i * n + 1] * x[i * n];
    }

    return y;
}

TEST_F(CppADCGPatternTest, DependentPatternMatcherDetached) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(model0, m, n, 6);
    testLibCreation("model0", model0, m, n, 6);
}

/**
 * @test All variables have a random index, no temporaries
 */
std::vector<ADCGD> modelRandom(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;
    vector<unsigned long> index0(repeat);
    vector<unsigned long> index1(repeat);
    for (size_t i = 0; i < repeat; i++) {
        if (i % 2 == 0) {
            index0[i] = i * n;
            index1[i] = i * n + 1;
        } else {
            index0[i] = (repeat - i) * n;
            index1[i] = (repeat - i) * n - 1;
        }
    }

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[index0[i]]);
        y[i * m + 1] = x[index1[i]] * x[index0[i]];
    }

    return y;
}

TEST_F(CppADCGPatternTest, random) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(modelRandom, m, n, 10);
    testLibCreation("modelRandom", modelRandom, m, n, 10);
}

/**
 * @test Some variables not indexed -> one constant temporary
 */
std::vector<ADCGD> model1(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * n]) + x[1] * x[2];
        y[i * m + 1] = x[i * n + 1] * x[i * n];
    }

    return y;
}

TEST_F(CppADCGPatternTest, DependentPatternMatcher) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(model1, m, n, 6);
    testLibCreation("model1", model1, m, n, 6);
}

/**
 * @test Some variables not indexed -> one constant temporary (defined outside 
 *       loop)
 */
std::vector<ADCGD> modelCommonTmp(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp = x[1] * x[2];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * n]) + tmp;
        y[i * m + 1] = x[i * n + 1] * x[i * n] + tmp;
    }

    return y;
}

TEST_F(CppADCGPatternTest, CommonTmp) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(modelCommonTmp, m, n, 6);
    testLibCreation("modelCommonTmp", modelCommonTmp, m, n, 6);
}

/**
 * @test Some variables not indexed -> 2 constant temporaries (defined outside
 *       loop) independents used by temporaries are not the same as the ones in
 *       the equation inside the loop
 */
std::vector<ADCGD> modelCommonTmp2(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    ADCGD tmp1 = x[1] * sin(x[1] + 4);
    ADCGD tmp2 = 2.5 * x[1];
    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = tmp1 * cos(x[i * n]) / tmp2;
        y[i * m + 1] = x[i * n + 1] * log(x[i * n + 1]);
    }

    return y;
}

TEST_F(CppADCGPatternTest, CommonTmp2) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(modelCommonTmp2, m, n, 6);
    testLibCreation("modelCommonTmp2", modelCommonTmp2, m, n, 6);
}

/**
 * @test All variables are indexed -> 3 temporaries (defined inside and outside
 *       the loop)
 */
std::vector<ADCGD> modelCommonTmp3(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 3;
    size_t n = 3;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);
    ADCGD tmp1 = sin(x[1] + 4);
    for (size_t i = 0; i < repeat; i++) {

        ADCGD tmp2 = 2.5 * x[i * n] + tmp1;
        ADCGD tmp3 = 5 * x[i * n];

        y[i * m] = cos(x[i * n]) / tmp2 + tmp3;
        y[i * m + 1] = x[i * n + 1] * log(x[i * n + 1]) * tmp3;
        y[i * m + 2] = x[i * n + 2] + tmp2;
    }

    return y;
}

TEST_F(CppADCGPatternTest, CommonTmp3) {
    size_t m = 3;
    size_t n = 3;

    testPatternDetection(modelCommonTmp3, m, n, 6);
    testLibCreation("modelCommonTmp3", modelCommonTmp3, m, n, 6);
}

/**
 * @test All variables are indexed but keep the same index after a given
 *       iteration
 */
std::vector<ADCGD> model2(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 1;
    size_t n = 2;

    // dependent variable vector 
    std::vector<ADCGD> y(repeat * m);

    ADCGD ax0;
    ADCGD ax1;
    for (size_t i = 0; i < repeat; i++) {
        if (i < repeat / 2) {
            ax0 = x[i * n];
            ax1 = x[i * n + 1];
        }
        y[i * m] = ax0 + ax0 * ax1;
    }

    return y;
}

TEST_F(CppADCGPatternTest, model2) {
    size_t m = 1;
    size_t n = 2;

    testPatternDetection(model2, m, n, 6);
    testLibCreation("model2", model2, m, n, 6);
}

/**
 * @test Some variables are not indexed, some indexed variables keep the same
 *       index after a given iteration
 */
std::vector<ADCGD> model3(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 1;
    size_t n = 2;

    // dependent variable vector 
    std::vector<ADCGD> y(repeat * m);

    ADCGD ax1;
    for (size_t i = 0; i < repeat; i++) {
        if (i < repeat / 2) {
            ax1 = x[i * n + 1];
        }
        y[i * m] = x[0] + x[0] * ax1 * ax1;
    }

    return y;
}

TEST_F(CppADCGPatternTest, model3) {
    size_t m = 1;
    size_t n = 2;

    testPatternDetection(model3, m, n, 6);
    testLibCreation("model3", model3, m, n, 6);
}

/**
 * @test 4 equations; 1 constant temporary; 1 equation with a constant value
 */
std::vector<ADCGD> model4Eq(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 4;
    size_t n = 4;
    size_t m2 = repeat * m;

    assert(x.size() == n * repeat);

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * n]) + 3 * x[1] * log(x[2]);
        y[i * m + 1] = x[i * n + 1] * x[i * n];
        y[i * m + 2] = x[i * n + 1] * x[i * n + 2];
        y[i * m + 3] = 5;
    }

    return y;
}

TEST_F(CppADCGPatternTest, Matcher4Eq) {
    size_t m = 4;
    size_t n = 4;

    testPatternDetection(model4Eq, m, n, 6);
    testLibCreation("model4Eq", model4Eq, m, n, 6);
}

/**
 * @test 
 */
std::vector<ADCGD> model4(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        ADCGD tmp = x[1] * x[i];
        y[i * m] = cos(x[i * n]) + tmp;
        y[i * m + 1] = x[i * n + 1] * x[i * n] + tmp;
    }

    return y;
}

TEST_F(CppADCGPatternTest, IndexedTmp) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(model4, m, n, 6);
    testLibCreation("indexedTmp", model4, m, n, 6);
}

/**
 * @test 2 loops required
 */
std::vector<ADCGD> model5(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2);

    for (size_t i = 0; i < repeat; i++) {
        ADCGD tmp = x[1] * x[i];
        y[i * m] = cos(x[i * n]) + tmp;

        if (i == 1) {
            for (size_t i2 = 0; i2 < repeat; i2++) {
                y[i2 * m + 1] = x[i2 * n + 1] * x[i2 * n] + tmp;
            }
        }
    }

    return y;
}

TEST_F(CppADCGPatternTest, DependentPatternMatcher5) {
    size_t m = 2;
    size_t n = 2;

    testPatternDetection(model5, m, n, 6, 2);
    testLibCreation("model5", model5, m, n, 6);
}

/**
 * @test using atomic functions
 */
std::vector<ADCGD> modelAtomic(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<double>*>& atoms) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    CGAbstractAtomicFun<Base>& atomic0 = *atoms[0];

    // dependent variable vector 
    std::vector<ADCGD> y(m2), ax(2), ay(1);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * n]);

        ax[0] = x[i * n];
        ax[1] = x[i * n + 1];
        atomic0(ax, ay);
        y[i * m + 1] = ay[0];
    }

    return y;
}

void atomicFunction(const std::vector<AD<double> >& x, std::vector<AD<double> >& y) {
    y[0] = x[1] * x[0];
}

TEST_F(CppADCGPatternTest, Atomic) {
    using namespace CppAD;

    size_t m = 2;
    size_t n = 2;

    // create atomic function
    std::vector<AD<double> > y(1), x(2);
    checkpoint<double> atomicfun("atomicFunc", atomicFunction, x, y);
    std::vector<atomic_base<double>*> atomics(1);
    atomics[0] = &atomicfun;


    testPatternDetectionWithAtomics(modelAtomic, atomics, m, n, 6);
    testLibCreationWithAtomics("modelAtomic", modelAtomic, atomics, m, n, 6);
}

/**
 * @using atomic functions; some indexed variables keep the same index after a
 *        given iteration
 */
std::vector<ADCGD> modelAtomic2(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<double>*>& atoms) {
    size_t m = 1;
    size_t n = 2;
    size_t m2 = repeat * m;

    CGAbstractAtomicFun<Base>& atomic0 = *atoms[0];

    // dependent variable vector 
    std::vector<ADCGD> y(m2), ax(2), ay(1);

    for (size_t i = 0; i < repeat; i++) {
        ax[0] = x[i * n];
        if (i < repeat / 2) {
            ax[1] = x[i * n + 1];
        }
        atomic0(ax, ay);
        y[i * m] = ay[0];
    }

    return y;
}

TEST_F(CppADCGPatternTest, Atomic2) {
    using namespace CppAD;

    size_t m = 1;
    size_t n = 2;

    // create atomic function
    std::vector<AD<double> > y(1), x(2);
    checkpoint<double> atomicfun("atomicFunc2", atomicFunction, x, y);
    std::vector<atomic_base<double>*> atomics(1);
    atomics[0] = &atomicfun;


    testPatternDetectionWithAtomics(modelAtomic2, atomics, m, n, 6);
    testLibCreationWithAtomics("modelAtomic2", modelAtomic2, atomics, m, n, 6);
}

/**
 * @test with an additional equation that does not belong to a loop
 */
std::vector<ADCGD> modelExtraFunc(std::vector<ADCGD>& x, size_t repeat) {
    size_t m = 2;
    size_t n = 2;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2 + 1);

    for (size_t i = 0; i < repeat; i++) {
        y[i * m] = cos(x[i * n]);
        y[i * m + 1] = 10 * x[i * n + 1] * x[i * n];
    }

    y.back() = log(x[0] * x[1]) + 3 * x[3];

    return y;
}

TEST_F(CppADCGPatternTest, modelExtraFunc) {
    size_t m = 2;
    size_t n = 2;
    size_t mExtra = 1; // equations not in loops

    testPatternDetection(modelExtraFunc, m, n, 6);
    testLibCreation("modelExtraFunc", modelExtraFunc, m, n, 6, mExtra);
}