/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2015 Ciengis
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

using namespace CppAD;
using namespace CppAD::cg;

typedef double Base;
typedef CppAD::cg::CG<Base> CGD;
typedef CppAD::AD<CGD> ADCGD;

namespace CppAD {
namespace cg {

class CppADCGPatternTestLoopDep : public CppADCGPatternTest {
protected:

    inline virtual void defineCustomSparsity(ADFun<CGD>& fun) {
        /**
         * only the lower left side (not currently required since the hessian is diagonal!)
         */
        std::vector<std::set<size_t> > hessSparAll = extra::hessianSparsitySet<std::vector<std::set<size_t> > >(fun);
        customHessSparsity_.resize(hessSparAll.size());
        for (size_t i = 0; i < hessSparAll.size(); i++) {
            std::set<size_t>::const_iterator it = hessSparAll[i].upper_bound(i); // only the lower left side
            if (it != hessSparAll[i].begin())
                customHessSparsity_[i].insert(hessSparAll[i].begin(), it);
        }
    }
};

}
}

/**
 * outer function
 */
std::vector<ADCGD> outerModel(const std::vector<ADCGD>& x,
                              size_t repeat,
                              atomic_base<CGD>& atomic) {
    size_t m = 4;
    size_t n = 4;
    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> y(m2), ax(4), ay(4);

    for (size_t i = 0; i < repeat; i++) {
        ax[0] = x[i * n];
        ax[1] = x[i * n + 1 ];
        ax[2] = x[i * n + 2 ];
        ax[3] = x[i * n + 3 ];
        atomic(ax, ay);
        y[i * m] = ay[0];
        y[i * m + 1] = ay[1];
        y[i * m + 2] = ay[2];
        y[i * m + 3] = ay[3];
    }

    return y;
}

void atomicFunction(const std::vector<AD<double> >& x, std::vector<AD<double> >& y) {
    y[0] = x[0] * x[0];
    y[1] = 1 * x[1] * x[1];
    y[2] = 2 * x[2] * x[2];
    y[3] = 3 * x[3] * x[3];
}

/**
 * @test Test the correct generation of loop inside the detected loops for 
 *       assignment of indexed dependent from arrays
 */
TEST_F(CppADCGPatternTestLoopDep, SimpleAtomic) {
    using namespace CppAD;


    size_t m = 4;
    size_t n = 4;
    size_t repeat = 6;

    // create atomic function
    std::vector<AD<double> > y(4), x(n);
    checkpoint<double> atomicfun("atomicFunc2", atomicFunction, x, y);
    CGAtomicFun<double> cgatomicfun(atomicfun, true);
    PatternTestModelWithAtom<CGD> model(outerModel, cgatomicfun);
    setModel(model);
    this->atoms_.push_back(&atomicfun);

    testPatternDetection(m, n, repeat);
    testLibCreation("modelSimpleAtomicDepLoop", m, n, repeat);
}
