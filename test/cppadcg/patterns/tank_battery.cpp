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
#include "../models/tank_battery.hpp"


namespace CppAD {

    class CppADCGPatternTankBatTest : public CppADCGPatternTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    protected:
        const size_t nTanks; // number of stages
        const size_t ns;
        const size_t nm;
        const size_t m; // total number of equations in the model 
        const size_t n; // number of independent variables
        std::vector<Base> xb; // values for the model
    public:

        inline CppADCGPatternTankBatTest(bool verbose = false, bool printValues = false) :
            CppADCGPatternTest(verbose, printValues),
            nTanks(6),
            ns(nTanks),
            nm(1),
            m(nTanks * 1),
            n(ns + nm + 1),
            xb(n, 1.0) {

            // this->testZeroOrder_ = false;
            // this->testJacobian_ = false;
            // this->testHessian_ = false;
            this->verbose_ = true;
        }

        inline void tape(std::auto_ptr<ADFun<CGD> >& fun) {
            /**
             * Tape model
             */
            std::vector<ADCGD> x(xb.size());
            for (size_t j = 0; j < xb.size(); j++)
                x[j] = xb[j];
            CppAD::Independent(x);
            if (xNorm_.size() > 0) {
                ASSERT_EQ(x.size(), xNorm_.size());
                for (size_t j = 0; j < x.size(); j++)
                    x[j] *= xNorm_[j];
            }

            std::vector<ADCGD> y = tankBatteryFunc(x);
            if (eqNorm_.size() > 0) {
                ASSERT_EQ(y.size(), eqNorm_.size());
                for (size_t i = 0; i < y.size(); i++)
                    y[i] /= eqNorm_[i];
            }

            fun.reset(new ADFun<CGD>());
            fun->Dependent(y);
        }

        /**
         * test
         */
        inline void test(ADFun<CGD>& fun) {

            std::string libName = "modelTankBattery";
            std::vector<atomic_base<Base>*> atoms;

            std::vector<std::set<size_t> > relatedDepCandidates(1);
            for (size_t i = 0; i < nTanks; i++) relatedDepCandidates[0].insert(i);

            testSourceCodeGen(fun, relatedDepCandidates, libName, atoms, xb, FORWARD, testJacobian_, testHessian_);
            if (testJacobian_) {
                testSourceCodeGen(fun, relatedDepCandidates, libName, atoms, xb, FORWARD, true, false, true);
                testSourceCodeGen(fun, relatedDepCandidates, libName, atoms, xb, REVERSE, true, false);
                testSourceCodeGen(fun, relatedDepCandidates, libName, atoms, xb, REVERSE, true, false, true);
            }

            if (testHessian_) {
                testSourceCodeGen(fun, relatedDepCandidates, libName, atoms, xb, FORWARD, false, true, false, true);
            }
        }

    };

}

using namespace CppAD;

/**
 * @test test the usage of loops for the generation of distillation model
 */
TEST_F(CppADCGPatternTankBatTest, tankBatteryAllVars) {
    using namespace CppAD;

    /**
     * Tape model
     */
    std::auto_ptr<ADFun<CGD> > fun;
    tape(fun);

    /**
     * test
     */
   // test(*fun.get());

}

TEST_F(CppADCGPatternTankBatTest, tankBattery) {
    using namespace CppAD::extra;

    std::auto_ptr<ADFun<CGD> > fun;
    tape(fun);

    /**
     * Determine the relevant elements
     */
    std::vector<std::set<size_t> > jacSparAll = jacobianSparsitySet<std::vector<std::set<size_t> > >(*fun.get());
    customJacSparsity_.resize(jacSparAll.size());
    for (size_t i = 0; i < jacSparAll.size(); i++) {
        // only differential information for states and controls
        std::set<size_t>::const_iterator itEnd = jacSparAll[i].upper_bound(ns + nm - 1);
        if (itEnd != jacSparAll[i].begin())
            customJacSparsity_[i].insert(jacSparAll[i].begin(), itEnd);
    }

    std::vector<std::set<size_t> > hessSparAll = hessianSparsitySet<std::vector<std::set<size_t> > >(*fun.get());
    customHessSparsity_.resize(hessSparAll.size());
    for (size_t i = 0; i < ns + nm; i++) {
        std::set<size_t>::const_iterator it = hessSparAll[i].upper_bound(i); // only the lower left side
        if (it != hessSparAll[i].begin())
            customHessSparsity_[i].insert(hessSparAll[i].begin(), it);
    }

    /**
     * test
     */
    test(*fun.get());

}