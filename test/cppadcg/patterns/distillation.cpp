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
#include "../models/distillation.hpp"


namespace CppAD {

    class CppADCGPatternDistillationTest : public CppADCGPatternTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    protected:
        const size_t nStage; // number of stages
        const size_t ns;
        const size_t nm;
        const size_t m; // total number of equations in the model 
        const size_t n; // number of independent variables
        std::vector<Base> xb; // values for the model
    public:

        inline CppADCGPatternDistillationTest(bool verbose = false, bool printValues = false) :
            CppADCGPatternTest(verbose, printValues),
            nStage(8),
            ns(nStage * 6),
            nm(4),
            m(nStage * 6),
            n(ns + nm + 4),
            xb(n, 1.0) {

            // this->testZeroOrder_ = false;
            // this->testJacobian_ = false;
            // this->testHessian_ = false;
            this->verbose_ = true;
            //this->epsilonA_ = std::numeric_limits<Base>::epsilon() * 4e3;

            xNorm_.resize(n, 1.0);

            size_t j = 0;
            for (size_t i = 0; i < nStage; i++, j++) xNorm_[j] = 12000; // mWater
            for (size_t i = 0; i < nStage; i++, j++) xNorm_[j] = 12000; // mEthanol[i]
            for (size_t i = 0; i < nStage; i++, j++) xNorm_[j] = 360; // T[i]
            for (size_t i = 0; i < nStage; i++, j++) xNorm_[j] = 0.5; // yWater[i]
            for (size_t i = 0; i < nStage; i++, j++) xNorm_[j] = 0.5; // yEthanol[i]
            for (size_t i = 0; i < nStage - 1; i++, j++) xNorm_[j] = 8; // V[i]
            xNorm_[j++] = 150e3; // Qc
            assert(j == ns);

            xNorm_[j++] = 250e3; // Qsteam 
            xNorm_[j++] = 0.1; // Fdistillate
            xNorm_[j++] = 2.5; // reflux
            xNorm_[j++] = 4; // Frectifier

            xNorm_[j++] = 30; // feed
            xNorm_[j++] = 1.01325e5; // P
            xNorm_[j++] = 0.7; // xFWater
            xNorm_[j++] = 366; // Tfeed

            j = 0;
            for (size_t i = 0; i < nStage; i++, j++) xb[j] = (12000 + 100 * i) / xNorm_[j]; // mWater
            for (size_t i = 0; i < nStage; i++, j++) xb[j] = (12000 - 100 * i) / xNorm_[j]; // mEthanol[i]
            for (size_t i = 0; i < nStage; i++, j++) xb[j] = (360 + i * 2) / xNorm_[j]; // T[i]
            for (size_t i = 0; i < nStage; i++, j++) xb[j] = (0.3 + 0.05 * i) / xNorm_[j]; // yWater[i]
            for (size_t i = 0; i < nStage; i++, j++) xb[j] = (0.7 - 0.05 * i) / xNorm_[j]; // yEthanol[i]
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

            std::vector<ADCGD> y = distillationFunc(x);
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

            std::string libName = "modelDistillation";
            std::vector<atomic_base<Base>*> atoms;

            std::vector<std::set<size_t> > relatedDepCandidates(6);
            size_t j = 0;
            for (size_t i = 0; i < nStage; i++, j++) relatedDepCandidates[0].insert(j); // mWater
            for (size_t i = 0; i < nStage; i++, j++) relatedDepCandidates[1].insert(j); // mEthanol
            for (size_t i = 0; i < nStage; i++, j++) relatedDepCandidates[4].insert(j); // T
            for (size_t i = 0; i < nStage; i++, j++) relatedDepCandidates[2].insert(j); // yWater
            for (size_t i = 0; i < nStage; i++, j++) relatedDepCandidates[3].insert(j); // yEthanol
            for (size_t i = 0; i < nStage - 1; i++, j++) relatedDepCandidates[5].insert(j); // V


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
TEST_F(CppADCGPatternDistillationTest, distillationAllVars) {
    /**
     * Tape model
     */
    std::auto_ptr<ADFun<CGD> > fun;
    tape(fun);

    /**
     * test
     */
    test(*fun.get());

}

TEST_F(CppADCGPatternDistillationTest, distillation) {
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