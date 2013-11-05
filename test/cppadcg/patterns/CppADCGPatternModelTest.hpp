#ifndef CPPAD_CG_TEST_CPPADCGPATTERNMODELTEST_INCLUDED
#define	CPPAD_CG_TEST_CPPADCGPATTERNMODELTEST_INCLUDED
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

namespace CppAD {

    class CppADCGPatternModelTest : public CppADCGPatternTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    protected:
        std::string modelName;
        const size_t ns;
        const size_t nm;
        const size_t m; // total number of equations in the model 
        const size_t n; // number of independent variables
        std::vector<Base> xb; // values for the model
    public:

        inline CppADCGPatternModelTest(const std::string& modelName_,
                                       const size_t ns_,
                                       const size_t nm_,
                                       const size_t npar_,
                                       const size_t m_,
                                       bool verbose = false,
                                       bool printValues = false) :
            CppADCGPatternTest(verbose, printValues),
            modelName(modelName_),
            ns(ns_),
            nm(nm_),
            m(m_),
            n(ns_ + nm_ + npar_),
            xb(n) {
            for (size_t j = 0; j < n; j++)
                xb[j] = 0.5 * (j + 1);
        }

        virtual std::vector<ADCGD> modelFunc(const std::vector<ADCGD>& x) = 0;

        virtual std::vector<std::set<size_t> > getRelatedCandidates() = 0;

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

            std::vector<ADCGD> y = modelFunc(x);
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

            std::string libName = "model" + modelName;
            std::vector<atomic_base<Base>*> atoms;

            std::vector<std::set<size_t> > relatedDepCandidates = getRelatedCandidates();

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

        inline void defineCustomSparsity(ADFun<CGD>& fun) {
            /**
             * Determine the relevant elements
             */
            std::vector<std::set<size_t> > jacSparAll = extra::jacobianSparsitySet<std::vector<std::set<size_t> > >(fun);
            customJacSparsity_.resize(jacSparAll.size());
            for (size_t i = 0; i < jacSparAll.size(); i++) {
                // only differential information for states and controls
                std::set<size_t>::const_iterator itEnd = jacSparAll[i].upper_bound(ns + nm - 1);
                if (itEnd != jacSparAll[i].begin())
                    customJacSparsity_[i].insert(jacSparAll[i].begin(), itEnd);
            }

            std::vector<std::set<size_t> > hessSparAll = extra::hessianSparsitySet<std::vector<std::set<size_t> > >(fun);
            customHessSparsity_.resize(hessSparAll.size());
            for (size_t i = 0; i < ns + nm; i++) {
                std::set<size_t>::const_iterator it = hessSparAll[i].upper_bound(i); // only the lower left side
                if (it != hessSparAll[i].begin())
                    customHessSparsity_[i].insert(hessSparAll[i].begin(), it);
            }
        }

    };

}

#endif