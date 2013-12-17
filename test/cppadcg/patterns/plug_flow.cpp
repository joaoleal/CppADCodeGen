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
#include "CppADCGPatternModelTest.hpp"
#include "../models/plug_flow.hpp"


namespace CppAD {

    const size_t nEls = 6; // number of discretization elements

    class CppADCGPatternPlugFlowTest : public CppADCGPatternModelTest {
    public:
        typedef double Base;
        typedef CppAD::CG<Base> CGD;
        typedef CppAD::AD<CGD> ADCGD;
    public:

        inline CppADCGPatternPlugFlowTest(bool verbose = false, bool printValues = false) :
            CppADCGPatternModelTest("PlugFlow",
                                    4 * nEls, // ns
                                    1, // nm
                                    4, // npar
                                    nEls * 4, // m
                                    verbose, printValues) {
            this->verbose_ = false;

            //this->epsilonA_ = std::numeric_limits<Base>::epsilon() * 1e2; 
            //this->hessianEpsilonA_ = std::numeric_limits<Base>::epsilon() * 1e2; 
            //this->hessianEpsilonR_ = std::numeric_limits<Base>::epsilon() * 1e2; 

            // states
            for (size_t j = 0; j < 6; j++) this->xb[j] = std::sqrt(14 / 1000.) - j * 0.01; // Ca (mol/l)
            for (size_t j = 0; j < 6; j++) this->xb[6 + j] = 1.2 - j * 0.01; // Cb (mol/l)
            for (size_t j = 0; j < 6; j++) this->xb[12 + j] = j * 0.01; // Cc (mol/l)
            for (size_t j = 0; j < 6; j++) this->xb[18 + j] = 50 + j * 0.1; // T (C)

            //controls
            this->xb[24] = 7; // Fin (l/min)

            //parameters
            this->xb[25] = std::sqrt(14 / 1000.); // Ca0 (mol/l)
            this->xb[26] = 0.1; // Cb0 (mol/l)
            this->xb[27] = 0.0; // Cc0 (mol/l)
            this->xb[28] = 20; // T0 (C)
        }

        virtual std::vector<ADCGD> modelFunc(const std::vector<ADCGD>& x) {
            PlugFlowModel<CGD> m;
            return m.model(x);
        }

        virtual std::vector<std::set<size_t> > getRelatedCandidates() {
            std::vector<std::set<size_t> > relatedDepCandidates(4);
            for (size_t i = 0; i < nEls; i++) {
                relatedDepCandidates[0].insert(0 * nEls + i);
                relatedDepCandidates[1].insert(1 * nEls + i);
                relatedDepCandidates[2].insert(2 * nEls + i);
                relatedDepCandidates[3].insert(3 * nEls + i);
            }
            return relatedDepCandidates;
        }

    };

}

using namespace CppAD;

/**
 * @test test the usage of loops for the generation of the plug flow model
 *       with the creation of differential information for all variables
 */
TEST_F(CppADCGPatternPlugFlowTest, plugflowAllVars) {
    modelName += "AllVars";

    /**
     * Tape model
     */
    std::auto_ptr<ADFun<CGD> > fun;
    this->tape(fun);

    //std::vector<std::set<size_t> > jacSparAll = extra::jacobianSparsitySet<std::vector<std::set<size_t> > >(*fun);
    //printSparsityPattern(jacSparAll, "jacobian", true);
    //this->customJacSparsity_.resize(fun->Range());
    //this->customJacSparsity_[2].insert(24);

    //std::vector<std::set<size_t> > hesSparAll = extra::hessianSparsitySet<std::vector<std::set<size_t> > >(*fun);
    //printSparsityPattern(hesSparAll, "hessian", true);
    //this->customHessSparsity_.resize(fun->Domain());
    //this->customHessSparsity_[0].insert(25);

    /**
     * test
     */
    this->test(*fun.get());

}

/**
 * @test test the usage of loops for the generation of the plug flow model
 *       with the creation of differential information only for states and
 *       controls
 */
TEST_F(CppADCGPatternPlugFlowTest, plugflow) {
    using namespace CppAD::extra;

    /**
     * Tape model
     */
    std::auto_ptr<ADFun<CGD> > fun;
    this->tape(fun);

    /**
     * Determine the relevant elements
     */
    this->defineCustomSparsity(*fun.get());

    /**
     * test
     */
    this->test(*fun.get());

}