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
#include "../models/cstr.hpp"
#include "../models/collocation.hpp"

namespace CppAD {
namespace cg {

/**
 * Collocation model using the CSTR model
 */
template<class T>
class CstrCollocationModel : public CollocationModel<T> {
public:

    CstrCollocationModel() :
        CollocationModel<T>(
        4, // ns
        2, // nm
        22) // npar
    {
    }

protected:

    virtual void atomicFunction(const std::vector<AD<CG<double> > >& x,
                                std::vector<AD<CG<double> > >& y) {
        y = CstrFunc(x);
    }

    virtual std::string getAtomicLibName() {
        return "cstrAtom";
    }
};

class CppADCGPatternCstrTest : public CppADCGPatternTest {
public:
    typedef double Base;
    typedef CppAD::cg::CG<Base> CGD;
    typedef CppAD::AD<CGD> ADCGD;
protected:
    static const size_t ns; // number of states in the CSTR model
    static const size_t nm; // number of controls in the CSTR model
    static const size_t npar; // number of parameters in the CSTR model

    static const size_t K; // order of collocation
    static const size_t m; // total number of equations in the collocation model 
    static const size_t na; // number of independent variables of the CSTR model
    static const size_t repeat; // number of time intervals
    std::unique_ptr<CstrCollocationModel<CGD> > colModel_;
    std::vector<Base> xx; // default CSTR model values
    std::vector<Base> x; // values for the collocation model
    std::unique_ptr<DynamicLib<double> > atomicDynamicLib_;
    std::unique_ptr<GenericModel<double> > atomicModel_;
    std::vector<std::string> flags_;
public:

    inline CppADCGPatternCstrTest(bool verbose = false, bool printValues = false) :
        CppADCGPatternTest(verbose, printValues),
        colModel_(new CstrCollocationModel<CGD>()),
        xx(na) {
        //this->verbose_ = true;
        this->hessianEpsilonA_ = std::numeric_limits<Base>::epsilon() * 4e6;
        this->hessianEpsilonR_ = std::numeric_limits<Base>::epsilon() * 1e3;

        flags_.push_back("-O0");
        flags_.push_back("-g");
        flags_.push_back("-ggdb");
        flags_.push_back("-D_FORTIFY_SOURCE=2");

        /**
         * CSTR model values
         */
        xx[0] = 0.3; // h
        xx[1] = 7.82e3; // Ca
        xx[2] = 304.65; // Tr
        xx[3] = 301.15; // Tj

        xx[4] = 2.3333e-04; // u1
        xx[5] = 6.6667e-05; // u2

        xx[6] = 6.2e14; // 
        xx[7] = 10080; //
        xx[8] = 2e3; //
        xx[9] = 10e3; //
        xx[10] = 1e-11; //
        xx[11] = 6.6667e-05; //
        xx[12] = 294.15; //
        xx[13] = 294.15; //
        xx[14] = 1000; //
        xx[15] = 4184; //Cp
        xx[16] = -33488; //deltaH
        xx[17] = 299.15; // Tj0
        xx[18] = 302.65; //   Tj2
        xx[19] = 7e5; // cwallj
        xx[20] = 1203; // csteam
        xx[21] = 3.22; //dsteam
        xx[22] = 950.0; //Ug
        xx[23] = 0.48649427192323; //vc6in
        xx[24] = 1000; //rhoj
        xx[25] = 4184; //Cpj
        xx[26] = 0.014; //Vj
        xx[27] = 1e-7; //cwallr

        colModel_->setTypicalAtomModelValues(xx);

        /**
         * collocation model values
         */
        xNorm_ = colModel_->getTypicalValues(repeat);
        x.resize(xNorm_.size(), 1.0);

#if 0
        x = xNorm_;
        xNorm_.clear();
#else
        for (size_t j = 0; j < xNorm_.size(); j++) {
            xNorm_[j] = xNorm_[j] != 0.0 ? xNorm_[j] : 1.0;
        }

        eqNorm_.resize(repeat * m);
        size_t e = 0;
        for (size_t i = 0; i < repeat; i++) {
            for (size_t k = 0; k < K; k++) {
                for (size_t j = 0; j < ns; j++) {
                    eqNorm_[e++] = xx[j];
                }
            }
        }
#endif

        this->setModel(*colModel_.get());
    }

    virtual void TearDown() {
        atomicDynamicLib_.reset(nullptr);
        atomicModel_.reset(nullptr);
        colModel_.reset(nullptr);

        CppADCGPatternTest::TearDown();
    }

    inline virtual void defineCustomSparsity(ADFun<CGD>& fun) {
        CppADCGPatternTest::defineCustomSparsity(fun);
        //std::vector<std::set<size_t> > hessSparAll = extra::hessianSparsitySet<std::vector<std::set<size_t> > >(fun);
        //printSparsityPattern(hessSparAll, "Full Hessian");
        //customHessSparsity_.resize(x.size());
        //customHessSparsity_[22].insert(22);
    }

};

const size_t CppADCGPatternCstrTest::ns = 4;
const size_t CppADCGPatternCstrTest::nm = 2;
const size_t CppADCGPatternCstrTest::npar = 22;

const size_t CppADCGPatternCstrTest::K = 3;
const size_t CppADCGPatternCstrTest::m = K * ns;
const size_t CppADCGPatternCstrTest::na = ns + nm + npar;

const size_t CppADCGPatternCstrTest::repeat = 6;

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;

/**
 * @test test the usage of loops for the generation of a orthogonal collocation
 *       model used to integrate a CSTR ODE model
 *       with differential information for all variables (including parameters)
 */
TEST_F(CppADCGPatternCstrTest, AtomicAllVars) {
    using namespace CppAD;

    /**
     * create atomic function for the ODE
     */
    colModel_->setIgnoreParameters(false);
    colModel_->setCompilerFlags(flags_);
    colModel_->createAtomicLib();
    atoms_.push_back(&colModel_->getDoubleAtomic());

    testPatternDetection(m, x, repeat);
    testLibCreation("modelCstrAtomicAllVars", m, repeat, x);
}

/**
 * @test test the usage of loops for the generation of a orthogonal collocation
 *       model used to integrate a CSTR ODE model,
 *       only with differential information for states and controls
 */
TEST_F(CppADCGPatternCstrTest, Atomic) {

    using namespace CppAD;

    /**
     * create atomic function for the ODE
     */
    colModel_->setIgnoreParameters(true);
    colModel_->setCompilerFlags(flags_);
    colModel_->createAtomicLib();
    atoms_.push_back(&colModel_->getDoubleAtomic());

    testPatternDetection(m, x, repeat);
    testLibCreation("modelCstrAtomic", m, repeat, x);

}