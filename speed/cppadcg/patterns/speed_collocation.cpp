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

#include "pattern_speed_test.hpp"
#include "../../../test/cppadcg/models/collocation.hpp"
#include "../../../test/cppadcg/models/plug_flow.hpp"

using namespace CppAD;
using namespace std;

typedef double Base;
typedef CppAD::CG<Base> CGD;

/**
 * Collocation model using the Plugflow model
 */
template<class T>
class PlugFlowCollocationModel : public CollocationModel<T> {
protected:
    size_t nEls_; // number of plugflow discretization elements
public:

    PlugFlowCollocationModel(size_t nEls) :
        CollocationModel<T>(5 * nEls, // ns
        2, // nm
        5), // npar
        nEls_(nEls) {
    }

protected:

    virtual void atomicFunction(const std::vector<AD<CG<double> > >& x,
                                std::vector<AD<CG<double> > >& y) {
        PlugFlowModel<CG<double> > m;
        y = m.model2(x, nEls_);
    }

    virtual std::string getAtomicLibName() {
        return "plugflow";
    }
};

/**
 * Speed test for the collocation model
 */
class CollocationPatternSpeedTest : public PatternSpeedTest {
protected:
    size_t nEls_;
    PlugFlowCollocationModel<double> modelCppAD_;
    PlugFlowCollocationModel<CG<double> > modelCppADCG_;
public:

    inline CollocationPatternSpeedTest(size_t nEls, bool verbose = false) :
        PatternSpeedTest("collocation", verbose),
        nEls_(nEls),
        modelCppAD_(nEls),
        modelCppADCG_(nEls) {

        /**
         * Prepare the atomic model
         */
        PlugFlowModel<CG<double> > m;
        modelCppAD_.setTypicalAtomModelValues(m.getTypicalValues(nEls_));
        modelCppADCG_.setTypicalAtomModelValues(m.getTypicalValues(nEls_));

        modelCppAD_.createAtomicLib();
        modelCppADCG_.createAtomicLib();

        CGAtomicGenericModel<double>& atom = modelCppAD_.getDoubleAtomic();
        std::vector<atomic_base<Base>*> atoms(1);
        atoms[0] = &atom;
        setAtomics(atoms);
    }

    inline std::vector<double> getTypicalValues(size_t repeat) {
        return modelCppAD_.getTypicalValues(repeat);
    }

    virtual std::vector<AD<CGD> > modelCppADCG(const std::vector<AD<CGD> >& x, size_t repeat) {
        return modelCppADCG_.evaluateModel(x, repeat);
    }

    virtual std::vector<AD<Base> > modelCppAD(const std::vector<AD<Base> >& x, size_t repeat) {
        return modelCppAD_.evaluateModel(x, repeat);
    }
};

int main(int argc, char **argv) {
    size_t repeat = PatternSpeedTest::parseProgramArguments(1, argc, argv, 10); // time intervals
    size_t nEls = PatternSpeedTest::parseProgramArguments(2, argc, argv, 10); // number of CSTR elements
    
    size_t K = 3;
    size_t ns = 5;
    CollocationPatternSpeedTest speed(nEls);
    speed.setNumberOfExecutions(30);
    speed.measureSpeed(K * ns * nEls, repeat, speed.getTypicalValues(repeat));
}
