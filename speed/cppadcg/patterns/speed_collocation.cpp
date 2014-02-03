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
#include "../../../test/cppad/cg/models/collocation.hpp"
#include "../../../test/cppad/cg/models/plug_flow.hpp"

namespace CppAD {
namespace cg {

typedef double Base;
typedef CppAD::cg::CG<Base> CGD;

using namespace std;

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
        y = m.model(x, nEls_);
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

        GenericModel<Base>* atom = modelCppAD_.getGenericModel();
        std::vector<GenericModel<Base>*> atoms(1);
        atoms[0] = atom;
        setExternalModels(atoms);
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

} // END cg namespace
} // END CppAD namespace

using namespace CppAD::cg;

int main(int argc, char **argv) {
    size_t repeat = PatternSpeedTest::parseProgramArguments(1, argc, argv, 10); // time intervals
    size_t nEls = PatternSpeedTest::parseProgramArguments(2, argc, argv, 10); // number of CSTR elements

    size_t K = 3;
    size_t ns = 5;
    CollocationPatternSpeedTest speed(nEls);
    speed.setNumberOfExecutions(30);
#if 0
    speed.preparation = false;
    speed.zeroOrder = false;
    speed.sparseJacobian = false;
    speed.cppADCG = false;
    speed.cppADCGLoops = true;
    speed.cppADCGLoopsLlvm = false;
    std::vector<std::string> compileFlags(3);
    compileFlags[0] = "-O2";
    compileFlags[1] = "-g";
    compileFlags[2] = "-ggdb";
    speed.setCompileFlags(compileFlags);
#endif
    speed.measureSpeed(K * ns * nEls, repeat, speed.getTypicalValues(repeat));
}
