/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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
#include "../../../../test/cppad/cg/models/plug_flow.hpp"

using namespace CppAD;
using namespace CppAD::cg;
using namespace std;

using Base = double;
using CGD = CppAD::cg::CG<Base>;

class PlugFlowPatternSpeedTest : public PatternSpeedTest {
public:

    explicit PlugFlowPatternSpeedTest(bool verbose = false) :
        PatternSpeedTest("plugflow", verbose) {
    }

    std::vector<AD<CGD> > modelCppADCG(const std::vector<AD<CGD> >& x,
                                       const std::vector<AD<CGD> >& p,
                                       size_t repeat) override {
        PlugFlowModel<CGD> m;
        return m.model2(x, repeat);
    }

    std::vector<AD<Base> > modelCppAD(const std::vector<AD<Base> >& x,
                                      const std::vector<AD<Base> >& p,
                                      size_t repeat) override {
        PlugFlowModel<Base> m;
        return m.model2(x, repeat);
    }
};

int main(int argc, char **argv) {
    size_t nEles = PatternSpeedTest::parseProgramArguments(1, argc, argv, 10);

    std::vector<Base> x = PlugFlowModel<Base>::getTypicalValues(nEles);
    std::vector<Base> p;
    std::vector<std::set<size_t> > relations = PlugFlowModel<Base>::getRelatedCandidates(nEles);

    std::vector<std::string> flags;
    //flags.push_back("-O2");

    PlugFlowPatternSpeedTest speed;
    //speed.cppAD = false;
    //speed.cppADCG = true;
    //speed.cppADCGLoops = false;
    //speed.cppADCGLoopsLlvm = false;
    //speed.zeroOrder = false;
    //speed.sparseJacobian = false;
    //speed.sparseHessian = false;
    speed.setNumberOfExecutions(30);
    speed.setCompileFlags(flags);
    speed.measureSpeed(relations, nEles, x, p);
}
