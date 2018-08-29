/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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
#include "CppADCGDynamicTest.hpp"

namespace CppAD {
namespace cg {

class CppADCGDynamicTest1 : public CppADCGDynamicTest {
public:

    inline CppADCGDynamicTest1(bool verbose = false, bool printValues = false) :
        CppADCGDynamicTest("dynamic_with_params", verbose, printValues) {
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ax) override {
        return std::vector<ADCGD>(); // not used
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ax,
                             const std::vector<ADCGD>& ap) override {
        std::vector<ADCGD> ay(2);

        ay[0] = cos(ax[0]) * ap[0];
        ay[1] = ax[1] * ax[2] + sin(ax[0]) + ap[1];

        return ay;
    }

};

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;
using namespace std;

TEST_F(CppADCGDynamicTest1, DynamicCustomElements) {
    // use a special object for source code generation
    using CGD = CG<double>;
    using ADCG = AD<CGD>;

    // independent variables
    std::vector<ADCG> ax(3);
    ax[0] = 1;
    ax[1] = 1;
    ax[2] = 1;

    std::vector<ADCG> ap(2);
    ap[0] = 2;
    ap[1] = 3;

    std::vector<double> x(ax.size());
    x[0] = 1;
    x[1] = 2;
    x[2] = 1;

    std::vector<double> p(ap.size());
    p[0] = 3;
    p[1] = 4;

    this->_forwardOne = false;
    this->_reverseOne = false;
    this->_reverseTwo = false;
    this->testDynamicFull(ax, ap, x, p);
}
