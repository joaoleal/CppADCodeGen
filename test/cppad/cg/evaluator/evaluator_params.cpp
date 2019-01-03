/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
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
#include "CppADCGEvaluatorTest.hpp"

using namespace CppAD;
using namespace CppAD::cg;

TEST_F(CppADCGEvaluatorTest, Params) {
    ModelType model = [](const std::vector<CGD>& x, const std::vector<CGD>& p) {
        // dependent variable vector
        std::vector<CGD> y(3);

        // model
        y[0] = p[0] * x[0] + x[1]; // CGD + CGD
        y[1] = y[0] + p[1]; // CGD + double
        y[2] = 1. + y[1]; // double + CGD
        return y;
    };

    this->test(model, std::vector<double>{0.5, 1.5}, std::vector<double>{2.5, 3.0});
}