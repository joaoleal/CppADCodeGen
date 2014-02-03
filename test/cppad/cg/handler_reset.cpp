/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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

#include <cppad/cg/cppadcg.hpp>

#include "gcc_load_dynamic.hpp"

namespace {

bool HandlerReset1() {
    using namespace CppAD;
    using namespace std;

    CodeHandler<double> handler(20);

    // independent variables of CppAD
    vector<CG<double> > u(3);
    handler.makeVariables(u);

    // independent variables of CppAD
    vector<AD<CG<double> > > U(3);
    U[0] = u[0];
    U[1] = u[1];
    U[2] = u[2];

    CppAD::Independent(U);

    // dependent variable of CppAD
    vector<AD<CG<double> > > w(3);
    w[0] = u[0] + 2.0;
    w[1] = u[1] + 3.0;
    w[2] = u[2] * 4.0;

    CppAD::ADFun<CG<double> > f(U, w);

    vector<CG<double> > dep = f.Forward(0, u);

    ostringstream code;
    handler.generateCode(code, dep);
    string code1 = code.str();

    code.str("");
    handler.generateCode(code, dep);
    string code2 = code.str();

    if (code1 == code2) {
        return true;
    } else {
        cerr << "\n ############################\n\n" << code1
                << "\n ############################\n\n" << code2
                << "\n ############################\n\n";
        return false;
    }

}

}

bool HandlerReset() {
    bool ok = true;
    ok &= HandlerReset1();
    return ok;
}
