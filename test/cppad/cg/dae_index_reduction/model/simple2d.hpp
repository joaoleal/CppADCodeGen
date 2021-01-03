#ifndef CPPAD_CG_TEST_SIMPLE2D_INCLUDED
#define CPPAD_CG_TEST_SIMPLE2D_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
 *    Copyright (C) 2016 Ciengis
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

namespace CppAD {
namespace cg {

template<class Base>
inline CppAD::ADFun<Base> Simple2D(std::vector<DaeVarInfo>& daeVar) {
    using namespace CppAD;
    using namespace std;
    using ADB = CppAD::AD<Base>;

    std::vector<ADB> U(5);
    CppAD::Independent(U);

    ADB x1 = U[0];
    ADB x2 = U[1];
    ADB t = U[2]; // time (not really used)

    ADB dx1dt = U[3];
    ADB dx2dt = U[4];

    daeVar.resize(U.size());
    daeVar[0] = DaeVarInfo("x1");
    daeVar[1] = DaeVarInfo("x2");
    daeVar[2].makeIntegratedVariable();
    daeVar[3] = DaeVarInfo(0);
    daeVar[4] = DaeVarInfo(1);

    // dependent variable vector
    std::vector<ADB> Z(2);
    Z[0] = dx1dt + dx2dt;
    Z[1] = x2 - 5;

    // create f: U -> Z and vectors used for derivative calculations
    return ADFun<Base> (U, Z);
}

template<class Base>
inline CppAD::ADFun<Base> Simple2DParam(std::vector<DaeVarInfo>& daeVar,
                                        const std::vector<Base>& p) {
    using namespace CppAD;
    using namespace std;
    using ADB = CppAD::AD<Base>;

    std::vector<ADB> U(5);
    std::vector<ADB> par(1);

    // use a special object for source code generation
    // declare independent variables, dynamic parameters, starting recording
    size_t abort_op_index = 0;
    bool record_compare = true;
    CppAD::Independent(U, abort_op_index, record_compare, par);

    ADB x1 = U[0];
    ADB x2 = U[1];
    ADB t = U[2]; // time (not really used)

    ADB dx1dt = U[3];
    ADB dx2dt = U[4];

    daeVar.resize(U.size());
    daeVar[0] = DaeVarInfo("x1");
    daeVar[1] = DaeVarInfo("x2");
    daeVar[2].makeIntegratedVariable();
    daeVar[3] = DaeVarInfo(0);
    daeVar[4] = DaeVarInfo(1);

    // dependent variable vector
    std::vector<ADB> Z(2);
    Z[0] = dx1dt + dx2dt;
    Z[1] = x2 - par[0];

    // create f: U -> Z and vectors used for derivative calculations
    return ADFun<Base> (U, Z);
}

} // END cg namespace
} // END CppAD namespace

#endif