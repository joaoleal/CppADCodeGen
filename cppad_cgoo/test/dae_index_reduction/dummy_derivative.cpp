/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_cgoo/cg.hpp>

#include "gcc_load_dynamic.hpp"
#include "pendulum.hpp"

inline bool DummyDerivPendulum2D() {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum2D<CGD > ();

    std::vector<bool> eqDifferentialInfo(5, true);
    eqDifferentialInfo[4] = false;
    std::vector<bool> varInfo(6, true);
    varInfo[5] = false;

    std::vector<double> x(6);
    std::vector<double> norm(6, 1.0);

    x[0] = -1.0; // x
    x[1] = 0.0; // y
    x[2] = 0.0; // vx
    x[3] = 0.0; // vy
    x[4] = 1.0; // Tension
    x[5] = 1.0; // length

    DummyDerivatives<double> dummyD(fun, eqDifferentialInfo, varInfo, x, norm);

    dummyD.reduceIndex();

    delete fun;

    return false;
}

inline bool DummyDerivPendulum3D() {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;
    typedef AD<CGD> ADCG;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum3D<CGD > ();

    std::vector<bool> eqDifferentialInfo(7, true);
    eqDifferentialInfo[6] = false;
    std::vector<bool> varInfo(7, true);

    //DummyDerivatives<double> dummyD(fun, eqDifferentialInfo, varInfo);

    //dummyD.reduceIndex();

    delete fun;

    return false;
}

bool DummyDeriv() {
    bool ok = true;
    ok &= DummyDerivPendulum2D();
    ok &= DummyDerivPendulum3D();
    return ok;
}
