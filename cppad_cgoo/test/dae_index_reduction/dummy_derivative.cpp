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

    std::vector<double> x(10);
    std::vector<double> normVar(10, 1.0);
    std::vector<double> normEq(5, 1.0);

    x[0] = -1.0; // x
    x[1] = 0.0; // y
    x[2] = 0.0; // vx
    x[3] = 0.0; // vy
    x[4] = 1.0; // Tension
    x[5] = 1.0; // length

    x[6] = 0.0; // dxdt
    x[7] = 0.0; // dydt
    x[8] = -1.0; // dvxdt
    x[9] = 9.80665; // dvydt

    std::vector<int> derivative(10, -1);
    derivative[0] = 6;
    derivative[1] = 7;
    derivative[2] = 8;
    derivative[3] = 9;

    std::vector<bool> timeDependent(10, true);
    timeDependent[5] = false;

    DummyDerivatives<double> dummyD(fun, derivative, timeDependent, x, normVar, normEq);

    dummyD.reduceIndex();
    
    dummyD.reduceEquations();

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

    std::vector<double> x(13);
    std::vector<double> normVar(13, 1.0);
    std::vector<double> normEq(7, 1.0);

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
