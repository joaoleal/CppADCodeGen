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

inline bool PantelidesPendulum2D() {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum2D<CGD > ();

    std::vector<int> derivative(10, -1);
    derivative[0] = 6;
    derivative[1] = 7;
    derivative[2] = 8;
    derivative[3] = 9;

    std::vector<bool> timeDependent(10, true);
    timeDependent[5] = false;

    Plantelides<double> pantelides(fun, derivative, timeDependent);

    pantelides.reduceIndex();

    delete fun;

    return false;
}

inline bool PantelidesPendulum3D() {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum3D<CGD > ();

    std::vector<int> derivative(13, -1);
    derivative[0] = 7;
    derivative[1] = 8;
    derivative[2] = 9;
    derivative[3] = 10;
    derivative[4] = 11;
    derivative[5] = 12;

    std::vector<bool> timeDependent(13, true);
    
    Plantelides<double> pantelides(fun, derivative, timeDependent);

    pantelides.reduceIndex();

    delete fun;

    return false;
}

bool Pantelides() {
    bool ok = true;
    ok &= PantelidesPendulum2D();
    ok &= PantelidesPendulum3D();
    return ok;
}
