/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_cgoo/dae_index_reduction/cg_pantelides.hpp>

#include "gcc_load_dynamic.hpp"
#include "pendulum.hpp"

inline bool PantelidesPendulum2D() {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum2D<CGD > ();

    std::vector<DaeVarInfo> daeVar(10);
    daeVar[5].makeTimeIndependent();
    daeVar[6] = 0;
    daeVar[7] = 1;
    daeVar[8] = 2;
    daeVar[9] = 3;

    Plantelides<double> pantelides(fun, daeVar);

    std::vector<DaeVarInfo> newDaeVar;
    ADFun<CGD>* reducedFun = pantelides.reduceIndex(newDaeVar);

    delete fun;
    delete reducedFun;

    return reducedFun != NULL;
}

inline bool PantelidesPendulum3D() {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum3D<CGD > ();

    std::vector<DaeVarInfo> daeVar(13);
    daeVar[7] = 0;
    daeVar[8] = 1;
    daeVar[9] = 2;
    daeVar[10] = 3;
    daeVar[11] = 4;
    daeVar[12] = 5;

    Plantelides<double> pantelides(fun, daeVar);

    std::vector<DaeVarInfo> newDaeVar;
    ADFun<CGD>* reducedFun = pantelides.reduceIndex(newDaeVar);

    delete fun;
    delete reducedFun;

    return reducedFun != NULL;
}

bool Pantelides() {
    bool ok = true;
    ok &= PantelidesPendulum2D();
    ok &= PantelidesPendulum3D();
    return ok;
}
