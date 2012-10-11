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

bool Pendulum2D() {
    using namespace CppAD;
    using namespace std;
    typedef CppAD::CG<double> CGD;
    typedef CppAD::AD<CGD> ADCG;

    std::vector<ADCG> U(5);
    Independent(U);

    ADCG x = U[0];
    ADCG y = U[1];
    ADCG vx = U[2]; // auxiliary variable
    ADCG vy = U[3]; // auxiliary variable
    ADCG T = U[4]; // tension

    double g = 9.80665; // gravity constant
    double L = 1.0; // fixed length

    // dependent variable vector 
    std::vector<ADCG> Z(5);
    Z[0] = vx; // dx/dt =
    Z[1] = vy; // dy/dt =
    Z[2] = T * x; // dvx/dt =
    Z[3] = T * y - g; // dvy/dt =
    Z[4] = x * x + y * y - L*L;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = new ADFun<CGD > (U, Z);

    std::vector<bool> eqDifferentialInfo(5, true);
    eqDifferentialInfo[4] = false;
    std::vector<bool> varInfo(5, true);

    Plantelides<double> pantelides(fun, eqDifferentialInfo, varInfo);

    pantelides.reduceIndex();

    delete fun;

    return false;

}

bool Pendulum3D() {
    using namespace CppAD;
    using namespace std;
    typedef CppAD::CG<double> CGD;
    typedef CppAD::AD<CGD> ADCG;

    std::vector<ADCG> U(7);
    Independent(U);

    ADCG x = U[0];
    ADCG y = U[1];
    ADCG z = U[2];
    ADCG vx = U[3]; // auxiliary variable
    ADCG vy = U[4]; // auxiliary variable
    ADCG vz = U[5]; // auxiliary variable
    ADCG T = U[6]; // tension

    double g = 9.80665; // gravity constant
    double L = 1.0; // fixed length

    // dependent variable vector 
    std::vector<ADCG> Z(7);
    Z[0] = vx; // dx/dt =
    Z[1] = vy; // dy/dt =
    Z[2] = vz; // dz/dt =
    Z[3] = T * x; // dvx/dt =
    Z[4] = T * y - g; // dvy/dt =
    Z[5] = T * z; // dvz/dt =
    Z[6] = x * x + y * y + z * z - L*L;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = new ADFun<CGD > (U, Z);

    std::vector<bool> eqDifferentialInfo(7, true);
    eqDifferentialInfo[6] = false;
    std::vector<bool> varInfo(7, true);

    Plantelides<double> pantelides(fun, eqDifferentialInfo, varInfo);

    pantelides.reduceIndex();

    delete fun;

    return false;

}

bool Pantelides() {
    bool ok = true;
    ok &= Pendulum2D();
    ok &= Pendulum3D();
    return ok;
}
