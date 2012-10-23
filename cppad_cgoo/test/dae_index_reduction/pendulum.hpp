#ifndef CPPADCGOO_TEST_PENDULUM_INCLUDED
#define	CPPADCGOO_TEST_PENDULUM_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

template<class Base>
inline CppAD::ADFun<Base>* Pendulum2D() {
    using namespace CppAD;
    using namespace std;
    typedef CppAD::AD<Base> ADB;

    std::vector<ADB> U(10);
    Independent(U);

    ADB x = U[0];
    ADB y = U[1];
    ADB vx = U[2]; // auxiliary variable
    ADB vy = U[3]; // auxiliary variable
    ADB T = U[4]; // tension
    ADB L = U[5]; // length  (parameter)

    ADB dxdt = U[6];
    ADB dydt = U[7];
    ADB dvxdt = U[8]; // auxiliary variable
    ADB dvydt = U[9]; // auxiliary variable

    double g = 9.80665; // gravity constant

    // dependent variable vector 
    std::vector<ADB> Z(5);
    Z[0] = dxdt - vx; // dx/dt =
    Z[1] = dydt - vy; // dy/dt =
    Z[2] = dvxdt - T * x; // dvx/dt =
    Z[3] = dvydt - (T * y - g); // dvy/dt =
    Z[4] = x * x + y * y - L*L;

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<Base > (U, Z);
}

template<class Base>
inline CppAD::ADFun<Base>* Pendulum3D() {
    using namespace CppAD;
    using namespace std;
    typedef CppAD::AD<Base> ADB;

    std::vector<ADB> U(13);
    Independent(U);

    ADB x = U[0];
    ADB y = U[1];
    ADB z = U[2];
    ADB vx = U[3]; // auxiliary variable
    ADB vy = U[4]; // auxiliary variable
    ADB vz = U[5]; // auxiliary variable
    ADB T = U[6]; // tension
    ADB dxdt = U[7];
    ADB dydt = U[8];
    ADB dzdt = U[9];
    ADB dvxdt = U[10]; // auxiliary variable
    ADB dvydt = U[11]; // auxiliary variable
    ADB dvzdt = U[12]; // auxiliary variable

    double g = 9.80665; // gravity constant
    double L = 1.0; // fixed length

    // dependent variable vector 
    std::vector<ADB> Z(7);
    Z[0] = dxdt - vx; // dx/dt =
    Z[1] = dydt - vy; // dy/dt =
    Z[2] = dzdt - vz; // dz/dt =
    Z[3] = dvxdt - T * x; // dvx/dt =
    Z[4] = dvydt - (T * y - g); // dvy/dt =
    Z[5] = dvzdt - T * z; // dvz/dt =
    Z[6] = x * x + y * y + z * z - L*L;

    // create f: U -> Z and vectors used for derivative calculations
    return new ADFun<Base > (U, Z);
}

#endif
