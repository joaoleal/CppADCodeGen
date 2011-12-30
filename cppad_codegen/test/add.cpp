/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_codegen/cppad_codegen.hpp>

#include "gcc_load_dynamic.hpp"

bool Add() {
    using namespace CppAD;
    using namespace std;

    bool ok = true;

    // independent variable vector, indices, values, and declaration
    std::vector<double> u(2);
    size_t s = 0;
    size_t t = 1;
    u[s] = 3.;
    u[t] = 2.;

    CPPAD_TEST_VECTOR< AD<double> > U(2);
    U[s] = u[s];
    U[t] = u[t];
    Independent(U);

    // dependent variable vector and indices
    CPPAD_TEST_VECTOR< AD<double> > Z(3);
    size_t x = 0;
    size_t y = 1;
    size_t z = 2;

    // dependent variable values
    Z[x] = U[s] + U[t]; // AD<double> + AD<double>
    Z[y] = Z[x] + 1.; // AD<double> + double
    Z[z] = 1. + Z[y]; // double + AD<double> 

    // create f: U -> Z and vectors used for derivative calculations
    ADFunCodeGen<double> f(U, Z);


    ok &= test0nJac("add", f, u, Z);
    
    return ok;
}
