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

namespace { // BEGIN empty namespace

    bool MulTestOne(void) {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector, indices, values, and declaration
        std::vector<double> uu(2);
        size_t s = 0;
        size_t t = 1;
        uu[s] = 3.;
        uu[t] = 2.;

        CPPAD_TEST_VECTOR< AD<double> > U(2);
        U[s] = uu[s];
        U[t] = uu[t];
        Independent(U);

        // assign some parameters
        AD<double> zero = 0.;
        AD<double> one = 1.;

        // dependent variable vector and indices
        CPPAD_TEST_VECTOR< AD<double> > Z(5);
        size_t x = 0;
        size_t y = 1;
        size_t z = 2;
        size_t u = 3;
        size_t v = 4;

        // assign the dependent variables
        Z[x] = U[s] * U[t]; // AD<double> * AD<double>
        Z[y] = Z[x] * 4.; // AD<double> *    double
        Z[z] = 4. * Z[y]; //    double  * AD<double> 
        Z[u] = one * Z[z]; // multiplication by parameter equal to one
        Z[v] = zero * Z[z]; // multiplication by parameter equal to zero

        // check multipilcation by zero results in a parameter
        ok &= Parameter(Z[v]);

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("MulTestOne", f, uu, Z);

        return ok;
    }

    bool MulTestTwo(void) {
        using namespace CppAD;

        bool ok = true;

        // independent variable vector
        std::vector<double> u(1);
        u[0] = .5;
        CPPAD_TEST_VECTOR< AD<double> > U(1);
        U[0] = u[0];
        Independent(U);

        AD<double> a = U[0] * 1.; // AD<double> * double
        AD<double> b = a * 2; // AD<double> * int
        AD<double> c = 3. * b; // double     * AD<double> 
        AD<double> d = 4 * c; // int        * AD<double> 

        // dependent variable vector 
        CPPAD_TEST_VECTOR< AD<double> > Z(1);
        Z[0] = U[0] * d; // AD<double> * AD<double>

        // create f: U -> Z and vectors used for derivative calculations
        ADFunCodeGen<double> f(U, Z);

        ok &= test0nJac("MulTestTwo", f, u, Z);

        return ok;
    }

} // END empty namespace

bool Mul(void) {
    bool ok = true;
    ok &= MulTestOne();
    ok &= MulTestTwo();
    return ok;
}
