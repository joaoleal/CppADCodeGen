/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad/cppad_code_gen.hpp>

#include "gcc_load_dynamic.hpp"

bool Parameter(void) {
    bool ok = true;
    using namespace CppAD;

    // number of different parameter values
    size_t n_parameter = 7;

    // number of parameter repeats
    size_t n_repeat = 5;

    // independent variable vector
    size_t j, n = n_parameter * n_repeat;
    std::vector<double> u(n);
    CPPAD_TEST_VECTOR< AD<double> > ax(n);
    for (j = 0; j < n; j++) {
        u[j] = double(j);
        ax[j] = u[j];
    }
    Independent(ax);

    // dependent variable vector and indices
    size_t i, m = n;
    CPPAD_TEST_VECTOR< AD<double> > ay(m);
    for (i = 0; i < m; i++) { // must avoid Float(k) = 0 because it would get optimized out	
        size_t k = (i % n_parameter);
        k = k * k * 10 + 1;
        j = i;
        ay[i] = ax[j] + double(k);
    }

    // create f: ax -> ay 
    ADFunCodeGen<double> f(ax, ay);

    ok = f.size_par() == n_parameter;

    ok &= test0nJac("parameter", f, u, ay);

    return ok;
}

