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

namespace {

    bool tan_case(bool tan_first) {
        using CppAD::AD;
        using CppAD::NearEqual;

        bool ok = true;
        double eps = 100. * std::numeric_limits<double>::epsilon();

        // independent variable vector, indices, values, and declaration
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = .7;

        CPPAD_TEST_VECTOR< AD<double> > ax(n);
        ax[0] = u[0];
        Independent(ax);

        // dependent variable vector and indices
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > ay(m);
        if (tan_first)
            ay[0] = atan(tan(ax[0]));
        else ay[0] = tan(atan(ax[0]));

        // create f: x -> y and vectors used for derivative calculations
        CppAD::ADFunCodeGen<double> f(ax, ay);

        ok &= test0nJac(std::string("tan_case_") + (tan_first ? "true" : "false"), f, u, ay, eps, eps);

        return ok;
    }

    bool tanh_case(bool tanh_first) {
        using CppAD::AD;
        using CppAD::NearEqual;

        bool ok = true;
        double eps = 100. * std::numeric_limits<double>::epsilon();

        // independent variable vector, indices, values, and declaration
        size_t n = 1;
        std::vector<double> u(n);
        u[0] = .5;

        CPPAD_TEST_VECTOR< AD<double> > ax(n);
        ax[0] = u[0];
        Independent(ax);

        // dependent variable vector and indices
        size_t m = 1;
        CPPAD_TEST_VECTOR< AD<double> > ay(m);
        AD<double> z;
        if (tanh_first) {
            z = tanh(ax[0]);
            ay[0] = .5 * log((1. + z) / (1. - z));
        } else {
            z = .5 * log((1. + ax[0]) / (1. - ax[0]));
            ay[0] = tanh(z);
        }

        // create f: x -> y and vectors used for derivative calculations
        CppAD::ADFunCodeGen<double> f(ax, ay);

        ok &= test0nJac(std::string("tanh_case_") + (tanh_first ? "true" : "false"), f, u, ay, eps, eps);

        return ok;
    }
}

bool Tan() {
    bool ok = true;
    ok &= tan_case(true);
    ok &= tan_case(false);
    ok &= tanh_case(true);
    ok &= tanh_case(false);
    return ok;
}
