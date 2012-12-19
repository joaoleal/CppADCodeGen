/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppadcg/cg.hpp>
#include "gcc_load_dynamic.hpp"
#include "add.hpp"

bool Add() {
    using namespace CppAD;
    using namespace std;

    std::vector<double> u(2);
    size_t s = 0;
    size_t t = 1;
    u[s] = 3.;
    u[t] = 2.;

    // create f: U -> Z and vectors used for derivative calculations   
    bool ok = test0nJac("add", &AddFunc<double >, &AddFunc<CG<double> >, u);

    return ok;
}
