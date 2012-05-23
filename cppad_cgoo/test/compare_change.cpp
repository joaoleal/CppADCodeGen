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
#include "compare_change.hpp"

using namespace CppAD;
using namespace std;



// ------------------------------- < ----------------------------

bool CompareChange1(const std::vector<std::vector<double> >& xV) {
    bool ok = true;

    int comparisons;
    ok &= test0("CompareChange<",
                &CompareChangeFunc1<double >,
                &CompareChangeFunc1<CG<double> >,
                xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed
    return ok;
}

// ------------------------------- > ----------------------------

bool CompareChange2(const std::vector<std::vector<double> >& xV) {
    bool ok = true;
    int comparisons;
    ok &= test0("CompareChange>",
                &CompareChangeFunc2<double >,
                &CompareChangeFunc2<CG<double> >,
                xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed
    return ok;
}

// ------------------------------- <= ----------------------------

bool CompareChange3(const std::vector<std::vector<double> >& xV) {
    bool ok = true;
    int comparisons;
    ok &= test0("CompareChange<=",
                &CompareChangeFunc3<double >,
                &CompareChangeFunc3<CG<double> >,
                xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed
    return ok;
}

// ------------------------------- >= ----------------------------

bool CompareChange4(const std::vector<std::vector<double> >& xV) {
    bool ok = true;
    int comparisons;
    ok &= test0("CompareChange>=",
                &CompareChangeFunc4<double >,
                &CompareChangeFunc4<CG<double> >,
                xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed
    return ok;
}

// ------------------------------- == ----------------------------

bool CompareChange5(const std::vector<std::vector<double> >& xV) {
    bool ok = true;
    int comparisons;
    ok &= test0("CompareChange>=",
                &CompareChangeFunc4<double >,
                &CompareChangeFunc4<CG<double> >,
                xV, comparisons);

    ok &= (comparisons == 4); // the first two comparisons do not change
    return ok;
}

// ------------------------------- != ----------------------------

bool CompareChange6(const std::vector<std::vector<double> >& xV) {
    bool ok = true;
    int comparisons;
    ok &= test0("CompareChange>=",
                &CompareChangeFunc4<double >,
                &CompareChangeFunc4<CG<double> >,
                xV, comparisons);

    ok &= (comparisons == 4); // the first two comparisons do not change
    return ok;
}

bool CompareChange() {
    std::vector<double> X(2);
    std::vector<std::vector<double> > xV;
    // create independent variables
    X[0] = 3.;
    X[1] = 4.;
    xV.push_back(X);
    X[0] = 4.;
    X[1] = 3.;
    xV.push_back(X);

    bool ok = true;

    ok &= !CompareChange1(xV);
    ok &= !CompareChange2(xV);
    ok &= !CompareChange3(xV);
    ok &= !CompareChange4(xV);
    ok &= !CompareChange5(xV);
    ok &= !CompareChange6(xV);

    return ok;
}