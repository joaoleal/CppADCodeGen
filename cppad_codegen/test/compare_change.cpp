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

bool CompareChange() {
    using namespace CppAD;

    bool ok = true;

    // ------------------------------- < ----------------------------

    // create independent variables
    CPPAD_TEST_VECTOR< AD<double> > X(2);
    X[0] = 3.;
    X[1] = 4.;
    Independent(X);

    // create dependent variables
    CPPAD_TEST_VECTOR< AD<double> > Y(6);

    // CondExp would never require retaping 
    if (X[0] < X[1]) // True variable < variable
        Y[0] = X[0];
    else Y[0] = X[1];
    if (X[1] < X[0]) // False variable < variable
        Y[1] = X[0];
    else Y[1] = X[1];
    if (3.5 < X[1]) // True parameter < variable
        Y[2] = X[0];
    else Y[2] = X[1];
    if (3.5 < X[0]) // False parameter < variable
        Y[3] = X[0];
    else Y[3] = X[1];
    if (X[0] < 4.) // True variable < parameter
        Y[4] = X[0];
    else Y[4] = X[1];
    if (X[1] < 4.) // False variable < parameter
        Y[5] = X[0];
    else Y[5] = X[1];

    // f : X -> Y
    ADFunCodeGen<double> *f = new ADFunCodeGen<double>(X, Y);

    // new argument value
    std::vector<double> x(X.size());
    x[0] = 4.;
    x[1] = 3.;

    std::vector<double> depCGen(Y.size());
    int comparisons;

    std::vector<std::vector<double> > xV;
    xV.push_back(x);

    ok &= run0(*f, "./tmp/test_CompareChange1.so", "CompareChange1", xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed

    // done with this function
    delete f;

    // ------------------------------- > ----------------------------
    // create independent variables
    Independent(X);

    if (X[0] > X[1]) // False variable > variable
        Y[0] = X[0];
    else Y[0] = X[1];
    if (X[1] > X[0]) // True variable > variable
        Y[1] = X[0];
    else Y[1] = X[1];
    if (3.5 > X[1]) // False parameter > variable
        Y[2] = X[0];
    else Y[2] = X[1];
    if (3.5 > X[0]) // True parameter > variable
        Y[3] = X[0];
    else Y[3] = X[1];
    if (X[0] > 3.) // False variable > parameter
        Y[4] = X[0];
    else Y[4] = X[1];
    if (X[1] > 3.) // True variable > parameter
        Y[5] = X[0];
    else Y[5] = X[1];

    // f : X -> Y
    f = new ADFunCodeGen<double> (X, Y);

    ok &= run0(*f, "./tmp/test_CompareChange2.so", "CompareChange2", xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed

    // done with this function
    delete f;

    // ------------------------------- <= ----------------------------
    // create independent variables
    Independent(X);

    if (X[0] <= X[1]) // True variable <= variable
        Y[0] = X[0];
    else Y[0] = X[1];
    if (X[1] <= X[0]) // False variable <= variable
        Y[1] = X[0];
    else Y[1] = X[1];
    if (4. <= X[1]) // True parameter <= variable
        Y[2] = X[0];
    else Y[2] = X[1];
    if (4. <= X[0]) // False parameter <= variable
        Y[3] = X[0];
    else Y[3] = X[1];
    if (X[0] <= 3.5) // True variable <= parameter
        Y[4] = X[0];
    else Y[4] = X[1];
    if (X[1] <= 3.5) // False variable <= parameter
        Y[5] = X[0];
    else Y[5] = X[1];

    // f : X -> Y
    f = new ADFunCodeGen<double> (X, Y);

    ok &= run0(*f, "./tmp/test_CompareChange3.so", "CompareChange3", xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed

    // done with this function
    delete f;


    // ------------------------------- >= ----------------------------
    // create independent variables
    Independent(X);

    if (X[0] >= X[1]) // False variable >= variable
        Y[0] = X[0];
    else Y[0] = X[1];
    if (X[1] >= X[0]) // True variable >= variable
        Y[1] = X[0];
    else Y[1] = X[1];
    if (3.5 >= X[1]) // False parameter >= variable
        Y[2] = X[0];
    else Y[2] = X[1];
    if (3.5 >= X[0]) // True parameter >= variable
        Y[3] = X[0];
    else Y[3] = X[1];
    if (X[0] >= 4.) // False variable >= parameter
        Y[4] = X[0];
    else Y[4] = X[1];
    if (X[1] >= 4.) // True variable >= parameter
        Y[5] = X[0];
    else Y[5] = X[1];

    // f : X -> Y
    f = new ADFunCodeGen<double> (X, Y);

    ok &= run0(*f, "./tmp/test_CompareChange4.so", "CompareChange4", xV, comparisons);

    ok &= (comparisons == 6); // all comparisons have changed

    // done with this function
    delete f;

    // ------------------------------- == ----------------------------
    // create independent variables
    Independent(X);

    if (X[0] == X[1]) // False variable == variable
        Y[0] = X[0];
    else Y[0] = X[1];
    if (X[0] == X[0]) // True variable == variable
        Y[1] = X[0];
    else Y[1] = X[1];
    if (3. == X[1]) // False parameter == variable
        Y[2] = X[0];
    else Y[2] = X[1];
    if (3. == X[0]) // True parameter == variable
        Y[3] = X[0];
    else Y[3] = X[1];
    if (X[0] == 4.) // False variable == parameter
        Y[4] = X[0];
    else Y[4] = X[1];
    if (X[1] == 4.) // True variable == parameter
        Y[5] = X[0];
    else Y[5] = X[1];


    // f : X -> Y
    f = new ADFunCodeGen<double> (X, Y);

    ok &= run0(*f, "./tmp/test_CompareChange4.so", "CompareChange4", xV, comparisons);

    // the first two comparisons do not change
    ok &= (comparisons == 4);

    // done with this function
    delete f;

    // ------------------------------- != ----------------------------
    // create independent variables
    Independent(X);

    if (X[0] != X[1]) // True variable != variable
        Y[0] = X[0];
    else Y[0] = X[1];
    if (X[0] != X[0]) // False variable != variable
        Y[1] = X[0];
    else Y[1] = X[1];
    if (3. != X[1]) // True parameter != variable
        Y[2] = X[0];
    else Y[2] = X[1];
    if (3. != X[0]) // False parameter != variable
        Y[3] = X[0];
    else Y[3] = X[1];
    if (X[0] != 4.) // True variable != parameter
        Y[4] = X[0];
    else Y[4] = X[1];
    if (X[1] != 4.) // False variable != parameter
        Y[5] = X[0];
    else Y[5] = X[1];

    // f : X -> Y
    f = new ADFunCodeGen<double> (X, Y);

    ok &= run0(*f, "./tmp/test_CompareChange5.so", "CompareChange5", xV, comparisons);

    // the first two comparisons do not change
    ok &= (comparisons == 4);

    // done with this function
    delete f;

    return ok;
}
