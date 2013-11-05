#ifndef CPPAD_CG_TEST_PLUG_FLOW_INCLUDED
#define CPPAD_CG_TEST_PLUG_FLOW_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

#include <assert.h>

template<class T>
std::vector<CppAD::AD<T> > plugFlowFunc(const std::vector<CppAD::AD<T> >& x) {
    using namespace CppAD;

    size_t j;
    
    // dependent variable vector 
    std::vector< AD<T> > y(6 * 4);

    // temporary variables
    std::vector< AD<T> > v(5);

    v[0] = 1.66666666666667e-05 * x[24];
    v[1] = 0.5 * 3.14159265358979 * x[25] * x[25];
    for (j = 0; j < 6; j++) {
        v[2] = 1000. * x[j];
        v[3] = 1000. * x[j + 6];
        v[4] = exp(35 - 83809.879272 / (8.31447215 * (x[j + 18] + 273.15))) * v[2] * v[3];
        y[j] = 0.001 * (1000. * x[(j < 2) ? j * -26 + 26 : j + -1] * v[0] - v[2] * v[0] + -1 * v[4] * v[1]) / v[1];
        y[j + 6] = 0.001 * (1000. * x[(j < 2) ? j * -21 + 27 : j + 5] * v[0] - v[3] * v[0] + -1 * v[4] * v[1]) / v[1];
        y[j + 12] = 0.001 * (1000. * x[(j < 2) ? j * -16 + 28 : j + 11] * v[0] - 1000. * x[j + 12] * v[0] + v[4] * v[1]) / v[1];
        y[j + 18] = 0;
    }

    return y;
}

#endif
