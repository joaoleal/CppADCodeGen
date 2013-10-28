#ifndef CPPAD_CG_TEST_TANK_BATTERY_INCLUDED
#define CPPAD_CG_TEST_TANK_BATTERY_INCLUDED
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
std::vector<CppAD::AD<T> > tankBatteryFunc(const std::vector<CppAD::AD<T> >& x) {
    using namespace CppAD;

    // dependent variable vector 
    std::vector< AD<T> > dhdt(6);

    // temporary variables
    std::vector< AD<T> > v(3);

    v[0] = 3.1415926535897931 * x[7] * x[7];
    v[1] = (1001.9158506700721 * 0.00050000000000000001 * sqrt(2. * 9.8066499999999994 * (x[0] / 1001.9158506700721) / v[0])) / 0.018015283300000001;
    dhdt[0] = 0.018015283300000001 * ((1001.9158506700721 * 0.016666666666666666 * x[6]) / 0.018015283300000001 - v[1]);

    v[2] = (1001.9158506700721 * 0.00050000000000000001 * sqrt(2. * 9.8066499999999994 * (x[1] / 1001.9158506700721) / v[0])) / 0.018015283300000001;
    dhdt[1] = 0.018015283300000001 * (v[1] - v[2]);

    v[1] = (1001.9158506700721 * 0.00050000000000000001 * sqrt(2. * 9.8066499999999994 * (x[2] / 1001.9158506700721) / v[0])) / 0.018015283300000001;
    dhdt[2] = 0.018015283300000001 * (v[2] - v[1]);

    v[2] = (1001.9158506700721 * 0.00050000000000000001 * sqrt(2. * 9.8066499999999994 * (x[3] / 1001.9158506700721) / v[0])) / 0.018015283300000001;
    dhdt[3] = 0.018015283300000001 * (v[1] - v[2]);

    v[1] = (1001.9158506700721 * 0.00050000000000000001 * sqrt(2. * 9.8066499999999994 * (x[4] / 1001.9158506700721) / v[0])) / 0.018015283300000001;
    dhdt[4] = 0.018015283300000001 * (v[2] - v[1]);

    dhdt[5] = 0.018015283300000001 * (v[1] - (1001.9158506700721 * 0.00050000000000000001 * sqrt(2. * 9.8066499999999994 * (x[5] / 1001.9158506700721) / v[0])) / 0.018015283300000001);


    return dhdt;
}

#endif
