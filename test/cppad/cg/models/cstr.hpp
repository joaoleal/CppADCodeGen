#ifndef CPPAD_CG_TEST_CSTR_INCLUDED
#define CPPAD_CG_TEST_CSTR_INCLUDED
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
std::vector<CppAD::AD<T> > CstrFunc(const std::vector<CppAD::AD<T> >& ind,
                                    const std::vector<CppAD::AD<T> >& par) {
    using namespace CppAD;

    // dependent variable vector 
    std::vector< AD<T> > dxdt(4);

    // temporary variables
    std::vector< AD<T> > v(17);

    v[0] = par[2] * par[4];
    v[1] = par[3] * par[5];
    v[2] = 0.1768325928384763 * ind[0];
    v[3] = 0.0180152833 * 1000 * (-13.851 + 0.64038 * ind[2] + -0.00191 * ind[2] * ind[2] + 1.8211e-06 * ind[2] * ind[2] * ind[2]);
    v[4] = v[2] * v[3];
    v[5] = ((v[4] - 1.0e-15 * ind[1] * v[2]) / 0.0180152833) / v[2];
    v[6] = (v[5] * v[2]) / (ind[1] * v[2] + v[5] * v[2]);
    v[7] = 1 - v[6];
    v[8] = (0.0180152833 * 1000 * (-13.851 + 0.64038 * ind[2] + -0.00191 * ind[2] * ind[2] + 1.8211e-06 * ind[2] * ind[2] * ind[2]) * ind[5]) / (1.0e-15 * v[7] + 0.0180152833 * v[6]);
    v[9] = exp(log(par[0]) - (8.31447215 * par[1]) / (8.31447215 * ind[2]));
    v[7] = v[0] + v[1] - v[7] * v[8] + -1 * v[9] * v[2];
    v[10] = 0.0180152833 * 1000 * (-13.851 + 0.64038 * par[6] + -0.00191 * par[6] * par[6] + 1.8211e-06 * par[6] * par[6] * par[6]) * par[4];
    v[11] = 0.0180152833 * 1000 * (-13.851 + 0.64038 * par[7] + -0.00191 * par[7] * par[7] + 1.8211e-06 * par[7] * par[7] * par[7]) * par[5];
    v[12] = par[6] * par[6];
    v[13] = ind[2] * ind[2];
    v[14] = par[7] * par[7];
    v[15] = ind[2] * ind[2];
    v[16] = par[16] * (0.1768325928384763 + 1.4906857141283567 * ind[0]) * (ind[2] - ind[3]);
    dxdt[2] = (1.0e-15 * v[0] * par[9] * (par[6] - ind[2]) + v[10] * (15340.863387921299 * (par[6] - ind[2]) + -58.009079435348092 * (par[6] * par[6] - ind[2] * ind[2]) + 0.1503353174209219 * (par[6] * par[6] * par[6] - ind[2] * ind[2] * ind[2]) + -0.00019588923145049848 * (par[6] * par[6] * par[6] * par[6] - ind[2] * ind[2] * ind[2] * ind[2]) + 1.0402389841962685e-07 * (v[12] * v[12] * par[6] - v[13] * v[13] * ind[2])) + 1.0e-15 * v[1] * par[9] * (par[7] - ind[2]) + v[11] * (15340.863387921299 * (par[7] - ind[2]) + -58.009079435348092 * (par[7] * par[7] - ind[2] * ind[2]) + 0.1503353174209219 * (par[7] * par[7] * par[7] - ind[2] * ind[2] * ind[2]) + -0.00019588923145049848 * (par[7] * par[7] * par[7] * par[7] - ind[2] * ind[2] * ind[2] * ind[2]) + 1.0402389841962685e-07 * (v[14] * v[14] * par[7] - v[15] * v[15] * ind[2])) - par[10] * v[9] * v[2] + 0 - v[16]) / (v[4] * (1.0e-15 * ind[1] * par[9] + 0.0180152833 * v[5] * (15340.863387921299 + -116.01815887069618 * ind[2] + 0.45100595226276569 * ind[2] * ind[2] + -0.00078355692580199391 * ind[2] * ind[2] * ind[2] + 5.2011949209813426e-07 * ind[2] * ind[2] * ind[2] * ind[2])) / (1.0e-15 * ind[1] + 0.0180152833 * v[5]) + par[21]);
    v[11] = (1.0e-15 * v[7] + 0.0180152833 * (v[10] / 0.0180152833 + v[11] / 0.0180152833 - v[6] * v[8]) - v[2] * 0.0180152833 * 1000 * (0.64038 + -0.00382 * ind[2] + 5.4633e-06 * ind[2] * ind[2]) * dxdt[2]) / v[3];
    dxdt[0] = v[11] / 0.1768325928384763;
    dxdt[1] = (v[7] - ind[1] * v[11]) / v[2];
    v[11] = par[11] * par[11];
    v[7] = ind[3] * ind[3];
    dxdt[3] = (0.0180152833 * 1000 * (-13.851 + 0.64038 * par[11] + -0.00191 * par[11] * par[11] + 1.8211e-06 * par[11] * par[11] * par[11]) * ind[4] * (15340.863387921299 * (par[11] - ind[3]) + -58.009079435348092 * (par[11] * par[11] - ind[3] * ind[3]) + 0.1503353174209219 * (par[11] * par[11] * par[11] - ind[3] * ind[3] * ind[3]) + -0.00019588923145049848 * (par[11] * par[11] * par[11] * par[11] - ind[3] * ind[3] * ind[3] * ind[3]) + 1.0402389841962685e-07 * (v[11] * v[11] * par[11] - v[7] * v[7] * ind[3])) + v[16]) / (0.0180152833 * (0.0180152833 * 1000 * (-13.851 + 0.64038 * ind[3] + -0.00191 * ind[3] * ind[3] + 1.8211e-06 * ind[3] * ind[3] * ind[3])) / 0.0180152833 * par[20] * (15340.863387921299 + -116.01815887069618 * ind[3] + 0.45100595226276569 * ind[3] * ind[3] + -0.00078355692580199391 * ind[3] * ind[3] * ind[3] + 5.2011949209813426e-07 * ind[3] * ind[3] * ind[3] * ind[3]) + par[13]);

    return dxdt;
}

template<class T>
std::vector<CppAD::AD<T> > CstrFunc(const std::vector<CppAD::AD<T> >& ind) {
    const std::vector<CppAD::AD<T> > ind2(ind.begin(), ind.begin() + 6);
    const std::vector<CppAD::AD<T> > par(ind.begin() + 6, ind.end());
    return CstrFunc(ind2, par);
}

#endif
