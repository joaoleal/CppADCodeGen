/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
#include "CppADCGPatternTest.hpp"
#include "../models/cstr.hpp"

typedef double Base;
typedef CppAD::CG<Base> CGD;
typedef CppAD::AD<CGD> ADCGD;

using namespace CppAD;

size_t ns = 4;
size_t nm = 2;
size_t npar = 23;

size_t K = 3;
size_t m = K * ns;
size_t na = ns + nm + npar;

std::vector<ADCGD> modelCollocation(std::vector<ADCGD>& x, size_t repeat, const std::vector<CGAbstractAtomicFun<double>*>& atoms) {
    CGAbstractAtomicFun<Base>& atomicCstr = *atoms[0];

    size_t m2 = repeat * m;

    // dependent variable vector 
    std::vector<ADCGD> dep(m2);

    std::vector<ADCGD> dxikdt(ns);
    std::vector<ADCGD> xik(na);

    // parameters
    for (size_t j = 0; j < npar; j++)
        xik[ns + nm + j] = x[j];

    size_t s = npar;
    size_t nvarsk = ns;
    size_t nMstart = npar + nvarsk * K * repeat + nvarsk;
    size_t eq = 0;

    for (size_t i = 0; i < repeat; i++) {
        size_t s0 = s;

        // controls
        for (size_t j = 0; j < nm; j++) {
            xik[ns + j] = x[nMstart + nm * i + j];
        }

        // K = 1
        s += nvarsk;

        for (size_t j = 0; j < ns; j++) {
            xik[j] = x[s + j]; // states
        }
        // xik[ns + nm + npar] = x[s + ns];// time

        atomicCstr(xik, dxikdt); // ODE
        for (size_t j = 0; j < ns; j++) {
            dep[eq + j] = dxikdt[j]
                    + 0.13797958971132715 * x[s0 + j]
                    + -0.10749149571305303 * x[s0 + nvarsk + j]
                    + -0.038928002823013501 * x[s0 + 2 * nvarsk + j]
                    + 0.008439908824739363 * x[s0 + 3 * nvarsk + j];
        }
        eq += ns;

        // K = 2
        s += nvarsk;

        for (size_t j = 0; j < ns; j++) {
            xik[j] = x[s + j]; // states
        }
        // xik[ns + nm + npar] = x[s + ns];// time

        atomicCstr(xik, dxikdt); // ODE
        for (size_t j = 0; j < ns; j++) {
            dep[eq + j] = dxikdt[j]
                    + -0.057979589711327127 * x[s0 + j]
                    + 0.11892800282301351 * x[s0 + nvarsk + j]
                    + -0.025841837620280327 * x[s0 + 2 * nvarsk + j]
                    + -0.035106575491406049 * x[s0 + 3 * nvarsk + j];
        }
        eq += ns;

        // K = 3

        s += nvarsk;
        for (size_t j = 0; j < ns; j++) {
            xik[j] = x[s + j]; // states
        }
        // xik[ns + nm + npar] = x[s + ns];// time

        atomicCstr(xik, dxikdt); // ODE
        for (size_t j = 0; j < ns; j++) {
            dep[eq + j] = dxikdt[j]
                    + 0.099999999999999978 * x[s0 + j]
                    + -0.18439908824739357 * x[s0 + nvarsk + j]
                    + 0.25106575491406025 * x[s0 + 2 * nvarsk + j]
                    + -0.16666666666666669 * x[s0 + 3 * nvarsk + j];
        }
        eq += ns;

    }

    return dep;
}

void atomicFunction(const std::vector<AD<double> >& x, std::vector<AD<double> >& y) {
    y = CstrFunc(x);
}

TEST_F(CppADCGPatternTest, Atomic) {
    using namespace CppAD;

    /**
     * 
     */
    std::vector<Base> xx(na);
    xx[0] = 0.3; // h
    xx[1] = 7.82e3; // Ca
    xx[2] = 304.65; // Tr
    xx[3] = 301.15; // Tj

    xx[4] = 2.3333e-04; // u1
    xx[5] = 6.6667e-05; // u2

    xx[6] = 6.2e14; // 
    xx[7] = 10080; //
    xx[8] = 2e3; //
    xx[9] = 10e3; //
    xx[10] = 1e-11; //
    xx[11] = 6.6667e-05; //
    xx[12] = 294.15; //
    xx[13] = 294.15; //
    xx[14] = 1000; //
    xx[15] = 4184; //Cp
    xx[16] = -33488; //deltaH
    xx[17] = 299.15; // Tj0
    xx[18] = 302.65; //   Tj2
    xx[19] = 7e5; // cwallj
    xx[20] = 1203; // csteam
    xx[21] = 3.22; //dsteam
    xx[22] = 950.0; //Ug
    xx[23] = 0.48649427192323; //vc6in
    xx[24] = 1000; //rhoj
    xx[25] = 4184; //Cpj
    xx[26] = 0.014; //Vj
    xx[27] = 1e-7; //cwallr

    /**
     * create atomic function for the ODE
     */
    std::vector<AD<double> > ay(ns), ax(na);
    for (size_t j = 0; j < na; j++)
        ax[j] = xx[j];

    checkpoint<double> atomicfun("atomicFunc", atomicFunction, ax, ay);
    std::vector<atomic_base<double>*> atomics(1);
    atomics[0] = &atomicfun;

    size_t repeat = 6;

    size_t nvarsk = ns;
    size_t nMstart = npar + nvarsk * K * repeat + nvarsk;

    std::vector<Base> x(nMstart + repeat * nm);
    // parameters
    for (size_t j = 0; j < npar; j++)
        x[j] = xx[ns + nm + j];

    size_t s = npar;

    for (size_t i = 0; i < repeat; i++) {
        // controls
        for (size_t j = 0; j < nm; j++) {
            x[nMstart + nm * i + j] = xx[ns + j];
        }

        // K = 1
        s += nvarsk;
        // states
        for (size_t j = 0; j < ns; j++) {
            x[s + j] = xx[j];
        }

        // K = 2
        s += nvarsk;
        // states
        for (size_t j = 0; j < ns; j++) {
            x[s + j] = xx[j];
        }

        // K = 3
        // states
        s += nvarsk;
        for (size_t j = 0; j < ns; j++) {
            x[s + j] = xx[j];
        }
    }


    testPatternDetectionWithAtomics(modelCollocation, atomics, m, x, repeat);

    testLibCreationWithAtomics(modelCollocation, atomics, m, x, repeat, "modelAtomicCstr", true, true);
}
