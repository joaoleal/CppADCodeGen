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
#include "../CppADCGDynamicTest.hpp"
#include "cstr.hpp"

using namespace CppAD;

namespace CppAD {

    class CstrDynamicTest : public CppADCGDynamicTest {
    public:

        inline CstrDynamicTest(bool verbose = false, bool printValues = false) :
            CppADCGDynamicTest("cstr", verbose, printValues) {
        }
    };
}

TEST_F(CstrDynamicTest, cstr) {
    // independent variable vector
    std::vector<double> x(28);
    x[0] = 0.3; // h
    x[1] = 7.82e3; // Ca
    x[2] = 304.65; // Tr
    x[3] = 301.15; // Tj

    x[4] = 2.3333e-04; // u1
    x[5] = 6.6667e-05; // u2

    x[6] = 6.2e14; // 
    x[7] = 10080; //
    x[8] = 2e3; //
    x[9] = 10e3; //
    x[10] = 1e-11; //
    x[11] = 6.6667e-05; //
    x[12] = 294.15; //
    x[13] = 294.15; //
    x[14] = 1000; //
    x[15] = 4184; //Cp
    x[16] = -33488; //deltaH
    x[17] = 299.15; // Tj0
    x[18] = 302.65; //   Tj2
    x[19] = 7e5; // cwallj
    x[20] = 1203; // csteam
    x[21] = 3.22; //dsteam
    x[22] = 950.0; //Ug
    x[23] = 0.48649427192323; //vc6in
    x[24] = 1000; //rhoj
    x[25] = 4184; //Cpj
    x[26] = 0.014; //Vj
    x[27] = 1e-7; //cwallr

    std::vector<AD<CG<double> > > u(28);
    for (size_t i = 0; i < x.size(); i++)
        u[i] = x[i];

    this->testDynamic1(u, x, &CstrFunc<CG<double> >, 1000, 1e-10, 1e-9);

}
