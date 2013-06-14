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
#include "CppADCGDynamicAtomicTest.hpp"
#include "cstr.hpp"

namespace CppAD {

    class CppADCGDynamicAtomicCstrTest : public CppADCGDynamicAtomicTest {
    protected:
        static const size_t n;
        static const size_t m;
    public:

        inline CppADCGDynamicAtomicCstrTest(bool verbose = false, bool printValues = false) :
            CppADCGDynamicAtomicTest("cstr_atomic", verbose, printValues) {
            this->verbose_ = false;
        }

        virtual std::vector<ADCGD> model(const std::vector<ADCGD>& u) {
            return CstrFunc<CGD>(u);
        }

    };

    const size_t CppADCGDynamicAtomicCstrTest::n = 28;
    const size_t CppADCGDynamicAtomicCstrTest::m = 4;
}

using namespace CppAD;
using namespace std;

TEST_F(CppADCGDynamicAtomicCstrTest, DynamicForRev) {
    using namespace std;
    using CppAD::vector;

    vector<Base> x(n);
    for (size_t i = 0; i < n; i++)
        x[i] = 1.0;
    
    vector<Base> xNorm(n);
    xNorm[0] = 0.3; // h
    xNorm[1] = 7.82e3; // Ca
    xNorm[2] = 304.65; // Tr
    xNorm[3] = 301.15; // Tj

    xNorm[4] = 2.3333e-04; // u1
    xNorm[5] = 6.6667e-05; // u2

    xNorm[6] = 6.2e14; // 
    xNorm[7] = 10080; //
    xNorm[8] = 2e3; //
    xNorm[9] = 10e3; //
    xNorm[10] = 1e-11; //
    xNorm[11] = 6.6667e-05; //
    xNorm[12] = 294.15; //
    xNorm[13] = 294.15; //
    xNorm[14] = 1000; //
    xNorm[15] = 4184; //Cp
    xNorm[16] = -33488; //deltaH
    xNorm[17] = 299.15; // Tj0
    xNorm[18] = 302.65; //   Tj2
    xNorm[19] = 7e5; // cwallj
    xNorm[20] = 1203; // csteam
    xNorm[21] = 3.22; //dsteam
    xNorm[22] = 950.0; //Ug
    xNorm[23] = 0.48649427192323; //vc6in
    xNorm[24] = 1000; //rhoj
    xNorm[25] = 4184; //Cpj
    xNorm[26] = 0.014; //Vj
    xNorm[27] = 1e-7; //cwallr

    vector<Base> eqNorm(m);
    eqNorm[0] = 0.3;
    eqNorm[1] = 7.82e3;
    eqNorm[2] = 304.65;
    eqNorm[3] = 301.15;

    this->testADFunAtomicLib(x, xNorm, eqNorm, 1e-14, 1e-13);
}

