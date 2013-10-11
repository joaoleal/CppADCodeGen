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
#include "../CppADCGDynamicTest.hpp"
#include "cstr.hpp"

using namespace CppAD;

namespace CppAD {

    class Cstr2DynamicTest : public CppADCGDynamicTest {
    public:

        inline Cstr2DynamicTest(bool verbose = false, bool printValues = false) :
            CppADCGDynamicTest("cstr2", verbose, printValues) {
        }

        virtual std::vector<ADCGD> model(const std::vector<ADCGD>& ind) {
            std::vector<ADCGD> eqs(8);
            std::vector<ADCGD> ind1(28);
            std::copy(ind.begin(), ind.begin() + 28, ind1.begin());
            std::vector<ADCGD> sys1 = CstrFunc<CG<double> >(ind1);
            std::copy(sys1.begin(), sys1.end(), eqs.begin());

            std::vector<ADCGD> ind2(28);
            std::copy(ind.begin() + 28, ind.end(), ind2.begin());
            std::vector<ADCGD> sys2 = CstrFunc<CG<double> >(ind2);
            std::copy(sys2.begin(), sys2.end(), eqs.begin() + sys1.size());
            return eqs;
        }

    };
}

TEST_F(Cstr2DynamicTest, cstr) {
    // independent variable vector
    std::vector<double> xNorm(28);
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


    std::vector<Base> eqNorm(4);
    eqNorm[0] = 0.3;
    eqNorm[1] = 7.82e3;
    eqNorm[2] = 304.65;
    eqNorm[3] = 301.15;

    
    std::vector<double> xNorm2(xNorm.size() * 2);
    std::copy(xNorm.begin(), xNorm.end(), xNorm2.begin());
    std::copy(xNorm.begin(), xNorm.end(), xNorm2.begin() + xNorm.size());

    std::vector<double> x2(xNorm2.size(), 1.0);
    std::vector<AD<CG<double> > > u2(x2.size());
    for (size_t i = 0; i < x2.size(); i++)
        u2[i] = x2[i];

    std::vector<Base> eqNorm2(eqNorm.size() * 2);
    std::copy(eqNorm.begin(), eqNorm.end(), eqNorm2.begin());
    std::copy(eqNorm.begin(), eqNorm.end(), eqNorm2.begin() + eqNorm.size());

    this->testDynamicFull(u2, x2, xNorm2, eqNorm2, 1000);

}
