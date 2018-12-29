/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
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



namespace CppAD {
namespace cg {

class CstrDynamicTest : public CppADCGDynamicTest {
public:

    inline CstrDynamicTest(bool verbose = false,
                           bool printValues = false) :
        CppADCGDynamicTest("cstr", verbose, printValues) {
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ind,
                             const std::vector<ADCGD>& par) override {
        assert(par.empty());
        return CstrFunc<CG<double> >(ind);
    }

};

class CstrDynamicWithParamsTest : public CppADCGDynamicTest {
public:

    inline CstrDynamicWithParamsTest(bool verbose = false,
                                     bool printValues = false) :
            CppADCGDynamicTest("cstr_params", verbose, printValues) {
    }

    std::vector<ADCGD> model(const std::vector<ADCGD>& ind,
                             const std::vector<ADCGD>& par) override {
        assert(!par.empty());
        return CstrFunc<CG<double> >(ind, par);
    }

};

} // END cg namespace
} // END CppAD namespace

using namespace CppAD;
using namespace CppAD::cg;

TEST_F(CstrDynamicTest, cstr) {
    // independent variable vector
    std::vector<double> x(28, 1.0);
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

    std::vector<AD<CG<double> > > u(28);
    for (size_t i = 0; i < x.size(); i++)
        u[i] = x[i];

    std::vector<Base> eqNorm(4);
    eqNorm[0] = 0.3;
    eqNorm[1] = 7.82e3;
    eqNorm[2] = 304.65;
    eqNorm[3] = 301.15;

    this->testDynamicFull(u, x, xNorm, eqNorm, 1000);

}

TEST_F(CstrDynamicWithParamsTest, cstr_params) {
    // independent variable vector
    std::vector<double> x(6, 1.0);
    std::vector<double> xNorm(6);
    xNorm[0] = 0.3; // h
    xNorm[1] = 7.82e3; // Ca
    xNorm[2] = 304.65; // Tr
    xNorm[3] = 301.15; // Tj

    xNorm[4] = 2.3333e-04; // u1
    xNorm[5] = 6.6667e-05; // u2

    std::vector<double> par(22, 1.0);
    par[0] = 6.2e14; //
    par[1] = 10080; //
    par[2] = 2e3; //
    par[3] = 10e3; //
    par[4] = 1e-11; //
    par[5] = 6.6667e-05; //
    par[6] = 294.15; //
    par[7] = 294.15; //
    par[8] = 1000; //
    par[9] = 4184; //Cp
    par[10] = -33488; //deltaH
    par[11] = 299.15; // Tj0
    par[12] = 302.65; //   Tj2
    par[13] = 7e5; // cwallj
    par[14] = 1203; // csteam
    par[15] = 3.22; //dsteam
    par[16] = 950.0; //Ug
    par[17] = 0.48649427192323; //vc6in
    par[18] = 1000; //rhoj
    par[19] = 4184; //Cpj
    par[20] = 0.014; //Vj
    par[21] = 1e-7; //cwallr

    std::vector<AD<CG<double> > > ax(x.size());
    for (size_t i = 0; i < x.size(); i++)
        ax[i] = x[i];

    std::vector<AD<CG<double> > > ap(par.size());
    for (size_t i = 0; i < par.size(); i++)
        ap[i] = par[i];

    std::vector<Base> eqNorm(4);
    eqNorm[0] = 0.3;
    eqNorm[1] = 7.82e3;
    eqNorm[2] = 304.65;
    eqNorm[3] = 301.15;

    this->testDynamicFull(ax, ap, x, par, xNorm, eqNorm, 1000);

}
