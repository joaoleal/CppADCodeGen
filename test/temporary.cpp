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
#include "CppADCGTest.hpp"

using namespace CppAD;

namespace CppAD {

    class CppADCGTempTest : public CppADCGTest {
    protected:
        typedef CppADCGTest::CGD CGD;
        typedef CppADCGTest::ADCGD ADCGD;
    public:

        inline CppADCGTempTest(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues) {
        }

        void testModel(ADFun<CGD>& f, size_t expectedTmp, size_t expectedArraySize) {
            size_t n = f.Domain();
            //size_t m = f.Range();

            CodeHandler<double> handler(10 + n * n);

            std::vector<CGD> indVars(n);
            handler.makeVariables(indVars);

            std::vector<CGD> dep = f.Forward(0, indVars);

            CLanguage<double> langC("double");
            CLangDefaultVariableNameGenerator<double> nameGen;

            handler.generateCode(std::cout, langC, dep, nameGen);

            ASSERT_EQ(handler.getTemporaryVariableCount(), expectedTmp);
            ASSERT_EQ(handler.getTemporaryArraySize(), expectedArraySize);
        }
    };
}

TEST_F(CppADCGTempTest, NoTemporary) {
    size_t n = 3;
    size_t m = 2;

    std::vector<ADCGD> u(n); // independent variable vector 
    u[0] = 1;
    u[0] = 2;
    u[0] = 3;
    Independent(u);

    std::vector<ADCGD> Z(m); // dependent variable vector 

    // model
    Z[0] = u[0] + u[1];
    ADCGD tmp = u[1] * u[2]; // this temporary variable should disapear
    Z[1] = tmp;

    ADFun<CGD> f(u, Z);
    testModel(f, 0, 0);
}

TEST_F(CppADCGTempTest, Temporary1) {
    size_t n = 3;
    size_t m = 2;

    std::vector<ADCGD> u(n); // independent variable vector 
    u[0] = 1;
    u[0] = 2;
    u[0] = 3;
    Independent(u);

    std::vector<ADCGD> Z(m); // dependent variable vector 

    // model
    Z[0] = u[0] + u[1];
    ADCGD tmp = u[1] * u[2]; // this temporary variable should NOT disapear
    Z[1] = tmp + tmp;

    ADFun<CGD> f(u, Z);
    testModel(f, 1, 0);
}

TEST_F(CppADCGTempTest, Temporary2) {
    size_t n = 3;
    size_t m = 2;

    std::vector<ADCGD> u(n); // independent variable vector 
    u[0] = 1;
    u[0] = 2;
    u[0] = 3;
    Independent(u);

    std::vector<ADCGD> Z(m); // dependent variable vector 

    // model
    Z[0] = u[0] + u[1];
    ADCGD tmp = u[1] * u[2];
    ADCGD tmp1 = tmp + 1;
    Z[1] = tmp1 + tmp1;

    ADFun<CGD> f(u, Z);
    testModel(f, 1, 0);
}