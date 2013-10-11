/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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

#include "CppADCGOperationTest.hpp"
#include "compare_change.hpp"

namespace CppAD {

    class CppADCGCompareTest : public CppADCGOperationTest {
    protected:
        std::vector<std::vector<double> > xV;
    public:

        inline CppADCGCompareTest(bool verbose = false, bool printValues = false) :
            CppADCGOperationTest(verbose, printValues) {
            std::vector<double> X(2);
            // create independent variables
            X[0] = 3.;
            X[1] = 4.;
            xV.push_back(X);
            X[0] = 4.;
            X[1] = 3.;
            xV.push_back(X);
        }
    };
}

using namespace CppAD;

// ------------------------------- < ----------------------------

TEST_F(CppADCGCompareTest, CompareChange1) {
    int comparisons;
    ASSERT_THROW(test0("CompareChange<",
                       &CompareChangeFunc1<double >,
                       &CompareChangeFunc1<CG<double> >,
                       xV, comparisons), CppAD::CGException);

    //ASSERT_EQ(comparisons, 6); // all comparisons have changed
}

// ------------------------------- > ----------------------------

TEST_F(CppADCGCompareTest, CompareChange2) {
    int comparisons;
    ASSERT_THROW(test0("CompareChange>",
                       &CompareChangeFunc2<double >,
                       &CompareChangeFunc2<CG<double> >,
                       xV, comparisons), CppAD::CGException);

    //ASSERT_EQ(comparisons, 6); // all comparisons have changed
}

// ------------------------------- <= ----------------------------

TEST_F(CppADCGCompareTest, CompareChange3) {
    int comparisons;
    ASSERT_THROW(test0("CompareChange<=",
                       &CompareChangeFunc3<double >,
                       &CompareChangeFunc3<CG<double> >,
                       xV, comparisons), CppAD::CGException);

    //ASSERT_EQ(comparisons, 6); // all comparisons have changed
}

// ------------------------------- >= ----------------------------

TEST_F(CppADCGCompareTest, CompareChange4) {
    int comparisons;
    ASSERT_THROW(test0("CompareChange>=",
                       &CompareChangeFunc4<double >,
                       &CompareChangeFunc4<CG<double> >,
                       xV, comparisons), CppAD::CGException);

    //ASSERT_EQ(comparisons, 6); // all comparisons have changed
}

// ------------------------------- == ----------------------------

TEST_F(CppADCGCompareTest, CompareChange5) {
    int comparisons;
    ASSERT_THROW(test0("CompareChange>=",
                       &CompareChangeFunc4<double >,
                       &CompareChangeFunc4<CG<double> >,
                       xV, comparisons), CppAD::CGException);

    //ASSERT_EQ(comparisons, 4); // the first two comparisons do not change
}

// ------------------------------- != ----------------------------

TEST_F(CppADCGCompareTest, CompareChange6) {
    int comparisons;
    ASSERT_THROW(test0("CompareChange>=",
                       &CompareChangeFunc4<double >,
                       &CompareChangeFunc4<CG<double> >,
                       xV, comparisons), CppAD::CGException);

    //ASSERT_EQ(comparisons, 4); // the first two comparisons do not change
}
