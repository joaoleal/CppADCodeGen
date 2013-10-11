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

#include <cppadcg/cg.hpp>
#include <gtest/gtest.h>

#include "CppADCGTest.hpp"

using namespace CppAD;
using namespace std;

TEST_F(CppADCGTest, multMatrixMatrixSparsityTrans) {
    std::vector<std::set<size_t> > aT(4); // a: 3 x 4
    aT[0].insert(0);
    aT[1].insert(1);
    aT[2].insert(1);
    aT[3].insert(0);
    aT[3].insert(2);
    std::vector<std::set<size_t> > b(4); // b: 4 x 2
    b[0].insert(0);
    b[2].insert(0);
    b[2].insert(1);
    b[3].insert(1);
    CppAD::vector<std::set<size_t> > rT(2); // r: 3 x 2 

    multMatrixMatrixSparsityTrans(aT, b, rT, 4, 2, 3);

    CppAD::vector<std::set<size_t> > rTExpected(2);
    rTExpected[0].insert(0);
    rTExpected[0].insert(1);
    rTExpected[1].insert(0);
    rTExpected[1].insert(1);
    rTExpected[1].insert(2);

    compareVectorSetValues(rT, rTExpected);
}

TEST_F(CppADCGTest, multMatrixMatrixSparsityTrans2) {
    size_t m = 2;
    size_t n = 7;
    size_t q = 33;
    std::vector<std::set<size_t> > sT(m);
    sT[0].insert(29);
    sT[1].insert(30);

    std::vector<std::set<size_t> > jac(m);
    jac[0].insert(0);
    jac[0].insert(2);
    jac[1].insert(1);
    jac[1].insert(3);

    CppAD::vector<std::set<size_t> > rT(n);

    multMatrixMatrixSparsityTrans(sT, jac, rT, m, n, q);

    CppAD::vector<std::set<size_t> > rTExpected(n);
    rTExpected[0].insert(29);
    rTExpected[1].insert(30);
    rTExpected[2].insert(29);
    rTExpected[3].insert(30);

    compareVectorSetValues(rT, rTExpected);
}

TEST_F(CppADCGTest, multMatrixMatrixSparsity) {
    size_t m = 4;
    size_t n = 3;
    size_t q = 2;
    std::vector<std::set<size_t> > a(m);
    a[0].insert(2);
    a[1].insert(1);
    a[2].insert(0);
    a[2].insert(1);
    a[2].insert(2);
    a[3].insert(0);

    std::vector<std::set<size_t> > b(n);
    b[0].insert(0);
    b[1].insert(1);
    b[2].insert(0);

    CppAD::vector<std::set<size_t> > r(m);

    multMatrixMatrixSparsity(a, b, r, m, n, q);

    CppAD::vector<std::set<size_t> > rExpected(m);
    rExpected[0].insert(0);
    rExpected[1].insert(1);
    rExpected[2].insert(0);
    rExpected[2].insert(1);
    rExpected[3].insert(0);

    compareVectorSetValues(r, rExpected);
}


TEST_F(CppADCGTest, multMatrixTransMatrixSparsity) {
    size_t m = 4;
    size_t n = 3;
    size_t q = 2;
    std::vector<std::set<size_t> > a(n);
    a[0].insert(2);
    a[0].insert(3);
    a[1].insert(1);
    a[1].insert(2);
    a[2].insert(0);
    a[2].insert(2);
    
    std::vector<std::set<size_t> > b(n);
    b[0].insert(0);
    b[1].insert(1);
    b[2].insert(0);

    CppAD::vector<std::set<size_t> > r(m);

    multMatrixTransMatrixSparsity(a, b, r, n, m, q);

    CppAD::vector<std::set<size_t> > rExpected(m);
    rExpected[0].insert(0);
    rExpected[1].insert(1);
    rExpected[2].insert(0);
    rExpected[2].insert(1);
    rExpected[3].insert(0);

    compareVectorSetValues(r, rExpected);
}
