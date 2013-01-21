/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
#include <cppadcg/dae_index_reduction/cg_pantelides.hpp>

#include "CppADCGIndexReductionTest.hpp"
#include "pendulum.hpp"

using namespace CppAD;

TEST_F(CppADCGIndexReductionTest, PantelidesPendulum2D) {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum2D<CGD > ();

    std::vector<DaeVarInfo> daeVar(10);
    daeVar[5].makeConstant();
    daeVar[6] = 0;
    daeVar[7] = 1;
    daeVar[8] = 2;
    daeVar[9] = 3;

    Plantelides<double> pantelides(fun, daeVar);

    std::vector<DaeVarInfo> newDaeVar;
    ADFun<CGD>* reducedFun = pantelides.reduceIndex(newDaeVar);

    delete fun;
    delete reducedFun;

    ASSERT_TRUE(reducedFun != NULL);
}

TEST_F(CppADCGIndexReductionTest, PantelidesPendulum3D) {
    using namespace CppAD;
    using namespace std;
    typedef CG<double> CGD;

    // create f: U -> Z and vectors used for derivative calculations
    ADFun<CGD>* fun = Pendulum3D<CGD > ();

    std::vector<DaeVarInfo> daeVar(13);
    daeVar[7] = 0;
    daeVar[8] = 1;
    daeVar[9] = 2;
    daeVar[10] = 3;
    daeVar[11] = 4;
    daeVar[12] = 5;

    Plantelides<double> pantelides(fun, daeVar);

    std::vector<DaeVarInfo> newDaeVar;
    ADFun<CGD>* reducedFun = pantelides.reduceIndex(newDaeVar);

    delete fun;
    delete reducedFun;

    ASSERT_TRUE(reducedFun != NULL);
}
