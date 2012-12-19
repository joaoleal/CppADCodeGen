#ifndef CPPADCGOO_TEST_COND_EXP_INCLUDED
#define	CPPADCGOO_TEST_COND_EXP_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <assert.h>

template<class T>
CppAD::ADFun<T>* CondExp_pvvvFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;
    using namespace std;

    assert(X.size() == 3);

    // parameter value
    AD<T> one = 1.;

    // dependent variable vector 
    std::vector< AD<T> > Y(5);

    // CondExp(parameter, variable, variable, variable)
    Y[0] = CondExpLt(one, X[0], X[1], X[2]);
    Y[1] = CondExpLe(one, X[0], X[1], X[2]);
    Y[2] = CondExpEq(one, X[0], X[1], X[2]);
    Y[3] = CondExpGe(one, X[0], X[1], X[2]);
    Y[4] = CondExpGt(one, X[0], X[1], X[2]);

    // create f: X -> Y 
    return new ADFun<T > (X, Y);
}

template<class T>
CppAD::ADFun<T>* CondExp_vpvvFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;
    using namespace std;

    assert(X.size() == 3);

    // parameter value
    AD<T> one = 1.;

    // dependent variable vector 
    std::vector< AD<T> > Y(5);

    // CondExp(variable, parameter, variable, variable)
    Y[0] = CondExpLt(X[0], one, X[1], X[2]);
    Y[1] = CondExpLe(X[0], one, X[1], X[2]);
    Y[2] = CondExpEq(X[0], one, X[1], X[2]);
    Y[3] = CondExpGe(X[0], one, X[1], X[2]);
    Y[4] = CondExpGt(X[0], one, X[1], X[2]);

    // create f: X -> Y 
    return new ADFun<T > (X, Y);
}

template<class T>
CppAD::ADFun<T>* CondExp_vvpvFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;
    using namespace std;

    assert(X.size() == 3);

    // parameter value
    AD<T> three = 3.;

    // dependent variable vector 
    std::vector< AD<T> > Y(5);

    // CondExp(variable, variable, parameter, variable)
    Y[0] = CondExpLt(X[0], X[1], three, X[2]);
    Y[1] = CondExpLe(X[0], X[1], three, X[2]);
    Y[2] = CondExpEq(X[0], X[1], three, X[2]);
    Y[3] = CondExpGe(X[0], X[1], three, X[2]);
    Y[4] = CondExpGt(X[0], X[1], three, X[2]);

    // create f: X -> Y 
    return new ADFun<T > (X, Y);
}

template<class T>
CppAD::ADFun<T>* CondExp_vvvpFunc(const std::vector<CppAD::AD<T> >& X) {
    using namespace CppAD;
    using namespace std;

    assert(X.size() == 3);

    // parameter value
    AD<T> three = 3.;

    // dependent variable vector 
    std::vector< AD<T> > Y(5);

    // CondExp(variable, variable, variable, parameter)
    Y[0] = CondExpLt(X[0], X[1], X[2], three);
    Y[1] = CondExpLe(X[0], X[1], X[2], three);
    Y[2] = CondExpEq(X[0], X[1], X[2], three);
    Y[3] = CondExpGe(X[0], X[1], X[2], three);
    Y[4] = CondExpGt(X[0], X[1], X[2], three);

    // create f: X -> Y 
    return new ADFun<T > (X, Y);
}

#endif