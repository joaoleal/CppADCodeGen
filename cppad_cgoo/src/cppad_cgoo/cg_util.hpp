#ifndef CPPAD_CG_UTIL_INCLUDED
#define	CPPAD_CG_UTIL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<class Base>
    inline std::vector<bool> jacobianForwardSparsity(ADFun<CppAD::CG<Base> >& fun) {
        size_t n = fun.Domain();

        std::vector<bool> r(n * n);
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++)
                r[j * n + k] = false;
            r[j * n + j] = true;
        }
        return fun.ForSparseJac(n, r);

    }

    template<class Base>
    inline std::vector<bool> jacobianReverseSparsity(ADFun<CppAD::CG<Base> >& fun) {
        size_t m = fun.Range();

        std::vector<bool> s(m * m);
        for (size_t i = 0; i < m; i++) {
            for (size_t k = 0; k < m; k++)
                s[i * m + k] = false;
            s[i * m + i] = true;
        }
        return fun.RevSparseJac(m, s);
    }

    template<class Base>
    inline std::vector<bool> jacobianSparsity(ADFun<CppAD::CG<Base> >& fun) {
        size_t m = fun.Range();
        size_t n = fun.Domain();

        if (n <= m) {
            // use forward mode 
            return jacobianForwardSparsity(fun);
        } else {
            // use reverse mode 
            return jacobianReverseSparsity(fun);
        }
    }

}

#endif

