#ifndef CPPAD_CG_BASE_DOUBLE_INCLUDED
#define	CPPAD_CG_BASE_DOUBLE_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<>
    inline CG<double> abs(const CG<double>& var) {
        if (var.isParameter()) {
            return CG<double> (fabs(var.getParameterValue()));
        } else {
            std::string operations = "fabs(" + var.operations() + ")";
            return CG<double>(*var.getCodeHandler(), operations, FUNCTION, &var);
        }
    }


}

#endif

