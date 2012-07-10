#ifndef CPPAD_CG_C_LANGUAGE_DOUBLE_INCLUDED
#define	CPPAD_CG_C_LANGUAGE_DOUBLE_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    /**
     * Specialization of the C language function names for doubles
     * 
     * \author Joao Leal
     */
    template<>
    inline const std::string& CLanguage<double>::absFuncName() {
        static const std::string name("fabs");
        return name;
    }
}

#endif