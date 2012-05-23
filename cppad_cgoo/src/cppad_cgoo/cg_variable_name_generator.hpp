#ifndef CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
#define	CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <string>

namespace CppAD {

    /**
     * Creates variables names for the source code.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class VariableNameGenerator {
    public:
        virtual std::string generateDependent(const CG<Base>& variable, size_t index) = 0;

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& variable) = 0;

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) = 0;
    };
}

#endif
