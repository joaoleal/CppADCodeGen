#ifndef CPPAD_CG_LANGUAGE_INCLUDED
#define	CPPAD_CG_LANGUAGE_INCLUDED
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
    class Language {
    protected:
        virtual void generateSourceCode(std::ostream& out,
                std::vector<SourceCodeFragment<Base> *>& independent,
                std::vector<CG<Base> >& dependent,
                const std::vector<SourceCodeFragment<Base>*>& variableOrder,
                VariableNameGenerator<Base>& nameGen) = 0;

        virtual size_t getMaximumCodeBlockVisit() = 0;

        virtual bool createsNewVariable(const SourceCodeFragment<Base>& op) = 0;

        virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) = 0;

        friend class CodeHandler<Base>;
    };

}

#endif

