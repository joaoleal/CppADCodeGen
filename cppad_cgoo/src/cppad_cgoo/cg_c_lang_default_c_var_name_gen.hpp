#ifndef CPPAD_CG_C_LANG_DEFAULT_VAR_NAME_GEN_INCLUDED
#define	CPPAD_CG_C_LANG_DEFAULT_VAR_NAME_GEN_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <iosfwd>

namespace CppAD {

    /**
     * Creates variables names for the source code.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class DefaultCVariableNameGenerator : public VariableNameGenerator<Base> {
    protected:
        // auxiliary string stream
        std::stringstream _ss;
    public:

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) {
            _ss.clear();
            _ss.str("");
            _ss << "dep[" << index << "]";

            return _ss.str();
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            _ss.clear();
            _ss.str("");
            _ss << "ind[" << (independent.variableID() - 1) << "]";

            return _ss.str();
        }

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) {
            _ss.clear();
            _ss.str("");
            _ss << "var" << variable.variableID();

            return _ss.str();
        }
    };
}

#endif
