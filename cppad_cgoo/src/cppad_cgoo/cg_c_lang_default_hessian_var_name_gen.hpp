#ifndef CPPAD_CG_C_LANG_DEFAULT_HESSIAN_VAR_NAME_GEN_INCLUDED
#define	CPPAD_CG_C_LANG_DEFAULT_HESSIAN_VAR_NAME_GEN_INCLUDED
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
     * Creates variables names for the source code generated for hessian
     * calculations.
     * The independet variables are considered to have been registered first as
     * variable in the code generation handler and then the multipliers.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangDefaultHessianVarNameGenerator : public CLangDefaultVariableNameGenerator<Base> {
    protected:
        // number of independent variables excluding the hessian multipliers
        size_t _indepCount;
        // array name of the independent variables
        std::string _multName;
    public:

        CLangDefaultHessianVarNameGenerator(size_t indepCount) :
            CLangDefaultVariableNameGenerator<Base>("hess", "ind", "var"),
            _indepCount(indepCount),
            _multName("mult") {
        }

        CLangDefaultHessianVarNameGenerator(const std::string& depName,
                                            const std::string& indepName,
                                            size_t indepCount,
                                            const std::string& multName,
                                            const std::string& tmpName) :
            CLangDefaultVariableNameGenerator<Base>(depName, indepName, tmpName),
            _indepCount(indepCount),
            _multName(multName) {
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            this->_ss.clear();
            this->_ss.str("");

            size_t index = independent.variableID() - 1;
            if (index < _indepCount) {
                this->_ss << this->_indepName << "[" << index << "]";
            } else {
                this->_ss << _multName << "[" << (index - _indepCount) << "]";
            }

            return this->_ss.str();
        }

    };
}

#endif
