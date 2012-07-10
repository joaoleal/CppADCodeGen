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
     * The independent variables are considered to have been registered first as
     * variable in the code generation handler and then the multipliers.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangDefaultHessianVarNameGenerator : public CLangDefaultVariableNameGenerator<Base> {
    protected:
        // the lowest variable ID used for the equation multipliers
        const size_t _minMultiplierID;
        // the lowest variable ID used for the dependent variables
        const size_t _minDependentID;
        // array name of the independent variables
        const std::string _multName;
    public:

        CLangDefaultHessianVarNameGenerator(size_t m, size_t n) :
            CLangDefaultVariableNameGenerator<Base>("hess", "ind", "var"),
            _minMultiplierID(n + 1),
            _minDependentID(_minMultiplierID + m),
            _multName("mult") {
        }

        CLangDefaultHessianVarNameGenerator(const std::string& depName,
                                            const std::string& indepName,
                                            const std::string& multName,
                                            const std::string& tmpName,
                                            size_t m, size_t n) :
            CLangDefaultVariableNameGenerator<Base>(depName, indepName, tmpName),
            _minMultiplierID(n + 1),
            _minDependentID(_minMultiplierID + m),
            _multName(multName) {
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            this->_ss.clear();
            this->_ss.str("");

            size_t id = independent.variableID();
            if (id < _minMultiplierID) {
                this->_ss << this->_indepName << "[" << (id - 1) << "]";
            } else {
                this->_ss << _multName << "[" << (id - _minMultiplierID) << "]";
            }

            return this->_ss.str();
        }

    };
}

#endif
