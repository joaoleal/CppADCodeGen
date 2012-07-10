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

namespace CppAD {

    /**
     * Creates variables names for the source code.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangDefaultVariableNameGenerator : public VariableNameGenerator<Base> {
    protected:
        // auxiliary string stream
        std::stringstream _ss;
        // array name of the dependent variables
        std::string _depName;
        // array name of the independent variables
        std::string _indepName;
        // array name of the temporary variables
        std::string _tmpName;
        // whether or not to use an array for the temporary variables
        bool _useArray4TempVars;
    public:

        CLangDefaultVariableNameGenerator() :
            _depName("dep"),
            _indepName("ind"),
            _tmpName("var"),
            _useArray4TempVars(true) {
        }

        CLangDefaultVariableNameGenerator(const std::string& depName,
                                          const std::string& indepName,
                                          const std::string& tmpName) :
            _depName(depName),
            _indepName(indepName),
            _tmpName(tmpName),
            _useArray4TempVars(true) {
        }

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) {
            _ss.clear();
            _ss.str("");

            _ss << _depName << "[" << index << "]";

            return _ss.str();
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            _ss.clear();
            _ss.str("");

            size_t id = independent.variableID();
            _ss << _indepName << "[" << (id - 1) << "]";

            return _ss.str();
        }

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) {
            _ss.clear();
            _ss.str("");

            size_t id = variable.variableID();
            if (_useArray4TempVars) {
                _ss << _tmpName << "[" << (id - this->_minTemporaryID) << "]";
            } else {
                _ss << _tmpName << id;
            }

            return _ss.str();
        }

        virtual bool isTemporaryVariablesArray() const {
            return _useArray4TempVars;
        }

        virtual const std::string& getTemporaryVariablesArrayName() const {
            return _tmpName;
        }
    };
}

#endif
