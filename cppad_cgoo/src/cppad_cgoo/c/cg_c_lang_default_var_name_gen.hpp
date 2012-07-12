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
        // the lowest variable ID used for the temporary variables
        size_t _minTemporaryID;
        // the highest variable ID used for the temporary variables
        size_t _maxTemporaryID;
    public:

        CLangDefaultVariableNameGenerator() :
            _depName("dep"),
            _indepName("ind"),
            _tmpName("var") {
            this->_independent.push_back(FuncArgument(_indepName));
            this->_dependent.push_back(FuncArgument(_depName));
            this->_temporary.push_back(FuncArgument(_tmpName));
        }

        CLangDefaultVariableNameGenerator(const std::string& depName,
                                          const std::string& indepName,
                                          const std::string& tmpName) :
            _depName(depName),
            _indepName(indepName),
            _tmpName(tmpName) {
            this->_independent.push_back(FuncArgument(_indepName));
            this->_dependent.push_back(FuncArgument(_depName));
            this->_temporary.push_back(FuncArgument(_tmpName));
        }

        virtual size_t getMinTemporaryVariableID() const {
            return _minTemporaryID;
        }

        virtual size_t getMaxTemporaryVariableID() const {
            return _maxTemporaryID;
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
            if (this->_temporary[0].array) {
                _ss << _tmpName << "[" << (id - this->_minTemporaryID) << "]";
            } else {
                _ss << _tmpName << id;
            }

            return _ss.str();
        }

        virtual void setTemporaryVariableID(size_t minTempID, size_t maxTempID) {
            _minTemporaryID = minTempID;
            _maxTemporaryID = maxTempID;

            // if
            //  _minTemporaryID == _maxTemporaryID + 1
            // then no temporary variables are being used
            assert(_minTemporaryID <= _maxTemporaryID + 1);
        }

    };
}

#endif
