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
     * Creates variables names for the source code generated for Hessian
     * calculations.
     * The independent variables are considered to have been registered first as
     * variable in the code generation handler and then the multipliers.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangDefaultHessianVarNameGenerator : public VariableNameGenerator<Base> {
    protected:
        VariableNameGenerator<Base>* _nameGen;
        // the lowest variable ID used for the equation multipliers
        const size_t _minMultiplierID;
        // array name of the independent variables
        const std::string _multName;
        // auxiliary string stream
        std::stringstream _ss;
    public:

        CLangDefaultHessianVarNameGenerator(VariableNameGenerator<Base>* nameGen,
                                            size_t n) :
            _nameGen(nameGen),
            _minMultiplierID(n + 1),
            _multName("mult") {

            CPPADCG_ASSERT_KNOWN(_nameGen != NULL, "The name generator must not be null");

            initialize();
        }

        CLangDefaultHessianVarNameGenerator(VariableNameGenerator<Base>* nameGen,
                                            const std::string& multName,
                                            size_t n) :
            _nameGen(nameGen),
            _minMultiplierID(n + 1),
            _multName(multName) {

            CPPADCG_ASSERT_KNOWN(_nameGen != NULL, "The name generator must not be null");
            CPPADCG_ASSERT_KNOWN(_multName.size() > 0, "The name for the multipliers must not be null");

            initialize();
        }

        virtual const std::vector<FuncArgument>& getDependent() const {
            return _nameGen->getDependent();
        }

        virtual const std::vector<FuncArgument>& getTemporaryVariables() const {
            return _nameGen->getTemporaryVariables();
        }

        virtual size_t getMinTemporaryVariableID() const {
            return _nameGen->getMinTemporaryVariableID();
        }

        virtual size_t getMaxTemporaryVariableID() const {
            return _nameGen->getMaxTemporaryVariableID();
        }

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) {
            return _nameGen->generateDependent(variable, index);
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            size_t id = independent.variableID();
            if (id < _minMultiplierID) {
                return _nameGen->generateIndependent(independent);
            }

            _ss.clear();
            _ss.str("");
            _ss << _multName << "[" << (id - _minMultiplierID) << "]";
            return _ss.str();
        }

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) {
            return _nameGen->generateTemporary(variable);
        }

        virtual void setTemporaryVariableID(size_t minTempID, size_t maxTempID) {
            _nameGen->setTemporaryVariableID(minTempID, maxTempID);
        }

    private:

        inline void initialize() {
            this->_independent = _nameGen->getIndependent(); // copy

            this->_independent.push_back(FuncArgument(_multName));
        }

    };
}

#endif
