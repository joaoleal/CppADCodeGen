#ifndef CPPAD_CG_C_LANG_DEFAULT_HESSIAN_VAR_NAME_GEN_INCLUDED
#define CPPAD_CG_C_LANG_DEFAULT_HESSIAN_VAR_NAME_GEN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    /**
     * Creates variables names for the source code generated for Hessian
     * calculations.
     * The independent variables are considered to have been registered first as
     * variable in the code generation handler and then the multipliers.
     * 
     * @author Joao Leal
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
            CPPADCG_ASSERT_KNOWN(_multName.size() > 0, "The name for the multipliers must not be empty");

            initialize();
        }

        virtual const std::vector<FuncArgument>& getDependent() const {
            return _nameGen->getDependent();
        }

        virtual const std::vector<FuncArgument>& getTemporary() const {
            return _nameGen->getTemporary();
        }

        virtual size_t getMinTemporaryVariableID() const {
            return _nameGen->getMinTemporaryVariableID();
        }

        virtual size_t getMaxTemporaryVariableID() const {
            return _nameGen->getMaxTemporaryVariableID();
        }

        virtual size_t getMaxTemporaryArrayVariableID() const {
            return _nameGen->getMaxTemporaryArrayVariableID();
        }

        virtual std::string generateDependent(size_t index) {
            return _nameGen->generateDependent(index);
        }

        virtual std::string generateIndependent(const OperationNode<Base>& independent) {
            size_t id = independent.getVariableID();
            if (id < _minMultiplierID) {
                return _nameGen->generateIndependent(independent);
            }

            _ss.clear();
            _ss.str("");
            _ss << _multName << "[" << (id - _minMultiplierID) << "]";
            return _ss.str();
        }

        virtual std::string generateTemporary(const OperationNode<Base>& variable) {
            return _nameGen->generateTemporary(variable);
        }

        virtual std::string generateTemporaryArray(const OperationNode<Base>& variable) {
            return _nameGen->generateTemporaryArray(variable);
        }

        virtual std::string generateIndexedDependent(const OperationNode<Base>& var,
                                                     const IndexPattern& ip) {
            return _nameGen->generateIndexedDependent(var, ip);
        }

        virtual std::string generateIndexedIndependent(const OperationNode<Base>& indexedIndep,
                                                       const IndexPattern& ip) {
            bool isX = indexedIndep.getInfo()[0] == 0;
            if (isX) {
                return _nameGen->generateIndexedIndependent(indexedIndep, ip);
            }

            CPPADCG_ASSERT_KNOWN(indexedIndep.getOperationType() == CGLoopIndexedIndepOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(indexedIndep.getArguments().size() > 0, "Invalid number of arguments");
            CPPADCG_ASSERT_KNOWN(indexedIndep.getArguments()[0].getOperation() != NULL, "Invalid argument");
            CPPADCG_ASSERT_KNOWN(indexedIndep.getArguments()[0].getOperation()->getOperationType() == CGIndexOp, "Invalid argument");
            const IndexOperationNode<Base>& index = static_cast<const IndexOperationNode<Base>&> (*indexedIndep.getArguments()[0].getOperation());

            _ss.clear();
            _ss.str("");

            _ss << _multName << "[" << CLanguage<Base>::indexPattern2String(ip, index.getIndex()) << "]";
            return _ss.str();
        }

        virtual void setTemporaryVariableID(size_t minTempID, size_t maxTempID, size_t maxTempArrayID) {
            _nameGen->setTemporaryVariableID(minTempID, maxTempID, maxTempArrayID);
        }

        inline virtual ~CLangDefaultHessianVarNameGenerator() {
        }

    private:

        inline void initialize() {
            this->_independent = _nameGen->getIndependent(); // copy

            this->_independent.push_back(FuncArgument(_multName));
        }

    };
}

#endif
