#ifndef CPPAD_CG_LANG_C_DEFAULT_HESSIAN_VAR_NAME_GEN_INCLUDED
#define CPPAD_CG_LANG_C_DEFAULT_HESSIAN_VAR_NAME_GEN_INCLUDED
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
namespace cg {

/**
 * Creates variables names for the source code generated for Hessian
 * calculations.
 * The independent variables are considered to have been registered first as
 * variable in the code generation handler and then the multipliers.
 * 
 * @author Joao Leal
 */
template<class Base>
class LangCDefaultHessianVarNameGenerator : public VariableNameGenerator<Base> {
protected:
    VariableNameGenerator<Base>* _nameGen;
    // the lowest variable ID used for the equation multipliers
    const size_t _minMultiplierID;
    // array name of the independent variables
    const std::string _multName;
    // auxiliary string stream
    std::stringstream _ss;
public:

    LangCDefaultHessianVarNameGenerator(VariableNameGenerator<Base>* nameGen,
                                        size_t n) :
        _nameGen(nameGen),
        _minMultiplierID(n + 1),
        _multName("mult") {

        CPPADCG_ASSERT_KNOWN(_nameGen != nullptr, "The name generator must not be NULL");

        initialize();
    }

    LangCDefaultHessianVarNameGenerator(VariableNameGenerator<Base>* nameGen,
                                        const std::string& multName,
                                        size_t n) :
        _nameGen(nameGen),
        _minMultiplierID(n + 1),
        _multName(multName) {

        CPPADCG_ASSERT_KNOWN(_nameGen != nullptr, "The name generator must not be null");
        CPPADCG_ASSERT_KNOWN(_multName.size() > 0, "The name for the multipliers must not be empty");

        initialize();
    }

    virtual const std::vector<FuncArgument>& getDependent() const override {
        return _nameGen->getDependent();
    }

    virtual const std::vector<FuncArgument>& getTemporary() const override {
        return _nameGen->getTemporary();
    }

    virtual size_t getMinTemporaryVariableID() const override {
        return _nameGen->getMinTemporaryVariableID();
    }

    virtual size_t getMaxTemporaryVariableID() const override {
        return _nameGen->getMaxTemporaryVariableID();
    }

    virtual size_t getMaxTemporaryArrayVariableID() const override {
        return _nameGen->getMaxTemporaryArrayVariableID();
    }

    virtual size_t getMaxTemporarySparseArrayVariableID() const override {
        return _nameGen->getMaxTemporarySparseArrayVariableID();
    }

    virtual std::string generateDependent(size_t index) override {
        return _nameGen->generateDependent(index);
    }

    virtual std::string generateIndependent(const OperationNode<Base>& independent) override {
        size_t id = independent.getVariableID();
        if (id < _minMultiplierID) {
            return _nameGen->generateIndependent(independent);
        }

        _ss.clear();
        _ss.str("");
        _ss << _multName << "[" << (id - _minMultiplierID) << "]";
        return _ss.str();
    }

    virtual std::string generateTemporary(const OperationNode<Base>& variable) override {
        return _nameGen->generateTemporary(variable);
    }

    virtual std::string generateTemporaryArray(const OperationNode<Base>& variable) override {
        return _nameGen->generateTemporaryArray(variable);
    }

    virtual std::string generateTemporarySparseArray(const OperationNode<Base>& variable) override {
        return _nameGen->generateTemporarySparseArray(variable);
    }

    virtual std::string generateIndexedDependent(const OperationNode<Base>& var,
                                                 const IndexPattern& ip) override {
        return _nameGen->generateIndexedDependent(var, ip);
    }

    virtual std::string generateIndexedIndependent(const OperationNode<Base>& indexedIndep,
                                                   const IndexPattern& ip) override {
        bool isX = indexedIndep.getInfo()[0] == 0;
        if (isX) {
            return _nameGen->generateIndexedIndependent(indexedIndep, ip);
        }

        CPPADCG_ASSERT_KNOWN(indexedIndep.getOperationType() == CGLoopIndexedIndepOp, "Invalid node type");
        CPPADCG_ASSERT_KNOWN(indexedIndep.getArguments().size() > 0, "Invalid number of arguments");
        CPPADCG_ASSERT_KNOWN(indexedIndep.getArguments()[0].getOperation() != nullptr, "Invalid argument");
        CPPADCG_ASSERT_KNOWN(indexedIndep.getArguments()[0].getOperation()->getOperationType() == CGIndexOp, "Invalid argument");
        const IndexOperationNode<Base>& index = static_cast<const IndexOperationNode<Base>&> (*indexedIndep.getArguments()[0].getOperation());

        _ss.clear();
        _ss.str("");

        _ss << _multName << "[" << LanguageC<Base>::indexPattern2String(ip, index.getIndex()) << "]";
        return _ss.str();
    }

    virtual const std::string& getIndependentArrayName(const OperationNode<Base>& indep) override {
        if (indep.getVariableID() < _minMultiplierID)
            return _nameGen->getIndependentArrayName(indep);
        else
            return _multName;
    }

    virtual size_t getIndependentArrayIndex(const OperationNode<Base>& indep) override {
        if (indep.getVariableID() < _minMultiplierID)
            return _nameGen->getIndependentArrayIndex(indep);
        else
            return indep.getVariableID() - _minMultiplierID;
    }

    virtual bool isConsecutiveInIndepArray(const OperationNode<Base>& indepFirst,
                                           const OperationNode<Base>& indepSecond) override {
        size_t id1 = indepFirst.getVariableID();
        size_t id2 = indepSecond.getVariableID();

        if ((id1 < _minMultiplierID) != (id2 < _minMultiplierID))
            return false;

        if (id1 < _minMultiplierID && id2 < _minMultiplierID)
            return _nameGen->isConsecutiveInIndepArray(indepFirst, indepSecond);
        else
            return id1 + 1 == id2;
    }

    virtual bool isInSameIndependentArray(const OperationNode<Base>& indep1,
                                          const OperationNode<Base>& indep2) override {
        size_t l1;
        if (indep1.getOperationType() == CGInvOp) {
            l1 = indep1.getVariableID() < _minMultiplierID ? 0 : 1;
        } else {
            l1 = indep1.getInfo()[0]; //CGLoopIndexedIndepOp
        }

        size_t l2;
        if (indep2.getOperationType() == CGInvOp) {
            l2 = indep2.getVariableID() < _minMultiplierID ? 0 : 1;
        } else {
            l2 = indep2.getInfo()[0]; //CGLoopIndexedIndepOp
        }

        return l1 == l2;
    }

    virtual void setTemporaryVariableID(size_t minTempID,
                                        size_t maxTempID,
                                        size_t maxTempArrayID,
                                        size_t maxTempSparseArrayID) override {
        _nameGen->setTemporaryVariableID(minTempID, maxTempID, maxTempArrayID, maxTempSparseArrayID);
    }

    inline virtual ~LangCDefaultHessianVarNameGenerator() {
    }

private:

    inline void initialize() {
        this->_independent = _nameGen->getIndependent(); // copy

        this->_independent.push_back(FuncArgument(_multName));
    }

};

} // END cg namespace
} // END CppAD namespace

#endif