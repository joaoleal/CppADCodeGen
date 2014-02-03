#ifndef CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
#define CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
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
 * Function arguments
 */
typedef struct FuncArgument {

    std::string name;
    bool array;

    inline FuncArgument() :
        array(false) {
    }

    inline FuncArgument(const std::string& nam, bool a = true) :
        name(nam),
        array(a) {
    }
} FuncArgument;

/**
 * Creates variables names for the source code.
 * 
 * @author Joao Leal
 */
template<class Base>
class VariableNameGenerator {

protected:
    std::vector<FuncArgument> _dependent;
    std::vector<FuncArgument> _independent;
    std::vector<FuncArgument> _temporary;
public:

    /**
     * Provides the dependent variable arguments used by a function.
     * 
     * @return the dependent variable arguments
     */
    virtual const std::vector<FuncArgument>& getDependent() const {
        return _dependent;
    }

    /**
     * Provides the independent variable arguments used by a function.
     * 
     * @return the independent variable arguments
     */
    virtual const std::vector<FuncArgument>& getIndependent() const {
        return _independent;
    }

    /**
     * Provides the temporary variable arguments used by a function.
     * 
     * @return the temporary variable arguments
     */
    virtual const std::vector<FuncArgument>& getTemporary() const {
        return _temporary;
    }

    /**
     * Provides the minimum variable ID of temporary variables.
     */
    virtual size_t getMinTemporaryVariableID() const = 0;

    /**
     * Provides the maximum used variable ID of temporary variables.
     */
    virtual size_t getMaxTemporaryVariableID() const = 0;

    /**
     * Provides the maximum variable ID of temporary dense array variables.
     */
    virtual size_t getMaxTemporaryArrayVariableID() const = 0;

    /**
     * Provides the maximum variable ID of temporary sparse array variables.
     */
    virtual size_t getMaxTemporarySparseArrayVariableID() const = 0;

    /**
     * Creates a name for a dependent variable.
     * 
     * @param index the dependent variable index
     * @return the generated name
     */
    virtual std::string generateDependent(size_t index) = 0;

    /**
     * Creates a name for a dependent variable.
     * 
     * @param variable the node representing the independent variable
     *                 (it should have an ID defined)
     * @return the generated name
     */
    virtual std::string generateIndependent(const OperationNode<Base>& variable) = 0;

    /**
     * Creates a name for a temporary variable.
     * 
     * @param variable the node representing the temporary variable
     *                 (it should have an ID defined)
     * @return the generated name
     */
    virtual std::string generateTemporary(const OperationNode<Base>& variable) = 0;

    /**
     * Creates a name for a temporary dense array variable.
     * 
     * @param variable the node representing the dense array variable
     *                 creation (it should have an ID defined)
     * @return the generated name
     */
    virtual std::string generateTemporaryArray(const OperationNode<Base>& variable) = 0;

    /**
     * Creates a name for a temporary sparse array variable.
     * 
     * @param variable the node representing the sparse array variable
     *                 creation (it should have an ID defined)
     * @return the generated name
     */
    virtual std::string generateTemporarySparseArray(const OperationNode<Base>& variable) = 0;

    /**
     * Creates a name for a reference to an indexed dependent variable 
     * expression.
     * 
     * @param var the node representing an indexed dependent variable
     * @param ip the index pattern
     * @return the generated name
     */
    virtual std::string generateIndexedDependent(const OperationNode<Base>& var,
                                                 const IndexPattern& ip) = 0;

    /**
     * Creates a name for a reference to an indexed independent variable 
     * expression.
     * 
     * @param var the node representing an indexed independent variable
     * @param ip the index pattern
     * @return the generated name
     */
    virtual std::string generateIndexedIndependent(const OperationNode<Base>& var,
                                                   const IndexPattern& ip) = 0;

    /**
     * Defines the ID ranges used by each variable type.
     * 
     * @param minTempID the lowest ID of temporary variables
     * @param maxTempID the highest used ID of temporary variables
     * @param maxTempArrayID the highest used ID of temporary dense array
     *                       variables
     * @param maxTempSparseArrayID the highest used ID of temporary sparse
     *                             array variables
     */
    virtual void setTemporaryVariableID(size_t minTempID,
                                        size_t maxTempID,
                                        size_t maxTempArrayID,
                                        size_t maxTempSparseArrayID) = 0;

    /**
     * Provides the name of the array where the provided independent.
     * It should only be called if independents are saved in an array.
     * 
     * @param indep the independent variable node (CGInvOp)
     * @return the array name
     */
    virtual const std::string& getIndependentArrayName(const OperationNode<Base>& indep) = 0;

    /**
     * Provides the index in the associated independent array of an 
     * independent variable.
     * It should only be called if independents are saved in an array.
     * 
     * @param indep the independent variable node (CGInvOp)
     * @return the index
     */
    virtual size_t getIndependentArrayIndex(const OperationNode<Base>& indep) = 0;

    /**
     * Whether or not two independent variables are considered to be part of
     * the same independent variable array at consecutive locations.
     * 
     * @param indepFirst the independent node (CGInvOp) with the lower index
     * @param indepSecond the independent node (CGInvOp) with the higher index
     * @return true if the independents are consecutive
     */
    virtual bool isConsecutiveInIndepArray(const OperationNode<Base>& indepFirst,
                                           const OperationNode<Base>& indepSecond) = 0;

    /**
     * Determines whether or not two independents are part of the same
     * independent variable array.
     * 
     * @param indep1 the first independent node (CGInvOp or CGLoopIndexedIndepOp)
     * @param indep2 the second independent node (CGInvOp or CGLoopIndexedIndepOp)
     * @return true if the independents are part of the same array
     */
    virtual bool isInSameIndependentArray(const OperationNode<Base>& indep1,
                                          const OperationNode<Base>& indep2) = 0;

    virtual void customFunctionVariableDeclarations(std::ostream& out) {
    }

    virtual void prepareCustomFunctionVariables(std::ostream& out) {
    }

    virtual void finalizeCustomFunctionVariables(std::ostream& out) {
    }

    inline virtual ~VariableNameGenerator() {
    }
};

} // END cg namespace
} // END CppAD namespace

#endif