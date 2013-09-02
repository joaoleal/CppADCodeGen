#ifndef CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
#define CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

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

        virtual const std::vector<FuncArgument>& getDependent() const {
            return _dependent;
        }

        virtual const std::vector<FuncArgument>& getIndependent() const {
            return _independent;
        }

        virtual const std::vector<FuncArgument>& getTemporary() const {
            return _temporary;
        }

        virtual size_t getMinTemporaryVariableID() const = 0;

        virtual size_t getMaxTemporaryVariableID() const = 0;

        virtual size_t getMaxTemporaryArrayVariableID() const = 0;

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) = 0;

        virtual std::string generateIndependent(const OperationNode<Base>& variable) = 0;

        virtual std::string generateTemporary(const OperationNode<Base>& variable) = 0;

        virtual std::string generateTemporaryArray(const OperationNode<Base>& variable) = 0;

        virtual std::string generateIndexedDependent(const OperationNode<Base>& var,
                                                     const IndexPattern& ip) = 0;

        virtual std::string generateIndexedIndependent(const OperationNode<Base>& var,
                                                       const IndexPattern& ip) = 0;

        virtual void setTemporaryVariableID(size_t minTempID, size_t maxTempID, size_t maxTempArrayID) = 0;

        virtual void customFunctionVariableDeclarations(std::ostream& out) {
        }

        virtual void prepareCustomFunctionVariables(std::ostream& out) {
        }

        virtual void finalizeCustomFunctionVariables(std::ostream& out) {
        }

        inline virtual ~VariableNameGenerator() {
        }
    };
}

#endif
