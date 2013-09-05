#ifndef CPPAD_CG_C_LANGUAGE_INCLUDED
#define CPPAD_CG_C_LANGUAGE_INCLUDED
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

#define CPPAD_CG_C_LANG_FUNCNAME(fn) \
inline virtual const std::string& fn ## FuncName() {\
    static const std::string name(#fn);\
    return name;\
}

namespace CppAD {

    /**
     * Generates code for the C language
     * 
     * @author Joao Leal
     */
    template<class Base>
    class CLanguage : public Language<Base> {
    public:
        static const std::string ATOMICFUN_STRUCT_DEFINITION;
    protected:
        static const std::string _C_COMP_OP_LT;
        static const std::string _C_COMP_OP_LE;
        static const std::string _C_COMP_OP_EQ;
        static const std::string _C_COMP_OP_GE;
        static const std::string _C_COMP_OP_GT;
        static const std::string _C_COMP_OP_NE;

    protected:
        // the type name of the Base class (e.g. "double")
        const std::string _baseTypeName;
        // spaces for 1 level indentation
        const std::string _spaces;
        // currrent identation
        std::string _indentation;
        // variable name used for the inlet variable
        std::string _inArgName;
        // variable name used for the oulet variable
        std::string _outArgName;
        // variable name used for the atomic functions array
        std::string _atomicArgName;
        // output stream for the generated source code
        std::ostringstream _code;
        // creates the variable names
        VariableNameGenerator<Base>* _nameGen;
        // auxiliary string stream
        std::stringstream _ss;
        //
        size_t _independentSize;
        //
        size_t _minTemporaryVarID;
        // maps atomic function IDs to their internal index
        std::map<size_t, size_t> _atomicFunctionId2Index;
        // maps atomic function IDs to their names
        std::map<size_t, std::string> _atomicFunctionId2Name;
        // maps the variable IDs to the their position in the dependent vector
        // (some IDs may be the same as the independent variables when dep = indep)
        std::map<size_t, size_t> _dependentIDs;
        // the dependent variable vector
        const std::vector<CG<Base> >* _dependent;
        // the temporary variables that may require a declaration
        std::map<size_t, OperationNode<Base>*> _temporary;
        // the operator used for assignment of dependent variables
        std::string _depAssignOperation;
        // whether or not to ignore assignment of constant zero values to dependent variables
        bool _ignoreZeroDepAssign;
        // the name of the function to be created (if the string is empty no function is created)
        std::string _functionName;
        // the arguments provided to local functions called by the main function
        std::string _localFunctionArguments;
        // the maximum number of assignment (~lines) per local function
        size_t _maxAssigmentsPerFunction;
        //
        std::map<std::string, std::string>* _sources;
        // the values in the temporary array
        std::vector<const Argument<Base>*> _tmpArrayValues;
        const std::set<const Index*>* _indexes;
        const std::vector<const IndexPattern*>* _loopDependentIndexPatterns;
        const std::vector<const IndexPattern*>* _loopIndependentIndexPatterns;
        std::vector<const LoopNodeInfo<Base>*> _currentLoopsUserDefined;
        std::vector<const LoopNodeInfo<Base>*> _currentLoops;
    private:
        std::string defaultFuncArgDcl_;
        std::string localFuncArgDcl_;
        std::string localFuncArgs_;
        std::string auxArrayName_;

    public:

        /**
         * Creates a C language source code generator
         * 
         * @param varTypeName variable data type (e.g. double)
         * @param spaces number of spaces for indentations
         */
        CLanguage(const std::string& varTypeName, size_t spaces = 3) :
            _baseTypeName(varTypeName),
            _spaces(spaces, ' '),
            _inArgName("in"),
            _outArgName("out"),
            _atomicArgName("atomicFun"),
            _nameGen(NULL),
            _independentSize(0),
            _dependent(NULL),
            _depAssignOperation("="),
            _ignoreZeroDepAssign(false),
            _maxAssigmentsPerFunction(0),
            _sources(NULL),
            _indexes(NULL),
            _loopDependentIndexPatterns(NULL),
            _loopIndependentIndexPatterns(NULL) {
        }

        inline const std::string& getArgumentIn() const {
            return _inArgName;
        }

        inline void setArgumentIn(const std::string& inArgName) {
            _inArgName = inArgName;
        }

        inline const std::string& getArgumentOut() const {
            return _outArgName;
        }

        inline void setArgumentOut(const std::string& outArgName) {
            _outArgName = outArgName;
        }

        inline const std::string& getArgumentAtomic() const {
            return _atomicArgName;
        }

        inline void setArgumentAtomic(const std::string& atomicArgName) {
            _atomicArgName = atomicArgName;
        }

        inline const std::string& getDependentAssignOperation() const {
            return _depAssignOperation;
        }

        inline void setDependentAssignOperation(const std::string& depAssignOperation) {
            _depAssignOperation = depAssignOperation;
        }

        inline bool isIgnoreZeroDepAssign() const {
            return _depAssignOperation;
        }

        inline void setIgnoreZeroDepAssign(bool ignore) {
            _ignoreZeroDepAssign = ignore;
        }

        virtual void setGenerateFunction(const std::string& functionName) {
            _functionName = functionName;
        }

        virtual void setMaxAssigmentsPerFunction(size_t maxAssigmentsPerFunction, std::map<std::string, std::string>* sources) {
            _maxAssigmentsPerFunction = maxAssigmentsPerFunction;
            _sources = sources;
        }

        virtual void setCurrentLoop(LoopAtomicFun<Base>* loop) {
            if (loop == NULL) {
                _currentLoopsUserDefined.clear();
            } else {
                _currentLoopsUserDefined.resize(1);
                _currentLoopsUserDefined[0] = loop;
            }
        }

        virtual std::string generateTemporaryVariableDeclaration(bool localOnly = false) {
            assert(_nameGen != NULL);

            // declare variables
            const std::vector<FuncArgument>& tmpArg = _nameGen->getTemporary();

            CPPADCG_ASSERT_KNOWN(tmpArg.size() == 2,
                                 "There must be two temporary variables");

            _ss << _spaces << "// auxiliary variables\n";
            /**
             * temporary variables
             */
            if (tmpArg[0].array) {
                size_t size = _nameGen->getMaxTemporaryVariableID() + 1 - _nameGen->getMinTemporaryVariableID();
                if (size > 0 || !localOnly) {
                    _ss << _spaces << _baseTypeName << " " << tmpArg[0].name << "[" << size << "];\n";
                }
            } else if (_temporary.size() > 0) {
                typename std::map<size_t, OperationNode<Base>*>::const_iterator it;

                for (it = _temporary.begin(); it != _temporary.end(); ++it) {
                    OperationNode<Base>* var = it->second;
                    if (var->getName() == NULL) {
                        var->setName(_nameGen->generateTemporary(*var));
                    }
                }

                OperationNode<Base>* var1 = _temporary.begin()->second;
                const std::string& varName1 = *var1->getName();
                _ss << _spaces << _baseTypeName << " " << varName1;

                it = _temporary.begin();
                for (it++; it != _temporary.end(); ++it) {
                    _ss << ", " << *it->second->getName();
                }
                _ss << ";\n";
            }

            /**
             * temporary array variables
             */
            size_t arraySize = _nameGen->getMaxTemporaryArrayVariableID();
            if (arraySize > 0 || !localOnly) {
                _ss << _spaces << _baseTypeName << " " << tmpArg[1].name << "[" << arraySize << "];\n";
                if (localOnly && arraySize > 0) {
                    _ss << _spaces << _baseTypeName << "* " << auxArrayName_ << ";\n";
                }
                if (arraySize > 0) {
                    _ss << _spaces << "unsigned long i;\n";
                }
            }

            // loop indexes
            createIndexDeclaration();

            // clean-up
            std::string code = _ss.str();
            _ss.str("");

            return code;
        }

        virtual std::string generateDependentVariableDeclaration() {
            const std::vector<FuncArgument>& depArg = _nameGen->getDependent();
            CPPADCG_ASSERT_KNOWN(depArg.size() > 0,
                                 "There must be at least one dependent argument");

            _ss << _spaces << "//dependent variables\n";
            for (size_t i = 0; i < depArg.size(); i++) {
                _ss << _spaces << argumentDeclaration(depArg[i]) << " = " << _outArgName << "[" << i << "];\n";
            }

            std::string code = _ss.str();
            _ss.str("");
            return code;
        }

        virtual std::string generateIndependentVariableDeclaration() {
            const std::vector<FuncArgument>& indArg = _nameGen->getIndependent();
            CPPADCG_ASSERT_KNOWN(indArg.size() > 0,
                                 "There must be at least one independent argument");

            _ss << _spaces << "//independent variables\n";
            for (size_t i = 0; i < indArg.size(); i++) {
                _ss << _spaces << "const " << argumentDeclaration(indArg[i]) << " = " << _inArgName << "[" << i << "];\n";
            }

            std::string code = _ss.str();
            _ss.str("");
            return code;
        }

        inline std::string generateArgumentAtomicDcl() const {
            return "struct CLangAtomicFun " + _atomicArgName;
        }

        virtual std::string generateDefaultFunctionArgumentsDcl() const {
            return _baseTypeName + " const *const * " + _inArgName +
                    ", " + _baseTypeName + "*const * " + _outArgName +
                    ", " + generateArgumentAtomicDcl();
        }

        virtual std::string generateDefaultFunctionArguments() const {
            return _inArgName + ", " + _outArgName + ", " + _atomicArgName;
        }

        inline void createIndexDeclaration();

        CPPAD_CG_C_LANG_FUNCNAME(abs)
        CPPAD_CG_C_LANG_FUNCNAME(acos)
        CPPAD_CG_C_LANG_FUNCNAME(asin)
        CPPAD_CG_C_LANG_FUNCNAME(atan)
        CPPAD_CG_C_LANG_FUNCNAME(cosh)
        CPPAD_CG_C_LANG_FUNCNAME(cos)
        CPPAD_CG_C_LANG_FUNCNAME(exp)
        CPPAD_CG_C_LANG_FUNCNAME(log)
        CPPAD_CG_C_LANG_FUNCNAME(sinh)
        CPPAD_CG_C_LANG_FUNCNAME(sin)
        CPPAD_CG_C_LANG_FUNCNAME(sqrt)
        CPPAD_CG_C_LANG_FUNCNAME(tanh)
        CPPAD_CG_C_LANG_FUNCNAME(tan)
        CPPAD_CG_C_LANG_FUNCNAME(pow)

        inline virtual ~CLanguage() {
        }

        /***********************************************************************
         * index patterns
         **********************************************************************/
        static inline std::string createIndexPattern(const IndexPattern& ip);

        static inline std::string createLinearIndexPattern(const LinearIndexPattern& lip);

    protected:

        virtual void generateSourceCode(std::ostream& out, LanguageGenerationData<Base>& info) {
            const bool createFunction = !_functionName.empty();
            const bool multiFunction = createFunction && _maxAssigmentsPerFunction > 0 && _sources != NULL;

            // clean up
            _code.str("");
            _ss.str("");
            _temporary.clear();
            _indentation = _spaces;
            defaultFuncArgDcl_ = "";
            localFuncArgDcl_ = "";
            localFuncArgs_ = "";
            auxArrayName_ = "";
            _currentLoops = _currentLoopsUserDefined;

            // save some info
            _independentSize = info.independent.size();
            _dependent = &info.dependent;
            _nameGen = &info.nameGen;
            _minTemporaryVarID = info.minTemporaryVarID;
            _atomicFunctionId2Index = info.atomicFunctionId2Index;
            _atomicFunctionId2Name = info.atomicFunctionId2Name;
            const std::vector<CG<Base> >& dependent = info.dependent;
            const std::vector<OperationNode<Base>*>& variableOrder = info.variableOrder;
            _tmpArrayValues.resize(_nameGen->getMaxTemporaryArrayVariableID());
            std::fill(_tmpArrayValues.begin(), _tmpArrayValues.end(), (Argument<Base>*) NULL);
            _indexes = &info.indexes;
            _loopDependentIndexPatterns = &info.loopDependentIndexPatterns;
            _loopIndependentIndexPatterns = &info.loopIndependentIndexPatterns;

            /**
             * generate variable names
             */
            //generate names for the independent variables
            typename std::vector<OperationNode<Base> *>::const_iterator it;
            for (it = info.independent.begin(); it != info.independent.end(); ++it) {
                OperationNode<Base>& op = **it;
                if (op.getName() == NULL) {
                    op.setName(_nameGen->generateIndependent(op));
                }
            }

            // generate names for the dependent variables (must be after naming the independent)
            for (size_t i = 0; i < dependent.size(); i++) {
                OperationNode<Base>* node = dependent[i].getOperationNode();
                if (node != NULL && node->getOperationType() != CGLoopEndOp && node->getName() == NULL) {
                    if (node->getOperationType() == CGLoopIndexedDepOp) {
                        assert(!_currentLoops.empty());
                        size_t pos = node->getInfo()[0];
                        const IndexPattern* ip = (*_loopDependentIndexPatterns)[pos];
                        node->setName(_nameGen->generateIndexedDependent(*node, *ip));

                    } else {
                        node->setName(_nameGen->generateDependent(dependent[i], i));
                    }
                }
            }

            /**
             * function variable declaration
             */
            const std::vector<FuncArgument>& indArg = _nameGen->getIndependent();
            const std::vector<FuncArgument>& depArg = _nameGen->getDependent();
            const std::vector<FuncArgument>& tmpArg = _nameGen->getTemporary();
            CPPADCG_ASSERT_KNOWN(indArg.size() > 0 && depArg.size() > 0,
                                 "There must be at least one dependent and one independent argument");
            CPPADCG_ASSERT_KNOWN(tmpArg.size() == 2,
                                 "There must be two temporary variables");

            if (createFunction) {
                defaultFuncArgDcl_ = generateDefaultFunctionArgumentsDcl();
                localFuncArgDcl_ = defaultFuncArgDcl_ + ", " + argumentDeclaration(tmpArg[0]) + ", " + argumentDeclaration(tmpArg[1]);
                localFuncArgs_ = generateDefaultFunctionArguments() + ", " + tmpArg[0].name + ", " + tmpArg[1].name;
            }

            auxArrayName_ = tmpArg[1].name + "p";

            /**
             * Determine the dependent variables that result from the same operations
             */
            // dependent variables indexes that are copies of other dependent variables
            std::set<size_t> dependentDuplicates;

            for (size_t i = 0; i < dependent.size(); i++) {
                OperationNode<Base>* node = dependent[i].getOperationNode();
                if (node != NULL) {
                    CGOpCode type = node->getOperationType();
                    if (type != CGInvOp && type != CGLoopEndOp) {
                        size_t varID = node->getVariableID();
                        if (varID > 0) {
                            std::map<size_t, size_t>::const_iterator it2 = _dependentIDs.find(varID);
                            if (it2 == _dependentIDs.end()) {
                                _dependentIDs[node->getVariableID()] = i;
                            } else {
                                // there can be several dependent variables with the same ID
                                dependentDuplicates.insert(i);
                            }
                        }
                    }
                }
            }

            // the names of local functions
            std::vector<std::string> localFuncNames;
            if (multiFunction) {
                localFuncNames.reserve(variableOrder.size() / _maxAssigmentsPerFunction);
            }

            /**
             * non-constant variables
             */
            if (variableOrder.size() > 0) {
                // generate names for temporary variables
                for (it = variableOrder.begin(); it != variableOrder.end(); ++it) {
                    OperationNode<Base>& op = **it;
                    if (op.getName() == NULL) {
                        if (!isDependent(op)) {
                            if (requiresVariableName(op) && op.getOperationType() != CGArrayCreationOp) {
                                op.setName(_nameGen->generateTemporary(op));
                            } else if (op.getOperationType() == CGArrayCreationOp) {
                                op.setName(_nameGen->generateTemporaryArray(op));
                            }
                        }
                    }
                }

                /**
                 * Source code generation magic!
                 */
                size_t assignCount = 0;
                for (it = variableOrder.begin(); it != variableOrder.end(); ++it) {
                    // check if a new function should start
                    if (assignCount >= _maxAssigmentsPerFunction && multiFunction && _currentLoops.empty()) {
                        assignCount = 0;
                        saveLocalFunction(localFuncNames);
                    }

                    OperationNode<Base>& op = **it;

                    // a dependent variable assigned by a loop does require any source code (its done inside the loop)
                    if (op.getOperationType() == CGDependentRefOp) {
                        continue; // nothing to do (this operation is right hand side only)
                    }

                    bool createsVar = directlyAssignsVariable(op); // do we need to do the assigment here?
                    if (!createsVar) {
                        printAssigmentStart(op);
                    }
                    printExpressionNoVarCheck(op);
                    if (!createsVar) {
                        printAssigmentEnd(op);
                    }

                    if (op.getOperationType() == CGArrayElementOp) {
                        size_t arrayId = op.getArguments()[0].getOperation()->getVariableID();
                        size_t pos = op.getInfo()[0];
                        _tmpArrayValues[arrayId - 1 + pos] = NULL; // this could probably be removed!
                    }

                    assignCount++;
                }

                if (localFuncNames.size() > 0 && assignCount > 0) {
                    assignCount = 0;
                    saveLocalFunction(localFuncNames);
                }
            }

            if (localFuncNames.size() > 0) {
                CPPADCG_ASSERT_KNOWN(tmpArg[0].array,
                                     "The temporary variables must be saved in an array in order to generate multiple functions");

                _code << ATOMICFUN_STRUCT_DEFINITION << "\n\n";
                // forward declarations
                for (size_t i = 0; i < localFuncNames.size(); i++) {
                    _code << "void " << localFuncNames[i] << "(" << localFuncArgDcl_ << ");\n";
                }
                _code << "\n"
                        << "void " << _functionName << "(" << defaultFuncArgDcl_ << ") {\n";
                _nameGen->customFunctionVariableDeclarations(_code);
                _code << generateIndependentVariableDeclaration() << "\n";
                _code << generateDependentVariableDeclaration() << "\n";
                _code << generateTemporaryVariableDeclaration(false) << "\n";
                _nameGen->prepareCustomFunctionVariables(_code);
                for (size_t i = 0; i < localFuncNames.size(); i++) {
                    _code << _spaces << localFuncNames[i] << "(" << localFuncArgs_ << ");\n";
                }
            }

            // dependent duplicates
            if (dependentDuplicates.size() > 0) {
                _code << _spaces << "// variable duplicates: " << dependentDuplicates.size() << "\n";
                for (std::set<size_t>::const_iterator it = dependentDuplicates.begin(); it != dependentDuplicates.end(); ++it) {
                    size_t index = *it;
                    const CG<Base>& dep = (*_dependent)[index];
                    std::string varName = _nameGen->generateDependent(dep, index);
                    const std::string& origVarName = *dep.getOperationNode()->getName();

                    _code << _spaces << varName << " " << _depAssignOperation << " " << origVarName << ";\n";
                }
            }

            // constant dependent variables 
            bool commentWritten = false;
            for (size_t i = 0; i < dependent.size(); i++) {
                if (dependent[i].isParameter()) {
                    if (!_ignoreZeroDepAssign || !dependent[i].IdenticalZero()) {
                        if (!commentWritten) {
                            _code << _spaces << "// dependent variables without operations\n";
                            commentWritten = true;
                        }
                        std::string varName = _nameGen->generateDependent(dependent[i], i);
                        _code << _spaces << varName << " " << _depAssignOperation << " ";
                        printParameter(dependent[i].getValue());
                        _code << ";\n";
                    }
                } else if (dependent[i].getOperationNode()->getOperationType() == CGInvOp) {
                    if (!commentWritten) {
                        _code << _spaces << "// dependent variables without operations\n";
                        commentWritten = true;
                    }
                    std::string varName = _nameGen->generateDependent(dependent[i], i);
                    const std::string& indepName = *dependent[i].getOperationNode()->getName();
                    _code << _spaces << varName << " " << _depAssignOperation << " " << indepName << ";\n";
                }
            }

            /**
             * encapsulate the code in a function
             */
            if (createFunction) {
                if (localFuncNames.empty()) {
                    _ss << "#include <math.h>\n\n"
                            << ATOMICFUN_STRUCT_DEFINITION << "\n\n"
                            << "void " << _functionName << "(" << defaultFuncArgDcl_ << ") {\n";
                    _nameGen->customFunctionVariableDeclarations(_ss);
                    _ss << generateIndependentVariableDeclaration() << "\n";
                    _ss << generateDependentVariableDeclaration() << "\n";
                    _ss << generateTemporaryVariableDeclaration(true) << "\n";
                    _nameGen->prepareCustomFunctionVariables(_ss);
                    _ss << _code.str();
                    _nameGen->finalizeCustomFunctionVariables(_ss);
                    _ss << "}\n\n";

                    out << _ss.str();

                    if (_sources != NULL) {
                        (*_sources)[_functionName + ".c"] = _ss.str();
                    }
                } else {
                    _nameGen->finalizeCustomFunctionVariables(_code);
                    _code << "}\n\n";

                    (*_sources)[_functionName + ".c"] = _code.str();
                }
            } else {
                out << _code.str();
            }
        }

        inline virtual void printAssigmentStart(OperationNode<Base>& op) {
            printAssigmentStart(op, createVariableName(op), isDependent(op));
        }

        inline virtual void printAssigmentStart(OperationNode<Base>& op, const std::string& varName, bool isDep) {
            if (!isDep) {
                _temporary[op.getVariableID()] = &op;
            }

            _code << _indentation << varName << " ";
            if (isDep) {
                if (op.getOperationType() == CGLoopIndexedDepOp) {
                    if (op.getInfo()[1] == 0) {
                        _code << "=";
                    } else {
                        _code << "+=";
                    }
                } else {
                    _code << _depAssignOperation;
                }
            } else {
                _code << "=";
            }
            _code << " ";
        }

        inline virtual void printAssigmentEnd(OperationNode<Base>& op) {
            _code << ";\n";
        }

        virtual std::string argumentDeclaration(const FuncArgument& funcArg) const {
            std::string dcl = _baseTypeName;
            if (funcArg.array) {
                dcl += "*";
            }
            return dcl + " " + funcArg.name;
        }

        virtual void saveLocalFunction(std::vector<std::string>& localFuncNames) {
            _ss << _functionName << "__" << (localFuncNames.size() + 1);
            std::string funcName = _ss.str();
            _ss.str("");

            _ss << "#include <math.h>\n\n"
                    << ATOMICFUN_STRUCT_DEFINITION << "\n\n"
                    << "void " << funcName << "(" << localFuncArgDcl_ << ") {\n";
            _nameGen->customFunctionVariableDeclarations(_ss);
            _ss << generateIndependentVariableDeclaration() << "\n";
            _ss << generateDependentVariableDeclaration() << "\n";
            size_t arraySize = _nameGen->getMaxTemporaryArrayVariableID() - 1;
            if (arraySize > 0) {
                _ss << _spaces << _baseTypeName << "* " << auxArrayName_ << ";\n";
                _ss << _spaces << "unsigned long i;\n";
            }

            // loop indexes
            createIndexDeclaration();

            _nameGen->prepareCustomFunctionVariables(_ss);
            _ss << _code.str();
            _nameGen->finalizeCustomFunctionVariables(_ss);
            _ss << "}\n\n";

            (*_sources)[funcName + ".c"] = _ss.str();
            localFuncNames.push_back(funcName);

            _code.str("");
            _ss.str("");
        }

        virtual bool createsNewVariable(const OperationNode<Base>& var) const {
            CGOpCode op = var.getOperationType();
            if (var.getTotalUsageCount() > 1) {
                if (op == CGLoopAtomicResultOp) {
                    const OperationNode<Base>* opArg = var.getArguments()[0].getOperation();
                    return opArg != NULL && opArg->getOperationType() != CGLoopIndexedIndepOp
                            && opArg->getOperationType() != CGInvOp;
                } else {
                    return op != CGArrayElementOp && op != CGIndexOp;
                }
            } else {
                return op == CGArrayCreationOp ||
                        op == CGAtomicForwardOp ||
                        op == CGAtomicReverseOp ||
                        op == CGComOpLt ||
                        op == CGComOpLe ||
                        op == CGComOpEq ||
                        op == CGComOpGe ||
                        op == CGComOpGt ||
                        op == CGComOpNe ||
                        op == CGLoopIndexedDepOp ||
                        op == CGIndexAssignOp;
            }
        }

        virtual bool requiresVariableName(const OperationNode<Base>& var) const {
            CGOpCode op = var.getOperationType();
            return (var.getTotalUsageCount() > 1 &&
                    op != CGAtomicForwardOp &&
                    op != CGAtomicReverseOp &&
                    op != CGLoopStartOp &&
                    op != CGLoopEndOp &&
                    op != CGIndexOp &&
                    op != CGIndexAssignOp);
        }

        /**
         * Whether or not this operation assign its expression to a variable by
         * itself.
         * 
         * @param var the operation node
         * @return 
         */
        virtual bool directlyAssignsVariable(const OperationNode<Base>& var) const {
            CGOpCode op = var.getOperationType();
            return isCondAssign(op) ||
                    op == CGArrayCreationOp ||
                    op == CGAtomicForwardOp ||
                    op == CGAtomicReverseOp ||
                    op == CGLoopForwardOp ||
                    op == CGLoopReverseOp ||
                    op == CGLoopStartOp ||
                    op == CGLoopEndOp ||
                    op == CGIndexAssignOp;
        }

        virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) const {
            return op == CGSignOp;
        }

        inline const std::string& createVariableName(OperationNode<Base>& var) {
            CGOpCode op = var.getOperationType();
            assert(var.getVariableID() > 0);
            assert(op != CGAtomicForwardOp);
            assert(op != CGAtomicReverseOp);
            assert(op != CGLoopForwardOp);
            assert(op != CGLoopReverseOp);
            assert(op != CGLoopStartOp);
            assert(op != CGLoopEndOp);
            assert(op != CGIndexOp);
            assert(op != CGIndexAssignOp);

            if (var.getName() == NULL) {
                if (op == CGArrayCreationOp) {
                    var.setName(_nameGen->generateTemporaryArray(var));
                } else if (op == CGLoopIndexedDepOp) {
                    assert(!_currentLoops.empty());
                    size_t pos = var.getInfo()[0];
                    const IndexPattern* ip = (*_loopDependentIndexPatterns)[pos];
                    var.setName(_nameGen->generateIndexedDependent(var, *ip));

                } else if (op == CGLoopIndexedIndepOp) {
                    assert(!_currentLoops.empty());
                    size_t pos = var.getInfo()[1];
                    const IndexPattern* ip = (*_loopIndependentIndexPatterns)[pos];
                    var.setName(_nameGen->generateIndexedIndependent(var, *ip));

                } else {
                    if (var.getVariableID() <= _independentSize) {
                        // independent variable
                        var.setName(_nameGen->generateIndependent(var));

                    } else if (var.getVariableID() < _minTemporaryVarID) {
                        // dependent variable
                        std::map<size_t, size_t>::const_iterator it = _dependentIDs.find(var.getVariableID());
                        assert(it != _dependentIDs.end());

                        size_t index = it->second;
                        const CG<Base>& dep = (*_dependent)[index];
                        var.setName(_nameGen->generateDependent(dep, index));

                    } else {
                        // temporary variable
                        if (requiresVariableName(var)) {
                            var.setName(_nameGen->generateTemporary(var));
                        }
                    }
                }
            }

            return *var.getName();
        }

        virtual void printIndependentVariableName(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 0, "Invalid number of arguments for independent variable");

            _code << _nameGen->generateIndependent(op);
        }

        virtual void print(const Argument<Base>& arg) {
            if (arg.getOperation() != NULL) {
                // expression
                printExpression(*arg.getOperation());
            } else {
                // parameter
                printParameter(*arg.getParameter());
            }
        }

        virtual void printExpression(OperationNode<Base>& op) throw (CGException) {
            if (op.getVariableID() > 0) {
                // use variable name
                _code << createVariableName(op);
            } else {
                // print expression code
                printExpressionNoVarCheck(op);
            }
        }

        virtual void printExpressionNoVarCheck(OperationNode<Base>& node) throw (CGException) {
            CGOpCode op = node.getOperationType();
            switch (op) {
                case CGArrayCreationOp:
                    printArrayCreationOp(node);
                    break;
                case CGArrayElementOp:
                    printArrayElementOp(node);
                    break;
                case CGAbsOp:
                case CGAcosOp:
                case CGAsinOp:
                case CGAtanOp:
                case CGCoshOp:
                case CGCosOp:
                case CGExpOp:
                case CGLogOp:
                case CGSinhOp:
                case CGSinOp:
                case CGSqrtOp:
                case CGTanhOp:
                case CGTanOp:
                    printUnaryFunction(node);
                    break;
                case CGAtomicForwardOp: // atomicFunction.forward(q, p, vx, vy, tx, ty)
                    printAtomicForwardOp(node);
                    break;
                case CGAtomicReverseOp: // atomicFunction.reverse(p, tx, ty, px, py)
                    printAtomicReverseOp(node);
                    break;
                case CGAddOp:
                    printOperationAdd(node);
                    break;
                case CGAliasOp:
                    printOperationAlias(node);
                    break;
                case CGComOpLt:
                case CGComOpLe:
                case CGComOpEq:
                case CGComOpGe:
                case CGComOpGt:
                case CGComOpNe:
                    printConditionalAssignment(node);
                    break;
                case CGDivOp:
                    printOperationDiv(node);
                    break;
                case CGInvOp:
                    printIndependentVariableName(node);
                    break;
                case CGMulOp:
                    printOperationMul(node);
                    break;
                case CGPowOp:
                    printPowFunction(node);
                    break;
                case CGSignOp:
                    printSignFunction(node);
                    break;
                case CGSubOp:
                    printOperationMinus(node);
                    break;

                case CGUnMinusOp:
                    printOperationUnaryMinus(node);
                    break;

                case CGIndexOp:
                    break; // nothing to do
                case CGIndexAssignOp:
                    printIndexAssign(node);
                    break;

                case CGLoopAtomicResultOp:
                    printOperationAlias(node); // just follow the argument
                    break;
                case CGLoopStartOp:
                    printLoopStart(node);
                    break;
                case CGLoopIndexedIndepOp:
                    printLoopIndexedIndep(node);
                    break;
                case CGLoopIndexedDepOp:
                    printLoopIndexedDep(node);
                    break;
                case CGLoopEndOp:
                    printLoopEnd(node);
                    break;
                default:
                    std::stringstream ss;
                    ss << "Unkown operation code '" << op << "'.";
                    throw CGException(ss.str());
            }
        }

        virtual void printUnaryFunction(OperationNode<Base>& op) throw (CGException) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 1, "Invalid number of arguments for unary function");

            switch (op.getOperationType()) {
                case CGAbsOp:
                    _code << absFuncName();
                    break;
                case CGAcosOp:
                    _code << acosFuncName();
                    break;
                case CGAsinOp:
                    _code << asinFuncName();
                    break;
                case CGAtanOp:
                    _code << atanFuncName();
                    break;
                case CGCoshOp:
                    _code << coshFuncName();
                    break;
                case CGCosOp:
                    _code << cosFuncName();
                    break;
                case CGExpOp:
                    _code << expFuncName();
                    break;
                case CGLogOp:
                    _code << logFuncName();
                    break;
                case CGSinhOp:
                    _code << sinhFuncName();
                    break;
                case CGSinOp:
                    _code << sinFuncName();
                    break;
                case CGSqrtOp:
                    _code << sqrtFuncName();
                    break;
                case CGTanhOp:
                    _code << tanhFuncName();
                    break;
                case CGTanOp:
                    _code << tanFuncName();
                    break;
                default:
                    std::stringstream ss;
                    ss << "Unkown function name for operation code '" << op.getOperationType() << "'.";
                    throw CGException(ss.str());
            }

            _code << "(";
            print(op.getArguments()[0]);
            _code << ")";
        }

        virtual void printPowFunction(OperationNode<Base>& op) throw (CGException) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for pow() function");

            _code << powFuncName() << "(";
            print(op.getArguments()[0]);
            _code << ", ";
            print(op.getArguments()[1]);
            _code << ")";
        }

        virtual void printSignFunction(OperationNode<Base>& op) throw (CGException) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 1, "Invalid number of arguments for sign() function");
            assert(op.getArguments()[0].getOperation() != NULL);
            assert(op.getArguments()[0].getOperation()->getVariableID() > 0);

            OperationNode<Base>& arg = *op.getArguments()[0].getOperation();

            const std::string& argName = createVariableName(arg);

            _code << "(" << argName << " " << _C_COMP_OP_GT << " ";
            printParameter(Base(0.0));
            _code << "?";
            printParameter(Base(1.0));
            _code << ":(" << argName << " " << _C_COMP_OP_LT << " ";
            printParameter(Base(0.0));
            _code << "?";
            printParameter(Base(-1.0));
            _code << ":";
            printParameter(Base(0.0));
            _code << "))";
        }

        virtual void printOperationAlias(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 1, "Invalid number of arguments for alias");
            print(op.getArguments()[0]);
        }

        virtual void printOperationAdd(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for addition");

            print(op.getArguments()[0]);
            _code << " + ";
            print(op.getArguments()[1]);
        }

        virtual void printOperationMinus(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for substraction");

            const Argument<Base>& left = op.getArguments()[0];
            const Argument<Base>& right = op.getArguments()[1];

            const OperationNode<Base>* opRight = right.getOperation();
            bool encloseRight = opRight != NULL &&
                    opRight->getVariableID() == 0 &&
                    opRight->getOperationType() != CGDivOp &&
                    opRight->getOperationType() != CGMulOp &&
                    !isFunction(opRight->getOperationType());

            print(left);
            _code << " - ";
            if (encloseRight) {
                _code << "(";
            }
            print(right);
            if (encloseRight) {
                _code << ")";
            }
        }

        virtual void printOperationDiv(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for division");

            const Argument<Base>& left = op.getArguments()[0];
            const Argument<Base>& right = op.getArguments()[1];

            const OperationNode<Base>* opLeft = left.getOperation();
            bool encloseLeft = opLeft != NULL &&
                    opLeft->getVariableID() == 0 &&
                    !isFunction(opLeft->getOperationType());
            const OperationNode<Base>* opRight = right.getOperation();
            bool encloseRight = opRight != NULL &&
                    opRight->getVariableID() == 0 &&
                    !isFunction(opRight->getOperationType());

            if (encloseLeft) {
                _code << "(";
            }
            print(left);
            if (encloseLeft) {
                _code << ")";
            }
            _code << " / ";
            if (encloseRight) {
                _code << "(";
            }
            print(right);
            if (encloseRight) {
                _code << ")";
            }
        }

        virtual void printOperationMul(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for multiplication");

            const Argument<Base>& left = op.getArguments()[0];
            const Argument<Base>& right = op.getArguments()[1];

            const OperationNode<Base>* opLeft = left.getOperation();
            bool encloseLeft = opLeft != NULL &&
                    opLeft->getVariableID() == 0 &&
                    opLeft->getOperationType() != CGDivOp &&
                    opLeft->getOperationType() != CGMulOp &&
                    !isFunction(opLeft->getOperationType());
            const OperationNode<Base>* opRight = right.getOperation();
            bool encloseRight = opRight != NULL &&
                    opRight->getVariableID() == 0 &&
                    opRight->getOperationType() != CGDivOp &&
                    opRight->getOperationType() != CGMulOp &&
                    !isFunction(opRight->getOperationType());

            if (encloseLeft) {
                _code << "(";
            }
            print(left);
            if (encloseLeft) {
                _code << ")";
            }
            _code << " * ";
            if (encloseRight) {
                _code << "(";
            }
            print(right);
            if (encloseRight) {
                _code << ")";
            }
        }

        virtual void printOperationUnaryMinus(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 1, "Invalid number of arguments for unary minus");

            const Argument<Base>& arg = op.getArguments()[0];

            const OperationNode<Base>* scf = arg.getOperation();
            bool enclose = scf != NULL &&
                    scf->getVariableID() == 0 &&
                    scf->getOperationType() != CGDivOp &&
                    scf->getOperationType() != CGMulOp &&
                    !isFunction(scf->getOperationType());

            _code << "-";
            if (enclose) {
                _code << "(";
            } else {
                _code << " "; // there may be several - together -> space required
            }
            print(arg);
            if (enclose) {
                _code << ")";
            }
        }

        virtual void printConditionalAssignment(OperationNode<Base>& op) {
            assert(op.getVariableID() > 0);

            const std::vector<Argument<Base> >& args = op.getArguments();
            const Argument<Base> &left = args[0];
            const Argument<Base> &right = args[1];
            const Argument<Base> &trueCase = args[2];
            const Argument<Base> &falseCase = args[3];

            bool isDep = isDependent(op);
            const std::string& varName = createVariableName(op);

            if ((trueCase.getParameter() != NULL && falseCase.getParameter() != NULL && *trueCase.getParameter() == *falseCase.getParameter()) ||
                    (trueCase.getOperation() != NULL && falseCase.getOperation() != NULL && trueCase.getOperation() == falseCase.getOperation())) {
                // true and false cases are the same
                printAssigmentStart(op, varName, isDep);
                print(trueCase);
                printAssigmentEnd(op);
            } else {
                _code << _indentation << "if( ";
                print(left);
                _code << " " << getComparison(op.getOperationType()) << " ";
                print(right);
                _code << " ) {\n";
                _code << _spaces;
                printAssigmentStart(op, varName, isDep);
                print(trueCase);
                printAssigmentEnd(op);
                _code << _indentation << "} else {\n";
                _code << _spaces;
                printAssigmentStart(op, varName, isDep);
                print(falseCase);
                printAssigmentEnd(op);
                _code << _indentation << "}\n";
            }
        }

        virtual void printArrayCreationOp(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() > 0, "Invalid number of arguments for array creation operation");

            const size_t id = op.getVariableID();
            const std::vector<Argument<Base> >& args = op.getArguments();
            const size_t argSize = args.size();

            if (argSize > 1 && args[0].getParameter() != NULL) {
                const Base& value = *args[0].getParameter();
                bool sameValue = true;
                for (size_t i = 1; i < argSize; i++) {
                    if (args[i].getParameter() == NULL || *args[i].getParameter() != value) {
                        sameValue = false;
                        break;
                    }
                }
                if (sameValue) {
                    bool assign = false;
                    for (size_t i = 0; i < argSize; i++) {
                        const Argument<Base>* oldArg = _tmpArrayValues[id - 1 + i];
                        if (oldArg == NULL || oldArg->getParameter() == NULL || *oldArg->getParameter() != value) {
                            assign = true;
                            break;
                        }
                    }
                    if (assign) {
                        _code << _indentation << "for(i = 0; i < " << argSize << "; i++) "
                                "(" << _nameGen->generateTemporaryArray(op) << ")[i] = " << value << ";\n";

                        for (size_t i = 0; i < argSize; i++) {
                            _tmpArrayValues[id - 1 + i] = &args[i];
                        }
                    }
                    return;
                }
            }

            bool firstElement = true;
            for (size_t i = 0; i < argSize; i++) {
                bool newValue = true;
                if (_tmpArrayValues[id - 1 + i] != NULL) {
                    const Argument<Base>& oldArg = *_tmpArrayValues[id - 1 + i];
                    if (oldArg.getParameter() != NULL) {
                        if (args[i].getParameter() != NULL) {
                            newValue = (*args[i].getParameter() != *oldArg.getParameter());
                        }
                    } else {
                        newValue = (args[i].getOperation() != oldArg.getOperation());
                    }
                }

                if (newValue) {
                    if (firstElement) {
                        _code << _indentation << auxArrayName_ << " = " << _nameGen->generateTemporaryArray(op) << "; // size: " << args.size() << "\n";
                        firstElement = false;
                    }
                    _code << _indentation << auxArrayName_ << "[" << i << "] = ";
                    print(args[i]);
                    _code << ";\n";

                    _tmpArrayValues[id - 1 + i] = &args[i];
                }
            }
        }

        virtual void printArrayElementOp(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for array element operation");
            CPPADCG_ASSERT_KNOWN(op.getArguments()[0].getOperation() != NULL, "Invalid argument for array element operation");
            CPPADCG_ASSERT_KNOWN(op.getInfo().size() == 1, "Invalid number of information indexes for array element operation");

            OperationNode<Base>& arrayOp = *op.getArguments()[0].getOperation();
            _code << "(" << _nameGen->generateTemporaryArray(arrayOp) << ")[" << op.getInfo()[0] << "]";
        }

        virtual void printAtomicForwardOp(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for atomic forward operation");
            CPPADCG_ASSERT_KNOWN(op.getInfo().size() == 3, "Invalid number of information elements for atomic forward operation");
            size_t id = op.getInfo()[0];
            size_t atomicIndex = _atomicFunctionId2Index.at(id);
            int q = op.getInfo()[1];
            int p = op.getInfo()[2];
            OperationNode<Base>& tx = *op.getArguments()[0].getOperation();
            OperationNode<Base>& ty = *op.getArguments()[1].getOperation();

            createVariableName(tx);
            createVariableName(ty);

            _code << _indentation << "atomicFun.forward(atomicFun.libModel, "
                    << atomicIndex << ", "
                    << q << ", "
                    << p << ", "
                    << *tx.getName() << ", " << tx.getArguments().size() << ", "
                    << *ty.getName() << ", " << ty.getArguments().size() << "); // "
                    << _atomicFunctionId2Name.at(id)
                    << "\n";

            /**
             * the values of ty are now changed
             */
            markArrayChanged(ty);
        }

        virtual void printAtomicReverseOp(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 4, "Invalid number of arguments for atomic forward operation");
            CPPADCG_ASSERT_KNOWN(op.getInfo().size() == 2, "Invalid number of information elements for atomic forward operation");
            size_t id = op.getInfo()[0];
            int p = op.getInfo()[1];
            size_t atomicIndex = _atomicFunctionId2Index.at(id);
            OperationNode<Base>& tx = *op.getArguments()[0].getOperation();
            OperationNode<Base>& ty = *op.getArguments()[1].getOperation();
            OperationNode<Base>& px = *op.getArguments()[2].getOperation();
            OperationNode<Base>& py = *op.getArguments()[3].getOperation();

            CPPADCG_ASSERT_KNOWN(tx.getArguments().size() == px.getArguments().size(), "Invalid array length");
            CPPADCG_ASSERT_KNOWN(ty.getArguments().size() == py.getArguments().size(), "Invalid array length");

            createVariableName(tx);
            createVariableName(ty);
            createVariableName(px);
            createVariableName(py);

            _code << _indentation << "atomicFun.reverse(atomicFun.libModel, "
                    << atomicIndex << ", "
                    << p << ", "
                    << *tx.getName() << ", " << *ty.getName() << ", "
                    << *px.getName() << ", " << *py.getName() << ", "
                    << tx.getArguments().size() << ", " << ty.getArguments().size() << "); // "
                    << _atomicFunctionId2Name.at(id)
                    << "\n";

            /**
             * the values of px are now changed
             */
            markArrayChanged(px);
        }

        inline void markArrayChanged(OperationNode<Base>& ty) {
            size_t id = ty.getVariableID();
            size_t tySize = ty.getArguments().size();
            for (size_t i = 0; i < tySize; i++) {
                _tmpArrayValues[id - 1 + i] = NULL;
            }
        }

        virtual void printLoopStart(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGLoopStartOp, "Invalid node type");

            LoopStartOperationNode<Base>& lnode = static_cast<LoopStartOperationNode<Base>&> (node);
            const LoopNodeInfo<Base>& loopInfo = lnode.getLoopInfo();

            _currentLoops.push_back(&loopInfo);

            const std::string& jj = loopInfo.getIndex().getName();
            std::string iterationCount;
            if (loopInfo.getIterationCountNode() != NULL) {
                iterationCount = loopInfo.getIterationCountNode()->getIndex().getName();
            } else {
                std::ostringstream oss;
                oss << loopInfo.getIterationCount();
                iterationCount = oss.str();
            }

            _code << _spaces << "for("
                    << jj << " = 0; "
                    << jj << " < " << iterationCount << "; "
                    << jj << "++) {\n";
            _indentation = _spaces + _spaces;
        }

        virtual void printLoopEnd(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGLoopEndOp, "Invalid node type");

            _code << _spaces << "}\n";
            _indentation = _spaces;
            _currentLoops.pop_back();
        }

        virtual void printLoopIndexedDep(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() >= 1, "Invalid number of arguments for loop indexed dependent operation");
            CPPADCG_ASSERT_KNOWN(!_currentLoops.empty(), "Not inside a loop");

            // CGLoopIndexedDepOp
            print(node.getArguments()[0]);
        }

        virtual void printLoopIndexedIndep(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGLoopIndexedIndepOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getInfo().size() == 1, "Invalid number of information elements for loop indexed independent operation");
            CPPADCG_ASSERT_KNOWN(!_currentLoops.empty(), "Not inside a loop");

            // CGLoopIndexedIndepOp
            size_t pos = node.getInfo()[1];
            const IndexPattern* ip = (*_loopIndependentIndexPatterns)[pos];
            _code << _nameGen->generateIndexedIndependent(node, *ip);
        }

        virtual void printIndexAssign(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGIndexAssignOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() > 0, "Invalid number of argumets for an index assigment operation");

            IndexAssignOperationNode<Base>& inode = static_cast<IndexAssignOperationNode<Base>&> (node);

            const IndexPattern& ip = inode.getIndexPattern();
            _code << _indentation << inode.getIndex().getName()
                    << " = " << createIndexPattern(ip) << ";\n";
        }

        inline bool isDependent(const OperationNode<Base>& arg) const {
            if (arg.getOperationType() == CGLoopIndexedDepOp) {
                return true;
            }
            size_t id = arg.getVariableID();
            return id > _independentSize && id < _minTemporaryVarID;
        }

        virtual void printParameter(const Base& value) {
            // make sure all digits of floating point values are printed
            _code.unsetf(std::ios::floatfield);
            //std::scientific 
            _code << std::setprecision(std::numeric_limits< Base >::digits10 + 2) << value;
        }

        virtual const std::string& getComparison(enum CGOpCode op) const {
            switch (op) {
                case CGComOpLt:
                    return _C_COMP_OP_LT;

                case CGComOpLe:
                    return _C_COMP_OP_LE;

                case CGComOpEq:
                    return _C_COMP_OP_EQ;

                case CGComOpGe:
                    return _C_COMP_OP_GE;

                case CGComOpGt:
                    return _C_COMP_OP_GT;

                case CGComOpNe:
                    return _C_COMP_OP_NE;

                default:
                    CPPAD_ASSERT_UNKNOWN(0);
            }
            throw CGException("Invalid comparison operator code"); // should never get here
        }

        static bool isFunction(enum CGOpCode op) {
            return isUnaryFunction(op) || op == CGPowOp;
        }

        static bool isUnaryFunction(enum CGOpCode op) {
            switch (op) {
                case CGAbsOp:
                case CGAcosOp:
                case CGAsinOp:
                case CGAtanOp:
                case CGCoshOp:
                case CGCosOp:
                case CGExpOp:
                case CGLogOp:
                case CGSinhOp:
                case CGSinOp:
                case CGSqrtOp:
                case CGTanhOp:
                case CGTanOp:
                    return true;
                default:
                    return false;
            }
        }

        static bool isCondAssign(enum CGOpCode op) {
            switch (op) {
                case CGComOpLt:
                case CGComOpLe:
                case CGComOpEq:
                case CGComOpGe:
                case CGComOpGt:
                case CGComOpNe:
                    return true;
                default:
                    return false;
            }
        }
    };

    template<class Base>
    const std::string CLanguage<Base>::_C_COMP_OP_LT = "<";
    template<class Base>
    const std::string CLanguage<Base>::_C_COMP_OP_LE = "<=";
    template<class Base>
    const std::string CLanguage<Base>::_C_COMP_OP_EQ = "==";
    template<class Base>
    const std::string CLanguage<Base>::_C_COMP_OP_GE = ">=";
    template<class Base>
    const std::string CLanguage<Base>::_C_COMP_OP_GT = ">";
    template<class Base>
    const std::string CLanguage<Base>::_C_COMP_OP_NE = "!=";

    template<class Base>
    const std::string CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION = "struct CLangAtomicFun {\n"
    "    void* libModel;\n"
    "    int (*forward)(void* libModel, int atomicIndex, int q, int p, const void* tx, unsigned long txSize, void* ty, unsigned long tySize);\n"
    "    int (*reverse)(void* libModel, int atomicIndex, int p, const void* tx, const void* ty, void* px, const void* py, unsigned long xSize, unsigned long ySize);\n"
    "};";

}

#endif