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
        const vector<CG<Base> >* _dependent;
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
        // indexes defined as function arguments
        std::vector<const IndexDclrOperationNode<Base>*> _funcArgIndexes;
        // all used indexes
        const std::set<const IndexDclrOperationNode<Base>*>* _indexes;
        const std::vector<const IndexPattern*>* _loopDependentIndexPatterns;
        const std::vector<const IndexPattern*>* _loopIndependentIndexPatterns;
        std::vector<const LoopStartOperationNode<Base>*> _currentLoops;
    private:
        std::string funcArgDcl_;
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

        virtual void setFunctionIndexArgument(const IndexDclrOperationNode<Base>& funcArgIndex) {
            _funcArgIndexes.resize(1);
            _funcArgIndexes[0] = &funcArgIndex;
        }

        virtual void setFunctionIndexArguments(const std::vector<const IndexDclrOperationNode<Base>*>& funcArgIndexes) {
            _funcArgIndexes = funcArgIndexes;
        }

        virtual const std::vector<const IndexDclrOperationNode<Base>*>& getFunctionIndexArguments() const {
            return _funcArgIndexes;
        }

        virtual void setMaxAssigmentsPerFunction(size_t maxAssigmentsPerFunction, std::map<std::string, std::string>* sources) {
            _maxAssigmentsPerFunction = maxAssigmentsPerFunction;
            _sources = sources;
        }

        virtual std::string generateTemporaryVariableDeclaration(bool localOnly = false,
                                                                 bool zeroArrayDependents = false) {
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
            }
            if (arraySize > 0 || zeroArrayDependents) {
                _ss << _spaces << "unsigned long i;\n";
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

        virtual std::string generateFunctionArgumentsDcl() const {
            std::string args = generateFunctionIndexArgumentsDcl();
            if (!args.empty())
                args += ", ";
            args += generateDefaultFunctionArgumentsDcl();

            return args;
        }

        virtual std::string generateDefaultFunctionArgumentsDcl() const {
            return _baseTypeName + " const *const * " + _inArgName +
                    ", " + _baseTypeName + "*const * " + _outArgName +
                    ", " + generateArgumentAtomicDcl();
        }

        virtual std::string generateFunctionIndexArgumentsDcl() const {
            std::string argtxt;
            for (size_t a = 0; a < _funcArgIndexes.size(); a++) {
                if (a > 0) argtxt += ", ";
                argtxt += "unsigned long " + *_funcArgIndexes[a]->getName();
            }
            return argtxt;
        }

        virtual std::string generateDefaultFunctionArguments() const {
            return _inArgName + ", " + _outArgName + ", " + _atomicArgName;
        }

        virtual std::string generateFunctionIndexArguments() const {
            std::string argtxt;
            for (size_t a = 0; a < _funcArgIndexes.size(); a++) {
                if (a > 0) argtxt += ", ";
                argtxt += *_funcArgIndexes[a]->getName();
            }
            return argtxt;
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

        static inline void printIndexCondExpr(std::ostringstream& out,
                                              const std::vector<size_t>& info,
                                              const std::string& index) {
            CPPADCG_ASSERT_KNOWN(info.size() > 1 && info.size() % 2 == 0, "Invalid number of information elements for an index condition expression operation");

            size_t infoSize = info.size();
            for (size_t e = 0; e < infoSize; e += 2) {
                if (e > 0) {
                    out << " || ";
                }
                size_t min = info[e];
                size_t max = info[e + 1];
                if (min == max) {
                    out << index << " == " << min;
                } else if (min == 0) {
                    out << index << " <= " << max;
                } else if (max == std::numeric_limits<size_t>::max()) {
                    out << min << " <= " << index;
                } else {
                    if (infoSize != 2)
                        out << "(";
                    out << min << " <= " << index << " && " << index << " <= " << max;
                    if (infoSize != 2)
                        out << ")";
                }
            }
        }

        /***********************************************************************
         * index patterns
         **********************************************************************/
        static inline std::string indexPattern2String(const IndexPattern& ip,
                                                      const IndexDclrOperationNode<Base>& index);

        static inline std::string indexPattern2String(const IndexPattern& ip,
                                                      const std::vector<const IndexDclrOperationNode<Base>*>& indexes);

        static inline std::string linearIndexPattern2String(const LinearIndexPattern& lip,
                                                            const IndexDclrOperationNode<Base>& index);

    protected:

        virtual void generateSourceCode(std::ostream& out, LanguageGenerationData<Base>& info) {
            const bool createFunction = !_functionName.empty();
            const bool multiFunction = createFunction && _maxAssigmentsPerFunction > 0 && _sources != NULL;

            // clean up
            _code.str("");
            _ss.str("");
            _temporary.clear();
            _indentation = _spaces;
            funcArgDcl_ = "";
            localFuncArgDcl_ = "";
            localFuncArgs_ = "";
            auxArrayName_ = "";
            _currentLoops.clear();

            // save some info
            _independentSize = info.independent.size();
            _dependent = &info.dependent;
            _nameGen = &info.nameGen;
            _minTemporaryVarID = info.minTemporaryVarID;
            _atomicFunctionId2Index = info.atomicFunctionId2Index;
            _atomicFunctionId2Name = info.atomicFunctionId2Name;
            const vector<CG<Base> >& dependent = info.dependent;
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
            for (size_t j = 0; j < _independentSize; j++) {
                OperationNode<Base>& op = *info.independent[j];
                if (op.getName() == NULL) {
                    op.setName(_nameGen->generateIndependent(op));
                }
            }

            // generate names for the dependent variables (must be after naming the independent)
            for (size_t i = 0; i < dependent.size(); i++) {
                OperationNode<Base>* node = dependent[i].getOperationNode();
                if (node != NULL && node->getOperationType() != CGLoopEndOp && node->getName() == NULL) {
                    if (node->getOperationType() == CGLoopIndexedDepOp) {
                        size_t pos = node->getInfo()[0];
                        const IndexPattern* ip = (*_loopDependentIndexPatterns)[pos];
                        node->setName(_nameGen->generateIndexedDependent(*node, *ip));

                    } else {
                        node->setName(_nameGen->generateDependent(i));
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
                funcArgDcl_ = generateFunctionArgumentsDcl();
                localFuncArgDcl_ = funcArgDcl_ + ", " + argumentDeclaration(tmpArg[0]) + ", " + argumentDeclaration(tmpArg[1]);
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
                    OperationNode<Base>& node = **it;
                    CGOpCode op = node.getOperationType();
                    if (!isDependent(node) && op != CGIndexDeclarationOp) {
                        // variable names for temporaries must always be created since they might have been used before with a different name/id
                        if (requiresVariableName(node) && op != CGArrayCreationOp) {
                            node.setName(_nameGen->generateTemporary(node));
                        } else if (op == CGArrayCreationOp) {
                            node.setName(_nameGen->generateTemporaryArray(node));
                        }
                    }
                }

                /**
                 * Source code generation magic!
                 */
                if (info.zeroDependents) {
                    // zero initial values
                    const std::vector<FuncArgument>& depArg = _nameGen->getDependent();
                    for (size_t i = 0; i < depArg.size(); i++) {
                        const FuncArgument& a = depArg[i];
                        if (a.array) {
                            _code << _indentation << "for(i = 0; i < " << _dependent->size() << "; i++) " << a.name << "[i]";
                        } else {
                            _code << _indentation << _nameGen->generateDependent(i);
                        }
                        _code << " = ";
                        printParameter(Base(0.0));
                        _code << ";\n";
                    }
                }

                size_t assignCount = 0;
                for (it = variableOrder.begin(); it != variableOrder.end(); ++it) {
                    // check if a new function should start
                    if (assignCount >= _maxAssigmentsPerFunction && multiFunction && _currentLoops.empty()) {
                        assignCount = 0;
                        saveLocalFunction(localFuncNames, localFuncNames.empty() && info.zeroDependents);
                    }

                    OperationNode<Base>& node = **it;

                    // a dependent variable assigned by a loop does require any source code (its done inside the loop)
                    if (node.getOperationType() == CGDependentRefRhsOp) {
                        continue; // nothing to do (this operation is right hand side only)
                    }

                    printAssigment(node);

                    assignCount++;
                }

                if (localFuncNames.size() > 0 && assignCount > 0) {
                    assignCount = 0;
                    saveLocalFunction(localFuncNames, false);
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
                        << "void " << _functionName << "(" << funcArgDcl_ << ") {\n";
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
                    std::string varName = _nameGen->generateDependent(index);
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
                        std::string varName = _nameGen->generateDependent(i);
                        _code << _spaces << varName << " " << _depAssignOperation << " ";
                        printParameter(dependent[i].getValue());
                        _code << ";\n";
                    }
                } else if (dependent[i].getOperationNode()->getOperationType() == CGInvOp) {
                    if (!commentWritten) {
                        _code << _spaces << "// dependent variables without operations\n";
                        commentWritten = true;
                    }
                    std::string varName = _nameGen->generateDependent(i);
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
                            << "void " << _functionName << "(" << funcArgDcl_ << ") {\n";
                    _nameGen->customFunctionVariableDeclarations(_ss);
                    _ss << generateIndependentVariableDeclaration() << "\n";
                    _ss << generateDependentVariableDeclaration() << "\n";
                    _ss << generateTemporaryVariableDeclaration(true, info.zeroDependents) << "\n";
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

        inline void printAssigment(OperationNode<Base>& node) {
            printAssigment(node, node);
        }

        inline void printAssigment(OperationNode<Base>& nodeName,
                                   const Argument<Base>& nodeRhs) {
            if (nodeRhs.getOperation() != NULL) {
                printAssigment(nodeName, *nodeRhs.getOperation());
            } else {
                printAssigmentStart(nodeName);
                printParameter(*nodeRhs.getParameter());
                printAssigmentEnd(nodeName);
            }
        }

        inline void printAssigment(OperationNode<Base>& nodeName,
                                   OperationNode<Base>& nodeRhs) {
            bool createsVar = directlyAssignsVariable(nodeRhs); // do we need to do the assigment here?
            if (!createsVar) {
                printAssigmentStart(nodeName);
            }
            printExpressionNoVarCheck(nodeRhs);
            if (!createsVar) {
                printAssigmentEnd(nodeRhs);
            }

            if (nodeRhs.getOperationType() == CGArrayElementOp) {
                size_t arrayId = nodeRhs.getArguments()[0].getOperation()->getVariableID();
                size_t pos = nodeRhs.getInfo()[0];
                _tmpArrayValues[arrayId - 1 + pos] = NULL; // this could probably be removed!
            }
        }

        inline virtual void printAssigmentStart(OperationNode<Base>& op) {
            printAssigmentStart(op, createVariableName(op), isDependent(op));
        }

        inline virtual void printAssigmentStart(OperationNode<Base>& node, const std::string& varName, bool isDep) {
            if (!isDep) {
                _temporary[node.getVariableID()] = &node;
            }

            _code << _indentation << varName << " ";
            if (isDep) {
                CGOpCode op = node.getOperationType();
                if (op == CGDependentMultiAssignOp || (op == CGLoopIndexedDepOp && node.getInfo()[1] == 1)) {
                    _code << "+=";
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

        virtual void saveLocalFunction(std::vector<std::string>& localFuncNames,
                                       bool zeroDependentArray) {
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
            }
            if (arraySize > 0 || zeroDependentArray) {
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
                return op != CGArrayElementOp && op != CGIndexOp && op != CGIndexDeclarationOp;
            } else {
                return (op == CGArrayCreationOp ||
                        op == CGAtomicForwardOp ||
                        op == CGAtomicReverseOp ||
                        op == CGComOpLt ||
                        op == CGComOpLe ||
                        op == CGComOpEq ||
                        op == CGComOpGe ||
                        op == CGComOpGt ||
                        op == CGComOpNe ||
                        op == CGLoopIndexedDepOp ||
                        op == CGIndexAssignOp ||
                        op == CGStartIfOp ||
                        op == CGElseIfOp ||
                        op == CGElseOp ||
                        op == CGEndIfOp) &&
                        op != CGCondResultOp;
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
                    op != CGIndexAssignOp &&
                    op != CGStartIfOp &&
                    op != CGElseIfOp &&
                    op != CGElseOp &&
                    op != CGEndIfOp &&
                    op != CGCondResultOp);
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
                    op == CGDependentMultiAssignOp ||
                    op == CGLoopStartOp ||
                    op == CGLoopEndOp ||
                    op == CGIndexAssignOp ||
                    op == CGStartIfOp ||
                    op == CGElseIfOp ||
                    op == CGElseOp ||
                    op == CGEndIfOp ||
                    op == CGCondResultOp ||
                    op == CGIndexDeclarationOp;
        }

        virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) const {
            return op == CGSignOp || op == CGCondResultOp;
        }

        inline const std::string& createVariableName(OperationNode<Base>& var) {
            CGOpCode op = var.getOperationType();
            assert(var.getVariableID() > 0);
            assert(op != CGAtomicForwardOp);
            assert(op != CGAtomicReverseOp);
            assert(op != CGLoopStartOp);
            assert(op != CGLoopEndOp);
            assert(op != CGIndexOp);
            assert(op != CGIndexAssignOp);
            assert(op != CGIndexDeclarationOp);

            if (var.getName() == NULL) {
                if (op == CGArrayCreationOp) {
                    var.setName(_nameGen->generateTemporaryArray(var));
                } else if (op == CGLoopIndexedDepOp) {
                    size_t pos = var.getInfo()[0];
                    const IndexPattern* ip = (*_loopDependentIndexPatterns)[pos];
                    var.setName(_nameGen->generateIndexedDependent(var, *ip));

                } else if (op == CGLoopIndexedIndepOp) {
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
                        var.setName(_nameGen->generateDependent(index));

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

                case CGDependentMultiAssignOp:
                    printDependentMultiAssign(node);
                    break;

                case CGIndexOp:
                    break; // nothing to do
                case CGIndexAssignOp:
                    printIndexAssign(node);
                    break;
                case CGIndexDeclarationOp:
                    break; // already done

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
                case CGIndexCondExprOp:
                    printIndexCondExprOp(node);
                    break;
                case CGStartIfOp:
                    printStartIf(node);
                    break;
                case CGElseIfOp:
                    printElseIf(node);
                    break;
                case CGElseOp:
                    printElse(node);
                    break;
                case CGEndIfOp:
                    printEndIf(node);
                    break;
                case CGCondResultOp:
                    printCondResult(node);
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

            bool encloseRight = encloseInParenthesesMul(right.getOperation());

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

        static inline bool encloseInParenthesesDiv(const OperationNode<Base>* node) {
            while (node != NULL) {
                if (node->getVariableID() != 0)
                    return false;
                if (node->getOperationType() == CGAliasOp)
                    node = node->getArguments()[0].getOperation();
                else
                    break;
            }
            return node != NULL &&
                    node->getVariableID() == 0 &&
                    !isFunction(node->getOperationType());
        }

        virtual void printOperationDiv(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for division");

            const Argument<Base>& left = op.getArguments()[0];
            const Argument<Base>& right = op.getArguments()[1];

            bool encloseLeft = encloseInParenthesesDiv(left.getOperation());
            bool encloseRight = encloseInParenthesesDiv(right.getOperation());

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

        static inline bool encloseInParenthesesMul(const OperationNode<Base>* node) {
            while (node != NULL) {
                if (node->getVariableID() != 0)
                    return false;
                else if (node->getOperationType() == CGAliasOp)
                    node = node->getArguments()[0].getOperation();
                else
                    break;
            }
            return node != NULL &&
                    node->getVariableID() == 0 &&
                    node->getOperationType() != CGDivOp &&
                    node->getOperationType() != CGMulOp &&
                    !isFunction(node->getOperationType());
        }

        virtual void printOperationMul(OperationNode<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.getArguments().size() == 2, "Invalid number of arguments for multiplication");

            const Argument<Base>& left = op.getArguments()[0];
            const Argument<Base>& right = op.getArguments()[1];

            bool encloseLeft = encloseInParenthesesMul(left.getOperation());
            bool encloseRight = encloseInParenthesesMul(right.getOperation());

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

            bool enclose = encloseInParenthesesMul(arg.getOperation());

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

        virtual void printDependentMultiAssign(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGDependentMultiAssignOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() > 0, "Invalid number of arguments");

            const std::vector<Argument<Base> >& args = node.getArguments();
            for (size_t a = 0; a < args.size(); a++) {
                bool useArg = false;
                const Argument<Base>& arg = args[a];
                if (arg.getParameter() != NULL) {
                    useArg = true;
                } else {
                    CGOpCode op = arg.getOperation()->getOperationType();
                    useArg = op != CGDependentRefRhsOp && op != CGLoopEndOp && op != CGEndIfOp;
                }

                if (useArg) {
                    printAssigment(node, arg); // ignore other arguments!
                    break;
                }
            }
        }

        virtual void printLoopStart(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGLoopStartOp, "Invalid node type");

            LoopStartOperationNode<Base>& lnode = static_cast<LoopStartOperationNode<Base>&> (node);
            _currentLoops.push_back(&lnode);

            const std::string& jj = *lnode.getIndex().getName();
            std::string iterationCount;
            if (lnode.getIterationCountNode() != NULL) {
                iterationCount = *lnode.getIterationCountNode()->getIndex().getName();
            } else {
                std::ostringstream oss;
                oss << lnode.getIterationCount();
                iterationCount = oss.str();
            }

            _code << _spaces << "for("
                    << jj << " = 0; "
                    << jj << " < " << iterationCount << "; "
                    << jj << "++) {\n";
            _indentation += _spaces;
        }

        virtual void printLoopEnd(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGLoopEndOp, "Invalid node type");

            _indentation.resize(_indentation.size() - _spaces.size());

            _code << _indentation << "}\n";

            _currentLoops.pop_back();
        }

        virtual void printLoopIndexedDep(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() >= 1, "Invalid number of arguments for loop indexed dependent operation");

            // CGLoopIndexedDepOp
            print(node.getArguments()[0]);
        }

        virtual void printLoopIndexedIndep(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGLoopIndexedIndepOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getInfo().size() == 1, "Invalid number of information elements for loop indexed independent operation");

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
            _code << _indentation << (*inode.getIndex().getName())
                    << " = " << indexPattern2String(ip, inode.getIndexPatternIndexes()) << ";\n";
        }

        virtual void printIndexCondExprOp(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGIndexCondExprOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() == 1, "Invalid number of argumets for an index condition expression operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[0].getOperation() != NULL, "Invalid argumet for an index condition expression operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[0].getOperation()->getOperationType() == CGIndexOp, "Invalid argumet for an index condition expression operation");

            const std::vector<size_t>& info = node.getInfo();

            IndexOperationNode<Base>& iterationIndexOp = static_cast<IndexOperationNode<Base>&> (*node.getArguments()[0].getOperation());
            const std::string& index = *iterationIndexOp.getIndex().getName();

            printIndexCondExpr(_code, info, index);
        }

        virtual void printStartIf(OperationNode<Base>& node) {
            /**
             * the first argument is the condition, following arguments are
             * just extra dependencies that must be defined outside the if
             */
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGStartIfOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() >= 1, "Invalid number of argumets for an 'if start' operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[0].getOperation() != NULL, "Invalid argumet for an 'if start' operation");

            _code << _indentation << "if(";
            printIndexCondExprOp(*node.getArguments()[0].getOperation());
            _code << ") {\n";

            _indentation += _spaces;
        }

        virtual void printElseIf(OperationNode<Base>& node) {
            /**
             * the first argument is the condition, the second argument is the 
             * if start node, the following arguments are assigments in the
             * previous if branch
             */
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGElseIfOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() >= 2, "Invalid number of argumets for an 'else if' operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[0].getOperation() != NULL, "Invalid argumet for an 'else if' operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[1].getOperation() != NULL, "Invalid argumet for an 'else if' operation");

            _indentation.resize(_indentation.size() - _spaces.size());

            _code << _indentation << "} else if(";
            printIndexCondExprOp(*node.getArguments()[1].getOperation());
            _code << ") {\n";

            _indentation += _spaces;
        }

        virtual void printElse(OperationNode<Base>& node) {
            /**
             * the first argument is the  if start node, the following arguments
             * are assigments in the previous if branch
             */
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGElseOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() >= 1, "Invalid number of argumets for an 'else' operation");

            _indentation.resize(_indentation.size() - _spaces.size());

            _code << _indentation << "} else {\n";

            _indentation += _spaces;
        }

        virtual void printEndIf(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGEndIfOp, "Invalid node type for an 'end if' operation");

            _indentation.resize(_indentation.size() - _spaces.size());

            _code << _indentation << "}\n";
        }

        virtual void printCondResult(OperationNode<Base>& node) {
            CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGCondResultOp, "Invalid node type");
            CPPADCG_ASSERT_KNOWN(node.getArguments().size() == 2, "Invalid number of argumets for an assigment inside an if/else operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[0].getOperation() != NULL, "Invalid argumet for an an assigment inside an if/else operation");
            CPPADCG_ASSERT_KNOWN(node.getArguments()[1].getOperation() != NULL, "Invalid argumet for an an assigment inside an if/else operation");

            // just follow the argument
            OperationNode<Base>& nodeArg = *node.getArguments()[1].getOperation();
            printAssigment(nodeArg);
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