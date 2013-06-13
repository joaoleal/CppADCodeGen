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
        // line indentation
        const std::string _spaces;
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
        std::map<size_t, SourceCodeFragment<Base>*> _temporary;
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
            _sources(NULL) {
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
                typename std::map<size_t, SourceCodeFragment<Base>*>::const_iterator it;

                for (it = _temporary.begin(); it != _temporary.end(); ++it) {
                    SourceCodeFragment<Base>* var = it->second;
                    if (var->getName() == NULL) {
                        var->setName(_nameGen->generateTemporary(*var));
                    }
                }

                SourceCodeFragment<Base>* var1 = _temporary.begin()->second;
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

    protected:

        virtual void generateSourceCode(std::ostream& out, LanguageGenerationData<Base>& info) {
            const bool createFunction = !_functionName.empty();
            const bool multiFunction = createFunction && _maxAssigmentsPerFunction > 0 && _sources != NULL;

            // clean up
            _code.str("");
            _ss.str("");
            _temporary.clear();
            defaultFuncArgDcl_ = "";
            localFuncArgDcl_ = "";
            localFuncArgs_ = "";
            auxArrayName_ = "";

            // save some info
            _independentSize = info.independent.size();
            _dependent = &info.dependent;
            _nameGen = &info.nameGen;
            _minTemporaryVarID = info.minTemporaryVarID;
            _atomicFunctionId2Index = info.atomicFunctionId2Index;
            _atomicFunctionId2Name = info.atomicFunctionId2Name;
            const std::vector<CG<Base> >& dependent = info.dependent;
            const std::vector<SourceCodeFragment<Base>*>& variableOrder = info.variableOrder;

            /**
             * generate variable names
             */
            //generate names for the independent variables
            typename std::vector<SourceCodeFragment<Base> *>::const_iterator it;
            for (it = info.independent.begin(); it != info.independent.end(); ++it) {
                SourceCodeFragment<Base>& op = **it;
                if (op.getName() == NULL) {
                    op.setName(_nameGen->generateIndependent(op));
                }
            }

            // generate names for the dependent variables (must be after naming the independent)
            for (size_t i = 0; i < dependent.size(); i++) {
                SourceCodeFragment<Base>* op = dependent[i].getSourceCodeFragment();
                if (op != NULL && op->getName() == NULL) {
                    op->setName(_nameGen->generateDependent(dependent[i], i));
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

            if (!_functionName.empty()) {
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
                SourceCodeFragment<Base>* op = dependent[i].getSourceCodeFragment();
                if (op != NULL && op->operation() != CGInvOp) {
                    size_t varID = op->variableID();
                    if (varID > 0) {
                        std::map<size_t, size_t>::const_iterator it2 = _dependentIDs.find(varID);
                        if (it2 == _dependentIDs.end()) {
                            _dependentIDs[op->variableID()] = i;
                        } else {
                            // there can be several dependent variables with the same ID
                            dependentDuplicates.insert(i);
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
                    SourceCodeFragment<Base>& op = **it;
                    if (op.getName() == NULL) {
                        if (!isDependent(op)) {
                            if (requiresVariableName(op) && op.operation() != CGArrayCreationOp) {
                                op.setName(_nameGen->generateTemporary(op));
                            } else if (op.operation() == CGArrayCreationOp) {
                                op.setName(_nameGen->generateTemporaryArray(op));
                            }
                        }
                    }
                }

                size_t assignCount = 0;
                for (it = variableOrder.begin(); it != variableOrder.end(); ++it) {
                    // check if a new function should start
                    if (assignCount >= _maxAssigmentsPerFunction && multiFunction) {
                        assignCount = 0;
                        saveLocalFunction(localFuncNames);
                    }

                    SourceCodeFragment<Base>& op = **it;

                    bool isDep = isDependent(op);
                    if (!isDep) {
                        _temporary[op.variableID()] = &op;
                    }

                    bool createsVar = directlyAssignsVariable(op);
                    if (createsVar) {
                        _code << _spaces << createVariableName(op) << " ";
                        _code << (isDep ? _depAssignOperation : "=") << " ";
                    }
                    printExpressionNoVarCheck(op);
                    if (createsVar) {
                        _code << ";\n";
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
                    const std::string& origVarName = *dep.getSourceCodeFragment()->getName();

                    _code << _spaces << varName << " " << _depAssignOperation << " " << origVarName << ";\n";
                }
            }

            // constant dependent variables 
            _code << _spaces << "// dependent variables without operations\n";
            for (size_t i = 0; i < dependent.size(); i++) {
                if (dependent[i].isParameter()) {
                    if (!_ignoreZeroDepAssign || !dependent[i].IdenticalZero()) {
                        std::string varName = _nameGen->generateDependent(dependent[i], i);
                        _code << _spaces << varName << " " << _depAssignOperation << " ";
                        printParameter(dependent[i].getValue());
                        _code << ";\n";
                    }
                } else if (dependent[i].getSourceCodeFragment()->operation() == CGInvOp) {
                    std::string varName = _nameGen->generateDependent(dependent[i], i);
                    const std::string& indepName = *dependent[i].getSourceCodeFragment()->getName();
                    _code << _spaces << varName << " " << _depAssignOperation << " " << indepName << ";\n";
                }
            }

            /**
             * encapsulate the code in a function
             */
            if (createFunction) {
                if (localFuncNames.empty()) {
                    _ss << "#include <math.h>\n"
                            "#include <string.h>\n\n"
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

            _ss << "#include <math.h>\n"
                    "#include <string.h>\n\n"
                    << ATOMICFUN_STRUCT_DEFINITION << "\n\n"
                    << "void " << funcName << "(" << localFuncArgDcl_ << ") {\n";
            _nameGen->customFunctionVariableDeclarations(_ss);
            _ss << generateIndependentVariableDeclaration() << "\n";
            _ss << generateDependentVariableDeclaration() << "\n";
            size_t arraySize = _nameGen->getMaxTemporaryArrayVariableID() - 1;
            if (arraySize > 0) {
                _ss << _spaces << _baseTypeName << "* " << auxArrayName_ << ";\n";
            }
            _nameGen->prepareCustomFunctionVariables(_ss);
            _ss << _code.str();
            _nameGen->finalizeCustomFunctionVariables(_ss);
            _ss << "}\n\n";

            (*_sources)[funcName + ".c"] = _ss.str();
            localFuncNames.push_back(funcName);

            _code.str("");
            _ss.str("");
        }

        virtual bool createsNewVariable(const SourceCodeFragment<Base>& var) const {
            CGOpCode op = var.operation();
            return (var.totalUsageCount() > 1 &&
                    op != CGArrayElementOp) ||
                    op == CGAtomicForwardOp ||
                    op == CGAtomicReverseOp ||
                    op == CGArrayCreationOp ||
                    op == CGComOpLt ||
                    op == CGComOpLe ||
                    op == CGComOpEq ||
                    op == CGComOpGe ||
                    op == CGComOpGt ||
                    op == CGComOpNe;
        }

        virtual bool requiresVariableName(const SourceCodeFragment<Base>& op) const {
            return (op.totalUsageCount() > 1 &&
                    op.operation() != CGAtomicForwardOp &&
                    op.operation() != CGAtomicReverseOp);
        }

        virtual bool directlyAssignsVariable(const SourceCodeFragment<Base>& var) const {
            CGOpCode op = var.operation();
            return !isCondAssign(op) &&
                    op != CGArrayCreationOp &&
                    op != CGAtomicForwardOp &&
                    op != CGAtomicReverseOp;
        }

        virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) const {
            return op == CGSignOp;
        }

        inline const std::string& createVariableName(SourceCodeFragment<Base>& var) {
            assert(var.variableID() > 0);
            assert(var.operation() != CGAtomicForwardOp && var.operation() != CGAtomicReverseOp);

            if (var.getName() == NULL) {
                if (var.operation() != CGArrayCreationOp) {
                    if (var.variableID() <= _independentSize) {
                        // independent variable
                        var.setName(_nameGen->generateIndependent(var));

                    } else if (var.variableID() < _minTemporaryVarID) {
                        // dependent variable
                        std::map<size_t, size_t>::const_iterator it = _dependentIDs.find(var.variableID());
                        assert(it != _dependentIDs.end());

                        size_t index = it->second;
                        const CG<Base>& dep = (*_dependent)[index];
                        var.setName(_nameGen->generateDependent(dep, index));

                    } else {
                        // temporary variable
                        if (requiresVariableName(var) && var.operation() != CGArrayCreationOp) {
                            var.setName(_nameGen->generateTemporary(var));
                        }
                    }
                } else {
                    var.setName(_nameGen->generateTemporaryArray(var));
                }
            }

            return *var.getName();
        }

        virtual void printIndependentVariableName(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 0, "Invalid number of arguments for independent variable");

            _code << _nameGen->generateIndependent(op);
        }

        virtual void print(const Argument<Base>& arg) {
            if (arg.operation() != NULL) {
                // expression
                printExpression(*arg.operation());
            } else {
                // parameter
                printParameter(*arg.parameter());
            }
        }

        virtual void printExpression(SourceCodeFragment<Base>& op) throw (CGException) {
            if (op.variableID() > 0) {
                // use variable name
                _code << createVariableName(op);
            } else {
                // print expression code
                printExpressionNoVarCheck(op);
            }
        }

        virtual void printExpressionNoVarCheck(SourceCodeFragment<Base>& op) throw (CGException) {
            switch (op.operation()) {
                case CGArrayCreationOp:
                    printArrayCreationOp(op);
                    break;
                case CGArrayElementOp:
                    printArrayElementOp(op);
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
                    printUnaryFunction(op);
                    break;
                case CGAtomicForwardOp: // atomicFunction.forward(q, p, vx, vy, tx, ty)
                    printAtomicForwardOp(op);
                    break;
                case CGAtomicReverseOp: // atomicFunction.reverse(p, tx, ty, px, py)
                    printAtomicReverseOp(op);
                    break;
                case CGAddOp:
                    printOperationAdd(op);
                    break;
                case CGAliasOp:
                    printOperationAlias(op);
                    break;
                case CGComOpLt:
                case CGComOpLe:
                case CGComOpEq:
                case CGComOpGe:
                case CGComOpGt:
                case CGComOpNe:
                    printConditionalAssignment(op);
                    break;
                case CGDivOp:
                    printOperationDiv(op);
                    break;
                case CGInvOp:
                    printIndependentVariableName(op);
                    break;
                case CGMulOp:
                    printOperationMul(op);
                    break;
                case CGPowOp:
                    printPowFunction(op);
                    break;
                case CGSignOp:
                    printSignFunction(op);
                    break;
                case CGSubOp:
                    printOperationMinus(op);
                    break;

                case CGUnMinusOp:
                    printOperationUnaryMinus(op);
                    break;
                default:
                    std::stringstream ss;
                    ss << "Unkown operation code '" << op.operation() << "'.";
                    throw CGException(ss.str());
            }
        }

        virtual void printUnaryFunction(SourceCodeFragment<Base>& op) throw (CGException) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 1, "Invalid number of arguments for unary function");

            switch (op.operation()) {
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
                    ss << "Unkown function name for operation code '" << op.operation() << "'.";
                    throw CGException(ss.str());
            }

            _code << "(";
            print(op.arguments()[0]);
            _code << ")";
        }

        virtual void printPowFunction(SourceCodeFragment<Base>& op) throw (CGException) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for pow() function");

            _code << powFuncName() << "(";
            print(op.arguments()[0]);
            _code << ", ";
            print(op.arguments()[1]);
            _code << ")";
        }

        virtual void printSignFunction(SourceCodeFragment<Base>& op) throw (CGException) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 1, "Invalid number of arguments for sign() function");
            assert(op.arguments()[0].operation() != NULL);
            assert(op.arguments()[0].operation()->variableID() > 0);

            SourceCodeFragment<Base>& arg = *op.arguments()[0].operation();

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

        virtual void printOperationAlias(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 1, "Invalid number of arguments for alias");
            print(op.arguments()[0]);
        }

        virtual void printOperationAdd(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for addition");

            print(op.arguments()[0]);
            _code << " + ";
            print(op.arguments()[1]);
        }

        virtual void printOperationMinus(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for substraction");

            const Argument<Base>& left = op.arguments()[0];
            const Argument<Base>& right = op.arguments()[1];

            const SourceCodeFragment<Base>* opRight = right.operation();
            bool encloseRight = opRight != NULL &&
                    opRight->variableID() == 0 &&
                    opRight->operation() != CGDivOp &&
                    opRight->operation() != CGMulOp &&
                    !isFunction(opRight->operation());

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

        virtual void printOperationDiv(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for division");

            const Argument<Base>& left = op.arguments()[0];
            const Argument<Base>& right = op.arguments()[1];

            const SourceCodeFragment<Base>* opLeft = left.operation();
            bool encloseLeft = opLeft != NULL &&
                    opLeft->variableID() == 0 &&
                    !isFunction(opLeft->operation());
            const SourceCodeFragment<Base>* opRight = right.operation();
            bool encloseRight = opRight != NULL &&
                    opRight->variableID() == 0 &&
                    !isFunction(opRight->operation());

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

        virtual void printOperationMul(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for multiplication");

            const Argument<Base>& left = op.arguments()[0];
            const Argument<Base>& right = op.arguments()[1];

            const SourceCodeFragment<Base>* opLeft = left.operation();
            bool encloseLeft = opLeft != NULL &&
                    opLeft->variableID() == 0 &&
                    opLeft->operation() != CGDivOp &&
                    opLeft->operation() != CGMulOp &&
                    !isFunction(opLeft->operation());
            const SourceCodeFragment<Base>* opRight = right.operation();
            bool encloseRight = opRight != NULL &&
                    opRight->variableID() == 0 &&
                    opRight->operation() != CGDivOp &&
                    opRight->operation() != CGMulOp &&
                    !isFunction(opRight->operation());

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

        virtual void printOperationUnaryMinus(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 1, "Invalid number of arguments for unary minus");

            const Argument<Base>& arg = op.arguments()[0];

            const SourceCodeFragment<Base>* scf = arg.operation();
            bool enclose = scf != NULL &&
                    scf->variableID() == 0 &&
                    scf->operation() != CGDivOp &&
                    scf->operation() != CGMulOp &&
                    !isFunction(scf->operation());

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

        virtual void printConditionalAssignment(SourceCodeFragment<Base>& op) {
            assert(op.variableID() > 0);

            const std::vector<Argument<Base> >& args = op.arguments();
            const Argument<Base> &left = args[0];
            const Argument<Base> &right = args[1];
            const Argument<Base> &trueCase = args[2];
            const Argument<Base> &falseCase = args[3];

            bool isDep = isDependent(op);
            const std::string& varName = createVariableName(op);

            if ((trueCase.parameter() != NULL && falseCase.parameter() != NULL && *trueCase.parameter() == *falseCase.parameter()) ||
                    (trueCase.operation() != NULL && falseCase.operation() != NULL && trueCase.operation() == falseCase.operation())) {
                // true and false cases are the same
                _code << _spaces << varName << " ";
                _code << (isDep ? _depAssignOperation : "=") << " ";
                print(trueCase);
                _code << ";\n";
            } else {
                _code << _spaces << "if( ";
                print(left);
                _code << " " << getComparison(op.operation()) << " ";
                print(right);
                _code << " ) {\n";
                _code << _spaces << _spaces << varName << " ";
                _code << (isDep ? _depAssignOperation : "=") << " ";
                print(trueCase);
                _code << ";\n";
                _code << _spaces << "} else {\n";
                _code << _spaces << _spaces << varName << " ";
                _code << (isDep ? _depAssignOperation : "=") << " ";
                print(falseCase);
                _code << ";\n";
                _code << _spaces << "}\n";
            }
        }

        virtual void printArrayCreationOp(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() > 0, "Invalid number of arguments for array creation operation");
            const std::vector<Argument<Base> >& args = op.arguments();

            if (args.size() > 1 && args[0].parameter() != NULL) {
                const Base& value = *args[0].parameter();
                bool sameValue = true;
                for (size_t i = 1; i < args.size(); i++) {
                    if (args[i].parameter() == NULL || *args[0].parameter() != value) {
                        sameValue = false;
                        break;
                    }
                }
                if (sameValue) {
                    _code << _spaces << "memset(" << _nameGen->generateTemporaryArray(op) << ", "
                            << value << ", "
                            << args.size() << " * sizeof(" << _baseTypeName << "));\n";
                    return;
                }
            }

            _code << _spaces << auxArrayName_ << " = " << _nameGen->generateTemporaryArray(op) << "; // size: " << args.size() << "\n";

            for (size_t i = 0; i < args.size(); i++) {
                _code << _spaces << auxArrayName_ << "[" << i << "] = ";
                print(args[i]);
                _code << ";\n";
            }
        }

        virtual void printArrayElementOp(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for array element operation");
            CPPADCG_ASSERT_KNOWN(op.arguments()[0].operation() != NULL, "Invalid argument for array element operation");
            CPPADCG_ASSERT_KNOWN(op.info().size() == 1, "Invalid number of information indexes for array element operation");

            SourceCodeFragment<Base>& arrayOp = *op.arguments()[0].operation();
            _code << "(" << _nameGen->generateTemporaryArray(arrayOp) << ")[" << op.info()[0] << "]";
        }

        virtual void printAtomicForwardOp(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 2, "Invalid number of arguments for atomic forward operation");
            CPPADCG_ASSERT_KNOWN(op.info().size() == 3, "Invalid number of information elements for atomic forward operation");
            size_t id = op.info()[0];
            size_t atomicIndex = _atomicFunctionId2Index.at(id);
            int q = op.info()[1];
            int p = op.info()[2];
            SourceCodeFragment<Base>& tx = *op.arguments()[0].operation();
            SourceCodeFragment<Base>& ty = *op.arguments()[1].operation();

            createVariableName(tx);
            createVariableName(ty);

            _code << _spaces << "atomicFun.forward(atomicFun.libModel, "
                    << atomicIndex << ", "
                    << q << ", "
                    << p << ", "
                    << *tx.getName() << ", " << tx.arguments().size() << ", "
                    << *ty.getName() << ", " << ty.arguments().size() << "); // "
                    << _atomicFunctionId2Name.at(id)
                    << "\n";
        }

        virtual void printAtomicReverseOp(SourceCodeFragment<Base>& op) {
            CPPADCG_ASSERT_KNOWN(op.arguments().size() == 4, "Invalid number of arguments for atomic forward operation");
            CPPADCG_ASSERT_KNOWN(op.info().size() == 2, "Invalid number of information elements for atomic forward operation");
            size_t id = op.info()[0];
            size_t atomicIndex = _atomicFunctionId2Index.at(id);
            int p = op.info()[1];
            SourceCodeFragment<Base>& tx = *op.arguments()[0].operation();
            SourceCodeFragment<Base>& ty = *op.arguments()[1].operation();
            SourceCodeFragment<Base>& px = *op.arguments()[2].operation();
            SourceCodeFragment<Base>& py = *op.arguments()[3].operation();

            CPPADCG_ASSERT_KNOWN(tx.arguments().size() == px.arguments().size(), "Invalid array length");
            CPPADCG_ASSERT_KNOWN(ty.arguments().size() == py.arguments().size(), "Invalid array length");

            createVariableName(tx);
            createVariableName(ty);
            createVariableName(px);
            createVariableName(py);

            _code << _spaces << "atomicFun.reverse(atomicFun.libModel, "
                    << atomicIndex << ", "
                    << p << ", "
                    << *tx.getName() << ", " << *ty.getName() << ", "
                    << *px.getName() << ", " << *py.getName() << ", "
                    << tx.arguments().size() << ", " << ty.arguments().size() << "); // "
                    << _atomicFunctionId2Name.at(id)
                    << "\n";
        }

        inline bool isDependent(const SourceCodeFragment<Base>& arg) const {
            size_t id = arg.variableID();
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