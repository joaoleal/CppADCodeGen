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
     * \author Joao Leal
     */
    template<class Base>
    class CLanguage : public Language<Base> {
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
        //
        std::string _inArgName;
        //
        std::string _outArgName;
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

    public:

        CLanguage(const std::string& varTypeName, size_t spaces = 3) :
            _baseTypeName(varTypeName),
            _spaces(spaces, ' '),
            _inArgName("in"),
            _outArgName("out"),
            _nameGen(NULL),
            _independentSize(0),
            _dependent(NULL),
            _depAssignOperation("="),
            _ignoreZeroDepAssign(false),
            _maxAssigmentsPerFunction(0),
            _sources(NULL) {
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

        virtual std::string generateTemporaryVariableDeclaration() {
            assert(_nameGen != NULL);

            // declare variables
            const std::vector<FuncArgument>& tmpArg = _nameGen->getTemporary();

            CPPADCG_ASSERT_KNOWN(tmpArg.size() == 1,
                                 "There must be one temporary variable");

            _ss << _spaces << "// auxiliary variables\n";
            if (tmpArg[0].array) {
                size_t size = _nameGen->getMaxTemporaryVariableID() + 1 - _nameGen->getMinTemporaryVariableID();
                _ss << _spaces << _baseTypeName << " " << tmpArg[0].name << "[" << size << "];\n";
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

            // save some info
            _independentSize = info.independent.size();
            _dependent = &info.dependent;
            _nameGen = &info.nameGen;
            _minTemporaryVarID = info.minTemporaryVarID;
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
            CPPADCG_ASSERT_KNOWN(tmpArg.size() == 1,
                                 "There must be one temporary variable");

            if (!_functionName.empty()) {
                defaultFuncArgDcl_ = _baseTypeName + " const *const * " + _inArgName +
                        ", " + _baseTypeName + " *const * " + _outArgName;

                localFuncArgDcl_ = defaultFuncArgDcl_ + ", " + argumentDeclaration(tmpArg[0]);
                localFuncArgs_ = _inArgName + ", " + _outArgName + ", " + tmpArg[0].name;
            }

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
                            op.setName(_nameGen->generateTemporary(op));
                        }
                    }
                }

                size_t assignCount = 0;
                for (it = variableOrder.begin(); it != variableOrder.end(); ++it) {
                    if (assignCount >= _maxAssigmentsPerFunction && multiFunction) {
                        assignCount = 0;
                        saveLocalFunction(localFuncNames);
                    }

                    SourceCodeFragment<Base>& op = **it;

                    bool isDep = isDependent(op);
                    if (!isDep) {
                        _temporary[op.variableID()] = &op;
                    }

                    bool condAssign = isCondAssign(op.operation());
                    if (!condAssign) {
                        _code << _spaces << createVariableName(op) << " ";
                        _code << (isDep ? _depAssignOperation : "=") << " ";
                    }
                    printExpressionNoVarCheck(op);
                    if (!condAssign) {
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

                // forward declarations
                for (size_t i = 0; i < localFuncNames.size(); i++) {
                    _code << "void " << localFuncNames[i] << "(" << localFuncArgDcl_ << ");\n";
                }
                _code << "\n"
                        << "void " << _functionName << "(" << defaultFuncArgDcl_ << ") {\n";
                _nameGen->customFunctionVariableDeclarations(_code);
                _code << generateIndependentVariableDeclaration() << "\n";
                _code << generateDependentVariableDeclaration() << "\n";
                _code << generateTemporaryVariableDeclaration() << "\n";
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
                        printParameter(dependent[i].getParameterValue());
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
                    _ss << "#include <math.h>\n\n"
                            << "void " << _functionName << "(" << defaultFuncArgDcl_ << ") {\n";
                    _nameGen->customFunctionVariableDeclarations(_ss);
                    _ss << generateIndependentVariableDeclaration() << "\n";
                    _ss << generateDependentVariableDeclaration() << "\n";
                    _ss << generateTemporaryVariableDeclaration() << "\n";
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

            _ss << "#include <math.h>\n\n"
                    << "void " << funcName << "(" << localFuncArgDcl_ << ") {\n";
            _nameGen->customFunctionVariableDeclarations(_ss);
            _ss << generateIndependentVariableDeclaration() << "\n";
            _ss << generateDependentVariableDeclaration() << "\n";
            _nameGen->prepareCustomFunctionVariables(_ss);
            _ss << _code.str();
            _nameGen->finalizeCustomFunctionVariables(_ss);
            _ss << "}\n\n";

            (*_sources)[funcName + ".c"] = _ss.str();
            localFuncNames.push_back(funcName);

            _code.str("");
            _ss.str("");
        }

        virtual bool createsNewVariable(const SourceCodeFragment<Base>& op) {
            return op.totalUsageCount() > 1 ||
                    op.operation() == CGComOpLt ||
                    op.operation() == CGComOpLe ||
                    op.operation() == CGComOpEq ||
                    op.operation() == CGComOpGe ||
                    op.operation() == CGComOpGt ||
                    op.operation() == CGComOpNe;
        }

        virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) {
            return op == CGSignOp;
        }

        inline const std::string& createVariableName(SourceCodeFragment<Base>& var) {
            assert(var.variableID() > 0);

            if (var.getName() == NULL) {
                if (var.variableID() <= _independentSize) {
                    var.setName(_nameGen->generateIndependent(var));

                } else if (var.variableID() < _minTemporaryVarID) {
                    std::map<size_t, size_t>::const_iterator it = _dependentIDs.find(var.variableID());
                    assert(it != _dependentIDs.end());

                    size_t index = it->second;
                    const CG<Base>& dep = (*_dependent)[index];
                    var.setName(_nameGen->generateDependent(dep, index));

                } else {
                    var.setName(_nameGen->generateTemporary(var));
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
                _code << createVariableName(op);
                return;
            }

            printExpressionNoVarCheck(op);
        }

        virtual void printExpressionNoVarCheck(SourceCodeFragment<Base>& op) throw (CGException) {
            switch (op.operation()) {
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
}

#endif