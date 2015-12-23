#ifndef CPPAD_CG_EVALUATOR_INCLUDED
#define CPPAD_CG_EVALUATOR_INCLUDED
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
 * Utility class used for some code transformation.
 * 
 * The color field of Operation nodes is used.
 *
 * @todo implement nonrecursive algorithm (so that there will never be any stack limit issues)
 */
template<class ScalarIn, class ScalarOut, class ActiveOut>
class EvaluatorBase {
protected:
    CodeHandler<ScalarIn>& handler_;
    const ActiveOut* indep_;
    std::vector<ActiveOut*> evals_;
    std::vector<CppAD::vector<ActiveOut>* > evalsArrays_;
    std::set<OperationNode<ScalarIn>*> evalsAtomic_;
    std::map<size_t, CppAD::atomic_base<ScalarOut>* > atomicFunctions_;
public:

    /**
     * @param handler The source code handler
     */
    EvaluatorBase(CodeHandler<ScalarIn>& handler) :
        handler_(handler),
        indep_(nullptr) {
    }

    /**
     * Provides an atomic function.
     * 
     * @param id The atomic function ID
     * @param atomic The atomic function
     * @return True if an atomic function with the same ID was already
     *         defined, false otherwise.
     */
    virtual bool addAtomicFunction(size_t id, atomic_base<ScalarOut>& atomic) {
        bool exists = atomicFunctions_.find(id) != atomicFunctions_.end();
        atomicFunctions_[id] = &atomic;
        return exists;
    }

    virtual void addAtomicFunctions(const std::map<size_t, atomic_base<ScalarOut>* >& atomics) {
        for (const auto& it : atomics) {
            atomic_base<ScalarOut>* atomic = it.second;
            if (atomic != nullptr) {
                atomicFunctions_[it.first] = atomic;
            }
        }
    }

    /**
     * Performs all the operations required to calculate the dependent 
     * variables with a (potentially) new data type
     * 
     * @param indepNew The new independent variables.
     * @param depOld Dependent variable vector (all variables must belong to
     *               the same code handler)
     * @return The dependent variable values
     */
    inline std::vector<ActiveOut> evaluate(const std::vector<ActiveOut>& indepNew,
                                           const std::vector<CG<ScalarIn> >& depOld) throw (CGException) {
        std::vector<ActiveOut> depNew(depOld.size());

        evaluate(indepNew.data(), indepNew.size(), depNew.data(), depOld.data(), depNew.size());

        return depNew;
    }

    /**
     * Performs all the operations required to calculate the dependent
     * variables with a (potentially) new data type
     *
     * @param indepNew The new independent variables.
     * @param indepSize The size of the array of independent variables.
     * @param depNew The new dependent variable vector that will be created.
     * @param depOld Dependent variable vector (all variables must belong to
     *               the same code handler)
     * @param depSize The size of the array of dependent variables.
     */
    inline void evaluate(const ActiveOut* indepNew,
                         size_t indepSize,
                         ActiveOut* depNew,
                         const CG<ScalarIn>* depOld,
                         size_t depSize) throw (CGException) {
        if (handler_.getIndependentVariableSize() != indepSize) {
            throw CGException("Invalid independent variable size. Expected ", handler_.getIndependentVariableSize(), " but got ", indepSize, ".");
        }

        CPPADCG_ASSERT_KNOWN(handler_.getIndependentVariableSize() == indepSize, "Invalid size the array of independent variables");

        indep_ = indepNew;

        clear(); // clean-up

        handler_.startNewOperationTreeVisit();

        for (size_t i = 0; i < depSize; i++) {
            depNew[i] = evalCG(depOld[i]);
        }

        clear(); // clean-up
    }

    inline virtual ~EvaluatorBase() {
        clear();
    }

protected:

    /**
     * clean-up
     */
    inline void clear() {
        for (const ActiveOut* it : evals_) {
            delete it;
        }
        evals_.clear();

        for (const CppAD::vector<ActiveOut>* it : evalsArrays_) {
            delete it;
        }
        evalsArrays_.clear();
    }

    inline ActiveOut evalCG(const CG<ScalarIn>& dep) throw (CGException) {
        if (dep.isParameter()) {
            // parameter
            return ActiveOut(dep.getValue());
        } else {
            return evalOperations(*dep.getOperationNode());
        }
    }

    inline ActiveOut evalArg(const Argument<ScalarIn>& arg) throw (CGException) {
        if (arg.getOperation() != nullptr) {
            return evalOperations(*arg.getOperation());
        } else {
            // parameter
            return ActiveOut(*arg.getParameter());
        }
    }

    inline const ActiveOut& evalOperations(OperationNode<ScalarIn>& node) throw (CGException) {
        using CppAD::vector;

        // check if this node was previously determined
        if (handler_.isVisited(node)) {
            return *evals_[node.getColor()];
        }

        // first evaluation of this node
        const std::vector<Argument<ScalarIn> >& args = node.getArguments();
        const CGOpCode code = node.getOperationType();
        ActiveOut result;
        switch (code) {
            case CGOpCode::Assign:
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for assign()");
                result = evalArg(args[0]);
                break;
            case CGOpCode::Abs: //  abs(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for abs()");
                result = abs(evalArg(args[0]));
                break;
            case CGOpCode::Acos: // acos(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for acos()");
                result = acos(evalArg(args[0]));
                break;
            case CGOpCode::Add: //  a + b
                CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for addition");
                result = evalArg(args[0]) + evalArg(args[1]);
                break;
            case CGOpCode::Alias:
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for alias");
                result = evalArg(args[0]);
                break;
                //case CGArrayCreationOp: // {a, b, c ...}
            case CGOpCode::ArrayElement: // x[i]
            {
                const std::vector<size_t>& info = node.getInfo();
                CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for array element");
                CPPADCG_ASSERT_KNOWN(args[0].getOperation() != nullptr, "Invalid argument for array element");
                CPPADCG_ASSERT_KNOWN(args[1].getOperation() != nullptr, "Invalid argument for array element");
                CPPADCG_ASSERT_KNOWN(info.size() == 1, "Invalid number of information data for array element");
                size_t index = info[0];
                vector<ActiveOut>& array = evalArrayCreationOperation(*args[0].getOperation()); // array creation
                evalAtomicOperation(*args[1].getOperation()); // atomic operation

                result = array[index];
                break;
            }
            case CGOpCode::Asin: // asin(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for asin()");
                result = asin(evalArg(args[0]));
                break;
            case CGOpCode::Atan: // atan(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for atan()");
                result = atan(evalArg(args[0]));
                break;
                //CGAtomicForwardOp
                //CGAtomicReverseOp
            case CGOpCode::ComLt: // result = left < right? trueCase: falseCase
                CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareLt, )");
                result = CondExpOp(CompareLt, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                break;
            case CGOpCode::ComLe: // result = left <= right? trueCase: falseCase
                CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareLe, )");
                result = CondExpOp(CompareLe, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                break;
            case CGOpCode::ComEq: // result = left == right? trueCase: falseCase
                CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareEq, )");
                result = CondExpOp(CompareEq, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                break;
            case CGOpCode::ComGe: // result = left >= right? trueCase: falseCase
                CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareGe, )");
                result = CondExpOp(CompareGe, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                break;
            case CGOpCode::ComGt: // result = left > right? trueCase: falseCase
                CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareGt, )");
                result = CondExpOp(CompareGt, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                break;
            case CGOpCode::ComNe: // result = left != right? trueCase: falseCase
                CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareNe, )");
                result = CondExpOp(CompareNe, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                break;
            case CGOpCode::Cosh: // cosh(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for cosh()");
                result = cosh(evalArg(args[0]));
                break;
            case CGOpCode::Cos: //  cos(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for cos()");
                result = cos(evalArg(args[0]));
                break;
            case CGOpCode::Div: // a / b
                CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for division");
                result = evalArg(args[0]) / evalArg(args[1]);
                break;
            case CGOpCode::Exp: //  exp(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for exp()");
                result = exp(evalArg(args[0]));
                break;
            case CGOpCode::Inv: //                             independent variable
                result = evalIndependent(node);
                break;
            case CGOpCode::Log: //  log(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for log()");
                result = log(evalArg(args[0]));
                break;
            case CGOpCode::Mul: // a * b
                CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for multiplication");
                result = evalArg(args[0]) * evalArg(args[1]);
                break;
            case CGOpCode::Pow: //  pow(a,   b)
                CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for pow()");
                result = pow(evalArg(args[0]), evalArg(args[1]));
                break;
                //case PriOp: //  PrintFor(text, parameter or variable, parameter or variable)
            case CGOpCode::Sign: // result = (x > 0)? 1.0:((x == 0)? 0.0:-1)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sign()");
                result = sign(evalArg(args[0]));
                break;
            case CGOpCode::Sinh: // sinh(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sinh()");
                result = sinh(evalArg(args[0]));
                break;
            case CGOpCode::Sin: //  sin(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sin()");
                result = sin(evalArg(args[0]));
                break;
            case CGOpCode::Sqrt: // sqrt(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sqrt()");
                result = sqrt(evalArg(args[0]));
                break;
            case CGOpCode::Sub: //  a - b
                CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for subtraction");
                result = evalArg(args[0]) - evalArg(args[1]);
                break;
            case CGOpCode::Tanh: //  tanh(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for tanh()");
                result = tanh(evalArg(args[0]));
                break;
            case CGOpCode::Tan: //  tan(variable)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for tan()");
                result = tan(evalArg(args[0]));
                break;
            case CGOpCode::UnMinus: // -(a)
                CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for unary minus");
                result = -evalArg(args[0]);
                break;
            default:
                result = evalUnsupportedOperation(node);
        }

        // save it for reuse
        node.setColor(evals_.size());
        handler_.markVisited(node);
        ActiveOut* resultPtr = new ActiveOut(result);
        if (evals_.size() == evals_.capacity())
            evals_.reserve(evals_.size() * 3 / 2 + 1);
        evals_.push_back(resultPtr);

        return *resultPtr;
    }

    virtual ActiveOut evalIndependent(OperationNode<ScalarIn>& node) {
        size_t index = handler_.getIndependentVariableIndex(node);
        return indep_[index];
    }

    inline CppAD::vector<ActiveOut>& evalArrayCreationOperation(OperationNode<ScalarIn>& node) throw (CGException) {
        using CppAD::vector;

        CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGOpCode::ArrayCreation, "Invalid array creation operation");

        // check if this array node was previously determined
        if (handler_.isVisited(node)) {
            return *evalsArrays_[node.getColor()];
        }

        const std::vector<Argument<ScalarIn> >& args = node.getArguments();
        vector<ActiveOut>* resultArray = new vector<ActiveOut>(args.size());

        // save it for reuse
        node.setColor(evalsArrays_.size());
        handler_.markVisited(node);
        evalsArrays_.push_back(resultArray);

        // define its elements
        for (size_t a = 0; a < args.size(); a++) {
            (*resultArray)[a] = evalArg(args[a]);
        }

        return *resultArray;
    }

    virtual void evalAtomicOperation(OperationNode<ScalarIn>& node) throw (CGException) {
        throw CGException("Evaluator is unable to handle atomic functions for these variable types");
    }

    virtual ActiveOut evalUnsupportedOperation(OperationNode<ScalarIn>& node) {
        throw CGException("Unknown operation code '", node.getOperationType(), "'");
    }

};

/**
 * 
 */
template<class ScalarIn, class ScalarOut, class ActiveOut = CppAD::AD<ScalarOut> >
class Evaluator : public EvaluatorBase<ScalarIn, ScalarOut, ActiveOut> {
public:

    inline Evaluator(CodeHandler<ScalarIn>& handler) :
        EvaluatorBase<ScalarIn, ScalarOut, ActiveOut>(handler) {
    }

    inline virtual ~Evaluator() {
    }
};

/**
 * 
 */
template<class ScalarIn, class ScalarOut>
class Evaluator<ScalarIn, ScalarOut, CppAD::AD<ScalarOut> > : public EvaluatorBase<ScalarIn, ScalarOut, CppAD::AD<ScalarOut> > {
public:
    typedef CppAD::AD<ScalarOut> ActiveOut;
protected:
    typedef EvaluatorBase<ScalarIn, ScalarOut, CppAD::AD<ScalarOut> > BaseClass;
    using BaseClass::evalsAtomic_;
    using BaseClass::atomicFunctions_;
    using BaseClass::handler_;
    using BaseClass::evalArrayCreationOperation;
public:

    inline Evaluator(CodeHandler<ScalarIn>& handler) :
        EvaluatorBase<ScalarIn, ScalarOut, ActiveOut>(handler) {
    }

    inline virtual ~Evaluator() {
    }

protected:

    virtual void evalAtomicOperation(OperationNode<ScalarIn>& node) throw (CGException) override {
        using CppAD::vector;

        if (evalsAtomic_.find(&node) != evalsAtomic_.end()) {
            return;
        }

        if (node.getOperationType() != CGOpCode::AtomicForward) {
            throw CGException("Evaluator can only handle zero forward mode for atomic functions");
        }

        const std::vector<size_t>& info = node.getInfo();
        const std::vector<Argument<ScalarIn> >& args = node.getArguments();
        CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for atomic forward mode");
        CPPADCG_ASSERT_KNOWN(info.size() == 3, "Invalid number of information data for atomic forward mode");

        // find the atomic function
        size_t id = info[0];
        typename std::map<size_t, atomic_base<ScalarOut>* >::const_iterator itaf = atomicFunctions_.find(id);
        atomic_base<ScalarOut>* atomicFunction = nullptr;
        if (itaf != atomicFunctions_.end()) {
            atomicFunction = itaf->second;
        }

        if (atomicFunction == nullptr) {
            std::stringstream ss;
            ss << "No atomic function defined in the evaluator for ";
            const std::string* atomName = handler_.getAtomicFunctionName(id);
            if (atomName != nullptr) {
                ss << "'" << *atomName << "'";
            } else
                ss << "id '" << id << "'";
            throw CGException(ss.str());
        }

        size_t p = info[2];
        if (p != 0) {
            throw CGException("Evaluator can only handle zero forward mode for atomic functions");
        }
        const vector<ActiveOut>& ax = evalArrayCreationOperation(*args[0].getOperation());
        vector<ActiveOut>& ay = evalArrayCreationOperation(*args[1].getOperation());

        (*atomicFunction)(ax, ay);

        evalsAtomic_.insert(&node);
    }
};

} // END cg namespace
} // END CppAD namespace

#endif