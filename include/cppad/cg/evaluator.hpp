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
 * @todo implement nonrecursive algorithm (so that there will never be any stack limit issues)
 */
template<class ScalarIn, class ScalarOut, class ActiveOut>
class EvaluatorBase {
protected:
    CodeHandler<ScalarIn>& handler_;
    const ActiveOut* indep_;
    CodeHandlerVector<ScalarIn, ActiveOut*> evals_;
    std::map<size_t, CppAD::vector<ActiveOut>* > evalsArrays_;
    bool underEval_;
public:

    /**
     * @param handler The source code handler
     */
    inline EvaluatorBase(CodeHandler<ScalarIn>& handler) :
        handler_(handler),
        indep_(nullptr),
        evals_(handler),
        underEval_(false) {
    }

    /**
     * @return true if this Evaluator is currently being used.
     */
    inline bool isUnderEvaluation() {
        return underEval_;
    }

    /**
     * Performs all the operations required to calculate the dependent 
     * variables with a (potentially) new data type
     * 
     * @param indepNew The new independent variables.
     * @param depOld Dependent variable vector (all variables must belong to
     *               the same code handler)
     * @return The dependent variable values
     * @throws CGException on error (such as an unhandled operation type)
     */
    inline std::vector<ActiveOut> evaluate(const std::vector<ActiveOut>& indepNew,
                                           const std::vector<CG<ScalarIn> >& depOld) {
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
     * @throws CGException on error (such as an unhandled operation type)
     */
    inline void evaluate(const ActiveOut* indepNew,
                         size_t indepSize,
                         ActiveOut* depNew,
                         const CG<ScalarIn>* depOld,
                         size_t depSize) {
        if (handler_.getIndependentVariableSize() != indepSize) {
            throw CGException("Invalid independent variable size. Expected ", handler_.getIndependentVariableSize(), " but got ", indepSize, ".");
        }

        CPPADCG_ASSERT_KNOWN(handler_.getIndependentVariableSize() == indepSize, "Invalid size the array of independent variables");

        if (underEval_) {
            throw CGException("The same evaluator cannot be used for simultaneous evaluations. "
                              "Either use a new one or wait for this one to finish its current evaluation.");
        }

        underEval_ = true;

        clear(); // clean-up from any previous call that might have failed
        evals_.fill(nullptr);
        evals_.adjustSize();

        try {

            indep_ = indepNew;

            for (size_t i = 0; i < depSize; i++) {
                depNew[i] = evalCG(depOld[i]);
            }

            clear(); // clean-up

        } catch (...) {
            underEval_ = false;
            throw;
        }

        underEval_ = false;
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

        for (const auto& p : evalsArrays_) {
            delete p.second;
        }
        evalsArrays_.clear();
    }

    inline ActiveOut evalCG(const CG<ScalarIn>& dep) {
        if (dep.isParameter()) {
            // parameter
            return ActiveOut(dep.getValue());
        } else {
            return evalOperations(*dep.getOperationNode());
        }
    }

    inline ActiveOut evalArg(const Argument<ScalarIn>& arg) {
        if (arg.getOperation() != nullptr) {
            return evalOperations(*arg.getOperation());
        } else {
            // parameter
            return ActiveOut(*arg.getParameter());
        }
    }

    inline const ActiveOut& evalOperations(OperationNode<ScalarIn>& node) {
        using CppAD::vector;

        CPPADCG_ASSERT_KNOWN(node.getHandlerPosition() < handler_.getManagedNodesCount(), "this node is not managed by the code handler");

        // check if this node was previously determined
        if (evals_[node] != nullptr) {
            return *evals_[node];
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
        CPPADCG_ASSERT_UNKNOWN(evals_[node] == nullptr);
        ActiveOut* resultPtr = new ActiveOut(result);
        evals_[node] = resultPtr;

        proccessActiveOut(node, *resultPtr);

        return *resultPtr;
    }

    virtual void proccessActiveOut(OperationNode<ScalarIn>& node,
                                   ActiveOut& a) {
    }

    virtual ActiveOut evalIndependent(OperationNode<ScalarIn>& node) {
        size_t index = handler_.getIndependentVariableIndex(node);
        return indep_[index];
    }

    inline CppAD::vector<ActiveOut>& evalArrayCreationOperation(OperationNode<ScalarIn>& node) {
        using CppAD::vector;

        CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGOpCode::ArrayCreation, "Invalid array creation operation");
        CPPADCG_ASSERT_KNOWN(node.getHandlerPosition() < handler_.getManagedNodesCount(), "this node is not managed by the code handler");

        // check if this node was previously determined
        auto it = evalsArrays_.find(node.getHandlerPosition());
        if (it != evalsArrays_.end()) {
            return *it->second;
        }

        const std::vector<Argument<ScalarIn> >& args = node.getArguments();
        vector<ActiveOut>* resultArray = new vector<ActiveOut>(args.size());

        // save it for reuse
        evalsArrays_[node.getHandlerPosition()] = resultArray;

        // define its elements
        for (size_t a = 0; a < args.size(); a++) {
            (*resultArray)[a] = evalArg(args[a]);
        }

        return *resultArray;
    }

    virtual void evalAtomicOperation(OperationNode<ScalarIn>& node) {
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

} // END cg namespace
} // END CppAD namespace

#endif