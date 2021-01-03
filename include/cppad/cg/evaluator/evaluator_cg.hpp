#ifndef CPPAD_CG_EVALUATOR_CG_INCLUDED
#define CPPAD_CG_EVALUATOR_CG_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
 *    Copyright (C) 2016 Ciengis
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
 * Specialization of class Evaluator for an output active type of CG<Base>.
 * This class should not be instantiated directly.
 */
template<class ScalarIn, class ScalarOut, class FinalEvaluatorType>
class EvaluatorCG : public EvaluatorOperations<ScalarIn, ScalarOut, CG<ScalarOut>, FinalEvaluatorType> {
    /**
     * must be friends with one of its super classes since there is a cast to
     * this type due to the curiously recurring template pattern (CRTP)
     */
    friend EvaluatorBase<ScalarIn, ScalarOut, CG<ScalarOut>, FinalEvaluatorType>;
    friend EvaluatorOperations<ScalarIn, ScalarOut, CG<ScalarOut>, FinalEvaluatorType>;
public:
    using ActiveIn = CG<ScalarIn>;
    using ActiveOut = CG<ScalarOut>;
    using NodeIn = OperationNode<ScalarIn>;
    using NodeOut = OperationNode<ScalarOut>;
    using ArgIn = Argument<ScalarIn>;
    using ArgOut = Argument<ScalarOut>;
protected:
    using Super = EvaluatorOperations<ScalarIn, ScalarOut, CG<ScalarOut>, FinalEvaluatorType>;
protected:
    /**
     * The source code handler used to create the evaluation results
     */
    CodeHandler<ScalarOut>* outHandler_;
    /**
     * Cache for the evaluation of atomic operations
     */
    std::map<const NodeIn*, std::vector<ScalarOut*>> atomicEvalResults_;
    /**
     * Whether or not the nodes with an operation type 'Pri' are printed out
     * during the evaluation.
     */
    bool printOutPriOperations_;
    using EvaluatorBase<ScalarIn, ScalarOut, CG<ScalarOut>, FinalEvaluatorType>::evals_;
public:

    inline EvaluatorCG(CodeHandler<ScalarIn>& handler) :
        Super(handler),
        outHandler_(nullptr),
        printOutPriOperations_(true) {
    }

    /**
     * Defines whether or not to print out the nodes with an operation type
     * 'Pri' during the evaluation.
     */
    inline void setPrintOutPrintOperations(bool print) {
        printOutPriOperations_ = print;
    }

    /**
     * Whether or not the nodes with an operation type 'Pri' are printed out
     * during the evaluation.
     */
    inline bool isPrintOutPrintOperations() const {
        return printOutPriOperations_;
    }

protected:

    /**
     * @note overrides the default analyzeOutIndeps() even though this method
     *        is not virtual (hides a method in EvaluatorBase)
     */
    inline void analyzeOutIndeps(const ActiveOut* indep,
                                 size_t n,
                                 const ActiveOut* par,
                                 size_t p) {
        CPPAD_ASSERT_KNOWN(indep != nullptr || n == 0, "null array with a non-zero size");
        outHandler_ = findHandler(ArrayView<const ActiveOut>(indep, n));

        if (par != nullptr && p > 0) {
            auto* outHandler = findHandler(ArrayView<const ActiveOut>(par, p));
            assert(outHandler == outHandler_);
        }
    }

    /**
     * @note overrides the default clear() even though this method
     *        is not virtual (hides a method in EvaluatorBase)
     */
    inline void clear() {
        Super::clear();

        for (const auto& it : atomicEvalResults_) {
            for (const ScalarOut* e : it.second) {
                delete e;
            }
        }
        atomicEvalResults_.clear();
    }


    /**
     * @note overrides the default processActiveOut() even though this method
     *        is not virtual (hides a method in EvaluatorOperations)
     */
    void processActiveOut(const NodeIn& node,
                          ActiveOut& a) {
        if (node.getName() != nullptr) {
            if (a.getOperationNode() != nullptr) {
                a.getOperationNode()->setName(*node.getName());
            }
        }
    }

    /**
     * @note overrides the default evalPrint() even though this method
     *        is not virtual (hides a method in EvaluatorOperations)
     */
    inline ActiveOut evalPrint(const NodeIn& node) {
        const std::vector<ArgIn>& args = node.getArguments();
        CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for print()")
        ActiveOut out(this->evalArg(args, 0));

        if (printOutPriOperations_) {
            const auto& nodePri = static_cast<const PrintOperationNode<ScalarIn>&>(node);
            std::cout << nodePri.getBeforeString() << out << nodePri.getAfterString();
        }

        if (out.getOperationNode() != nullptr) {
            const auto& nodePri = static_cast<const PrintOperationNode<ScalarIn>&>(node);
            ActiveOut out2(*outHandler_->makePrintNode(nodePri.getBeforeString(), *out.getOperationNode(), nodePri.getAfterString()));
            if (out.isValueDefined())
                out2.setValue(out.getValue());
            return out2;
        }

        return out;
    }

    /**
     * @note overrides the default evalAtomicOperation() even though this method
     *        is not virtual (hides a method in EvaluatorOperations)
     */
    void evalAtomicOperation(const NodeIn& node) {
        CGOpCode op = node.getOperationType();
        CPPADCG_ASSERT_KNOWN(op == CGOpCode::AtomicForward || op == CGOpCode::AtomicReverse,
                             "Invalid operation type")

        // check if this node was previously determined
        if (evals_[node] != nullptr) {
            return; // *evals_[node];
        }

        const std::vector<size_t>& info = node.getInfo();
        const std::vector<Argument<ScalarIn> >& inArgs = node.getArguments();

        CPPADCG_ASSERT_KNOWN(info.size() == 3, "Invalid number of information data for atomic operation")
        size_t order_up = info[2];
        size_t p1 = order_up + 1;

        CPPADCG_ASSERT_KNOWN(inArgs.size() == 2 * p1 + 1, "Invalid number of information data for atomic operation")

        if (outHandler_ == nullptr) {
            throw CGException("Evaluator is unable to determine the new CodeHandler for an atomic operation");
        }

        std::vector<Argument<ScalarOut> > outArgs(inArgs.size());

        std::vector<std::vector<ScalarOut>> outVals(inArgs.size());
        bool valuesDefined = true;
        bool allParameters = true;

        for (size_t i = 0; i < inArgs.size(); i++) {
            auto* a = inArgs[i].getOperation();
            CPPADCG_ASSERT_KNOWN(a != nullptr, "Invalid argument for atomic operation")

            outArgs[i] = asArgument(makeArray(*a, outVals[i], valuesDefined, allParameters));
        }

        this->saveEvaluation(node, ActiveOut(*outHandler_->makeNode(op, info, outArgs)));

        if (valuesDefined) {
            const std::map<size_t, CGAbstractAtomicFun<ScalarIn>*>& afun = this->handler_.getAtomicFunctions();
            size_t id = info[0];
            size_t order_low = info[1];
            if (op == CGOpCode::AtomicForward) {
                auto itAFun = afun.find(id);
                if (itAFun == afun.end()) {
                    if (allParameters)
                        throw CGException("Atomic function for ID ", id, " is not defined in evaluator");
                    else
                        return; // the atomic function is not available but it is not essential
                }
                auto& atomic = *itAFun->second;

                //size_t m = outVals[1].size() / (order_up + 1);
                size_t n = outVals[0].size() / (order_up + 1);
                size_t npar = outVals[2].size();
                size_t ndyn = n + npar;

                CppAD::vector<ActiveOut> taylor_x(ndyn * (order_up + 1));
                CppAD::vector<ActiveOut> taylor_y(outVals[1].size());
                CppAD::vector<ActiveOut> par_x(ndyn);

                for (size_t i = 0; i < npar; ++i)
                    par_x[n + i] = ActiveIn(outVals[2][i]);

                CppAD::vector<ad_type_enum> type_x(n);
                std::fill(type_x.data(), type_x.data() + n, ad_type_enum::variable_enum);
                std::fill(type_x.data() + n, type_x.data() + ndyn, ad_type_enum::dynamic_enum);

                auto need_y = size_t(ad_type_enum::variable_enum);

                for (size_t j = 0; j < outVals[0].size(); ++j)
                    taylor_x[j] = ActiveIn(outVals[0][j]);
                for (size_t j = 0; j < npar; ++j)
                    taylor_x[(n + j) * (order_up + 1)] = par_x[n+j];

                for (size_t i = 0; i < taylor_y.size(); ++i)
                    taylor_y[i] = ActiveIn(outVals[1][i]);

                atomic.forward(par_x, type_x, need_y, order_low, order_up, taylor_x, taylor_y);

                std::vector<ScalarOut*>& yOut = atomicEvalResults_[&node];
                assert(yOut.empty());
                yOut.resize(taylor_y.size());
                for (size_t i = 0; i < taylor_y.size(); ++i) {
                    if (taylor_y[i].isValueDefined())
                        yOut[i] = new ScalarOut(taylor_y[i].getValue());
                }
            }
        }

    }

    /**
     * @note overrides the default evalArrayElement() even though this method
     *        is not virtual (hides a method in EvaluatorOperations)
     */
    inline ActiveOut evalArrayElement(const NodeIn& node) {
        // check if this node was previously determined
        if (evals_[node] != nullptr) {
            return *evals_[node];
        }

        const std::vector<ArgIn>& args = node.getArguments();
        const std::vector<size_t>& info = node.getInfo();
        CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for array element")
        CPPADCG_ASSERT_KNOWN(args[0].getOperation() != nullptr, "Invalid argument for array element")
        CPPADCG_ASSERT_KNOWN(args[1].getOperation() != nullptr, "Invalid argument for array element")
        CPPADCG_ASSERT_KNOWN(info.size() == 1, "Invalid number of information data for array element")
        size_t index = info[0];

        ArgOut arrayArg = asArgument(makeArray(*args[0].getOperation()));

        auto& thisOps = static_cast<FinalEvaluatorType&>(*this);
        const NodeIn& atomicNode = *args[1].getOperation();
        thisOps.evalAtomicOperation(atomicNode); // atomic operation
        ArgOut atomicArg = *evals_[atomicNode]->getOperationNode();

        ActiveOut out(*outHandler_->makeNode(CGOpCode::ArrayElement, {index}, {arrayArg, atomicArg}));

        auto it = atomicEvalResults_.find(&atomicNode);
        if (it != atomicEvalResults_.end()) {
            const std::vector<ScalarOut*>& yOut = it->second;
            if (index < yOut.size() && yOut[index] != nullptr)
                out.setValue(*yOut[index]);
        }

        return out;
    }

    inline ActiveOut makeArray(const NodeIn& node) {
        if (node.getOperationType() == CGOpCode::ArrayCreation) {
            return makeDenseArray(node);
        } else {
            return makeSparseArray(node);
        }
    }

    inline ActiveOut makeArray(const NodeIn& node,
                               std::vector<ScalarOut>& values,
                               bool& valuesDefined,
                               bool& allParameters) {
        const CppAD::vector<ActiveOut>* arrayActiveOut;
        ActiveOut result;

        if (node.getOperationType() == CGOpCode::ArrayCreation) {
            result = makeDenseArray(node);
            arrayActiveOut = this->evalsArrays_[node.getHandlerPosition()];
        } else {
            result = makeSparseArray(node);
            arrayActiveOut = this->evalsSparseArrays_[node.getHandlerPosition()];
        }

        processArray(*arrayActiveOut, values, valuesDefined, allParameters);

        return result;
    }

    inline ActiveOut makeDenseArray(const NodeIn& node) {
        CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGOpCode::ArrayCreation, "Invalid array creation operation")
        CPPADCG_ASSERT_KNOWN(node.getHandlerPosition() < this->handler_.getManagedNodesCount(), "this node is not managed by the code handler")

        // check if this node was previously determined
        if (evals_[node] != nullptr) {
            return *evals_[node];
        }

        if (outHandler_ == nullptr) {
            throw CGException("Evaluator is unable to determine the new CodeHandler for an array creation operation");
        }

        // values
        const CppAD::vector<ActiveOut>& array = this->evalArrayCreationOperation(node);

        // makeDenseArray() never called directly by EvaluatorOperations
        return *this->saveEvaluation(node, ActiveOut(*outHandler_->makeNode(CGOpCode::ArrayCreation, {}, asArguments(array))));
    }

    inline ActiveOut makeSparseArray(const NodeIn& node) {
        CPPADCG_ASSERT_KNOWN(node.getOperationType() == CGOpCode::SparseArrayCreation, "Invalid sparse array creation operation")
        CPPADCG_ASSERT_KNOWN(node.getHandlerPosition() < this->handler_.getManagedNodesCount(), "this node is not managed by the code handler")

        // check if this node was previously determined
        if (evals_[node] != nullptr) {
            return *evals_[node];
        }

        if (outHandler_ == nullptr) {
            throw CGException("Evaluator is unable to determine the new CodeHandler for a sparse array creation operation");
        }

        // values
        const CppAD::vector<ActiveOut>& array = this->evalSparseArrayCreationOperation(node);

        // makeSparseArray() never called directly by EvaluatorOperations
        return *this->saveEvaluation(node, ActiveOut(*outHandler_->makeNode(CGOpCode::SparseArrayCreation, node.getInfo(), asArguments(array))));
    }

    static inline void processArray(const CppAD::vector<ActiveOut>& array,
                                    std::vector<ScalarOut>& values,
                                    bool& valuesDefined,
                                    bool& allParameters) {
        values.resize(array.size());
        for (size_t i = 0; i < array.size(); i++) {
            if (!array[i].isValueDefined()) {
                valuesDefined = false;
                allParameters = false;
                break;
            } else {
                values[i] = array[i].getValue();
                if (!array[i].isParameter())
                    allParameters = false;
            }
        }
    }

    static inline bool isParameters(const CppAD::vector<ActiveOut>& tx) {
        for (size_t i = 0; i < tx.size(); i++) {
            if (!tx[i].isParameter()) {
                return false;
            }
        }
        return true;
    }

    static inline bool isValuesDefined(const std::vector<ArgOut>& tx) {
        for (size_t i = 0; i < tx.size(); i++) {
            if (tx[i].getOperationNode() != nullptr) {
                return false;
            }
        }
        return true;
    }

};

/**
 * Specialization of Evaluator for an output active type of CG<Base>
 */
template<class ScalarIn, class ScalarOut>
class Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > : public EvaluatorCG<ScalarIn, ScalarOut, Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > > {
protected:
    using Super = EvaluatorCG<ScalarIn, ScalarOut, Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > >;
public:

    inline Evaluator(CodeHandler<ScalarIn>& handler) :
        Super(handler) {
    }
};

} // END cg namespace
} // END CppAD namespace

#endif
