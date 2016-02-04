#ifndef CPPAD_CG_EVALUATOR_ADCG_INCLUDED
#define CPPAD_CG_EVALUATOR_ADCG_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
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
 * Specialization for an output active type of AD<CG<Base>>
 */
template<class ScalarIn, class BaseOut>
class Evaluator<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> > > : public EvaluatorBase<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> > > {
public:
    typedef CG<BaseOut> ScalarOut;
    typedef CppAD::AD<ScalarOut> ActiveOut;
protected:
    typedef EvaluatorBase<ScalarIn, ScalarOut, ActiveOut> Super;
    using Super::evalsAtomic_;
    using Super::atomicFunctions_;
    using Super::handler_;
    using Super::evalArrayCreationOperation;
public:

    inline Evaluator(CodeHandler<ScalarIn>& handler) :
        Super(handler) {
    }

    inline virtual ~Evaluator() {
    }

protected:

    virtual void proccessActiveOut(OperationNode<ScalarIn>& node,
                                   ActiveOut& a) override {
        if (node.getName() != nullptr) {
            ScalarOut a2(CppAD::Value(a));
            if (a2.getOperationNode() != nullptr) {
                a2.getOperationNode()->setName(*node.getName());
            }
        }
    }

    /**
     * @throws CGException on an internal evaluation error
     */
    virtual void evalAtomicOperation(OperationNode<ScalarIn>& node) override {
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