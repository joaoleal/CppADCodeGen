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
class Evaluator<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> > > : public EvaluatorAD<ScalarIn, CG<BaseOut>, Evaluator<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> > > > {
    /**
     * must be friends with one of its super classes since there is a cast to
     * this type due to the  curiously recurring template pattern (CRTP)
     */
    friend EvaluatorBase<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> >, Evaluator<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> > > >;
public:
    typedef CG<BaseOut> ScalarOut;
    typedef CppAD::AD<ScalarOut> ActiveOut;
    typedef EvaluatorAD<ScalarIn, ScalarOut, Evaluator<ScalarIn, CG<BaseOut>, CppAD::AD<CG<BaseOut> > > > Super;
protected:
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

    /**
     * @note overrides the default proccessActiveOut() even though this method
     *        is not virtual (hides a method in EvaluatorOperations)
     */
    void proccessActiveOut(const OperationNode<ScalarIn>& node,
                                   ActiveOut& a) {
        if (node.getName() != nullptr) {
            ScalarOut a2(CppAD::Value(a));
            if (a2.getOperationNode() != nullptr) {
                a2.getOperationNode()->setName(*node.getName());
            }
        }
    }

};

} // END cg namespace
} // END CppAD namespace

#endif