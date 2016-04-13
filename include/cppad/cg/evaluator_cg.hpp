#ifndef CPPAD_CG_EVALUATOR_CG_INCLUDED
#define CPPAD_CG_EVALUATOR_CG_INCLUDED
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
 * Specialization of Evaluator for an output active type of CG<Base>
 */
template<class ScalarIn, class ScalarOut>
class Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > : public EvaluatorOperations<ScalarIn, ScalarOut, CG<ScalarOut>, Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > > {
    /**
     * must be friends with one of its super classes since there is a cast to
     * this type due to the curiously recurring template pattern (CRTP)
     */
    friend EvaluatorBase<ScalarIn, ScalarOut, CG<ScalarOut>, Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > >;
public:
    typedef CG<ScalarOut> ActiveOut;
protected:
    typedef EvaluatorOperations<ScalarIn, ScalarOut, CG<ScalarOut>, Evaluator<ScalarIn, ScalarOut, CG<ScalarOut> > > Super;
public:

    inline Evaluator(CodeHandler<ScalarIn>& handler) :
        Super(handler) {
    }

protected:

    /**
     * @note overrides the default proccessActiveOut() even though this method
     *        is not virtual (hides a method in EvaluatorOperations)
     */
    void proccessActiveOut(const OperationNode<ScalarIn>& node,
                           ActiveOut& a) {
        if (node.getName() != nullptr) {
            if (a.getOperationNode() != nullptr) {
                a.getOperationNode()->setName(*node.getName());
            }
        }
    }

};

} // END cg namespace
} // END CppAD namespace

#endif