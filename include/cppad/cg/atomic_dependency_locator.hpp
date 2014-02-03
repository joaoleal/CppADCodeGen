#ifndef CPPAD_CG_ATOMIC_DEPENDENCY_LOCATOR_INCLUDED
#define CPPAD_CG_ATOMIC_DEPENDENCY_LOCATOR_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

template<class Base>
class AtomicDependencyLocator {
private:
    ADFun<CG<Base> >& fun_;
    std::map<size_t, std::set<size_t> > atomicIndeps_;
    std::map<OperationNode<Base>*, std::set<size_t> > indeps_;
public:

    inline AtomicDependencyLocator(ADFun<CG<Base> >& fun) :
        fun_(fun) {
    }

    inline std::map<size_t, std::set<size_t> > findAtomicsUsage() {
        using CppAD::vector;

        if (!atomicIndeps_.empty()) {
            return atomicIndeps_;
        }

        size_t m = fun_.Range();
        size_t n = fun_.Domain();

        CodeHandler<Base> handler;

        vector<CG<Base> > x(n);
        handler.makeVariables(x);
        for (size_t j = 0; j < n; j++)
            x[j].getOperationNode()->setColor(j);

        vector<CG<Base> > dep = fun_.Forward(0, x);

        for (size_t i = 0; i < m; i++) {
            findAtomicsUsage(dep[i].getOperationNode());
        }

        return atomicIndeps_;
    }

private:

    inline std::set<size_t> findAtomicsUsage(OperationNode<Base>* node) {
        if (node == nullptr)
            return std::set<size_t>();

        CGOpCode op = node->getOperationType();
        if (op == CGInvOp) {
            std::set<size_t> indeps;
            indeps.insert(node->getColor());
            return indeps;
        }

        if (node->getUsageCount() > 0) {
            // been here before
            return indeps_.at(node);
        }

        node->increaseUsageCount();

        std::set<size_t> indeps;
        const std::vector<Argument<Base> >& args = node->getArguments();
        for (size_t a = 0; a < args.size(); a++) {
            std::set<size_t> aindeps = findAtomicsUsage(args[a].getOperation());
            indeps.insert(aindeps.begin(), aindeps.end());
        }
        indeps_[node] = indeps;

        if (op == CGAtomicForwardOp) {
            CPPADCG_ASSERT_UNKNOWN(node->getInfo().size() > 1);
            size_t id = node->getInfo()[0];
            atomicIndeps_[id].insert(indeps.begin(), indeps.end());
        }

        return indeps;
    }
};

} // END cg namespace
} // END CppAD namespace

#endif