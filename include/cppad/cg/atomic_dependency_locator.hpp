#ifndef CPPAD_CG_ATOMIC_DEPENDENCY_LOCATOR_INCLUDED
#define CPPAD_CG_ATOMIC_DEPENDENCY_LOCATOR_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2019 Joao Leal
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

class AtomicInputVars {
public:
    /**
     * the outer independent variable indexes which affect a call of that atomic function
     */
    std::set<size_t> outerIndeps;
    /**
     * the outer parameter indexes which affect a call of that atomic function
     */
    std::set<size_t> outerParams;
public:
    inline AtomicInputVars() = default;

    inline AtomicInputVars(AtomicInputVars&&) noexcept = default;

    inline AtomicInputVars(const AtomicInputVars&) = default;

    inline AtomicInputVars(std::set<size_t> outerIndeps,
                           std::set<size_t> outerParams) :
            outerIndeps(std::move(outerIndeps)),
            outerParams(std::move(outerParams)) {
    }

    inline AtomicInputVars& operator=(AtomicInputVars&&) = default;

    inline AtomicInputVars& operator=(const AtomicInputVars&) = default;

    virtual ~AtomicInputVars() = default;
};

/**
 * Utility class that holds information on how atomic functions are used
 */
template<class Base>
class AtomicUseInfo {
public:
    /**
     * the atomic function
     */
    CGAbstractAtomicFun<Base>* atom;
    /**
     * pairs of independent and dependent sizes
     */
    std::set<std::pair<size_t, size_t>> sizes;
    /**
     * the outer variable indexes which affect a call of that atomic function
     */
    AtomicInputVars outerVars;
public:
    inline AtomicUseInfo() :
            atom(nullptr) {
    }
};

/**
 * Finds atomic functions in a CppAD tape and collects some information on how
 * they are used.
 */
template<class Base>
class AtomicDependencyLocator {
private:
    ADFun<CG<Base> >& fun_;
    std::map<size_t, AtomicUseInfo<Base>> atomicInfo_;
    std::map<OperationNode<Base>*, AtomicInputVars> indeps_;
    CodeHandler<Base> handler_;
public:

    explicit inline AtomicDependencyLocator(ADFun<CG<Base> >& fun) :
        fun_(fun) {
    }

    inline const std::map<size_t, AtomicUseInfo<Base>>& findAtomicsUsage() {
        if (!atomicInfo_.empty()) {
            return atomicInfo_;
        }

        size_t m = fun_.Range();
        size_t n = fun_.Domain();

        std::vector<CG<Base> > x(n);
        handler_.makeVariables(x);

        std::vector<CG<Base> > p(fun_.size_dyn_ind());
        handler_.makeParameters(p);

        fun_.new_dynamic(p);

        // make sure the position in the code handler is the same as the independent index
        assert(x.size() == 0 || (x[0].getOperationNode()->getHandlerPosition() == 0 && x[x.size() - 1].getOperationNode()->getHandlerPosition() == x.size() - 1));

        std::vector<CG<Base> > dep = fun_.Forward(0, x);

        for (size_t i = 0; i < m; i++) {
            findAtomicsUsage(dep[i].getOperationNode());
        }

        const auto& regAtomics = handler_.getAtomicFunctions();
        for (auto& pair: atomicInfo_) {
            size_t id = pair.first;

            pair.second.atom = regAtomics.at(id);
        }

        return atomicInfo_;
    }

private:

    inline AtomicInputVars findAtomicsUsage(OperationNode<Base>* node) {
        if (node == nullptr)
            return AtomicInputVars({}, {});

        CGOpCode op = node->getOperationType();
        if (op == CGOpCode::Inv) {
            // particular case where the position in the code handler is the same as the independent index
            return AtomicInputVars({node->getHandlerPosition()}, {});
        } else if (op == CGOpCode::InvPar) {
            // particular case where the position in the code handler is the same as (parameter index + number of independents)
            return AtomicInputVars({}, {node->getHandlerPosition() - fun_.Domain()});
        }

        if (handler_.isVisited(*node)) {
            // been here before
            return indeps_.at(node);
        }

        handler_.markVisited(*node);

        AtomicInputVars atomIn({}, {});
        const std::vector<Argument<Base> >& args = node->getArguments();
        for (size_t a = 0; a < args.size(); a++) {
            auto* argNode = args[a].getOperation();
            AtomicInputVars aindeps = findAtomicsUsage(argNode);
            atomIn.outerIndeps.insert(aindeps.outerIndeps.begin(), aindeps.outerIndeps.end());
            atomIn.outerParams.insert(aindeps.outerParams.begin(), aindeps.outerParams.end());
        }
        indeps_[node] = atomIn;

        if (op == CGOpCode::AtomicForward) {
            CPPADCG_ASSERT_UNKNOWN(node->getInfo().size() > 1)
            CPPADCG_ASSERT_UNKNOWN(node->getArguments().size() > 1)
            size_t id = node->getInfo()[0];

#ifndef NDEBUG
            size_t p = node->getInfo()[2];
            CPPADCG_ASSERT_UNKNOWN(p == 0)
#endif

            OperationNode<Base>* tx = node->getArguments()[0].getOperation();
            OperationNode<Base>* ty = node->getArguments()[1].getOperation();

            CPPADCG_ASSERT_UNKNOWN(tx != nullptr && tx->getOperationType() == CGOpCode::ArrayCreation)
            CPPADCG_ASSERT_UNKNOWN(ty != nullptr && ty->getOperationType() == CGOpCode::ArrayCreation)

            auto& info = atomicInfo_[id];
            info.outerVars.outerIndeps.insert(atomIn.outerIndeps.begin(), atomIn.outerIndeps.end());
            info.outerVars.outerParams.insert(atomIn.outerParams.begin(), atomIn.outerParams.end());
            info.sizes.insert(std::pair<size_t, size_t>(tx->getArguments().size(),
                                                        ty->getArguments().size()));
        }

        return atomIn;
    }
};

} // END cg namespace
} // END CppAD namespace

#endif