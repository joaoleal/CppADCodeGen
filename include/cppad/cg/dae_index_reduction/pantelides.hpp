#ifndef CPPAD_CG_PANTELIDES_INCLUDED
#define CPPAD_CG_PANTELIDES_INCLUDED
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

#include <cppad/cg/dae_index_reduction/dae_index_reduction.hpp>
#include <cppad/cg/dae_index_reduction/augment_path_depth_lookahead.hpp>
#include <cppad/cg/dae_index_reduction/bipartite_graph.hpp>

namespace CppAD {
namespace cg {

/**
 * Pantelides DAE index reduction algorithm
 */
template<class Base>
class Pantelides : public DaeIndexReduction<Base> {
protected:
    typedef CppAD::cg::CG<Base> CGBase;
    typedef CppAD::AD<CGBase> ADCG;
protected:
    //
    BipartiteGraph<Base> graph_;
    // typical values;
    std::vector<Base> x_;
    // new index reduced model
    ADFun<CGBase>* reducedFun_;
    AugmentPathDepthLookahead<Base> defaultAugmentPath_;
    AugmentPath<Base>* augmentPath_;
public:

    /**
     * Creates the DAE index reduction algorithm that implements the 
     * Pantelides method.
     * 
     * @param fun The DAE model
     * @param varInfo DAE model variable classification
     * @param eqName Equation names (it can be an empty vector)
     * @param x Typical variable values (used to avoid NaNs in CppAD checks)
     */
    Pantelides(ADFun<CG<Base> >* fun,
               const std::vector<DaeVarInfo>& varInfo,
               const std::vector<std::string>& eqName,
               const std::vector<Base>& x) :
        DaeIndexReduction<Base>(fun),
        graph_(fun, varInfo, eqName, *this),
        x_(x),
        reducedFun_(nullptr),
        augmentPath_(&defaultAugmentPath_) {

    }

    Pantelides(const Pantelides& p) = delete;

    Pantelides& operator=(const Pantelides& p) = delete;

    virtual ~Pantelides() {
    }

    AugmentPath<Base>& getAugmentPath() const {
        return *augmentPath_;
    }

    void setAugmentPath(AugmentPath<Base>& a) const {
        augmentPath_ = &a;
    }

    /**
     * Defines whether or not original names saved by using
     * CppAD::PrintFor(0, "", val, name)
     * should be kept by also adding PrintFor operations in the reduced model.
     */
    inline void setPreserveNames(bool p) {
        graph_.setPreserveNames(p);
    }

    /**
     * Whether or not original names saved by using
     * CppAD::PrintFor(0, "", val, name)
     * should be kept by also adding PrintFor operations in the reduced model.
     */
    inline bool isPreserveNames() const {
        return graph_.isPreserveNames();
    }

    /**
     * Performs the DAE differentiation index reductions
     * 
     * @param newVarInfo Variable related information of the reduced index
     *                   model
     * @param equationInfo Equation related information of the reduced index
     *                     model
     * @return the reduced index model (must be deleted by user)
     * @throws CGException on failure
     */
    virtual inline ADFun<CG<Base> >* reduceIndex(std::vector<DaeVarInfo>& newVarInfo,
                                                 std::vector<DaeEquationInfo>& equationInfo) override {
        if (reducedFun_ != nullptr)
            throw CGException("reduceIndex() can only be called once!");

        if (this->verbosity_ >= Verbosity::High)
            log() << "########  Pantelides method  ########\n";

        augmentPath_->setLogger(*this);

        detectSubset2Dif();

        if (this->verbosity_ >= Verbosity::High)
            graph_.printResultInfo("Pantelides");

        reducedFun_ = graph_.generateNewModel(newVarInfo, equationInfo, x_);

        return reducedFun_;
    }

    /**
     * Provides the differentiation index. It can only be called after
     * reduceIndex().
     *
     * @return the DAE differentiation index.
     * @throws CGException
     */
    inline size_t getDifferentiationIndex() const {
        return graph_.getDifferentiationIndex();
    }

protected:
    using DaeIndexReduction<Base>::log;

    /**
     * 
     */
    inline void detectSubset2Dif() {
        auto& vnodes = graph_.variables();
        auto& enodes = graph_.equations();

        Enode<Base>* ll;

        size_t Ndash = enodes.size();
        for (size_t k = 0; k < Ndash; k++) {
            Enode<Base>* i = enodes[k];

            if (this->verbosity_ >= Verbosity::High)
                log() << "Outer loop: equation k = " << *i << "\n";

            bool pathfound = false;
            while (!pathfound) {

                /**
                 * delete all V-nodes with A!=0 and their incident edges
                 * from the graph
                 */
                for (Vnode<Base>* jj : vnodes) {
                    if (!jj->isDeleted() && jj->derivative() != nullptr) {
                        jj->deleteNode(log(), this->verbosity_);
                    }
                }

                graph_.uncolorAll();

                pathfound = augmentPath_->augmentPath(*i);

                if (!pathfound) {
                    const size_t vsize = vnodes.size(); // the size might change
                    for (size_t l = 0; l < vsize; ++l) {
                        Vnode<Base>* jj = vnodes[l];
                        if (jj->isColored() && !jj->isDeleted()) {
                            // add new variable derivatives of colored variables
                            graph_.createDerivate(*jj);
                        }
                    }

                    const size_t esize = enodes.size(); // the size might change
                    for (size_t l = 0; l < esize; l++) {
                        ll = enodes[l];
                        if (ll->isColored()) {
                            // add new derivative equations for colored equations and create edges
                            graph_.createDerivate(*ll);
                        }
                    }

                    // structural check to avoid infinite recursion
                    for (size_t l = esize; l < enodes.size(); l++) {
                        ll = enodes[l];
                        const std::vector<Vnode<Base>*>& nvars = ll->originalVariables();
                        bool ok = false;
                        for (Vnode<Base>* js : nvars) {
                            if (js->equations().size() > 1) {
                                ok = true;
                                break;
                            }
                        }
                        if (!ok)
                            throw CGException("Invalid equation structure. The model appears to be over-defined.");
                    }

                    for (Vnode<Base>* jj : vnodes) {
                        if (jj->isColored() && !jj->isDeleted()) {
                            Vnode<Base>* jDiff = jj->derivative();
                            jDiff->setAssignmentEquation(*jj->assigmentEquation()->derivative(), log(), this->verbosity_);
                        }
                    }

                    i = i->derivative();

                    if (this->verbosity_ >= Verbosity::High)
                        log() << "Set current equation to (i=" << i->index() << ") " << *i << "\n";

                }
            }

        }

    }

};

} // END cg namespace
} // END CppAD namespace

#endif