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
#include <cppad/cg/dae_index_reduction/bipartite.hpp>
#include <cppad/cg/dae_index_reduction/dae_equation_info.hpp>
#include <cppad/cg/dae_index_reduction/time_diff.hpp>

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
    // typical values;
    std::vector<Base> x_;
    // original sparsity
    std::vector<bool> sparsity_;
    // Bipartite graph ([equation i][variable j])
    std::vector<Vnode<Base>*> vnodes_;
    std::vector<Enode<Base>*> enodes_;
    // new index reduced model
    ADFun<CGBase>* reducedFun_;
    // the maximum order of the time derivatives in the original model
    int origMaxTimeDivOrder_;
    // the number of time dependent variables in the original model
    size_t origTimeDependentCount_;
private:
    int timeOrigVarIndex_; // time index in the original user model (may not exist)
public:

    /**
     * Creates the DAE index reduction algorithm that implements the 
     * Pantelides method.
     * 
     * @param fun The DAE model
     * @param varInfo DAE model variable classification
     * @param x typical variable values (used to avoid NaNs in CppAD checks)
     */
    Pantelides(ADFun<CG<Base> >* fun,
               const std::vector<DaeVarInfo>& varInfo,
               const std::vector<Base>& x) :
        DaeIndexReduction<Base>(fun, varInfo),
        x_(x),
        reducedFun_(nullptr),
        origMaxTimeDivOrder_(0),
        origTimeDependentCount_(0),
        timeOrigVarIndex_(-1) {

        using namespace CppAD::extra;
        using namespace std;
        using std::vector;

        const size_t m = fun->Range(); // equation count
        const size_t n = fun->Domain(); // total variable count

        // create equation nodes
        enodes_.reserve(1.2 * m + 1);
        enodes_.resize(m);
        for (size_t i = 0; i < m; i++) {
            enodes_[i] = new Enode<Base > (i);
        }

        // locate the time variable (if present)
        for (size_t dj = 0; dj < n; dj++) {
            if (this->varInfo_[dj].isIntegratedVariable()) {
                if (timeOrigVarIndex_ >= 0) {
                    throw CGException("More than one time variable (integrated variable) defined");
                }
                timeOrigVarIndex_ = dj;
            }
        }

        // determine the order of each time derivative
        vector<int> derivOrder = determineVariableDiffOrder(this->varInfo_);
        map<int, vector<size_t> > order2Tape;
        for (size_t tape = 0; tape < derivOrder.size(); ++tape) {
            order2Tape[derivOrder[tape]].push_back(tape);
        }
        origMaxTimeDivOrder_ = *std::max_element(derivOrder.begin(), derivOrder.end());

        /**
         * generate names for the variables
         */
        std::string timeVarName;
        if (timeOrigVarIndex_ < 0) {
            timeVarName = "t";
        } else {
            if (this->varInfo_[timeOrigVarIndex_].getName().empty()) {
                this->varInfo_[timeOrigVarIndex_].setName("t");
            }
            timeVarName = this->varInfo_[timeOrigVarIndex_].getName();
        }

        stringstream ss;
        for (int order = 0; order <= origMaxTimeDivOrder_; order++) {
            //size_t j = 0; j < this->varInfo_.size(); j++
            const vector<size_t>& tapeIndexes = order2Tape[order];
            if (order < 0) {
                for (size_t p = 0; p < tapeIndexes.size(); ++p) {
                    DaeVarInfo& var = this->varInfo_[tapeIndexes[p]];
                    if (var.getName().empty()) {
                        ss << "p" << p;
                        var.setName(ss.str());
                        ss.str("");
                        ss.clear();
                    }
                }

            } else if (order == 0) {
                for (size_t p = 0; p < tapeIndexes.size(); ++p) {
                    DaeVarInfo& var = this->varInfo_[tapeIndexes[p]];
                    if (var.getName().empty()) {
                        ss << "x" << p;
                        var.setName(ss.str());
                        ss.str("");
                        ss.clear();
                    }
                }
            } else if (order > 0) {
                for (size_t p = 0; p < tapeIndexes.size(); ++p) {
                    DaeVarInfo& var = this->varInfo_[tapeIndexes[p]];
                    if (var.getName().empty()) {
                        const DaeVarInfo& deriv = this->varInfo_[var.getAntiDerivative()];
                        var.setName("d" + deriv.getName() + "d" + timeVarName);
                    }
                }
            }
        }

        // sort the variables according to the time derivative order (constants are kept out)
        vector<size_t> new2Tape;
        vector<int> tape2New(n, -1);
        new2Tape.reserve(n);
        for (int order = 0; order <= origMaxTimeDivOrder_; order++) {
            const vector<size_t>& tapeIndexes = order2Tape[order];
            for (size_t p = 0; p < tapeIndexes.size(); ++p) {
                size_t tapeIndex = tapeIndexes[p];
                tape2New[tapeIndex] = new2Tape.size();
                new2Tape.push_back(tapeIndex);
            }
        }

        // create the variable nodes
        origTimeDependentCount_ = new2Tape.size();
        vnodes_.resize(origTimeDependentCount_);
        for (size_t j = 0; j < vnodes_.size(); j++) {
            size_t tapeIndex = new2Tape[j];
            int tapeIndex0 = this->varInfo_[tapeIndex].getAntiDerivative();
            const std::string& name = this->varInfo_[tapeIndex].getName();

            CPPADCG_ASSERT_UNKNOWN(this->varInfo_[tapeIndex].isFunctionOfIntegrated());

            if (tapeIndex0 < 0) {
                // generate the variable name
                vnodes_[j] = new Vnode<Base > (j, tapeIndex, name);
            } else {
                Vnode<Base>* derivativeOf = vnodes_[tape2New[tapeIndex0]];
                vnodes_[j] = new Vnode<Base > (j, tapeIndex, derivativeOf, name);
            }
        }

        // create the edges
        sparsity_ = jacobianSparsity < vector<bool>, CGBase > (*fun);

        for (size_t i = 0; i < m; i++) {
            for (size_t p = 0; p < n; p++) {
                int j = tape2New[p];
                if (j >= 0 && sparsity_[i * n + p]) {
                    enodes_[i]->addVariable(vnodes_[j]);
                }
            }
        }

        // make sure the system is not under or over determined
        size_t nvar = 0;
        for (size_t j = 0; j < vnodes_.size(); j++) {
            const Vnode<Base>* jj = vnodes_[j];
            if (!jj->isParameter() && // exclude constants
                    (jj->antiDerivative() != nullptr || // derivatives
                    jj->derivative() == nullptr) // algebraic variables
                    ) {
                nvar++;
            }
        }

        if (nvar != m) {
            stringstream ss;
            ss << "The system is not well determined. "
                    "The of number of equations (" << enodes_.size() << ") does not match the number of unknown variables (" << nvar << ").";
            throw CGException(ss.str());
        }
    }

    Pantelides(const Pantelides& p) = delete;
    Pantelides& operator=(const Pantelides& p) = delete;

    /**
     * Performs the DAE differentiation index reductions
     * 
     * @param newVarInfo Variable related information of the reduced index
     *                   model
     * @param equationInfo Equation related information of the reduced index
     *                     model
     * @return the reduced index model (must be deleted by user)
     */
    virtual inline ADFun<CG<Base> >* reduceIndex(std::vector<DaeVarInfo>& newVarInfo,
                                                 std::vector<DaeEquationInfo>& equationInfo) throw (CGException) {
        if (reducedFun_ != nullptr)
            throw CGException("reduceIndex() can only be called once!");

        detectSubset2Dif();

        if (this->verbosity_ >= VERBOSITY_HIGH)
            printResultInfo();

        generateNewModel(newVarInfo, equationInfo);

        return reducedFun_;
    }

    /**
     * Provides the differentiation index. It can only be called after
     * reduceIndex().
     * 
     * @return the DAE differentiation index.
     */
    inline size_t getDifferentiationIndex() const throw (CGException) {
        size_t origM = this->fun_->Range();
        if (origM == enodes_.size()) {
            // no index reduction performed: it is either an index 1 DAE or an ODE
            bool isDAE = false;
            for (size_t j = 0; j < this->varInfo_.size(); j++) {
                const DaeVarInfo& jj = this->varInfo_[j];
                if (jj.getDerivative() < 0 && !jj.isIntegratedVariable() && jj.isFunctionOfIntegrated()) {
                    isDAE = true; // found algebraic variable
                    break;
                }
            }
            if (!isDAE) {
                return 0;
            } else {
                return 1;
            }
        }

        size_t index = 0;
        for (size_t i = origM; i < enodes_.size(); i++) {
            Enode<Base>* ii = enodes_[i];
            size_t eqOrder = 0;
            if (ii->derivative() == nullptr) {
                Enode<Base>* eq = ii;
                while (eq->derivativeOf() != nullptr) {
                    eq = eq->derivativeOf();
                    eqOrder++;
                }
                if (eqOrder > index)
                    index = eqOrder;
            }
        }
        return index + 1; // one extra differentiation to get an ODE
    }

    inline void printResultInfo() {
        std::cout << "\nPantelides DAE differentiation index reduction:\n\n"
                "   Equations count: " << enodes_.size() << "\n";
        for (Enode<Base>* ii : enodes_) {
            std::cout << "      " << ii->index() << " - " << *ii << "\n";
        }

        std::cout << "\n   Variable count: " << vnodes_.size() << "\n";

        for (const Vnode<Base>* jj : vnodes_) {
            std::cout << "      " << jj->index() << " - " << *jj;
            if (jj->assigmentEquation() != nullptr) {
                std::cout << " assigned to " << *jj->assigmentEquation() << "\n";
            } else if (jj->isParameter()) {
                std::cout << " is a parameter (time independent)\n";
            } else {
                std::cout << " NOT assigned to any equation\n";
            }
        }

        std::cout << "\n   Degrees of freedom: " << vnodes_.size() - enodes_.size() << "\n";
    }

    virtual ~Pantelides() {
        for (size_t i = 0; i < enodes_.size(); i++)
            delete enodes_[i];

        for (size_t j = 0; j < vnodes_.size(); j++)
            delete vnodes_[j];
    }

protected:

    static inline std::vector<int> determineVariableDiffOrder(const std::vector<DaeVarInfo>& varInfo) {
        size_t n = varInfo.size();
        // determine the order of each time derivative
        std::vector<int> derivOrder(n, 0);
        for (size_t dj = 0; dj < n; dj++) {
            size_t j0;
            derivOrder[dj] = determineVariableDiffOrder(varInfo, dj, j0);
        }

        return derivOrder;
    }

    static inline int determineVariableDiffOrder(const std::vector<DaeVarInfo>& varInfo, size_t index, size_t& j0) {
        int derivOrder = -1;
        j0 = index;
        if (varInfo[index].isFunctionOfIntegrated()) {
            derivOrder = 0;
            while (varInfo[j0].getAntiDerivative() >= 0) {
                CPPADCG_ASSERT_UNKNOWN(j0 < varInfo.size());
                CPPADCG_ASSERT_UNKNOWN(varInfo[j0].isFunctionOfIntegrated());
                derivOrder++;
                j0 = varInfo[j0].getAntiDerivative();
            }
        }

        return derivOrder;
    }

    /**
     * 
     */
    inline void detectSubset2Dif() {
        Enode<Base>* ll;

        size_t Ndash = enodes_.size();
        for (size_t k = 0; k < Ndash; k++) {
            Enode<Base>* i = enodes_[k];

            if (this->verbosity_ >= VERBOSITY_HIGH)
                std::cout << "Outer loop: equation k = " << *i << "\n";

            bool pathfound = false;
            while (!pathfound) {

                /**
                 * delete all V-nodes with A!=0 and their incident edges
                 * from the graph
                 */
                for (Vnode<Base>* jj : vnodes_) {
                    if (!jj->isDeleted() && jj->derivative() != nullptr) {
                        jj->deleteNode(this->verbosity_);
                    }
                }

                uncolorAll();

                pathfound = augmentPath(*i);

                if (!pathfound) {
                    const size_t vsize = vnodes_.size(); // the size might change
                    for (size_t l = 0; l < vsize; ++l) {
                        Vnode<Base>* jj = vnodes_[l];
                        if (jj->isColored() && !jj->isDeleted()) {
                            // add new variable derivatives of colored variables
                            size_t newVarCount = vnodes_.size() - origTimeDependentCount_;
                            size_t tapeIndex = this->varInfo_.size() + newVarCount;

                            Vnode<Base>* jDiff = new Vnode<Base > (vnodes_.size(), tapeIndex, jj);
                            vnodes_.push_back(jDiff);

                            if (this->verbosity_ >= VERBOSITY_HIGH)
                                std::cout << "Created " << *jDiff << "\n";

                        }
                    }

                    const size_t esize = enodes_.size(); // the size might change
                    for (size_t l = 0; l < esize; l++) {
                        ll = enodes_[l];
                        if (ll->isColored()) {
                            // add new derivative equations for colored equations
                            Enode<Base>* lDiff = new Enode<Base > (enodes_.size(), ll);
                            enodes_.push_back(lDiff);

                            // differentiate newI and create edges!!!
                            dirtyDifferentiateEq(*ll, *lDiff);

                            if (this->verbosity_ >= VERBOSITY_HIGH)
                                std::cout << "Created " << *lDiff << "\n";

                        }
                    }

                    // structural check to avoid infinite recursion
                    for (size_t l = esize; l < enodes_.size(); l++) {
                        ll = enodes_[l];
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

                    for (Vnode<Base>* jj : vnodes_) {
                        if (jj->isColored() && !jj->isDeleted()) {
                            Vnode<Base>* jDiff = jj->derivative();
                            jDiff->setAssigmentEquation(*jj->assigmentEquation()->derivative(), this->verbosity_);
                        }
                    }

                    i = i->derivative();

                    if (this->verbosity_ >= VERBOSITY_HIGH)
                        std::cout << "Set current equation to (i=" << i->index() << ") " << *i << "\n";

                }
            }

        }

    }

    /**
     * 
     * @param i The equation node
     * @return true if an augmented path was found
     */
    bool augmentPath(Enode<Base>& i) {
        i.color(this->verbosity_);

        const std::vector<Vnode<Base>*>& vars = i.variables();

        // first look for derivative variables
        for (Vnode<Base>* jj : vars) {
            if (jj->antiDerivative() != nullptr && jj->assigmentEquation() == nullptr) {
                jj->setAssigmentEquation(i, this->verbosity_);
                return true;
            }
        }

        // look for algebraic variables
        for (Vnode<Base>* jj : vars) {
            if (jj->antiDerivative() == nullptr && jj->assigmentEquation() == nullptr) {
                jj->setAssigmentEquation(i, this->verbosity_);
                return true;
            }
        }


        for (Vnode<Base>* jj : vars) {
            if (!jj->isColored()) {
                jj->color(this->verbosity_);

                Enode<Base>& k = *jj->assigmentEquation();
                bool pathFound = augmentPath(k);
                if (pathFound) {
                    jj->setAssigmentEquation(i, this->verbosity_);
                    return true;
                }
            }
        }

        return false;
    }

    inline void uncolorAll() {
        for (Vnode<Base>* j : vnodes_) {
            j->uncolor();
        }

        for (Enode<Base>* i : enodes_) {
            i->uncolor();
        }
    }

    /**
     * Adds a new equation assuming the new equation differential contains 
     * all variables present in the original equation and their time
     * derivatives (not exactly correct but it works because the 
     * potentially extra variables are removed later)
     * 
     * @param i equation node to differentiate
     */
    inline void dirtyDifferentiateEq(Enode<Base>& i, Enode<Base>& newI) throw (CGException) {
        for (Vnode<Base>* jj : i.originalVariables()) {
            newI.addVariable(jj);
            if (jj->derivative() != nullptr) {
                newI.addVariable(jj->derivative());
            }
        }
    }

    /**
     * Creates a new tape for the index 1 model
     */
    inline void generateNewModel(std::vector<DaeVarInfo>& newVarInfo,
                                 std::vector<DaeEquationInfo>& equationInfo) {
        using std::vector;

        vector<vector<Enode<Base>*> > newEquations;

        // find new equations that must be generated by differentiation
        vector<Enode<Base>* > newEqs;
        size_t origM = this->fun_->Range();
        for (size_t i = 0; i < origM; i++) {
            if (enodes_[i]->derivative() != nullptr) {
                CPPADCG_ASSERT_UNKNOWN(enodes_[i]->derivativeOf() == nullptr);
                newEqs.push_back(enodes_[i]->derivative());
            }
        }

        while (newEqs.size() > 0) {
            newEquations.push_back(newEqs);
            newEqs.clear();
            vector<Enode<Base>*>& eqs = newEquations.back();
            for (size_t i = 0; i < eqs.size(); i++) {
                if (eqs[i]->derivative() != nullptr) {
                    newEqs.push_back(eqs[i]->derivative());
                }
            }
        }

        if (newEquations.empty()) {
            // nothing to do
            return;
        }

        /**
         * Add the relationship between variables and derivatives
         */
        CPPADCG_ASSERT_UNKNOWN(reducedFun_ == nullptr);

        /**
         * Prepare the output information
         */
        newVarInfo = this->varInfo_; // copy
        size_t newVars = vnodes_.size() - origTimeDependentCount_;
        newVarInfo.reserve(this->varInfo_.size() + newVars);
        for (size_t j = origTimeDependentCount_; j < vnodes_.size(); j++) {
            // new variable derivative added by the Pantelides method
            Vnode<Base>* jj = vnodes_[j];
            CPPADCG_ASSERT_UNKNOWN(jj->antiDerivative() != nullptr);
            size_t antiDeriv = jj->antiDerivative()->tapeIndex();
            size_t id = newVarInfo.size();
            newVarInfo.push_back(DaeVarInfo(antiDeriv, jj->name(), id)); // create the new variable
            DaeVarInfo& newVar = newVarInfo.back();
            DaeVarInfo& newAntiDeriv = newVarInfo[antiDeriv];

            newAntiDeriv.setDerivative(jj->tapeIndex()); // update the antiderivative
            newVar.setOrder(newAntiDeriv.getOrder() + 1);
            newVar.setOriginalAntiDerivative(newVar.getOrder() == 1 ? newAntiDeriv.getOriginalIndex() : newAntiDeriv.getOriginalAntiDerivative());
            if (jj->derivative() != nullptr) {
                newVar.setDerivative(jj->derivative()->tapeIndex());
            }
        }

        std::map<Enode<Base>*, Vnode<Base>*> assigments;
        for (size_t j = 0; j < vnodes_.size(); j++) {
            Vnode<Base>* jj = vnodes_[j];
            if (jj->assigmentEquation() != nullptr) {
                assigments[jj->assigmentEquation()] = jj;
            }
        }

        equationInfo.resize(enodes_.size());
        for (size_t i = 0; i < enodes_.size(); i++) {
            Enode<Base>* ii = enodes_[i];
            int derivativeOf = ii->derivativeOf() != nullptr ? ii->derivativeOf()->index() : -1;
            int origIndex = ii->derivativeOf() == nullptr ? i : -1;
            int assignedVarIndex = assigments.count(ii) > 0 ? assigments[ii]->tapeIndex() : -1;

            equationInfo[i] = DaeEquationInfo(i, origIndex, derivativeOf, assignedVarIndex);
        }

        size_t timeTapeIndex;
        {
            CodeHandler<Base> handler;

            vector<CGBase> indep0(this->fun_->Domain());
            handler.makeVariables(indep0);

            const vector<CGBase> dep0 = this->fun_->Forward(0, indep0);

            /**
             * generate a new tape
             */

            vector<ADCG> indepNew;
            if (timeOrigVarIndex_ >= 0) {
                indepNew = vector<ADCG > (newVarInfo.size()); // variables + time (vnodes include time)
                timeTapeIndex = timeOrigVarIndex_;
            } else {
                indepNew = vector<ADCG > (newVarInfo.size() + 1); // variables + time (new time variable added)
                timeTapeIndex = indepNew.size() - 1;
            }

            // initialize with the user provided values
            for (size_t j = 0; j < x_.size(); j++) {
                indepNew[j] = x_[j];
            }
            Independent(indepNew);

            // variables with the relationship between x dxdt and t
            vector<ADCG> indep2 = prepareTimeDependentVariables(indepNew, newVarInfo, timeTapeIndex);
            indep2.resize(indep0.size());

            Evaluator<Base, CGBase> evaluator0(handler, dep0);
            vector<ADCG> depNew = evaluator0.evaluate(indep2);
            depNew.resize(enodes_.size());

            try {
                reducedFun_ = new ADFun<CGBase > (indepNew, depNew);
            } catch (const std::exception& ex) {
                throw CGException(std::string("Failed to create ADFun: ") + ex.what());
            }

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << "Original model:\n";
                printModel(reducedFun_, newVarInfo);
            }
        }


        /**
         * generate the system of equations by repeatedly differentiating
         * and adding equations to the DAE system
         */
        for (size_t d = 0; d < newEquations.size(); d++) {
            vector<Enode<Base>*>& equations = newEquations[d];

            size_t m = reducedFun_->Domain(); // total variable count
            //size_t n = reducedFun_->Range(); // equation count

            /**
             * register operations from the other equations
             */
            CodeHandler<Base> handler0;

            vector<CGBase> indep0(m);
            handler0.makeVariables(indep0);

            vector<CGBase> dep = reducedFun_->Forward(0, indep0);

            /**
             * register operations used to differentiate the equations
             */
            //forwardTimeDiff(equations, dep, timeTapeIndex);
            reverseTimeDiff(equations, dep, timeTapeIndex);

            delete reducedFun_; // not needed anymore
            reducedFun_ = nullptr;

            /**
             * reconstruct the new system of equations 
             */
            vector<ADCG> indep2;
            vector<ADCG> indepNew;

            if (d < newEquations.size() - 1) {
                indepNew.resize(m);
            } else if (timeOrigVarIndex_ < 0) {
                // the very last model creation
                indepNew.resize(m - 1); // take out time (it was added by this function and not the user)
            } else {
                // the very last model creation
                indepNew.resize(m);
            }

            for (size_t j = 0; j < x_.size(); j++) {
                indepNew[j] = x_[j];
            }
            Independent(indepNew);

            if (d < newEquations.size() - 1) {
                // variables with the relationship between x, dxdt and t
                indep2 = prepareTimeDependentVariables(indepNew, newVarInfo, timeTapeIndex);
            } else {
                indep2 = indepNew;
                indep2.resize(m);
            }

            Evaluator<Base, CGBase> evaluator(handler0, dep);
            vector<ADCG> depNew = evaluator.evaluate(indep2);

            try {
                reducedFun_ = new ADFun<CGBase > (indepNew, depNew);
            } catch (const std::exception& ex) {
                throw CGException(std::string("Failed to create ADFun: ") + ex.what());
            }

            if (this->verbosity_ >= VERBOSITY_HIGH) {
                std::cout << equations.size() << " new equations:\n";
                printModel(reducedFun_, newVarInfo);
            }
        }

    }

    inline void forwardTimeDiff(const std::vector<Enode<Base>*>& equations,
                                std::vector<CG<Base> >& dep,
                                size_t tapeTimeIndex) const {

        size_t m = reducedFun_->Domain();

        std::vector<CGBase> u(m, CGBase(0));
        u[tapeTimeIndex] = CGBase(1);
        std::vector<CGBase> v;
        try {
            v = reducedFun_->Forward(1, u);
        } catch (const std::exception& ex) {
            throw CGException(std::string("Failed to determine model Jacobian (forward mode): ") + ex.what());
        }

        for (size_t e = 0; e < equations.size(); e++) {
            dep[equations[e]->index()] = v[equations[e]->derivativeOf()->index()];
        }
    }

    inline void reverseTimeDiff(const std::vector<Enode<Base>*>& equations,
                                std::vector<CG<Base> >& dep,
                                size_t tapeTimeIndex) const {
        size_t m = reducedFun_->Domain();
        size_t n = reducedFun_->Range();
        std::vector<CGBase> u(m);
        std::vector<CGBase> v(n);

        for (size_t e = 0; e < equations.size(); e++) {
            size_t i = equations[e]->derivativeOf()->index();
            if (reducedFun_->Parameter(i)) { // return zero for this component of f
                dep[equations[e]->index()] = 0;
            } else {
                // set v to the i-th coordinate direction
                v[i] = 1;

                // compute the derivative of this component of f
                try {
                    u = reducedFun_->Reverse(1, v);
                } catch (const std::exception& ex) {
                    throw CGException(std::string("Failed to determine model Jacobian (reverse mode): ") + ex.what());
                }

                // reset v to vector of all zeros
                v[i] = 0;

                // return the result
                dep[equations[e]->index()] = u[tapeTimeIndex];
            }
        }
    }

    /**
     * Introduces a dependency with respect to time in the provided
     * variables.
     * 
     * @param indepOrig  The variables without time dependency 
     *                    (in the original variable order).
     * @return The new variables with the time dependency 
     *          (in the original variable order).
     */
    inline std::vector<CppAD::AD<CG<Base> > > prepareTimeDependentVariables(const std::vector<ADCG>& indepOrig,
                                                                            const std::vector<DaeVarInfo>& newVarInfo,
                                                                            size_t timeTapeIndex) const {
        CPPADCG_ASSERT_UNKNOWN(timeTapeIndex < indepOrig.size());

        using std::vector;
        typedef CppAD::AD<CGBase> ADCGBase;

        vector<ADCGBase> indepOut(indepOrig.size());
        vector<ADCGBase> ax(3); // function inputs
        vector<ADCGBase> ay(1); // function output

        ax[2] = indepOrig[timeTapeIndex]; // time

        for (size_t j = 0; j < newVarInfo.size(); j++) {
            const DaeVarInfo& jj = newVarInfo[j];
            if (jj.getDerivative() >= 0) {
                ax[0] = indepOrig[j]; // x
                ax[1] = indepOrig[jj.getDerivative()]; // dxdt
                time_var(0, ax, ay);
                indepOut[j] = ay[0];
            } else {
                indepOut[j] = indepOrig[j];
            }
        }

        if (newVarInfo.size() < indepOrig.size()) {
            indepOut[indepOut.size() - 1] = indepOrig[timeTapeIndex];
        }

        return indepOut;
    }

    inline void printModel(ADFun<CG<Base> >* fun) {
        printModel(fun, this->varInfo_);
    }

    /**
     * Prints out a DAE model to the standard output.
     * 
     * @param fun  The taped model
     */
    inline static void printModel(ADFun<CG<Base> >* fun, const std::vector<DaeVarInfo>& varInfo) {
        std::vector<std::string> vnames(varInfo.size());
        for (size_t p = 0; p < varInfo.size(); p++) {
            vnames[p] = varInfo[p].getName();
        }
        printModel(fun, vnames);
    }

    /**
     * Prints out a DAE model to the standard output.
     * 
     * @param fun  The taped model
     * @param vnodes  The independent variables
     */
    inline static void printModel(ADFun<CG<Base> >* fun,
                                  const std::vector<std::string>& indepNames) {
        using CppAD::vector;

        CPPADCG_ASSERT_UNKNOWN(fun != nullptr);
        CPPADCG_ASSERT_UNKNOWN(fun->Domain() == indepNames.size() || fun->Domain() == indepNames.size() + 1); // with or without time

        CodeHandler<Base> handler;

        vector<CGBase> indep0(fun->Domain());
        handler.makeVariables(indep0);

        vector<CGBase> dep0 = fun->Forward(0, indep0);

        LanguageC<double> langC("double");

        /**
         * create variable names
         */
        std::vector<std::string> depNames;

        /**
         * generate the source code
         */
        LangCCustomVariableNameGenerator<double> nameGen(depNames, indepNames, "res");

        std::ostringstream code;
        handler.generateCode(code, langC, dep0, nameGen);
        std::cout << "\n" << code.str() << std::endl;
    }

};

} // END cg namespace
} // END CppAD namespace

#endif