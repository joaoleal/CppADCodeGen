#ifndef CPPAD_CG_LOOP_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_LOOP_ATOMIC_FUN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    /**
     * Independent variable positions
     */
    class LoopPosition {
    public:
        size_t tape;
        size_t atomic;
        size_t original;

        inline LoopPosition() :
            tape(-1),
            atomic(-1),
            original(-1) {
        }

        inline LoopPosition(size_t t, size_t a, size_t o) :
            tape(t),
            atomic(a),
            original(o) {
        }
    };

    /**
     * Temporary variable positions
     */
    class LoopPositionTmp {
    public:
        size_t tape;
        size_t atomic;
        // the independent variables that this temporary variable depends on
        std::set<size_t> originalIndeps;

        inline LoopPositionTmp() :
            tape(-1),
            atomic(-1) {
        }

        inline LoopPositionTmp(size_t t, size_t a, const std::set<size_t>& o) :
            tape(t),
            atomic(a),
            originalIndeps(o) {
        }
    };

    template<class Base>
    class LoopOperationGraph {
    public:
        OperationNode<Base>* loopStart;
        OperationNode<Base>* loopEnd;
        std::vector<OperationNode<Base>*> indexedResults; // y, jac, hess...
        std::vector<OperationNode<Base>*> indexedIndependents;

        LoopOperationGraph() :
            loopStart(NULL),
            loopEnd(NULL) {

        }
    };

    /**
     * An atomic function for source code generation within loops
     * 
     * @author Joao Leal
     */
    template <class Base>
    class LoopAtomicFun : public BaseAbstractAtomicFun<Base> {
    public:
        typedef CppAD::CG<Base> CGB;
        typedef Argument<Base> Arg;
    protected:
        const size_t loopId_;
        /**
         * The tape for a single loop iteration
         */
        ADFun<CGB> * const fun_;
        /**
         * Number of loop iterations
         */
        const size_t iterationCount_;
        /**
         * the total number of independent variables of the atomic function
         * (NOT the tape)
         */
        size_t nFull_;
        /**
         * loop tape dependent variable count (number of equation patterns)
         */
        const size_t m_;
        /**
         * The dependent variables ([tape equation][iteration])
         */
        std::vector<std::vector<LoopPosition> > dependentIndexes_;
        /**
         * The indexed independent variables ([tape variable][iteration])
         */
        std::vector<std::vector<LoopPosition> > indexedIndepIndexes_;
        /**
         * The non-indexed independent variables ([tape variable])
         */
        std::vector<LoopPosition> nonIndexedIndepIndexes_;
        /**
         * The independent variables related with temporary variables of the
         * original model.
         */
        std::vector<LoopPositionTmp> temporaryIndependents_;
        /**
         * Maps the independent variable indexes in the original model to this
         * atomic function
         */
        std::map<size_t, size_t> indepOrig2full_;
        /**
         * 
         */
        std::map<size_t, LoopPosition*> depOrigIndexes_;
        /**
         *
         */
        vector<IndexPattern*> indepIndexPatterns_;
        /**
         * 
         */
        vector<IndexPattern*> depIndexPatterns_;

        vector<std::set<size_t> > jacTape_;
        CustomPosition custom_jac_;
        vector<std::set<size_t> > hessTape_;
        std::map<size_t, vector<std::set<size_t> > > eqHessTape_;
        std::map<size_t, vector<std::set<size_t> > > hess_;
        CustomPosition custom_hess_;
        LoopOperationGraph<Base>* graphForward_;
    public:

        /**
         * Creates a new atomic function that is responsible for defining the
         * dependencies to calls of a user atomic function.
         * 
         * @param name The atomic function name.
         */
        LoopAtomicFun(const std::string& name,
                      ADFun<CGB>* fun,
                      size_t iterationCount,
                      const std::vector<std::vector<size_t> >& dependentOrigIndexes,
                      const std::vector<std::vector<size_t> >& indexedIndepOrigIndexes,
                      const std::vector<size_t>& nonIndexedIndepOrigIndexes,
                      const std::vector<std::set<size_t> >& temporaryIndependents) :
            BaseAbstractAtomicFun<Base>(name),
            loopId_(createNewLoopId()),
            fun_(fun),
            iterationCount_(iterationCount),
            nFull_(0),
            m_(dependentOrigIndexes.size()),
            dependentIndexes_(m_, std::vector<LoopPosition>(iterationCount)),
            indexedIndepIndexes_(indexedIndepOrigIndexes.size(), std::vector<LoopPosition>(iterationCount)),
            nonIndexedIndepIndexes_(nonIndexedIndepOrigIndexes.size()),
            temporaryIndependents_(temporaryIndependents.size()),
            graphForward_(NULL) {
            CPPADCG_ASSERT_KNOWN(fun != NULL, "fun cannot be NULL");

            /**
             * dependents
             */
            for (size_t i = 0; i < m_; i++) {
                for (size_t it = 0; it < iterationCount_; it++) {
                    size_t orig = dependentOrigIndexes[i][it];
                    dependentIndexes_[i][it] = LoopPosition(i, it * m_ + i, orig);
                    depOrigIndexes_[orig] = &dependentIndexes_[i][it];
                }
            }

            /**
             * independents
             */


            size_t nIndexed = indexedIndepOrigIndexes.size();

            // indexed
            for (size_t it = 0; it < iterationCount_; it++) {
                for (size_t j = 0; j < nIndexed; j++) {
                    size_t orig = indexedIndepOrigIndexes[j][it];
                    size_t full;
                    std::map<size_t, size_t>::const_iterator iter = indepOrig2full_.find(orig);
                    if (iter != indepOrig2full_.end()) {
                        full = iter->second;
                    } else {
                        full = indepOrig2full_.size();
                        indepOrig2full_[orig] = full;
                    }

                    indexedIndepIndexes_[j][it] = LoopPosition(j, full, orig);
                }
            }

            // nonindexed
            size_t nNonIndexed = nonIndexedIndepOrigIndexes.size();
            for (size_t j = 0; j < nNonIndexed; j++) {
                size_t orig = nonIndexedIndepOrigIndexes[j];
                size_t full;
                std::map<size_t, size_t>::const_iterator iter = indepOrig2full_.find(orig);
                if (iter != indepOrig2full_.end()) {
                    full = iter->second;
                } else {
                    full = indepOrig2full_.size();
                    indepOrig2full_[orig] = full;
                }

                nonIndexedIndepIndexes_[j] = LoopPosition(nIndexed + j, full, orig);
            }

            // temporary
            size_t fullStart = indepOrig2full_.size();
            for (size_t j = 0; j < temporaryIndependents.size(); j++) {
                temporaryIndependents_[j] = LoopPositionTmp(nIndexed + nNonIndexed + j, fullStart + j, temporaryIndependents[j]);
            }

            nFull_ = fullStart + temporaryIndependents_.size();

        }

        template <class ADVector>
        inline void operator()(const ADVector& ax, ADVector& ay, size_t id = 0) {
            this->atomic_base<CGB>::operator()(ax, ay, id);
        }

        /**
         * Provides a unique identifier for this loop.
         * 
         * @return a unique identifier ID
         */
        inline size_t getLoopId() const {
            return loopId_;
        }

        inline size_t getIterationCount() const {
            return iterationCount_;
        }

        inline size_t getLoopDependentCount() const {
            return m_ * iterationCount_;
        }

        inline size_t getLoopIndependentCount() const {
            return nFull_;
        }

        /**
         * Provides the dependent variables indexes ([tape equation][iteration])
         */
        inline const std::vector<std::vector<LoopPosition> >& getDependentIndexes() const {
            return dependentIndexes_;
        }

        /**
         * Provides the indexed independent variables ([tape variable][iteration])
         */
        inline const std::vector<std::vector<LoopPosition> >& getIndexedIndepIndexes() const {
            return indexedIndepIndexes_;
        }

        /**
         * Provides the non-indexed independent variables ([tape variable])
         */
        inline const std::vector<LoopPosition>& getNonIndexedIndepIndexes() const {
            return nonIndexedIndepIndexes_;
        }

        /**
         * Provides the independent variables related with temporary variables of the
         * original model.
         */
        inline const std::vector<LoopPositionTmp>& getTemporaryIndependents() const {
            return temporaryIndependents_;
        }

        /**
         * Maps the independent variable indexes in the original model to this
         * atomic function
         */
        inline const std::map<size_t, size_t>& getOriginalIndepIndex2Atomic() const {
            return indepOrig2full_;
        }

        inline size_t getTapeDependentIndex(size_t originalDepIndex) const {
            return depOrigIndexes_.at(originalDepIndex)->tape;
        }

        virtual vector<CGB> insertIntoModel(const std::vector<CGB>& independents,
                                            const std::vector<OperationNode<Base>* >& temporaryOrigVarOrder) {
            CodeHandler<Base>* handler = independents[0].handler_;
            assert(handler != NULL);
            assert(temporaryOrigVarOrder.size() == temporaryIndependents_.size());

            vector<CGB> x(nFull_);

            // place independents
            std::map<size_t, size_t>::const_iterator ito2a;
            for (ito2a = indepOrig2full_.begin(); ito2a != indepOrig2full_.end(); ++ito2a) {
                x[ito2a->second] = independents[ito2a->first];
            }
            //place temporaries
            size_t nTmps = temporaryOrigVarOrder.size();
            for (size_t j = 0; j < nTmps; j++) {
                OperationNode<Base>* tmp = temporaryOrigVarOrder[j];
                x[indepOrig2full_.size() + j] = CG<Base> (*handler, Argument<Base>(*tmp));
            }

            std::vector<size_t> opInfo(4);
            opInfo[0] = loopId_;
            opInfo[1] = nFull_;
            opInfo[2] = getLoopDependentCount();
            opInfo[3] = 0; // zero-order
            std::vector<Argument<Base> > args = BaseAbstractAtomicFun<Base>::asArguments(x);

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGLoopForwardOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerLoop(*this); // is this really required?

            size_t mFull = getLoopDependentCount();
            vector<CGB> y(mFull);

            opInfo.resize(1);
            args.resize(1);
            for (size_t i = 0; i < mFull; i++) {
                opInfo[0] = i;
                args[0] = Argument<Base>(*atomicOp);

                y[i] = CGB(*handler, new OperationNode<Base>(CGLoopResultOp, opInfo, args));
            }

            return y;
        }

        inline void detectIndexPatterns() {
            if (indepIndexPatterns_.size() > 0)
                return; // already done

            indepIndexPatterns_.resize(indexedIndepIndexes_.size());
            for (size_t j = 0; j < indepIndexPatterns_.size(); j++) {
                vector<size_t> indexes(iterationCount_);
                for (size_t it = 0; it < iterationCount_; it++) {
                    indexes[it] = indexedIndepIndexes_[j][it].original;
                }
                indepIndexPatterns_[j] = IndexPattern::detect(indexes);
            }

            depIndexPatterns_.resize(dependentIndexes_.size());
            for (size_t j = 0; j < depIndexPatterns_.size(); j++) {
                vector<size_t> indexes(iterationCount_);
                for (size_t it = 0; it < iterationCount_; it++) {
                    indexes[it] = dependentIndexes_[j][it].original;
                }
                depIndexPatterns_[j] = IndexPattern::detect(indexes);
            }
        }

        inline const vector<IndexPattern*>& getDependentIndexPatterns() const {
            return depIndexPatterns_;
        }

        inline const vector<IndexPattern*>& getIndependentIndexPatterns() const {
            return indepIndexPatterns_;
        }

        inline LoopOperationGraph<Base>* getForwardOperationGraph() const {
            return graphForward_;
        }

        virtual LoopOperationGraph<Base>* generateForwardGraph(CodeHandler<Base>& handler,
                                                               size_t p,
                                                               const std::vector<Argument<Base> >& args) {
            assert(graphForward_ == NULL);

            //size_t mTape = fun_->Range();
            size_t nTape = fun_->Domain();
            std::vector<size_t> startEndInfo(2);
            startEndInfo[0] = loopId_;
            startEndInfo[1] = iterationCount_;

            std::vector<size_t> info(1);

            graphForward_ = new LoopOperationGraph<Base>();
            std::vector<Argument<Base> > startArgs(temporaryIndependents_.size());
            for (size_t j = 0; j < temporaryIndependents_.size(); j++) {
                const LoopPositionTmp& pos = temporaryIndependents_[j];
                startArgs[j] = args[pos.atomic];
            }

            graphForward_->loopStart = new OperationNode<Base>(CGLoopStartOp, startEndInfo, startArgs);
            handler.manageOperationNode(graphForward_->loopStart);

            // indexed independents
            size_t nIndexed = indexedIndepIndexes_.size();
            graphForward_->indexedIndependents.resize(nIndexed);
            vector<CGB> tx(nTape);
            std::vector<Argument<Base> > indexedArgs(iterationCount_ + 1);
            indexedArgs[0] = Argument<Base>(*graphForward_->loopStart);
            for (size_t j = 0; j < nIndexed; j++) {
                for (size_t it = 0; it < iterationCount_; it++) {
                    //assert(args[j].getOperation() != NULL);
                    indexedArgs[it + 1] = args[indexedIndepIndexes_[j][it].atomic];
                }
                info[0] = j;
                graphForward_->indexedIndependents[j] = new OperationNode<Base>(CGLoopIndexedIndepOp, info, indexedArgs);
                tx[j] = CGB(handler, graphForward_->indexedIndependents[j]);
            }
            // non indexed
            size_t nNonIndexed = nonIndexedIndepIndexes_.size();
            for (size_t j = 0; j < nNonIndexed; j++) {
                tx[nIndexed + j] = CGB(handler, args[nonIndexedIndepIndexes_[j].atomic]);
            }
            // temporaries
            for (size_t j = 0; j < temporaryIndependents_.size(); j++) {
                tx[nIndexed + nNonIndexed + j] = CGB(handler, args[temporaryIndependents_[j].atomic]);
            }

            vector<CGB> ty = fun_->Forward(p, tx);

            size_t ty_size = ty.size();
            graphForward_->indexedResults.resize(ty_size);
            std::vector<Argument<Base> > endArgs(ty_size);
            indexedArgs.resize(1);
            for (size_t i = 0; i < ty_size; i++) {
                indexedArgs[0] = Argument<Base>(*ty[i].getOperationNode());
                info[0] = i;
                OperationNode<Base>* yIndexed = new OperationNode<Base>(CGLoopIndexedDepOp, info, indexedArgs);
                graphForward_->indexedResults[i] = yIndexed;
                handler.manageOperationNode(yIndexed);
                endArgs[i] = Argument<Base>(*yIndexed);
            }

            graphForward_->loopEnd = new OperationNode<Base>(CGLoopEndOp, startEndInfo, endArgs);
            handler.manageOperationNode(graphForward_->loopEnd);

            return graphForward_;
        }

        virtual bool forward(size_t q,
                             size_t p,
                             const vector<bool>& vx,
                             vector<bool>& vy,
                             const vector<CGB>& tx,
                             vector<CGB>& ty) {

            bool valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(tx);

            if (vx.size() > 0) {
                size_t mFull = this->getLoopDependentCount();
                CPPADCG_ASSERT_KNOWN(vx.size() >= nFull_, "Invalid vx size");
                CPPADCG_ASSERT_KNOWN(vy.size() >= mFull, "Invalid vy size");
                determineFullJacobianSparsity(true);

                const vector<std::set<size_t> >& jacSparsity = custom_jac_.getFullElements();
                for (size_t i = 0; i < mFull; i++) {
                    std::set<size_t>::const_iterator it;
                    for (it = jacSparsity[i].begin(); it != jacSparsity[i].end(); ++it) {
                        size_t j = *it;
                        if (vx[j]) {
                            vy[i] = true;
                            break;
                        }
                    }
                }
            }

            bool allParameters = BaseAbstractAtomicFun<Base>::isParameters(tx);
            if (allParameters) {
                vector<Base> tyb;
                if (!evalForwardValues(q, p, tx, tyb, ty.size()))
                    return false;

                assert(tyb.size() == ty.size());
                for (size_t i = 0; i < ty.size(); i++) {
                    ty[i] = tyb[i];
                }
                return true;
            }

            vector<bool> vyLocal;
            if (p == 0) {
                vyLocal = vy;
            } else if (p >= 1) {
                /**
                 * Use the jacobian sparsity to determine which elements
                 * will always be zero
                 */
                size_t m = ty.size() / (p + 1);
                size_t n = tx.size() / (p + 1);

                vector< std::set<size_t> > r(n);
                for (size_t j = 0; j < n; j++) {
                    if (!tx[j * (p + 1) + 1].isParameter() || !tx[j * (p + 1) + 1].IdenticalZero())
                        r[j].insert(0);
                }
                vector< std::set<size_t> > s(m);
                this->for_sparse_jac(1, r, s);

                vyLocal.resize(ty.size());
                for (size_t i = 0; i < vyLocal.size(); i++) {
                    vyLocal[i] = true;
                }

                for (size_t i = 0; i < m; i++) {
                    vyLocal[i * (p + 1) + 1] = s[i].size() > 0;
                }

                if (p == 1) {
                    bool allZero = true;
                    for (size_t i = 0; i < vyLocal.size(); i++) {
                        if (vyLocal[i]) {
                            allZero = false;
                            break;
                        }
                    }

                    if (allZero) {
                        for (size_t i = 0; i < ty.size(); i++) {
                            ty[i] = Base(0.0);
                        }
                        return true;
                    }
                }
            }

            vector<Base> tyb;
            if (valuesDefined) {
                if (!evalForwardValues(q, p, tx, tyb, ty.size()))
                    return false;
            }

            CodeHandler<Base>* handler = BaseAbstractAtomicFun<Base>::findHandler(tx);
            assert(handler != NULL);

            std::vector<Argument<Base> > args = BaseAbstractAtomicFun<Base>::asArguments(tx);

            std::vector<size_t> opInfo(4);
            opInfo[0] = loopId_;
            opInfo[1] = nFull_;
            opInfo[2] = getLoopDependentCount();
            opInfo[3] = p;

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGLoopForwardOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerLoop(*this);

            opInfo.resize(1);
            args.resize(1);
            args[0] = Argument<Base>(*atomicOp);

            for (size_t i = 0; i < ty.size(); i++) {
                if (vyLocal.size() == 0 || vyLocal[i]) {
                    opInfo[0] = i;
                    ty[i] = CGB(*handler, new OperationNode<Base>(CGLoopResultOp, opInfo, args));
                    if (valuesDefined) {
                        ty[i].setValue(tyb[i]);
                    }
                } else {
                    CPPADCG_ASSERT_KNOWN(tyb.size() == 0 || IdenticalZero(tyb[i]), "Invalid value");
                    ty[i] = 0; // not a variable (zero)
                }
            }

            return true;
        }

        virtual bool reverse(size_t p,
                             const vector<CGB>& tx,
                             const vector<CGB>& ty,
                             vector<CGB>& px,
                             const vector<CGB>& py) {

            bool allParameters = BaseAbstractAtomicFun<Base>::isParameters(tx);
            if (allParameters) {
                allParameters = BaseAbstractAtomicFun<Base>::isParameters(ty);
                if (allParameters) {
                    allParameters = BaseAbstractAtomicFun<Base>::isParameters(py);
                }
            }

            if (allParameters) {
                vector<Base> pxb;

                if (!evalReverseValues(p, tx, ty, pxb, py))
                    return false;

                assert(pxb.size() == px.size());

                for (size_t i = 0; i < px.size(); i++) {
                    px[i] = pxb[i];
                }
                return true;
            }

            /**
             * Use the jacobian sparsity to determine which elements
             * will always be zero
             */
            vector<bool> vxLocal(px.size());
            for (size_t j = 0; j < vxLocal.size(); j++) {
                vxLocal[j] = true;
            }

            // k == 0
            size_t m = ty.size() / (p + 1);
            size_t n = tx.size() / (p + 1);

            vector< std::set<size_t> > rt(m);
            for (size_t i = 0; i < m; i++) {
                if (!py[i * (p + 1)].isParameter() || !py[i * (p + 1)].IdenticalZero()) {
                    rt[i].insert(0);
                }
            }
            vector< std::set<size_t> > st(n);
            this->rev_sparse_jac(1, rt, st);

            for (size_t j = 0; j < n; j++) {
                vxLocal[j * (p + 1) + p] = st[j].size() > 0;
            }

            if (p >= 1) {
                /**
                 * Use the hessian sparsity to determine which elements
                 * will always be zero
                 */
                vector<bool> vx(n);
                vector<bool> s(m);
                vector<bool> t(n);
                vector< std::set<size_t> > r(n);
                vector< std::set<size_t> > u(m);
                vector< std::set<size_t> > v(n);

                for (size_t j = 0; j < n; j++) {
                    vx[j] = !tx[j * (p + 1)].isParameter();
                    if (!tx[j * (p + 1) + 1].isParameter() || !tx[j * (p + 1) + 1].IdenticalZero()) {
                        r[j].insert(0);
                    }
                }
                for (size_t i = 0; i < m; i++) {
                    s[i] = !py[i * (p + 1) + 1].isParameter() || !py[i * (p + 1) + 1].IdenticalZero();
                }

                this->rev_sparse_hes(vx, s, t, 1, r, u, v);

                for (size_t j = 0; j < n; j++) {
                    vxLocal[j * (p + 1) + p - 1] = v[j].size() > 0;
                }
            }

            bool allZero = true;
            for (size_t j = 0; j < vxLocal.size(); j++) {
                if (vxLocal[j]) {
                    allZero = false;
                    break;
                }
            }

            if (allZero) {
                for (size_t j = 0; j < px.size(); j++) {
                    px[j] = Base(0.0);
                }
                return true;
            }

            bool valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(tx);
            if (valuesDefined) {
                valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(ty);
                if (valuesDefined) {
                    valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(py);
                }
            }

            vector<Base> pxb;
            if (valuesDefined) {
                if (!evalReverseValues(p, tx, ty, pxb, py))
                    return false;
            }

            CodeHandler<Base>* handler = BaseAbstractAtomicFun<Base>::findHandler(tx);
            if (handler == NULL) {
                handler = BaseAbstractAtomicFun<Base>::findHandler(ty);
                if (handler == NULL) {
                    handler = BaseAbstractAtomicFun<Base>::findHandler(py);
                }
            }
            assert(handler != NULL);

            std::vector<Argument<Base> > args(tx.size() + py.size());
            BaseAbstractAtomicFun<Base>::appendAsArguments(args.begin(), tx);
            BaseAbstractAtomicFun<Base>::appendAsArguments(args.begin() + tx.size(), py);

            std::vector<size_t> opInfo(4);
            opInfo[0] = loopId_;
            opInfo[1] = nFull_;
            opInfo[2] = getLoopDependentCount();
            opInfo[3] = p;

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGLoopReverseOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerLoop(*this);

            opInfo.resize(1);
            args.resize(1);
            args[0] = Argument<Base>(*atomicOp);
            for (size_t j = 0; j < px.size(); j++) {
                if (vxLocal[j]) {
                    opInfo[0] = j;
                    px[j] = CGB(*handler, new OperationNode<Base>(CGLoopResultOp, opInfo, args));
                    if (valuesDefined) {
                        px[j].setValue(pxb[j]);
                    }
                } else {
                    // CPPADCG_ASSERT_KNOWN(pxb.size() == 0 || IdenticalZero(pxb[j]), "Invalid value");
                    // pxb[j] might be non-zero but it is not required (it might have been used to determine other pxbs)
                    px[j] = Base(0); // not a variable (zero)
                }
            }

            return true;
        }

        virtual bool for_sparse_jac(size_t q,
                                    const vector<std::set<size_t> >& r,
                                    vector<std::set<size_t> >& s) {
            size_t m = iterationCount_ * m_;
            if (!custom_jac_.isFullDefined()) {
                determineFullJacobianSparsity(true);
            }

            for (size_t i = 0; i < s.size(); i++) {
                s[i].clear();
            }
            CppAD::multMatrixMatrixSparsity(custom_jac_.getFullElements(), r, s, m, nFull_, q);

            return true;
        }

        virtual bool rev_sparse_jac(size_t q,
                                    const vector<std::set<size_t> >& rt,
                                    vector<std::set<size_t> >& st) {
            size_t m = iterationCount_ * m_;
            if (!custom_jac_.isFullDefined()) {
                determineFullJacobianSparsity(false);
            }

            for (size_t i = 0; i < st.size(); i++) {
                st[i].clear();
            }
            CppAD::multMatrixMatrixSparsityTrans(rt, custom_jac_.getFullElements(), st, m, nFull_, q);

            return true;
        }

        virtual bool rev_sparse_hes(const vector<bool>& vx,
                                    const vector<bool>& s,
                                    vector<bool>& t,
                                    size_t q,
                                    const vector<std::set<size_t> >& r,
                                    const vector<std::set<size_t> >& u,
                                    vector<std::set<size_t> >& v) {
            size_t m = iterationCount_ * m_;

            for (size_t i = 0; i < nFull_; i++) {
                v[i].clear();
            }

            if (!custom_jac_.isFullDefined()) {
                determineFullJacobianSparsity(false); //improve this by determining the best choice: forward/reverse
            }
            const vector<std::set<size_t> >& jacSparsity = custom_jac_.getFullElements();

            /**
             *  V(x)  =  f'^T(x) U(x)  +  Sum(  s(x)i  f''(x)  R(x)   )
             */
            // f'^T(x) U(x)
            CppAD::multMatrixTransMatrixSparsity(jacSparsity, u, v, m, nFull_, q);

            // Sum(  s(x)i  f''(x)  R(x)   )
            bool allSelected = true;
            for (size_t i = 0; i < m; i++) {
                if (!s[i]) {
                    allSelected = false;
                    break;
                }
            }

            if (allSelected) {
                if (!custom_hess_.isFullDefined()) {
                    determineFullHessianSparsity();
                }
                const vector<std::set<size_t> >& sF2 = custom_hess_.getFullElements();
                CppAD::multMatrixTransMatrixSparsity(sF2, r, v, nFull_, nFull_, q); // f''^T * R
            } else {
                vector<std::set<size_t> > sparsitySF2R(nFull_);
                for (size_t i = 0; i < m; i++) {
                    if (s[i]) {
                        std::map<size_t, vector<std::set<size_t> > >::const_iterator itH = hess_.find(i);
                        const vector<std::set<size_t> >* spari;
                        if (itH == hess_.end()) {
                            vector<std::set<size_t> >& hi = determineFullHessianSparsity(i);
                            spari = &hi;
                            custom_hess_.filter(hi);
                        } else {
                            spari = &itH->second;
                        }
                        CppAD::addMatrixSparsity(*spari, sparsitySF2R);
                    }
                }
                CppAD::multMatrixTransMatrixSparsity(sparsitySF2R, r, v, nFull_, nFull_, q); // f''^T * R
            }

            /**
             * S(x) * f'(x)
             */
            std::set<size_t>::const_iterator it;
            for (size_t i = 0; i < m; i++) {
                if (s[i]) {
                    for (it = jacSparsity[i].begin(); it != jacSparsity[i].end(); ++it) {
                        size_t j = *it;
                        t[j] = true;
                    }
                }
            }

            return true;
        }

        virtual ~LoopAtomicFun() {
            delete fun_;
            for (size_t i = 0; i < indepIndexPatterns_.size(); i++) {
                delete indepIndexPatterns_[i];
            }
            for (size_t i = 0; i < depIndexPatterns_.size(); i++) {
                delete depIndexPatterns_[i];
            }
            delete graphForward_;
        }

    protected:

        void determineFullJacobianSparsity(bool forward) {
            if (forward) {
                jacTape_ = CppAD::jacobianForwardSparsitySet<vector<std::set<size_t> > >(*fun_);
                fun_->size_forward_set(0);
            } else {
                jacTape_ = CppAD::jacobianReverseSparsitySet<vector<std::set<size_t> > >(*fun_);
            }
            size_t mFull = iterationCount_ * m_;
            vector<std::set<size_t> > fullJac(mFull);

            for (size_t it = 0; it < iterationCount_; it++) {
                for (size_t d = 0; d < m_; d++) {
                    size_t i = it * m_ + d;
                    const std::set<size_t>& firstitd = jacTape_[d];

                    std::set<size_t>::const_iterator itj1;
                    for (itj1 = firstitd.begin(); itj1 != firstitd.end(); ++itj1) {
                        size_t j1 = *itj1; // tape independent index
                        size_t jFull = getAtomicLoopIndex(it, j1);
                        fullJac[i].insert(jFull);
                    }
                }
            }

            custom_jac_.setFullElements(fullJac);
        }

        inline void determineFullHessianSparsity() {
            size_t mTape = fun_->Range();
            size_t nTape = fun_->Domain();
            assert(m_ == mTape);

            /**
             * Make sure that an independent of the original model is not used
             * by more than one indexed tape independent and/or a nonindexed
             * independent.
             */
            for (size_t d = 0; d < mTape; d++) {
                bool detachedVars = isIndependentVariablesFullyDetached(d);
                if (!detachedVars) {
                    throw CGException("Unable exploit the hessian structure of a loop where the independent variables are not fully detached");
                }
            }


            hessTape_ = CppAD::hessianSparsitySet<vector<std::set<size_t> > >(*fun_); // f''(x)

            vector<std::set<size_t> > fullHess(nFull_);

            for (size_t j1 = 0; j1 < nTape; j1++) {
                const std::set<size_t>& row = hessTape_[j1];
                if (row.empty())
                    continue;

                for (size_t it = 0; it < iterationCount_; it++) {
                    size_t jFull1 = getAtomicLoopIndex(it, j1);
                    std::set<size_t>& rowFull = fullHess[jFull1];

                    std::set<size_t>::const_iterator itj2;
                    for (itj2 = row.begin(); itj2 != row.end(); ++itj2) {
                        size_t j2 = *itj2; // tape independent index
                        size_t jFull2 = getAtomicLoopIndex(it, j2);
                        rowFull.insert(jFull2);
                    }
                }
            }

            custom_hess_.setFullElements(fullHess);
        }

        inline vector<std::set<size_t> >& determineFullHessianSparsity(size_t i) {
            size_t mTape = fun_->Range();
            size_t nTape = fun_->Domain();
            assert(m_ == mTape);
            size_t iteration = i / iterationCount_;
            size_t d = iteration * mTape + (i % iterationCount_); // tape equation/dependent

            /**
             * Make sure that an independent of the original model is not used
             * by more than one indexed tape independent and/or a nonindexed
             * independent.
             */
            bool detachedVars = isIndependentVariablesFullyDetached(d);
            if (!detachedVars) {
                throw CGException("Unable to exploit the hessian structure of a loop where the independent variables are not fully detached");
            }

            const vector<std::set<size_t> >& hessd = eqHessTape_[d] = CppAD::hessianSparsitySet<vector<std::set<size_t> > >(*fun_, d); // f''_i(x)

            vector<std::set<size_t> >& fullHessi = hess_[i];
            fullHessi.resize(nFull_);

            for (size_t j1 = 0; j1 < nTape; j1++) {
                const std::set<size_t>& row = hessd[j1];
                if (row.empty())
                    continue;

                size_t jFull1 = getAtomicLoopIndex(iteration, j1);
                std::set<size_t>& rowFull = fullHessi[jFull1];

                std::set<size_t>::const_iterator itj2;
                for (itj2 = row.begin(); itj2 != row.end(); ++itj2) {
                    size_t j2 = *itj2; // tape independent index
                    size_t jFull2 = getAtomicLoopIndex(iteration, j2);
                    rowFull.insert(jFull2);
                }
            }

            return fullHessi;
        }

        /**
         * Checks if an independent of the original model is NOT used by more
         * than one indexed tape independent and/or a nonindexed independent.
         */
        inline bool isIndependentVariablesFullyDetached(size_t d) const {
            std::set<size_t> origIndexedIndependents;
            const std::set<size_t>& eqJac = jacTape_[d];

            /**
             ****************************************************************** 
             ******************************************************************
             *          TODO : temporaries also depend on independents
             * 
             ******************************************************************* 
             ******************************************************************* 
             */
            size_t nIndexed = indexedIndepIndexes_.size();
            size_t nNonIndexed = nonIndexedIndepIndexes_.size();

            std::set<size_t>::const_iterator itj;
            for (itj = eqJac.begin(); itj != eqJac.end(); ++itj) {
                size_t j1 = *itj; // tape independent index
                if (j1 < nIndexed) {
                    // indexed 
                    for (size_t it = 0; it < iterationCount_; it++) {
                        size_t origj = indexedIndepIndexes_[j1][it].original;
                        if (origIndexedIndependents.find(origj) != origIndexedIndependents.end()) {
                            return false;
                        }
                        origIndexedIndependents.insert(origj);
                    }
                } else {

                    if (j1 < nIndexed + nNonIndexed) {
                        // the index does not change
                        size_t origj = nonIndexedIndepIndexes_[j1 - nIndexed].original;
                        if (origIndexedIndependents.find(origj) != origIndexedIndependents.end()) {
                            return false;
                        }
                    } else {
                        // temporary variables (the index does not change)
                        const std::set<size_t>& origs = temporaryIndependents_[j1 - (nIndexed + nNonIndexed)].originalIndeps;
                        std::set<size_t>::const_iterator itoj;
                        for (itoj = origs.begin(); itoj != origs.end(); ++itoj) {
                            size_t origj = *itoj;
                            if (origIndexedIndependents.find(origj) != origIndexedIndependents.end()) {
                                return false;
                            }
                        }
                    }
                }
            }

            return true;
        }

        /**
         * Determines the index of an independent variable of the loop tape in 
         * the atomic loop function
         * 
         * @param iteration The iteration number
         * @param j The independent variable index in the loop tape
         * @return the index in the atomic loop function
         */
        inline size_t getAtomicLoopIndex(size_t iteration, size_t j) const {
            size_t nIndexed = indexedIndepIndexes_.size();
            size_t nNonIndexed = nonIndexedIndepIndexes_.size();

            if (j < nIndexed) {
                return indexedIndepIndexes_[j][iteration].atomic;
            } else if (j < nIndexed + nNonIndexed) {
                return nonIndexedIndepIndexes_[j - nIndexed].atomic;
            } else {
                return temporaryIndependents_[j - nIndexed - nNonIndexed].atomic;
            }
        }

    private:

        static size_t createNewLoopId() {
            CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;
            static size_t count = 0;
            count++;
            return count;
        }

        inline bool evalForwardValues(size_t q,
                                      size_t p,
                                      const vector<CGB>& tx,
                                      vector<Base>& tyb,
                                      size_t ty_size) {
            vector<Base> txb(tx.size());
            tyb.resize(ty_size);

            for (size_t i = 0; i < txb.size(); i++) {
                txb[i] = tx[i].getValue();
            }

            return atomicForward(q, p, txb, tyb);
        }

        /**
         * Used to evaluate function values and forward mode function values and
         * derivatives.
         * 
         * @param q Lowerest order for this forward mode calculation.
         * @param p Highest order for this forward mode calculation.
         * @param tx Taylor coefficients corresponding to \c x for this
         *           calculation
         * @param ty Taylor coefficient corresponding to \c y for this 
         *           calculation
         * @return true on success, false otherwise
         */
        bool atomicForward(size_t q,
                           size_t p,
                           const vector<Base>& tx,
                           vector<Base>& ty) {

            if (p == 0) {
                //size_t mTape = fun_->Range();
                size_t nTape = fun_->Domain();

                vector<CGB> xTape(nTape);

                for (size_t j = 0; j < nonIndexedIndepIndexes_.size(); j++) {
                    const LoopPosition& pos = nonIndexedIndepIndexes_[j];
                    xTape[pos.tape] = tx[pos.atomic];
                }

                for (size_t j = 0; j < temporaryIndependents_.size(); j++) {
                    const LoopPositionTmp& pos = temporaryIndependents_[j];
                    xTape[pos.tape] = tx[pos.atomic];
                }

                /**
                 * loop...
                 */
                for (size_t it = 0; it < iterationCount_; it++) {
                    // place indexed independents
                    for (size_t j = 0; j < indexedIndepIndexes_.size(); j++) {
                        const LoopPosition& pos = indexedIndepIndexes_[j][it];
                        xTape[pos.tape] = tx[pos.atomic];
                    }

                    vector<CGB> yTape = fun_->Forward(0, xTape);

                    // place dependents
                    for (size_t i = 0; i < dependentIndexes_.size(); i++) {
                        const LoopPosition& pos = dependentIndexes_[i][it];
                        ty[pos.atomic] = yTape[pos.tape].getValue();
                    }
                }

                return true;
            }

            return false;
        }

        inline bool evalReverseValues(size_t p,
                                      const vector<CGB>& tx,
                                      const vector<CGB>& ty,
                                      vector<Base>& pxb,
                                      const vector<CGB>& py) {
            vector<Base> txb(tx.size());
            vector<Base> tyb(ty.size());
            pxb.resize(tx.size());
            vector<Base> pyb(py.size());

            for (size_t i = 0; i < txb.size(); i++) {
                txb[i] = tx[i].getValue();
            }
            for (size_t i = 0; i < tyb.size(); i++) {
                tyb[i] = ty[i].getValue();
            }
            for (size_t i = 0; i < pyb.size(); i++) {
                pyb[i] = py[i].getValue();
            }

            return atomicReverse(p, txb, tyb, pxb, pyb);
        }

        /**
         * Used to evaluate reverse mode function derivatives.
         * 
         * @param p Highest order for this forward mode calculation.
         * @param tx Taylor coefficients corresponding to \c x for this
         *           calculation
         * @param ty Taylor coefficient corresponding to \c y for this 
         *           calculation
         * @param px Partials w.r.t. the \c x Taylor coefficients.
         * @param py Partials w.r.t. the \c y Taylor coefficients
         * @return true on success, false otherwise
         */
        bool atomicReverse(size_t p,
                           const vector<Base>& tx,
                           const vector<Base>& ty,
                           vector<Base>& px,
                           const vector<Base>& py) {
            return false;
        }

    };

}

#endif