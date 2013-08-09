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
         * the total number of dependent variables of the atomic function
         * (NOT the tape)
         */
        const size_t mFull_;
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
         * Maps the original dependent variable indexes to their positions in 
         * the loop
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
        /**
         * Jacobian sparsity pattern of the tape
         */
        vector<std::set<size_t> > jacTapeSparsity_;
        CustomPosition custom_jac_;
        vector<std::set<size_t> > hessTape_;
        std::map<size_t, vector<std::set<size_t> > > eqHessTape_;
        std::map<size_t, vector<std::set<size_t> > > hess_;
        CustomPosition custom_hess_;
        /**
         * 
         */
        std::map<size_t, bool> equationFullyDetached_;
        /**
         * Maps original model independent index to the possitions it may be
         * associated with in the tape (also includes temporary variables)
         * for a given tape equation (equation -> orig index -> tape indexes)
         */
        std::vector<std::map<size_t, std::set<size_t> > > orig2Tape_;
    public:

        /**
         * Creates a new atomic function that is responsible for defining the
         * dependencies to calls of a user atomic function.
         * 
         * @param name The atomic function name.
         * @param fun
         * @param iterationCount
         * @param dependentOrigIndexes
         * @param indexedIndepOrigIndexes
         * @param nonIndexedIndepOrigIndexes
         * @param temporaryIndependents
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
            mFull_(iterationCount * dependentOrigIndexes.size()),
            nFull_(0),
            m_(dependentOrigIndexes.size()),
            dependentIndexes_(m_, std::vector<LoopPosition>(iterationCount)),
            indexedIndepIndexes_(indexedIndepOrigIndexes.size(), std::vector<LoopPosition>(iterationCount)),
            nonIndexedIndepIndexes_(nonIndexedIndepOrigIndexes.size()),
            temporaryIndependents_(temporaryIndependents.size()) {
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
                assert(!temporaryIndependents[j].empty());
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

        inline ADFun<CGB>* getTape() const {
            return fun_;
        }

        inline size_t getLoopDependentCount() const {
            return mFull_;
        }

        inline size_t getLoopIndependentCount() const {
            return nFull_;
        }

        inline size_t getTapeDependentCount() const {
            return m_;
        }

        inline size_t getTapeIndependentCount() const {
            return fun_->Domain();
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

        inline const LoopPosition& getTapeDependentIndex(size_t originalDepIndex) const {
            return *depOrigIndexes_.at(originalDepIndex);
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
                x[indepOrig2full_.size() + j] = handler->createCG(Argument<Base>(*tmp));
            }

            std::vector<size_t> opInfo(6);
            opInfo[0] = loopId_;
            opInfo[1] = nFull_;
            opInfo[2] = mFull_;
            opInfo[3] = 0; // zero-order
            opInfo[4] = mFull_; // dimension of results of this atomic function call
            opInfo[5] = m_; // dimension of results relative to the tape
            std::vector<Argument<Base> > args = asArguments(x);

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGLoopForwardOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerLoop(*this); // is this really required?

            vector<CGB> y(mFull_);

            opInfo.resize(1);
            args.resize(1);
            for (size_t i = 0; i < mFull_; i++) {
                opInfo[0] = i;
                args[0] = Argument<Base>(*atomicOp);

                y[i] = handler->createCG(new OperationNode<Base>(CGLoopAtomicResultOp, opInfo, args));
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

        inline bool isTemporary(size_t tapeIndex) const {
            size_t nIndexed = indexedIndepIndexes_.size();
            size_t nNonIndexed = nonIndexedIndepIndexes_.size();

            return nIndexed + nNonIndexed <= tapeIndex;
        }

        /**
         * Checks if all variables of an equation pattern are completelly
         * unrelated, that is, if all variables in the original are only
         * associated with indexed tape independents or non-indexed tape 
         * independents.
         * 
         * @param the index of the equation pattern in the tape.
         * @return true if all variables for the equation pattern are 
         *         completelly unrelated
         */
        inline bool isIndependentVariablesFullyDetached(size_t d) {
            assert(d < m_);

            std::map<size_t, bool>::const_iterator iteqd = equationFullyDetached_.find(d);
            if (iteqd != equationFullyDetached_.end()) {
                return iteqd->second;
            }

            equationFullyDetached_[d] = false;

            // iteration -> orignal index -> indexed
            std::vector < std::map<size_t, bool> > origIndexedIndependents(iterationCount_);
            const std::set<size_t>& eqJac = jacTapeSparsity_[d];

            size_t nIndexed = indexedIndepIndexes_.size();
            size_t nNonIndexed = nonIndexedIndepIndexes_.size();

            std::map<size_t, bool>::const_iterator itIndexed;
            std::set<size_t>::const_iterator itj;
            for (itj = eqJac.begin(); itj != eqJac.end(); ++itj) {
                size_t j1 = *itj; // tape independent index
                if (j1 < nIndexed) {
                    // indexed independents
                    for (size_t it = 0; it < iterationCount_; it++) {
                        size_t origj = indexedIndepIndexes_[j1][it].original;
                        itIndexed = origIndexedIndependents[it].find(origj);
                        if (itIndexed != origIndexedIndependents[it].end()) {
                            return false;
                        }
                        origIndexedIndependents[it][origj] = true; // indexed
                    }

                } else if (j1 < nIndexed + nNonIndexed) {
                    // the index does not change
                    size_t origj = nonIndexedIndepIndexes_[j1 - nIndexed].original;
                    for (size_t it = 0; it < iterationCount_; it++) {
                        itIndexed = origIndexedIndependents[it].find(origj);
                        if (itIndexed != origIndexedIndependents[it].end() && itIndexed->second) {
                            return false;
                        }
                        origIndexedIndependents[it][origj] = false; // not indexed
                    }
                } else {
                    // temporary variables (the index does not change)
                    const std::set<size_t>& origs = temporaryIndependents_[j1 - (nIndexed + nNonIndexed)].originalIndeps;
                    std::set<size_t>::const_iterator itoj;
                    for (itoj = origs.begin(); itoj != origs.end(); ++itoj) {
                        size_t origj = *itoj;
                        for (size_t it = 0; it < iterationCount_; it++) {
                            itIndexed = origIndexedIndependents[it].find(origj);
                            if (itIndexed != origIndexedIndependents[it].end() && itIndexed->second) {
                                return false;
                            }
                            origIndexedIndependents[it][origj] = false; // not indexed
                        }
                    }
                }

            }

            equationFullyDetached_[d] = true;
            return true;
        }

        std::set<size_t> getIndependentTapeIndexes(size_t tapeI, size_t origIndep) {
            assert(tapeI < m_);

            if (orig2Tape_.empty()) {
                orig2Tape_.resize(m_);

                size_t nIndexed = indexedIndepIndexes_.size();
                size_t nNonIndexed = nonIndexedIndepIndexes_.size();

                if (jacTapeSparsity_.size() == 0) {
                    determineFullJacobianSparsity();
                }

                for (size_t i = 0; i < m_; i++) {
                    const std::set<size_t>& row = jacTapeSparsity_[i];
                    std::set<size_t>::const_iterator itr;
                    for (itr = row.begin(); itr != row.end(); ++itr) {
                        size_t j = *itr;

                        if (j < nIndexed) {
                            // indexed
                            for (size_t it = 0; it < iterationCount_; it++) {
                                const LoopPosition& pos = indexedIndepIndexes_[j][it];
                                orig2Tape_[i][pos.original].insert(pos.tape);
                            }
                        } else if (j < nIndexed + nNonIndexed) {
                            // non-indexed
                            const LoopPosition& pos = nonIndexedIndepIndexes_[j - nIndexed];
                            orig2Tape_[i][pos.original].insert(pos.tape);
                        } else {
                            // temporaries
                            const LoopPositionTmp& pos = temporaryIndependents_[j - (nIndexed + nNonIndexed)];
                            const std::set<size_t>& origs = pos.originalIndeps;
                            std::set<size_t>::const_iterator ito;
                            for (ito = origs.begin(); ito != origs.end(); ++ito) {
                                orig2Tape_[i][*ito].insert(pos.tape);
                            }
                        }
                    }
                }

            }

            return orig2Tape_[tapeI].at(origIndep);
        }

        virtual bool forward(size_t q,
                             size_t p,
                             const vector<bool>& vx,
                             vector<bool>& vy,
                             const vector<CGB>& tx,
                             vector<CGB>& ty) {

            bool valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(tx);

            if (vx.size() > 0) {
                CPPADCG_ASSERT_KNOWN(vx.size() >= nFull_, "Invalid vx size");
                CPPADCG_ASSERT_KNOWN(vy.size() >= mFull_, "Invalid vy size");

                for (size_t i = 0; i < mFull_; i++) {
                    /**
                     * all dependents must be considered a variable
                     * (otherwise CppAD will remember constants and will not 
                     *  make the association with the loop results)
                     */
                    vy[i] = true;
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

            CodeHandler<Base>* handler = findHandler(tx);
            assert(handler != NULL);

            std::vector<Argument<Base> > args = asArguments(tx);

            std::vector<size_t> opInfo(6);
            opInfo[0] = loopId_;
            opInfo[1] = mFull_;
            opInfo[2] = nFull_;
            opInfo[3] = p;
            opInfo[4] = ty.size(); // dimension of results of this atomic function call
            opInfo[5] = m_ * (p + 1); // dimension of results relative to the tape

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGLoopForwardOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerLoop(*this);

            vector<size_t> tapeIndex(ty.size());
            for (size_t k = 0; k <= p; k++) {
                for (size_t j = 0; j < dependentIndexes_.size(); j++) {
                    for (size_t it = 0; it < iterationCount_; it++) {
                        const LoopPosition& pos = dependentIndexes_[j][it];
                        tapeIndex[pos.atomic * (p + 1) + k] = pos.tape * (p + 1) + k;
                    }
                }
            }

            opInfo.resize(2);
            args.resize(1);
            args[0] = Argument<Base>(*atomicOp);

            for (size_t i = 0; i < ty.size(); i++) {
                if (vyLocal.size() == 0 || vyLocal[i]) {
                    opInfo[0] = i; // atomic index
                    opInfo[1] = tapeIndex[i]; // tape index
                    ty[i] = handler->createCG(new OperationNode<Base>(CGLoopAtomicResultOp, opInfo, args));
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

            CodeHandler<Base>* handler = findHandler(tx);
            if (handler == NULL) {
                handler = findHandler(ty);
                if (handler == NULL) {
                    handler = findHandler(py);
                }
            }
            assert(handler != NULL);

            std::vector<Argument<Base> > args(tx.size() + py.size());
            BaseAbstractAtomicFun<Base>::appendAsArguments(args.begin(), tx);
            BaseAbstractAtomicFun<Base>::appendAsArguments(args.begin() + tx.size(), py);

            std::vector<size_t> opInfo(6);
            opInfo[0] = loopId_;
            opInfo[1] = mFull_;
            opInfo[2] = nFull_;
            opInfo[3] = p;
            opInfo[4] = px.size(); // dimension of results of this atomic function call
            opInfo[5] = getTapeIndependentCount()*(p + 1); // dimension of results of this atomic function call

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGLoopReverseOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerLoop(*this);

            vector<size_t> tapeIndex(px.size());
            for (size_t k = 0; k <= p; k++) {
                for (size_t j = 0; j < indexedIndepIndexes_.size(); j++) {
                    for (size_t it = 0; it < iterationCount_; it++) {
                        const LoopPosition& pos = indexedIndepIndexes_[j][it];
                        tapeIndex[pos.atomic * (p + 1) + k] = pos.tape * (p + 1) + k;
                    }
                }
                for (size_t j = 0; j < nonIndexedIndepIndexes_.size(); j++) {
                    const LoopPosition& pos = nonIndexedIndepIndexes_[j];
                    tapeIndex[pos.atomic * (p + 1) + k] = pos.tape * (p + 1) + k;
                }
                for (size_t j = 0; j < temporaryIndependents_.size(); j++) {
                    const LoopPositionTmp& pos = temporaryIndependents_[j];
                    tapeIndex[pos.atomic * (p + 1) + k] = pos.tape * (p + 1) + k;
                }
            }

            opInfo.resize(2);
            args.resize(1);
            args[0] = Argument<Base>(*atomicOp);
            for (size_t j = 0; j < px.size(); j++) {
                if (vxLocal[j]) {
                    opInfo[0] = j; //atomic index
                    opInfo[1] = tapeIndex[j]; //tape index
                    px[j] = handler->createCG(new OperationNode<Base>(CGLoopAtomicResultOp, opInfo, args));
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
            CppAD::multMatrixTransMatrixSparsity(jacSparsity, u, v, mFull_, nFull_, q);

            // Sum(  s(x)i  f''(x)  R(x)   )
            bool allSelected = true;
            for (size_t i = 0; i < mFull_; i++) {
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
                for (size_t i = 0; i < mFull_; i++) {
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
            for (size_t i = 0; i < mFull_; i++) {
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
        }

    protected:

        void determineFullJacobianSparsity() {
            determineFullJacobianSparsity(m_ < fun_->Domain());
        }

        void determineFullJacobianSparsity(bool forward) {
            if (forward) {
                jacTapeSparsity_ = CppAD::jacobianForwardSparsitySet<vector<std::set<size_t> > >(*fun_);
                fun_->size_forward_set(0);
            } else {
                jacTapeSparsity_ = CppAD::jacobianReverseSparsitySet<vector<std::set<size_t> > >(*fun_);
            }
            size_t mFull = iterationCount_ * m_;
            vector<std::set<size_t> > fullJac(mFull);

            for (size_t it = 0; it < iterationCount_; it++) {
                for (size_t d = 0; d < m_; d++) {
                    size_t i = it * m_ + d;
                    const std::set<size_t>& firstitd = jacTapeSparsity_[d];

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
            size_t nTape = fun_->Domain();
            assert(m_ == fun_->Range());

            /**
             * Make sure that an independent of the original model is not used
             * by more than one indexed tape independent and/or a nonindexed
             * independent.
             */
            for (size_t d = 0; d < m_; d++) {
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
            size_t nTape = fun_->Domain();
            assert(m_ == fun_->Range());
            size_t iteration = i / iterationCount_;
            size_t d = iteration * m_ + (i % iterationCount_); // tape equation/dependent

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
            } else if (p == 1) {
                size_t nTape = fun_->Domain();

                vector<CGB> xTape(nTape);
                vector<CGB> txTape(nTape * 2);

                for (size_t j = 0; j < nonIndexedIndepIndexes_.size(); j++) {
                    const LoopPosition& pos = nonIndexedIndepIndexes_[j];
                    xTape[pos.tape] = tx[pos.atomic * 2];
                    txTape[pos.tape * 2] = tx[pos.atomic * 2];
                    txTape[pos.tape * 2 + 1] = tx[pos.atomic * 2 + 1];
                }

                for (size_t j = 0; j < temporaryIndependents_.size(); j++) {
                    const LoopPositionTmp& pos = temporaryIndependents_[j];
                    xTape[pos.tape] = tx[pos.atomic * 2];
                    txTape[pos.tape * 2] = tx[pos.atomic * 2];
                    txTape[pos.tape * 2 + 1] = tx[pos.atomic * 2 + 1];
                }

                /**
                 * loop...
                 */
                for (size_t it = 0; it < iterationCount_; it++) {
                    // place indexed independents
                    for (size_t j = 0; j < indexedIndepIndexes_.size(); j++) {
                        const LoopPosition& pos = indexedIndepIndexes_[j][it];
                        xTape[pos.tape] = tx[pos.atomic * 2];
                        txTape[pos.tape * 2] = tx[pos.atomic * 2];
                        txTape[pos.tape * 2 + 1] = tx[pos.atomic * 2 + 1];
                    }

                    fun_->Forward(0, xTape);
                    vector<CGB> tyTape = fun_->Forward(1, txTape);

                    // place dependents
                    for (size_t i = 0; i < dependentIndexes_.size(); i++) {
                        const LoopPosition& pos = dependentIndexes_[i][it];
                        ty[pos.atomic * 2] = tyTape[pos.tape * 2].getValue();
                        ty[pos.atomic * 2 + 1] = tyTape[pos.tape * 2 + 1].getValue();
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
            if (p == 0) {
                size_t nTape = fun_->Domain();

                vector<CGB> xTape(nTape);
                vector<CGB> pyTape(m_);

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

                    for (size_t i = 0; i < m_; i++) {
                        const LoopPosition& pos = dependentIndexes_[i][it];
                        pyTape[pos.tape] = py[pos.atomic];
                    }

                    fun_->Forward(0, xTape);
                    vector<CGB> pxTape = fun_->Reverse(1, pyTape);

                    // place dependents
                    for (size_t j = 0; j < indexedIndepIndexes_.size(); j++) {
                        const LoopPosition& pos = indexedIndepIndexes_[j][it];
                        px[pos.atomic] = pxTape[pos.tape].getValue(); // should it be += ?
                    }
                }

                return true;
            }

            return false;
        }

    };

}

#endif