#ifndef CPPAD_CG_ABSTRACT_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_ABSTRACT_ATOMIC_FUN_INCLUDED
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
     * An atomic function for source code generation
     * 
     * @author Joao Leal
     */
    template <class Base>
    class CGAbstractAtomicFun : public BaseAbstractAtomicFun<Base> {
    public:
        typedef CppAD::CG<Base> CGB;
        typedef Argument<Base> Arg;
    protected:
        const size_t id_;
        bool standAlone_;

    protected:

        /**
         * Creates a new atomic function that is responsible for defining the
         * dependencies to calls of a user atomic function.
         * 
         * @param name The atomic function name.
         * @param standAlone Whether or not forward and reverse function calls
         *                   do not require the Taylor coefficients for the 
         *                   dependent variables (ty) and any previous
         *                   evaluation of other forward/reverse modes. 
         */
        CGAbstractAtomicFun(const std::string& name, bool standAlone = false) :
            BaseAbstractAtomicFun<Base>(name),
            id_(createNewId()),
            standAlone_(standAlone) {
            CPPADCG_ASSERT_KNOWN(!name.empty(), "The atomic function name cannot be empty");
        }

    public:

        template <class ADVector>
        void operator()(const ADVector& ax, ADVector& ay, size_t id = 0) {
            this->BaseAbstractAtomicFun<Base>::operator()(ax, ay, id);
        }

        /**
         * Provides a unique identifier for this atomic function type.
         * 
         * @return a unique identifier ID
         */
        size_t getId() const {
            return id_;
        }

        virtual bool forward(size_t q,
                             size_t p,
                             const vector<bool>& vx,
                             vector<bool>& vy,
                             const vector<CGB>& tx,
                             vector<CGB>& ty) {

            bool valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(tx);
            CPPADCG_ASSERT_KNOWN(valuesDefined || vx.size() == 0,
                                 "Values must be defined in order to call atomic function");
            if (!valuesDefined && vx.size() > 0)
                return false;

            bool allParameters = BaseAbstractAtomicFun<Base>::isParameters(tx);
            if (allParameters) {
                vector<Base> tyb;
                if (!evalForwardValues(q, p, vx, vy, tx, tyb, ty.size()))
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
                if (!evalForwardValues(q, p, vx, vy, tx, tyb, ty.size()))
                    return false;
            }

            CodeHandler<Base>* handler = findHandler(tx);
            assert(handler != NULL);

            OperationNode<Base>* txArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, tx);
            OperationNode<Base>* tyArray;

            if (standAlone_ && p > 0) {
                tyArray = BaseAbstractAtomicFun<Base>::makeZeroArray(*handler, ty);
            } else {
                tyArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, ty);
            }

            std::vector<size_t> opInfo(3);
            opInfo[0] = id_;
            opInfo[1] = q;
            opInfo[2] = p;
            std::vector<Argument<Base> > args(2);
            args[0] = Argument<Base>(*txArray);
            args[1] = Argument<Base>(*tyArray);

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGAtomicForwardOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerAtomicFunction(*this);

            opInfo.resize(1);
            args.resize(2);
            for (size_t i = 0; i < ty.size(); i++) {
                if (vyLocal.size() == 0 || vyLocal[i]) {
                    opInfo[0] = i;
                    args[0] = Argument<Base>(*tyArray);
                    args[1] = Argument<Base>(*atomicOp);

                    ty[i] = CGB(*handler, new OperationNode<Base>(CGArrayElementOp, opInfo, args));
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

            OperationNode<Base>* txArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, tx);
            OperationNode<Base>* tyArray;
            OperationNode<Base>* pxArray = BaseAbstractAtomicFun<Base>::makeZeroArray(*handler, px);
            OperationNode<Base>* pyArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, py);

            if (standAlone_) {
                tyArray = BaseAbstractAtomicFun<Base>::makeZeroArray(*handler, ty);
            } else {
                tyArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, ty);
            }

            std::vector<size_t> opInfo(2);
            opInfo[0] = id_;
            opInfo[1] = p;
            std::vector<Argument<Base> > args(4);
            args[0] = Argument<Base>(*txArray);
            args[1] = Argument<Base>(*tyArray);
            args[2] = Argument<Base>(*pxArray);
            args[3] = Argument<Base>(*pyArray);

            OperationNode<Base>* atomicOp = new OperationNode<Base>(CGAtomicReverseOp, opInfo, args);
            handler->manageOperationNode(atomicOp);
            handler->registerAtomicFunction(*this);

            opInfo.resize(1);
            args.resize(2);
            for (size_t j = 0; j < px.size(); j++) {
                if (vxLocal[j]) {
                    opInfo[0] = j;
                    args[0] = Argument<Base>(*pxArray);
                    args[1] = Argument<Base>(*atomicOp);

                    px[j] = CGB(*handler, new OperationNode<Base>(CGArrayElementOp, opInfo, args));
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

        virtual ~CGAbstractAtomicFun() {
        }

    protected:

        /**
         * Used to evaluate function values and forward mode function values and
         * derivatives.
         * 
         * @param q Lowerest order for this forward mode calculation.
         * @param p Highest order for this forward mode calculation.
         * @param vx If size not zero, which components of \c x are variables
         * @param vy If size not zero, which components of \c y are variables
         * @param tx Taylor coefficients corresponding to \c x for this
         *           calculation
         * @param ty Taylor coefficient corresponding to \c y for this 
         *           calculation
         * @return true on success, false otherwise
         */
        virtual bool atomicForward(size_t q,
                                   size_t p,
                                   const vector<bool>& vx,
                                   vector<bool>& vy,
                                   const vector<Base>& tx,
                                   vector<Base>& ty) = 0;
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
        virtual bool atomicReverse(size_t p,
                                   const vector<Base>& tx,
                                   const vector<Base>& ty,
                                   vector<Base>& px,
                                   const vector<Base>& py) = 0;

    private:

        inline bool evalForwardValues(size_t q,
                                      size_t p,
                                      const vector<bool>& vx,
                                      vector<bool>& vy,
                                      const vector<CGB>& tx,
                                      vector<Base>& tyb,
                                      size_t ty_size) {
            vector<Base> txb(tx.size());
            tyb.resize(ty_size);

            for (size_t i = 0; i < txb.size(); i++) {
                txb[i] = tx[i].getValue();
            }

            return atomicForward(q, p, vx, vy, txb, tyb);
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

        static size_t createNewId() {
            CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;
            static size_t count = 0;
            count++;
            return count;
        }

    };

}

#endif