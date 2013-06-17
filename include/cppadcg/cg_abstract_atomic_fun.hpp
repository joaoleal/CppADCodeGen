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
    class CGAbstractAtomicFun : public atomic_base<CppAD::CG<Base> > {
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
         *                   dependent variables (ty) and the previous
         *                   evaluation of other forward/reverse modes.
         */
        CGAbstractAtomicFun(const std::string& name, bool standAlone = false) :
            atomic_base<CGB>(name.c_str()),
            id_(createNewId()),
            standAlone_(standAlone) {

        }

        template <class ADVector>
        void operator()(const ADVector& ax, ADVector& ay, size_t id = 0) {
            this->atomic_base<CGB>::operator()(ax, ay, id);
        }

    public:

        virtual bool forward(size_t q,
                             size_t p,
                             const vector<bool>& vx,
                             vector<bool>& vy,
                             const vector<CGB>& tx,
                             vector<CGB>& ty) {

            bool valuesDefined = true;
            for (size_t i = 0; i < tx.size(); i++) {
                if (!tx[i].isValueDefined()) {
                    valuesDefined = false;
                    break;
                }
            }

            CPPADCG_ASSERT_KNOWN(valuesDefined || vx.size() == 0,
                                 "Values must be defined in order to call atomic function");
            if (!valuesDefined && vx.size() > 0)
                return false;

            vector<Base> txb;
            vector<Base> tyb;
            if (valuesDefined) {
                txb.resize(tx.size());
                tyb.resize(ty.size());
                for (size_t i = 0; i < txb.size(); i++) {
                    txb[i] = tx[i].getValue();
                }

                if (!atomicForward(q, p, vx, vy, txb, tyb))
                    return false;
            }

            bool allParameters = isParameters(tx);
            if (allParameters) {
                assert(tyb.size() == ty.size());

                for (size_t i = 0; i < ty.size(); i++) {
                    ty[i] = tyb[i];
                }
            } else {
                CodeHandler<Base>* handler = findHandler(tx);
                assert(handler != NULL);

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
                }


                bool allZero = false;
                if (vyLocal.size() > 0) {
                    allZero = true;
                    for (size_t i = 0; i < vyLocal.size(); i++) {
                        if (vyLocal[i]) {
                            allZero = false;
                            break;
                        }
                    }
                }

                if (allZero) {
                    for (size_t i = 0; i < ty.size(); i++) {
                        ty[i] = Base(0.0);
                    }
                } else {
                    SourceCodeFragment<Base>* txArray = makeArray(*handler, tx);
                    SourceCodeFragment<Base>* tyArray;

                    if (standAlone_ && p > 0) {
                        tyArray = makeZeroArray(*handler, ty);
                    } else {
                        tyArray = makeArray(*handler, ty);
                    }

                    std::vector<size_t> opInfo(3);
                    opInfo[0] = id_;
                    opInfo[1] = q;
                    opInfo[2] = p;
                    std::vector<Argument<Base> > args(2);
                    args[0] = Argument<Base>(*txArray);
                    args[1] = Argument<Base>(*tyArray);

                    SourceCodeFragment<Base>* atomicOp = new SourceCodeFragment<Base>(CGAtomicForwardOp, opInfo, args);
                    handler->manageSourceCodeBlock(atomicOp);
                    handler->registerAtomicFunction(id_, this->afun_name());

                    opInfo.resize(1);
                    args.resize(2);
                    for (size_t i = 0; i < ty.size(); i++) {
                        if (tyb.size() == 0 || vyLocal.size() == 0 || vyLocal[i]) {
                            opInfo[0] = i;
                            args[0] = Argument<Base>(*tyArray);
                            args[1] = Argument<Base>(*atomicOp);

                            ty[i] = CGB(*handler, new SourceCodeFragment<Base>(CGArrayElementOp, opInfo, args));
                            if (valuesDefined) {
                                ty[i].setValue(tyb[i]);
                            }
                        } else {
                            ty[i] = tyb[i]; // not a variable (zero)
                        }
                    }
                }
            }

            return true;
        }

        virtual bool reverse(size_t p,
                             const vector<CGB>& tx,
                             const vector<CGB>& ty,
                             vector<CGB>& px,
                             const vector<CGB>& py) {

            bool valuesDefined = true;
            for (size_t i = 0; i < tx.size(); i++) {
                if (!tx[i].isValueDefined()) {
                    valuesDefined = false;
                    break;
                }
            }
            if (valuesDefined) {
                for (size_t i = 0; i < ty.size(); i++) {
                    if (!ty[i].isValueDefined()) {
                        valuesDefined = false;
                        break;
                    }
                }
                if (valuesDefined) {
                    for (size_t i = 0; i < py.size(); i++) {
                        if (!py[i].isValueDefined()) {
                            valuesDefined = false;
                            break;
                        }
                    }
                }
            }

            vector<Base> txb, tyb, pxb, pyb;
            if (valuesDefined) {
                txb.resize(tx.size());
                tyb.resize(ty.size());
                pxb.resize(px.size());
                pyb.resize(py.size());

                for (size_t i = 0; i < txb.size(); i++) {
                    txb[i] = tx[i].getValue();
                }
                for (size_t i = 0; i < tyb.size(); i++) {
                    tyb[i] = ty[i].getValue();
                }
                for (size_t i = 0; i < pyb.size(); i++) {
                    pyb[i] = py[i].getValue();
                }

                if (!atomicReverse(p, txb, tyb, pxb, pyb))
                    return false;
            }

            bool allParameters = isParameters(tx);
            if (allParameters) {
                allParameters = isParameters(ty);
                if (allParameters) {
                    allParameters = isParameters(py);
                }
            }

            if (allParameters) {
                assert(pxb.size() == px.size());

                for (size_t i = 0; i < px.size(); i++) {
                    px[i] = pxb[i];
                }
            } else {
                CodeHandler<Base>* handler = findHandler(tx);
                if (handler == NULL) {
                    handler = findHandler(ty);
                    if (handler == NULL) {
                        handler = findHandler(py);
                    }
                }
                assert(handler != NULL);

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
                    vxLocal[j * (p + 1)] = st[j].size() > 0;
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
                        r[j].insert(0);
                    }
                    for (size_t i = 0; i < m; i++) {
                        s[i] = !py[i * (p + 1) + 1].isParameter() || !py[i * (p + 1) + 1].IdenticalZero();
                    }

                    this->rev_sparse_hes(vx, s, t, 1, r, u, v);

                    for (size_t j = 0; j < n; j++) {
                        vxLocal[j * (p + 1) + 1] = v[j].size() > 0;
                    }
                }

                bool allZero = false;
                if (vxLocal.size() > 0) {
                    allZero = true;
                    for (size_t j = 0; j < vxLocal.size(); j++) {
                        if (vxLocal[j]) {
                            allZero = false;
                            break;
                        }
                    }
                }

                if (allZero) {
                    for (size_t j = 0; j < px.size(); j++) {
                        px[j] = Base(0.0);
                    }
                } else {

                    SourceCodeFragment<Base>* txArray = makeArray(*handler, tx);
                    SourceCodeFragment<Base>* tyArray;
                    SourceCodeFragment<Base>* pxArray = makeZeroArray(*handler, px);
                    SourceCodeFragment<Base>* pyArray = makeArray(*handler, py);

                    if (standAlone_) {
                        tyArray = makeZeroArray(*handler, ty);
                    } else {
                        tyArray = makeArray(*handler, ty);
                    }

                    std::vector<size_t> opInfo(2);
                    opInfo[0] = id_;
                    opInfo[1] = p;
                    std::vector<Argument<Base> > args(4);
                    args[0] = Argument<Base>(*txArray);
                    args[1] = Argument<Base>(*tyArray);
                    args[2] = Argument<Base>(*pxArray);
                    args[3] = Argument<Base>(*pyArray);

                    SourceCodeFragment<Base>* atomicOp = new SourceCodeFragment<Base>(CGAtomicReverseOp, opInfo, args);
                    handler->manageSourceCodeBlock(atomicOp);
                    handler->registerAtomicFunction(id_, this->afun_name());

                    opInfo.resize(1);
                    args.resize(2);
                    for (size_t j = 0; j < px.size(); j++) {
                        if (pxb.size() == 0 || vxLocal.size() == 0 || vxLocal[j]) {
                            opInfo[0] = j;
                            args[0] = Argument<Base>(*pxArray);
                            args[1] = Argument<Base>(*atomicOp);
                            px[j] = CGB(*handler, new SourceCodeFragment<Base>(CGArrayElementOp, opInfo, args));
                            if (valuesDefined) {
                                px[j].setValue(pxb[j]);
                            }
                        } else {
                            px[j] = pxb[j]; // not a variable (zero)
                        }
                    }
                }
            }

            return true;
        }

        virtual ~CGAbstractAtomicFun() {
        }

    protected:

        /**
         * Used to evaluate function values and forward mode function vales and
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

        /**
         * Provides a unique identifier for this atomic function type.
         * 
         * @return a unique identifier ID
         */
        size_t getId() const {
            return id_;
        }

    private:

        static size_t createNewId() {
            CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL;
            static size_t count = 0;
            count++;
            return count;
        }

        static inline CodeHandler<Base>* findHandler(const vector<CGB>& ty) {
            for (size_t i = 0; i < ty.size(); i++) {
                if (ty[i].getCodeHandler() != NULL) {
                    return ty[i].getCodeHandler();
                }
            }
            return NULL;
        }

        static inline std::vector<Arg> asArguments(const vector<CGB>& tx) {
            std::vector<Arg> arguments(tx.size());
            for (size_t i = 0; i < arguments.size(); i++) {
                if (tx[i].isParameter()) {
                    arguments[i] = Arg(tx[i].getValue());
                } else {
                    arguments[i] = Arg(*tx[i].getSourceCodeFragment());
                }
            }
            return arguments;
        }

        static inline SourceCodeFragment<Base>* makeArray(CodeHandler<Base>& handler,
                                                          const vector<CGB>& tx) {
            if (tx.size() > 0) {
                SourceCodeFragment<Base>* op = tx[0].getSourceCodeFragment();
                if (op != NULL && op->operation() == CGArrayElementOp) {
                    SourceCodeFragment<Base>* otherArray = op->arguments()[0].operation();
                    bool reuseArray = true;
                    for (size_t i = 0; i < tx.size(); i++) {
                        op = tx[i].getSourceCodeFragment();
                        if (op == NULL ||
                                op->operation() != CGArrayElementOp ||
                                op->arguments()[0].operation() != otherArray ||
                                op->info()[0] != i) {
                            reuseArray = false;
                            break;
                        }
                    }
                    if (reuseArray) {
                        return otherArray;
                    }
                }
            }
            std::vector<Arg> arrayArgs = asArguments(tx);
            std::vector<size_t> info; // empty
            SourceCodeFragment<Base>* array = new SourceCodeFragment<Base>(CGArrayCreationOp, info, arrayArgs);
            handler.manageSourceCodeBlock(array);
            return array;
        }

        static inline SourceCodeFragment<Base>* makeZeroArray(CodeHandler<Base>& handler,
                                                              const vector<CGB>& tx) {
            vector<CGB> tx2(tx.size());
            std::vector<Arg> arrayArgs = asArguments(tx2);
            std::vector<size_t> info; // empty
            SourceCodeFragment<Base>* array = new SourceCodeFragment<Base>(CGArrayCreationOp, info, arrayArgs);
            handler.manageSourceCodeBlock(array);
            return array;
        }

        static inline bool isParameters(const vector<CGB>& tx) {
            for (size_t i = 0; i < tx.size(); i++) {
                if (!tx[i].isParameter()) {
                    return false;
                }
            }
            return true;
        }

    };

}

#endif