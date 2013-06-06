#ifndef CPPAD_CG_ATOMIC_BASE_INCLUDED
#define CPPAD_CG_ATOMIC_BASE_INCLUDED

#include "cg_cg.hpp"
#include "cg_source_code_fragment.hpp"

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
    class CGAtomicFun : public atomic_base<CppAD::CG<Base> > {
    public:
        typedef CppAD::CG<Base> CGB;
        typedef Argument<Base> Arg;
    protected:
        atomic_base<Base>& atomicFun_;
        const size_t id_;
    public:

        CGAtomicFun(atomic_base<Base>& atomicFun) :
            atomic_base<CGB>(atomicFun.afun_name().c_str()),
            atomicFun_(atomicFun),
            id_(createNewId()) {

        }

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

                if (!atomicFun_.forward(q, p, vx, vy, txb, tyb))
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

                SourceCodeFragment<Base>* txArray = makeArray(*handler, tx);
                SourceCodeFragment<Base>* tyArray = makeArray(*handler, ty);
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

                vector<bool> vyLocal;
                if (p == 0) {
                    vyLocal = vy;
                } else {
                    if (p >= 1) {
                        /**
                         * Use the jacobian sparsity to determine which elements
                         * will always be zero
                         */
                        size_t m = ty.size() / (p + 1);
                        size_t n = tx.size() / (p + 1);

                        vector< std::set<size_t> > r(n);
                        for (size_t j = 0; j < n; j++) {
                            r[j].insert(0);
                        }
                        vector< std::set<size_t> > s(m);
                        for_sparse_jac(1, r, s);

                        vyLocal.resize(ty.size());
                        for (size_t i = 0; i < vyLocal.size(); i++) {
                            vyLocal[i] = true;
                        }

                        for (size_t i = 0; i < m; i++) {
                            vyLocal[i * (p + 1) + 1] = s[i].size() > 0;
                        }
                    }
                }

                opInfo.resize(1);
                args.resize(2);
                for (size_t i = 0; i < ty.size(); i++) {
                    if (vyLocal.size() == 0 || vyLocal[i]) {
                        opInfo[0] = i;
                        args[0] = Argument<Base>(*tyArray);
                        args[1] = Argument<Base>(*atomicOp);

                        ty[i] = CGB(*handler, new SourceCodeFragment<Base>(CGArrayElementOp, opInfo, args));
                        if (valuesDefined) {
                            ty[i].setValue(tyb[i]);
                        }
                    } else {
                        ty[i] = tyb[i]; // not a variable
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

                if (!atomicFun_.reverse(p, txb, tyb, pxb, pyb))
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
                    rt[i].insert(0);
                }
                vector< std::set<size_t> > st(n);
                rev_sparse_jac(1, rt, st);

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

                    rev_sparse_hes(vx, s, t, 1, r, u, v);

                    for (size_t j = 0; j < n; j++) {
                        vxLocal[j * (p + 1) + 1] = v[j].size() > 0;
                    }
                }

                SourceCodeFragment<Base>* txArray = makeArray(*handler, tx);
                SourceCodeFragment<Base>* tyArray = makeArray(*handler, ty);
                SourceCodeFragment<Base>* pxArray = makeArray(*handler, px);
                SourceCodeFragment<Base>* pyArray = makeArray(*handler, py);

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
                    if (vxLocal.size() == 0 || vxLocal[j]) {
                        opInfo[0] = j;
                        args[0] = Argument<Base>(*pxArray);
                        args[1] = Argument<Base>(*atomicOp);
                        px[j] = CGB(*handler, new SourceCodeFragment<Base>(CGArrayElementOp, opInfo, args));
                        if (valuesDefined) {
                            px[j].setValue(pxb[j]);
                        }
                    } else {
                        px[j] = pxb[j]; // not a variable
                    }
                }

            }

            return true;
        }

        virtual bool for_sparse_jac(size_t q,
                                    const vector< std::set<size_t> >& r,
                                    vector< std::set<size_t> >& s) {
            return atomicFun_.for_sparse_jac(q, r, s);
        }

        virtual bool for_sparse_jac(size_t q,
                                    const vector<bool>& r,
                                    vector<bool>& s) {
            return atomicFun_.for_sparse_jac(q, r, s);
        }

        virtual bool rev_sparse_jac(size_t q,
                                    const vector< std::set<size_t> >& rt,
                                    vector< std::set<size_t> >& st) {
            return atomicFun_.rev_sparse_jac(q, rt, st);
        }

        virtual bool rev_sparse_jac(size_t q,
                                    const vector<bool>& rt,
                                    vector<bool>& st) {
            return atomicFun_.rev_sparse_jac(q, rt, st);
        }

        virtual bool rev_sparse_hes(const vector<bool>& vx,
                                    const vector<bool>& s,
                                    vector<bool>& t,
                                    size_t q,
                                    const vector< std::set<size_t> >& r,
                                    const vector< std::set<size_t> >& u,
                                    vector< std::set<size_t> >& v) {
            return atomicFun_.rev_sparse_hes(vx, s, t, q, r, u, v);
        }

        virtual bool rev_sparse_hes(const vector<bool>& vx,
                                    const vector<bool>& s,
                                    vector<bool>& t,
                                    size_t q,
                                    const vector<bool>& r,
                                    const vector<bool>& u,
                                    vector<bool>& v) {
            return atomicFun_.rev_sparse_hes(vx, s, t, q, r, u, v);
        }

        virtual ~CGAtomicFun() {
        }

    protected:

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