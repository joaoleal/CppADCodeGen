#ifndef CPPAD_CG_ATOMIC_FUN_BRIDGE_INCLUDED
#define CPPAD_CG_ATOMIC_FUN_BRIDGE_INCLUDED
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

    /**
     * An atomic function wrapper for atomic functions using the CppAD::CG
     * type.
     * This is class can be useful when a CppAD::ADFun<CppAD::CG > is going to
     * be used to create a compiled model library but has not been compiled yet.
     * 
     * @author Joao Leal
     */
    template <class Base>
    class CGAtomicFunBridge : public CGAbstractAtomicFun<Base> {
    public:
        typedef CppAD::CG<Base> CGB;
        typedef CppAD::AD<CGB> ADCGD;
    protected:
        ADFun<CGB>& fun_;
        bool cacheSparsities_;
        CustomPosition custom_jac_;
        CustomPosition custom_hess_;
        std::map<size_t, vector<std::set<size_t> > > hess_;
    public:

        /**
         * Creates a new atomic function wrapper.
         * 
         * @param name The atomic function name
         * @param fun The atomic function to be wrapped
         * @param standAlone Whether or not forward and reverse function calls
         *                   do not require the Taylor coefficients for the 
         *                   dependent variables (ty) and the previous
         *                   evaluation of other forward/reverse modes.
         * @param cacheSparsities Whether or not to cache information related 
         *                        with sparsity evaluation.
         */
        CGAtomicFunBridge(const std::string& name,
                          CppAD::ADFun<CGB>& fun,
                          bool standAlone = false,
                          bool cacheSparsities = true) :
            CGAbstractAtomicFun<Base>(name, standAlone),
            fun_(fun),
            cacheSparsities_(cacheSparsities) {
            this->option(CppAD::atomic_base<CGB>::set_sparsity_enum);
        }

        template <class ADVector>
        void operator()(const ADVector& ax, ADVector& ay, size_t id = 0) {
            this->CGAbstractAtomicFun<Base>::operator()(ax, ay, id);
        }

        template<class VectorSize>
        inline void setCustomSparseJacobianElements(const VectorSize& row,
                                                    const VectorSize& col) {
            custom_jac_ = CustomPosition(fun_.Range(), fun_.Domain(), row, col);
        }

        template<class VectorSet>
        inline void setCustomSparseJacobianElements(const VectorSet& elements) {
            custom_jac_ = CustomPosition(fun_.Range(), fun_.Domain(), elements);
        }

        template<class VectorSize>
        inline void setCustomSparseHessianElements(const VectorSize& row,
                                                   const VectorSize& col) {
            size_t n = fun_.Domain();
            custom_hess_ = CustomPosition(n, n, row, col);
        }

        template<class VectorSet>
        inline void setCustomSparseHessianElements(const VectorSet& elements) {
            size_t n = fun_.Domain();
            custom_hess_ = CustomPosition(n, n, elements);
        }

        virtual bool for_sparse_jac(size_t q,
                                    const vector<std::set<size_t> >& r,
                                    vector<std::set<size_t> >& s) {
            if (cacheSparsities_ || custom_jac_.isFilterDefined()) {
                size_t n = fun_.Domain();
                size_t m = fun_.Range();
                if (!custom_jac_.isFullDefined()) {
                    custom_jac_.setFullElements(extra::jacobianForwardSparsitySet<vector<std::set<size_t> > >(fun_));
                    fun_.size_forward_set(0);
                }

                for (size_t i = 0; i < s.size(); i++) {
                    s[i].clear();
                }
                CppAD::multMatrixMatrixSparsity(custom_jac_.getFullElements(), r, s, m, n, q);
            } else {
                s = fun_.ForSparseJac(q, r);
                fun_.size_forward_set(0);
            }

            return true;
        }

        virtual bool rev_sparse_jac(size_t q,
                                    const vector<std::set<size_t> >& rt,
                                    vector<std::set<size_t> >& st) {

            if (cacheSparsities_ || custom_jac_.isFilterDefined()) {
                size_t n = fun_.Domain();
                size_t m = fun_.Range();
                if (!custom_jac_.isFullDefined()) {
                    custom_jac_.setFullElements(extra::jacobianReverseSparsitySet<vector<std::set<size_t> > >(fun_));
                }

                for (size_t i = 0; i < st.size(); i++) {
                    st[i].clear();
                }
                CppAD::multMatrixMatrixSparsityTrans(rt, custom_jac_.getFullElements(), st, m, n, q);
            } else {
                st = fun_.RevSparseJac(q, rt, true);
            }

            return true;
        }

        virtual bool rev_sparse_hes(const vector<bool>& vx,
                                    const vector<bool>& s,
                                    vector<bool>& t,
                                    size_t q,
                                    const vector<std::set<size_t> >& r,
                                    const vector<std::set<size_t> >& u,
                                    vector<std::set<size_t> >& v) {
            using namespace CppAD::extra;

            if (cacheSparsities_ || custom_jac_.isFilterDefined() || custom_hess_.isFilterDefined()) {
                size_t n = fun_.Domain();
                size_t m = fun_.Range();

                for (size_t i = 0; i < n; i++) {
                    v[i].clear();
                }

                if (!custom_jac_.isFullDefined()) {
                    custom_jac_.setFullElements(jacobianSparsitySet<vector<std::set<size_t> > >(fun_));
                }
                const vector<std::set<size_t> >& jacSparsity = custom_jac_.getFullElements();

                /**
                 *  V(x)  =  f'^T(x) U(x)  +  Sum(  s(x)i  f''(x)  R(x)   )
                 */
                // f'^T(x) U(x)
                CppAD::multMatrixTransMatrixSparsity(jacSparsity, u, v, m, n, q);

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
                        custom_hess_.setFullElements(hessianSparsitySet<vector<std::set<size_t> > >(fun_)); // f''(x)
                    }
                    const vector<std::set<size_t> >& sF2 = custom_hess_.getFullElements();
                    CppAD::multMatrixTransMatrixSparsity(sF2, r, v, n, n, q); // f''^T * R
                } else {
                    vector<std::set<size_t> > sparsitySF2R(n);
                    for (size_t i = 0; i < m; i++) {
                        if (s[i]) {
                            std::map<size_t, vector<std::set<size_t> > >::const_iterator itH = hess_.find(i);
                            const vector<std::set<size_t> >* spari;
                            if (itH == hess_.end()) {
                                vector<std::set<size_t> >& hi = hess_[i] = hessianSparsitySet<vector<std::set<size_t> > >(fun_, i); // f''_i(x)
                                spari = &hi;
                                custom_hess_.filter(hi);
                            } else {
                                spari = &itH->second;
                            }
                            CppAD::addMatrixSparsity(*spari, sparsitySF2R);
                        }
                    }
                    CppAD::multMatrixTransMatrixSparsity(sparsitySF2R, r, v, n, n, q); // f''^T * R
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
            } else {
                size_t m = fun_.Range();
                size_t n = fun_.Domain();

                t = fun_.RevSparseJac(1, s);
                vector<std::set<size_t> > a = fun_.RevSparseJac(q, u, true);

                // set version of s
                vector<std::set<size_t> > set_s(1);
                for (size_t i = 0; i < m; i++) {
                    if (s[i])
                        set_s[0].insert(i);
                }

                fun_.ForSparseJac(q, r);
                v = fun_.RevSparseHes(q, set_s, true);

                std::set<size_t>::const_iterator itr;
                for (size_t i = 0; i < n; i++) {
                    for (itr = a[i].begin(); itr != a[i].end(); itr++) {
                        size_t j = *itr;
                        CPPAD_ASSERT_UNKNOWN(j < q);
                        v[i].insert(j);
                    }
                }

                fun_.size_forward_set(0);
            }

            return true;
        }

        virtual ~CGAtomicFunBridge() {
        }

    protected:

        virtual void zeroOrderDependency(const vector<bool>& vx,
                                         vector<bool>& vy) {
            CppAD::zeroOrderDependency(fun_, vx, vy);
        }

        virtual bool atomicForward(size_t q,
                                   size_t p,
                                   const vector<Base>& tx,
                                   vector<Base>& ty) {
            vector<CGB> txcg(tx.size());
            toCG(tx, txcg);

            vector<CGB> tycg = fun_.Forward(p, txcg);
            fromCG(tycg, ty);

            fun_.capacity_taylor(0);

            return true;
        }

        virtual bool atomicReverse(size_t p,
                                   const vector<Base>& tx,
                                   const vector<Base>& ty,
                                   vector<Base>& px,
                                   const vector<Base>& py) {
            vector<CGB> txcg(tx.size());
            vector<CGB> pycg(py.size());

            toCG(tx, txcg);
            toCG(py, pycg);

            fun_.Forward(p, txcg);

            vector<CGB> pxcg = fun_.Reverse(p + 1, pycg);
            fromCG(pxcg, px);

            fun_.capacity_taylor(0);
            return true;
        }

    private:

        static void toCG(const vector<Base>& from, vector<CGB>& to) {
            CPPAD_ASSERT_UNKNOWN(from.size() == to.size());

            for (size_t i = 0; i < from.size(); i++) {
                to[i] = from[i];
            }
        }

        static void fromCG(const vector<CGB>& from, vector<Base>& to) {
            CPPAD_ASSERT_UNKNOWN(from.size() == to.size());

            for (size_t i = 0; i < from.size(); i++) {
                CPPADCG_ASSERT_KNOWN(from[i].isValueDefined(), "No value defined")
                to[i] = from[i].getValue();
            }
        }

    private:
        CGAtomicFunBridge(const CGAtomicFunBridge& orig); // not implemented
        CGAtomicFunBridge& operator=(const CGAtomicFunBridge& rhs); // not implemented

    };

}

#endif