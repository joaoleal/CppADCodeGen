#ifndef CPPAD_CG_ATOMIC_FUN_BRIDGE_INCLUDED
#define CPPAD_CG_ATOMIC_FUN_BRIDGE_INCLUDED
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
         */
        CGAtomicFunBridge(const std::string& name,
                          CppAD::ADFun<CGB>& fun,
                          bool standAlone = false) :
            CGAbstractAtomicFun<Base>(name, standAlone),
            fun_(fun) {

        }

        virtual bool for_sparse_jac(size_t q,
                                    const vector< std::set<size_t> >& r,
                                    vector< std::set<size_t> >& s) {
            s = fun_.ForSparseJac(q, r);
            fun_.size_forward_set(0);
            return true;
        }

        virtual bool for_sparse_jac(size_t q,
                                    const vector<bool>& r,
                                    vector<bool>& s) {
            s = fun_.ForSparseJac(q, r);
            fun_.size_forward_bool(0);
            return true;
        }

        virtual bool rev_sparse_jac(size_t q,
                                    const vector< std::set<size_t> >& rt,
                                    vector< std::set<size_t> >& st) {
            st = fun_.RevSparseJac(q, rt, true);
            return true;
        }

        virtual bool rev_sparse_jac(size_t q,
                                    const vector<bool>& rt,
                                    vector<bool>& st) {
            st = fun_.RevSparseJac(q, rt, true);
            return true;
        }

        virtual bool rev_sparse_hes(const vector<bool>& vx,
                                    const vector<bool>& s,
                                    vector<bool>& t,
                                    size_t q,
                                    const vector< std::set<size_t> >& r,
                                    const vector< std::set<size_t> >& u,
                                    vector< std::set<size_t> >& v) {
            size_t m = fun_.Range();
            size_t n = fun_.Domain();

            t = fun_.RevSparseJac(1, s);
            vector< std::set<size_t> > a(n);
            a = fun_.RevSparseJac(q, u, true);

            // set version of s
            vector< std::set<size_t> > set_s(1);
            CPPAD_ASSERT_UNKNOWN(set_s[0].empty());

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
            return true;
        }

        virtual bool rev_sparse_hes(const vector<bool>& vx,
                                    const vector<bool>& s,
                                    vector<bool>& t,
                                    size_t q,
                                    const vector<bool>& r,
                                    const vector<bool>& u,
                                    vector<bool>& v) {
            size_t n = fun_.Domain();

            t = fun_.RevSparseJac(1, s);

            vector<bool> a(n * q);
            a = fun_.RevSparseJac(q, u, true);

            fun_.ForSparseJac(q, r);
            v = fun_.RevSparseHes(q, s, true);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < q; j++)
                    v[ i * q + j ] |= a[ i * q + j];
            }

            fun_.size_forward_set(0);

            return true;
        }

        virtual ~CGAtomicFunBridge() {
        }

    protected:

        virtual bool atomicForward(size_t q,
                                   size_t p,
                                   const vector<bool>& vx,
                                   vector<bool>& vy,
                                   const vector<Base>& tx,
                                   vector<Base>& ty) {

            if (vx.size() > 0) {
                zeroOrderDependency(fun_, vx, vy);
            }
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