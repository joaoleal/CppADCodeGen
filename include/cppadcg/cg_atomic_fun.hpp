#ifndef CPPAD_CG_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_ATOMIC_FUN_INCLUDED
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
    class CGAtomicFun : public CGAbstractAtomicFun<Base> {
    protected:
        atomic_base<Base>& atomicFun_;
    public:

        /**
         * Creates a new atomic function wrapper that is responsible for 
         * defining the dependencies to calls of a user atomic function.
         * 
         * @param atomicFun The atomic function to the called by the compiled
         *                  source.
         * @param standAlone Whether or not forward and reverse function calls
         *                   do not require the Taylor coefficients for the 
         *                   dependent variables (ty) and the previous
         *                   evaluation of other forward/reverse modes.
         */
        CGAtomicFun(atomic_base<Base>& atomicFun, bool standAlone = false) :
            CGAbstractAtomicFun<Base>(atomicFun.afun_name(), standAlone),
            atomicFun_(atomicFun) {

        }

        template <class ADVector>
        void operator()(const ADVector& ax, ADVector& ay, size_t id = 0) {
            this->CGAbstractAtomicFun<Base>::operator()(ax, ay, id);
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

        virtual bool atomicForward(size_t q,
                                   size_t p,
                                   const vector<bool>& vx,
                                   vector<bool>& vy,
                                   const vector<Base>& tx,
                                   vector<Base>& ty) {
            return atomicFun_.forward(q, p, vx, vy, tx, ty);
        }

        virtual bool atomicReverse(size_t p,
                                   const vector<Base>& tx,
                                   const vector<Base>& ty,
                                   vector<Base>& px,
                                   const vector<Base>& py) {
            return atomicFun_.reverse(p, tx, ty, px, py);
        }
    };

}

#endif