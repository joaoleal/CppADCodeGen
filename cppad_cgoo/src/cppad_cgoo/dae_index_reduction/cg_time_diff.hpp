#ifndef CPPAD_CG_TIME_DIFF_INCLUDED
#define	CPPAD_CG_TIME_DIFF_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    // ----------------------------------------------------------------------
    // forward mode routine called by CppAD

    template<class Base>
    bool time_diff_forward(size_t id,
                           size_t order,
                           size_t n,
                           size_t m,
                           const vector<bool>& vx,
                           vector<bool>& vzy,
                           const vector<CG<Base> >& tx,
                           vector<CG<Base> >& tzy) {
        assert(n == 3); // [x, dxdt, t]
        assert(m == 1);
        assert(tx.size() >= (order + 1) * n);
        assert(tzy.size() >= (order + 1) * m);

        size_t n_order = order + 1;
        const size_t xIndex = 0; // index of the variable in the argument list
        const size_t dxdtIndex = 1; // index of the time derivative variable in the argument list
        const size_t timeIndex = 2; // index of the time variable in the argument list

        // check if this is during the call to time_var(id, ax, ay)
        if (vx.size() > 0) {
            assert(vx.size() >= n);
            assert(vzy.size() >= m);

            vzy[0] = vx[0] || vx[1] || vx[2];
        }

        if (order == 0) {
            tzy[0] = tx[0];
        } else if (order == 1) {
            const CG<Base>& ttime = tx[timeIndex * n_order + order]; //
            const CG<Base>& txx = tx[xIndex * n_order + order]; //
            assert(ttime.isParameter());
            assert(txx.isParameter());
            if (ttime.getParameterValue() > 0) {
                assert(txx.getParameterValue() == 0);
                tzy[1] = ttime * tx[dxdtIndex * n_order + 0]; // transform x(t) into dx(t)/dt
            } else {
                tzy[1] = txx; // do nothing
            }

        } else {
            return false; // not implemented
        }

        // All orders are implemented and there are no possible errors
        return true;
    }
    // ----------------------------------------------------------------------
    // reverse mode routine called by CppAD

    template<class Base>
    bool time_diff_reverse(size_t id,
                           size_t order,
                           size_t n,
                           size_t m,
                           const vector<CG<Base> >& tx,
                           const vector<CG<Base> >& tzy,
                           vector<CG<Base> >& px,
                           const vector<CG<Base> >& pzy) {
        return false;
    }
    // ----------------------------------------------------------------------
    // forward Jacobian sparsity routine called by CppAD

    template<class Base>
    bool time_diff_for_jac_sparse(size_t id,
                                  size_t n,
                                  size_t m,
                                  size_t q,
                                  const vector< std::set<size_t> >& r,
                                  vector< std::set<size_t> >& s) {
        assert(n == 3);
        assert(m == 1);
        assert(r.size() >= n);
        assert(s.size() >= m);

        // sparsity for z and y are the same as for x
        s[0] = r[0]; // x
        s[0].insert(r[1].begin(), r[1].end()); // dxdt
        s[0].insert(r[2].begin(), r[2].end()); // t

        return true;
    }
    // ----------------------------------------------------------------------
    // reverse Jacobian sparsity routine called by CppAD

    template<class Base>
    bool time_diff_rev_jac_sparse(size_t id,
                                  size_t n,
                                  size_t m,
                                  size_t q,
                                  vector< std::set<size_t> >& r,
                                  const vector< std::set<size_t> >& s) {
        return false;
    }
    // ----------------------------------------------------------------------
    // reverse Hessian sparsity routine called by CppAD

    template<class Base>
    bool time_diff_rev_hes_sparse(size_t id,
                                  size_t n,
                                  size_t m,
                                  size_t q,
                                  const vector< std::set<size_t> >& r,
                                  const vector<bool>& s,
                                  vector<bool>& t,
                                  const vector< std::set<size_t> >& u,
                                  vector< std::set<size_t> >& v) {
        return false;
    }
    // ---------------------------------------------------------------------
    // Declare the AD<CG<Base> > routine time_var(id, ax, ay)
    template<class Base>
    CPPAD_USER_ATOMIC(time_var,
                      std::vector,
                      CG<Base>,
                      time_diff_forward<Base>,
                      time_diff_reverse<Base>,
                      time_diff_for_jac_sparse<Base>,
                      time_diff_rev_jac_sparse<Base>,
                      time_diff_rev_hes_sparse<Base>)
}

#endif
