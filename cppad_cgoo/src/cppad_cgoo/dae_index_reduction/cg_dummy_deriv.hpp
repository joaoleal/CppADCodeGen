#ifndef CPPAD_CG_DUMMY_DERIV_INCLUDED
#define	CPPAD_CG_DUMMY_DERIV_INCLUDED

#include "cg_bipartite.hpp"

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
namespace CppAD {

    /**
     * Dummy derivatives DAE index reduction algorithm
     */
    template<class Base>
    class DummyDerivatives : public Plantelides<Base> {
    public:

        DummyDerivatives(ADFun<CG<Base> >* fun,
                         const std::vector<bool>& eqDifferentialInfo,
                         const std::vector<bool>& varInfo) :
            Plantelides<Base>(fun, eqDifferentialInfo, varInfo) {
        }

        virtual inline void reduceIndex() {
            Plantelides<Base>::reduceIndex();

            std::set<Vnode<Base>* > vars; // variables of interest
            for (size_t j = 0; j < this->vnodes_.size(); j++) {
                Vnode<Base>* jj = this->vnodes_[j];
                if (jj->derivativeOf() != NULL && jj->derivative() == NULL) {
                    vars.insert(jj);
                }
            }

            std::set<Enode<Base>* > eqs; // equations of interest
            eqs.insert(this->enodes_.begin(), this->enodes_.end());

            while (true) {

                /**
                 * Consider all of the current equations that are
                 * differentiated versions of the original ones.
                 * Collect their predecessors and let them be the
                 * current equations.
                 */
                std::set<Enode<Base>* > newEqs;
                typename std::set<Enode<Base>* >::const_iterator i;
                for (i = eqs.begin(); i != eqs.end(); ++i) {
                    Enode<Base>* ii = *i;
                    if (ii->derivativeOf() != NULL) {
                        newEqs.insert(ii->derivativeOf());
                    }
                }
                eqs = newEqs;

                if (eqs.empty()) {
                    break;
                }

                /**
                 * Consider all current unknowns that are at least of
                 * order one. Collect their predecessors of one order
                 * less and let them be the current candidates for
                 * elimination.
                 */
                std::set<Vnode<Base>* > varsNew;
                std::map<Vnode<Base>*, int> order;
                for (size_t j = 0; j < varsNew.size(); j++) {
                    Vnode<Base>* jj = varsNew[j];
                    Vnode<Base>* v = jj->derivativeOf();
                    if (v != NULL) {
                        varsNew.insert(v);
                        int o = 1;
                        while (v != NULL) {
                            v = v->derivativeOf();
                            o++;
                        }
                        order[jj->derivativeOf()] = o;
                    }
                }
                vars = varsNew;

                // Exploit the current equations for elimination of candidates


            }
        }

    };
}


#endif

