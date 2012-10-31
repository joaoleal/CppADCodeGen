#ifndef CPPAD_CG_GRAPH_MOD_INCLUDED
#define	CPPAD_CG_GRAPH_MOD_INCLUDED

#include "cg_cg.hpp"

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<class Base>
    inline void CodeHandler<Base>::substituteIndependent(const CG<Base>& indep,
    const CG<Base>& dep) throw (CGException) {
        substituteIndependent(indep.getSourceCodeFragment(), dep.getSourceCodeFragment());
    }

    template<class Base>
    inline void CodeHandler<Base>::substituteIndependent(SourceCodeFragment<Base>* indep,
    SourceCodeFragment<Base>* dep) throw (CGException) {
        using std::vector;
        typedef CG<Base> CGBase;
        typedef AD<CGBase> ADCG;

        assert(indep != NULL);
        assert(dep != NULL);

        //check if the independent variable belongs to this handler
        size_t indepIndex = getIndependentVariableIndex(*indep);

        //check if the dependent variable belongs to this handler
        typename std::vector<SourceCodeFragment<Base> *>::const_iterator it =
                std::find(_codeBlocks.begin(), _codeBlocks.end(), dep);
        if (it == _codeBlocks.end()) {
            throw CGException("The dependent variable does not belong to this handler");
        }

        // determine the expression for the independent variable
        CGBase dummyExp = solveFor(dep, indep);

        Argument<Base> arg;
        // change the independent variable
        if (dummyExp.isVariable()) {
            arg = Argument<Base > (*dummyExp.getSourceCodeFragment());
        } else {
            // create a bogus variable to avoid searching for all occurrences of the independent variable
            arg = Argument<Base > (dummyExp.getParameterValue());
        }

        indep->makeAlias(arg);

        // remove the substituted variable from the independent variable vector
        _independentVariables.erase(_independentVariables.begin() + indepIndex);
    }

}

#endif
