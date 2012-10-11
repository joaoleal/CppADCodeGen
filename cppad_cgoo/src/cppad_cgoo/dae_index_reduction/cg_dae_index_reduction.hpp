#ifndef CPPAD_CG_DAE_INDEX_REDUCTION_INCLUDED
#define	CPPAD_CG_DAE_INDEX_REDUCTION_INCLUDED
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
     * Base class for algorithms that perform automatic (differentiation) index
     * reduction of semi-implicit DAEs.
     */
    template<class Base>
    class DaeIndexReduction {
    protected:
        ADFun<CG<Base> >* fun_;
        // defines the equations that are differential: dxdt = ...
        const std::vector<bool> eqDifferentialInfo_;
        // defines the independent variables that dependent on time
        const std::vector<bool> varInfo_;
    public:

        /**
         * Creates a new DAE model index reduction algorithm.
         * 
         * \param fun  The original (high index) model
         * \param eqDifferentialInfo  Defines the equations that are
         *                             differential
         * \param varInfo  Defines the independent variables that dependent on
         *                  time
         */
        DaeIndexReduction(ADFun<CG<Base> >* fun,
                          const std::vector<bool>& eqDifferentialInfo,
                          const std::vector<bool>& varInfo) :
            fun_(fun),
            eqDifferentialInfo_(eqDifferentialInfo),
            varInfo_(varInfo) {
            assert(fun_ != NULL);
            assert(eqDifferentialInfo_.size() == fun->Range());
            assert(varInfo_.size() == fun->Domain());
        }
    };
}

#endif	

