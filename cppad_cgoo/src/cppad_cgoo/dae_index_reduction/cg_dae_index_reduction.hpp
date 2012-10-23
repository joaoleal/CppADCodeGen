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
     * reduction of implicit DAEs.
     */
    template<class Base>
    class DaeIndexReduction {
    protected:
        /**
         * The original model
         */
        ADFun<CG<Base> >* fun_;
        /** Defines the variable index that represents the time derivative. 
         * A negative value means that there isn't any time derivative for the
         * current variable.
         */
        const std::vector<int> derivative_;
        // defines the independent variables that dependent on time
        const std::vector<bool> timeDependent_;
    public:

        /**
         * Creates a new DAE model index reduction algorithm.
         * 
         * \param fun  The original (high index) model
         * \param derivative  Defines the variable index that represents the 
         *                    time derivative.  A negative value means that 
         *                    there isn't any time derivative for the current 
         *                    variable.
         * \param timeDepebdent Defines the time dependent variables.
         */
        DaeIndexReduction(ADFun<CG<Base> >* fun,
                          const std::vector<int>& derivative,
                          const std::vector<bool>& timeDependent) :
            fun_(fun),
            derivative_(derivative),
            timeDependent_(timeDependent) {
            assert(fun_ != NULL);
            assert(derivative_.size() == fun->Domain());
            assert(timeDependent_.size() == fun->Domain());
        }
    };
}

#endif	

