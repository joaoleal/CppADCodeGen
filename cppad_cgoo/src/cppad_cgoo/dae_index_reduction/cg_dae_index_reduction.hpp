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
     * DAE variable information
     */
    class DaeVarInfo {
    private:
        /** The index of the variable for which this variable is the time
         * derivative. A negative value means that the current variable isn't 
         * a time derivative.
         */
        int derivativeOf_;
        // defines the independent variables that dependent on time
        bool timeDependent_;
    public:

        inline DaeVarInfo() :
            derivativeOf_(-1),
            timeDependent_(true) {
        }

        inline DaeVarInfo(int derivativeOf) :
            derivativeOf_(derivativeOf),
            timeDependent_(true) {
        }

        inline int getDerivativeOf() const {
            return derivativeOf_;
        }

        inline bool isTimeDependent() const {
            return timeDependent_;
        }

        inline void makeTimeIndependent() {
            timeDependent_ = false;
            derivativeOf_ = -1;
        }
    };

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
        // DAE variable information
        const std::vector<DaeVarInfo> varInfo_;
    public:

        /**
         * Creates a new DAE model index reduction algorithm.
         * 
         * \param fun  The original (high index) model
         * \param varInfo  DAE  system variable information (in the same order 
         *                 as in the tape)
         */
        DaeIndexReduction(ADFun<CG<Base> >* fun,
                          const std::vector<DaeVarInfo>& varInfo) :
            fun_(fun),
            varInfo_(varInfo) {
            assert(fun_ != NULL);
            assert(varInfo_.size() == fun->Domain());
        }
    };
}

#endif	

