#ifndef CPPAD_CG_DAE_VAR_INFO_INCLUDED
#define CPPAD_CG_DAE_VAR_INFO_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
        
        inline virtual ~DaeVarInfo() {
        }
    };

}

#endif
