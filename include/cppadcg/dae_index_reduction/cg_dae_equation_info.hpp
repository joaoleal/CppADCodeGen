#ifndef CPPAD_CG_DAE_EQUATION_INFO_HPP
#define	CPPAD_CG_DAE_EQUATION_INFO_HPP
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
     * DAE equation information
     */
    class DaeEquationInfo {
    private:
        /**
         * The equation index in the original user model
         */
        int originalIndex_;
        /**
         * The index of the equation that was differentiated to obtain this
         * equation. A negative value means that the current equation isn't a
         * differentiation of an existing equation.
         */
        int derivativeOf_;
        /**
         * The variable index associated with this equation. A negative value is
         * used if this equation does not have an assigned variable.
         */
        int assignedVarIndex_;
        /**
         * Whether or not if it is an explicit differential equation
         */
        bool explicit_;
    public:

        inline DaeEquationInfo() :
            originalIndex_(-1),
            derivativeOf_(-1),
            assignedVarIndex_(-1),
            explicit_(false) {
        }

        inline DaeEquationInfo(int originalIndex, int derivativeOf, int assignedVarIndex, bool explicitEq = false) :
            originalIndex_(originalIndex),
            derivativeOf_(derivativeOf),
            assignedVarIndex_(assignedVarIndex),
            explicit_(explicitEq) {
        }

        inline int getDerivativeOf() const {
            return derivativeOf_;
        }

        inline void setDerivativeOf(int derivativeOf) {
            derivativeOf_ = derivativeOf;
        }

        inline int getAssignedVarIndex() const {
            return assignedVarIndex_;
        }

        inline void setAssignedVarIndex(int assignedVarIndex) {
            assignedVarIndex_ = assignedVarIndex;
        }

        inline int getOriginalIndex() const {
            return originalIndex_;
        }

        inline bool isExplicit() const {
            return explicit_;
        }

        inline void setExplicit(bool explicitEq) {
            explicit_ = explicitEq;
        }

        inline virtual ~DaeEquationInfo() {
        }
    };

}

#endif
