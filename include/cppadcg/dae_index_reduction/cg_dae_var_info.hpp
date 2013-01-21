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
        /** Defines the independent variables that dependent on the integrated
         *  variable
         */
        bool integratedDependent_;
        /** Whether or not this variable is an integrated variable (usually 
         * associated with time)
         */
        bool integratedVariable_;
        /**
         * A custom variable name (keep it empty to use an automatically 
         * generated name)
         */
        std::string name_;
    public:

        inline DaeVarInfo(const std::string& name = "") :
            derivativeOf_(-1),
            integratedDependent_(true),
            integratedVariable_(false),
            name_(name) {
        }

        inline DaeVarInfo(int derivativeOf, const std::string& name = "") :
            derivativeOf_(derivativeOf),
            integratedDependent_(true),
            integratedVariable_(false),
            name_(name) {
        }

        /**
         * The index of the variable that the current variable is the derivative
         * of. A negative value means that the current variable isn't a 
         * derivative.
         * 
         * \return The index of the variable for which this variable is the 
         *         derivative of.
         */
        inline int getDerivativeOf() const {
            return derivativeOf_;
        }

        /**
         * Determines whether or not this variable depends on the 
         * independent/integrated variables.
         * 
         * \return true if it is a parameter that does not depend on the
         *              integrated variables
         */
        inline bool isFunctionOfIntegrated() const {
            return integratedDependent_;
        }

        /**
         * Defines this variable as a parameter/constant that does not depend on
         * the independent/integrated variables
         */
        inline void makeConstant() {
            integratedVariable_ = false;
            integratedDependent_ = false;
            derivativeOf_ = -1;
        }

        /**
         * Defines this variable as an integrated variable, also known
         * as the independent variable of the DAE system (usually associated
         * with time)
         */
        inline void makeIntegratedVariable() {
            integratedVariable_ = true;
            integratedDependent_ = false;
            derivativeOf_ = -1;
        }

        /**
         * Determines whether or not this is an integrated variable, also known
         * as the independent variable of the DAE system (typically time).
         * 
         * \return true if it is the integrated variable
         */
        inline bool isIntegratedVariable() const {
            return integratedVariable_;
        }

        /**
         * Returns the custom variable name. If the string is empty an 
         * automatically generated name will be used.
         * 
         * \return the custom variable name
         */
        inline const std::string& getName() const {
            return name_;
        }

        /**
         * Defines a custom variable name. If the string is empty an 
         * automatically generated name will be used.
         * 
         * \param name the custom variable name
         */
        inline void setName(const std::string& name) {
            name_ = name;
        }

        inline virtual ~DaeVarInfo() {
        }
    };

}

#endif
