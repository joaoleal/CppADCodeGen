#ifndef CPPAD_CG_BASE_DOUBLE_INCLUDED
#define CPPAD_CG_BASE_DOUBLE_INCLUDED
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
     * Specialization of the abs operation for doubles
     * 
     * @author Joao Leal
     */
    template<>
    inline CG<double> abs(const CG<double>& var) {
        if (var.isParameter()) {
            return CG<double> (fabs(var.getValue()));
        } else {
            CG<double> result(*var.getCodeHandler(), new OperationNode<double>(CGAbsOp, var.argument()));
            if(var.isValueDefined()) {
                result.setValue(fabs(var.getValue()));
            }
            return result;
        }
    }


}

#endif

