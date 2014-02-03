#ifndef CPPAD_CG_BASE_DOUBLE_INCLUDED
#define CPPAD_CG_BASE_DOUBLE_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
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
inline cg::CG<double> abs(const cg::CG<double>& var) {
    using namespace CppAD::cg;

    if (var.isParameter()) {
        return CG<double> (fabs(var.getValue()));
    } else {
        CG<double> result(*var.getCodeHandler(), new OperationNode<double>(CGAbsOp, var.argument()));
        if (var.isValueDefined()) {
            result.setValue(fabs(var.getValue()));
        }
        return result;
    }
}

/**
 * Specialization of the numeric_limits for doubles
 */
template <>
class numeric_limits<cg::CG<double> > {
public:

    static cg::CG<double> epsilon() {
        return std::numeric_limits<double>::epsilon();
    }

    static cg::CG<double> min() {
        return std::numeric_limits<double>::min();
    }

    static cg::CG<double> max() {
        return std::numeric_limits<double>::max();
    }
};

/**
 * Specialization of the machine epsilon for CG<double>
 */
template <>
inline cg::CG<double> epsilon<cg::CG<double> >() {
    return std::numeric_limits<double>::epsilon();
}

} // END CppAD namespace

#endif

