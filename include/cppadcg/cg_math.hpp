#ifndef CPPAD_CG_MATH_INCLUDED
#define CPPAD_CG_MATH_INCLUDED
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

#define CPPAD_CG_CREATE_OPERATION(OpName, OpCode)                              \
    template<class Base>                                                       \
    inline CG<Base> OpName(const CG<Base>& var) {                              \
        if (var.isParameter()) {                                               \
            return CG<Base> (OpName(var.getValue()));                          \
        } else {                                                               \
            CG<Base> result(*var.getCodeHandler(), new OperationNode<Base>(OpCode, var.argument()));\
            if(var.isValueDefined())                                           \
                result.setValue(OpName(var.getValue()));                       \
            return result;                                                     \
        }                                                                      \
    }

    CPPAD_CG_CREATE_OPERATION(abs, CGAbsOp)
    CPPAD_CG_CREATE_OPERATION(acos, CGAcosOp)
    CPPAD_CG_CREATE_OPERATION(asin, CGAsinOp)
    CPPAD_CG_CREATE_OPERATION(atan, CGAtanOp)
    CPPAD_CG_CREATE_OPERATION(cos, CGCosOp)
    CPPAD_CG_CREATE_OPERATION(cosh, CGCoshOp)
    CPPAD_CG_CREATE_OPERATION(exp, CGExpOp)
    CPPAD_CG_CREATE_OPERATION(log, CGLogOp)
    CPPAD_CG_CREATE_OPERATION(sin, CGSinOp)
    CPPAD_CG_CREATE_OPERATION(sinh, CGSinhOp)
    CPPAD_CG_CREATE_OPERATION(sqrt, CGSqrtOp)
    CPPAD_CG_CREATE_OPERATION(tan, CGTanOp)
    CPPAD_CG_CREATE_OPERATION(tanh, CGTanhOp)
}
#endif

