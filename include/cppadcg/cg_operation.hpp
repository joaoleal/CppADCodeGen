#ifndef CPPAD_CG_OPERATION_INCLUDED
#define CPPAD_CG_OPERATION_INCLUDED
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
     * Possible operations
     * 
     * @author Joao Leal
     */
    enum CGOpCode {
        CGAbsOp, //  abs(variable)
        CGAcosOp, // asin(variable)
        CGAddOp, //  a + b
        CGAliasOp, //  alias (reference to another operation)
        CGArrayCreationOp, // {a, b, c ...}
        CGArrayElementOp, // x[i]
        CGAsinOp, // asin(variable)
        CGAtanOp, // atan(variable)
        CGAtomicForwardOp, // atomicFunction.forward(q, p, vx, vy, tx, ty)
        CGAtomicReverseOp, // atomicFunction.reverse(p, tx, ty, px, py)
        CGComOpLt, // result = left < right? trueCase: falseCase
        CGComOpLe, // result = left <= right? trueCase: falseCase
        CGComOpEq, // result = left == right? trueCase: falseCase
        CGComOpGe, // result = left >= right? trueCase: falseCase
        CGComOpGt, // result = left > right? trueCase: falseCase
        CGComOpNe, // result = left != right? trueCase: falseCase
        CGCoshOp, // cosh(variable)
        CGCosOp, //  cos(variable)
        CGDivOp, // a / b
        CGExpOp, //  exp(variable)
        CGInvOp, //                             independent variable
        CGLogOp, //  log(variable)
        CGMulOp, // a * b
        CGPowOp, //  pow(a,   b)
        //        PriOp, //  PrintFor(text, parameter or variable, parameter or variable)
        CGSignOp, // result = (x > 0)? 1.0:((x == 0)? 0.0:-1)
        CGSinhOp, // sinh(variable)
        CGSinOp, //  sin(variable)
        CGSqrtOp, // sqrt(variable)
        CGSubOp, //  a - b
        CGTanhOp, //  tanh(variable)
        CGTanOp, //  tan(variable)
        CGUnMinusOp, // -(a)
        CGDependentMultiAssignOp, // operation which associates a dependent variables with loops and regular operations
        CGDependentRefRhsOp, // operation referencing a dependent variable (right hand side only)
        CGIndexDeclarationOp, // an integer index declaration
        CGIndexOp, // an integer index
        CGIndexAssignOp, // assigment of an integer index to an index pattern expression
        CGLoopStartOp, // for() {}
        CGLoopIndexedIndepOp, // indexed independent used by a loop
        CGLoopIndexedDepOp, // indexed output from a loop
        CGLoopEndOp, // endfor
        CGIndexCondExprOp, // a condition expression return a boolean
        CGStartIfOp, // the start of an if statement
        CGElseIfOp, // else if()
        CGElseOp, // else
        CGEndIfOp, // end of if
        CGCondResultOp // assigment inside an if branch
    };

    inline std::ostream& operator <<(std::ostream& os, const CGOpCode& op) {
        switch (op) {
            case CGAbsOp:
                os << "abs()";
                break;
            case CGAcosOp:
                os << "acos()";
                break;
            case CGAddOp:
                os << "a+b";
                break;
            case CGArrayCreationOp:
                os << "new array[size]";
                break;
            case CGArrayElementOp:
                os << "array[i]";
                break;
            case CGAsinOp:
                os << "asin()";
                break;
            case CGAtanOp:
                os << "atan()";
                break;
            case CGAtomicForwardOp:
                os << "atomicFunction.forward(q, p, vx, vy, tx, ty)";
                break;
            case CGAtomicReverseOp:
                os << "atomicFunction.reverse(p, tx, ty, px, py)";
                break;
            case CGComOpLt:
                os << "result = left < right? trueCase: falseCase";
                break;
            case CGComOpLe:
                os << "result = left <= right? trueCase: falseCase";
                break;
            case CGComOpEq:
                os << "result = left == right? trueCase: falseCase";
                break;
            case CGComOpGe:
                os << "result = left >= right? trueCase: falseCase";
                break;
            case CGComOpGt:
                os << "result = left > right? trueCase: falseCase";
                break;
            case CGComOpNe:
                os << "result = left != right? trueCase: falseCase";
                break;
            case CGCoshOp:
                os << "cosh()";
                break;
            case CGCosOp:
                os << "cos()";
                break;
            case CGDivOp:
                os << "a / b";
                break;
            case CGExpOp:
                os << "exp()";
                break;
            case CGInvOp:
                os << "independent()";
                break;
            case CGLogOp:
                os << "log()";
                break;
            case CGMulOp:
                os << "a * b";
                break;
            case CGPowOp:
                os << "pow(a, b)";
                break;
                //        PriOp: //  PrintFor(text, parameter or variable, parameter or variable)
            case CGSignOp:
                os << "sign()";
                break;
            case CGSinhOp:
                os << "sinh()";
                break;
            case CGSinOp:
                os << "sin()";
                break;
            case CGSqrtOp:
                os << "sqrt()";
                break;
            case CGSubOp:
                os << "a - b";
                break;
            case CGTanhOp:
                os << "tanh()";
                break;
            case CGTanOp:
                os << "tan()";
                break;
            case CGUnMinusOp:
                os << "-(a)";
                break;
            case CGIndexDeclarationOp:
                os << "index declaration";
                break;
            case CGIndexOp:
                os << "index";
                break;
            case CGIndexAssignOp:
                os << "index = expression()";
                break;
            case CGLoopStartOp:
                os << "for";
                break;
            case CGLoopIndexedIndepOp:
                os << "loopIndexedDep";
                break;
            case CGLoopIndexedDepOp:
                os << "loopIndexedIndep";
                break;
            case CGLoopEndOp:
                os << "endfor";
                break;
            case CGIndexCondExprOp:
                os << "";
                break;
            case CGStartIfOp:
                os << "if()";
                break;
            case CGElseIfOp:
                os << "else if()";
                break;
            case CGElseOp:
                os << "else";
                break;
            case CGEndIfOp:
                os << "endif";
                break;
            case CGCondResultOp:
                os << "ifResult =";
                break;
            default:
                os << "\?\?\?\?()";
                break;
        }
        return os;
    }

}

#endif
