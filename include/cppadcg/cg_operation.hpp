#ifndef CPPAD_CG_OPERATION_INCLUDED
#define CPPAD_CG_OPERATION_INCLUDED
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
     * Possible operations
     * 
     * @author Joao Leal
     */
    enum CGOpCode {
        CGAssignOp, // a = b
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
        CGPriOp, //  PrintFor(text, parameter or variable, parameter or variable)
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
        CGIndexAssignOp, // assignment of an integer index to an index pattern expression
        CGLoopStartOp, // for() {}
        CGLoopIndexedIndepOp, // indexed independent used by a loop
        CGLoopIndexedDepOp, // indexed output for a dependent variable from a loop
        CGLoopIndexedTmpOp, // indexed output for a temporary variable from a loop
        CGLoopEndOp, // endfor
        CGTmpDclOp, // marks the beginning of the use of a temporary variable across several scopes (used by CGLoopIndexedTmpOp)
        CGTmpOp, // reference to a temporary variable defined by CGTmpDclOp
        CGIndexCondExprOp, // a condition expression which returns a boolean
        CGStartIfOp, // the start of an if statement
        CGElseIfOp, // else if()
        CGElseOp, // else
        CGEndIfOp, // end of if
        CGCondResultOp // assignment inside an if branch
    };

    inline std::ostream& operator <<(std::ostream& os, const CGOpCode& op) {
        switch (op) {
            case CGAssignOp:
                os << "$1 = $2";
                break;
            case CGAbsOp:
                os << "abs($1)";
                break;
            case CGAcosOp:
                os << "acos($1)";
                break;
            case CGAddOp:
                os << "$1 + $2";
                break;
            case CGAliasOp:
                os << "alias($1)";
                break;
            case CGArrayCreationOp:
                os << "new array[size]";
                break;
            case CGArrayElementOp:
                os << "array[i]";
                break;
            case CGAsinOp:
                os << "asin($1)";
                break;
            case CGAtanOp:
                os << "atan($1)";
                break;
            case CGAtomicForwardOp:
                os << "atomicFunction.forward(q, p, vx, vy, tx, ty)";
                break;
            case CGAtomicReverseOp:
                os << "atomicFunction.reverse(p, tx, ty, px, py)";
                break;
            case CGComOpLt:
                os << "result = ($1 < $2)? $3 : $4";
                break;
            case CGComOpLe:
                os << "result = ($1 <= $2)? $3 : $4";
                break;
            case CGComOpEq:
                os << "result = ($1 == $2)? $3 : $4";
                break;
            case CGComOpGe:
                os << "result = ($1 >= $2)? $3 : $4";
                break;
            case CGComOpGt:
                os << "result = ($1 > $2)? $3 : $4";
                break;
            case CGComOpNe:
                os << "result = ($1 != $2)? $3 : $4";
                break;
            case CGCoshOp:
                os << "cosh($1)";
                break;
            case CGCosOp:
                os << "cos($1)";
                break;
            case CGDivOp:
                os << "$1 / $2";
                break;
            case CGExpOp:
                os << "exp($1)";
                break;
            case CGInvOp:
                os << "independent()";
                break;
            case CGLogOp:
                os << "log($1)";
                break;
            case CGMulOp:
                os << "$1 * $2";
                break;
            case CGPowOp:
                os << "pow($1, $2)";
                break;
            case CGPriOp:
                os << "print($1)";
                break;
            case CGSignOp:
                os << "sign($1)";
                break;
            case CGSinhOp:
                os << "sinh($1)";
                break;
            case CGSinOp:
                os << "sin($1)";
                break;
            case CGSqrtOp:
                os << "sqrt($1)";
                break;
            case CGSubOp:
                os << "$1 - $2";
                break;
            case CGTanhOp:
                os << "tanh($1)";
                break;
            case CGTanOp:
                os << "tan($1)";
                break;
            case CGUnMinusOp:
                os << "-($1)";
                break;
            case CGDependentMultiAssignOp:
                os << "dep($1) = ($2) + ...";
                break;
            case CGDependentRefRhsOp:
                os << "depref($1)";
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
                os << "loopIndexedIndep";
                break;
            case CGLoopIndexedDepOp:
                os << "loopIndexedDep";
                break;
            case CGLoopIndexedTmpOp:
                os << "loopIndexedTmp";
                break;
            case CGLoopEndOp:
                os << "endfor";
                break;
            case CGTmpDclOp:
                os << "declare tempVar";
                break;
            case CGTmpOp:
                os << "tempVar";
                break;
            case CGIndexCondExprOp:
                os << "bool(index expression)";
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
