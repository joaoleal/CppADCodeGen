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
namespace cg {

/**
 * Possible operations
 * 
 * @author Joao Leal
 */
enum class CGOpCode {
    Assign, // a = b
    Abs, //  abs(variable)
    Acos, // acos(variable)
    Acosh, // acosh(variable)
    Add, //  a + b
    Alias, //  alias (reference to another operation)
    ArrayCreation, // dense array {a, b, c ...}
    SparseArrayCreation, // {a, b, c ...}; {index1, index2, index3, ...};
    ArrayElement, // x[i]
    Asin, // asin(variable)
    Asinh, // asinh(variable)
    Atan, // atan(variable)
    Atanh, // atanh(variable)
    AtomicForward, // atomicFunction.forward(q, p, vx, vy, tx, ty)
    AtomicReverse, // atomicFunction.reverse(p, tx, ty, px, py)
    ComLt, // result = left < right? trueCase: falseCase
    ComLe, // result = left <= right? trueCase: falseCase
    ComEq, // result = left == right? trueCase: falseCase
    ComGe, // result = left >= right? trueCase: falseCase
    ComGt, // result = left > right? trueCase: falseCase
    ComNe, // result = left != right? trueCase: falseCase
    Cosh, // cosh(variable)
    Cos, //  cos(variable)
    Div, // a / b
    Erf, // erf(variable)
    Exp, //  exp(variable)
    Expm1, //  expm1(variable)
    Inv, //                             independent variable
    Log, //  log(variable)
    Log1p, //  log1p(variable)
    Mul, // a * b
    Pow, //  pow(a,   b)
    Pri, //  PrintFor(text, parameter or variable, parameter or variable)
    Sign, // result = (x > 0)? 1.0:((x == 0)? 0.0:-1)
    Sinh, // sinh(variable)
    Sin, //  sin(variable)
    Sqrt, // sqrt(variable)
    Sub, //  a - b
    Tanh, //  tanh(variable)
    Tan, //  tan(variable)
    UnMinus, // -(a)
    DependentMultiAssign, // operation which associates a dependent variables with loops and regular operations
    DependentRefRhs, // operation referencing a dependent variable (right hand side only)
    IndexDeclaration, // an integer index declaration
    Index, // an integer index
    IndexAssign, // assignment of an integer index to an index pattern expression
    LoopStart, // for() {}
    LoopIndexedIndep, // indexed independent used by a loop
    LoopIndexedDep, // indexed output for a dependent variable from a loop
    LoopIndexedTmp, // indexed output for a temporary variable from a loop
    LoopEnd, // endfor
    TmpDcl, // marks the beginning of the use of a temporary variable across several scopes (used by LoopIndexedTmp)
    Tmp, // reference to a temporary variable defined by TmpDcl
    IndexCondExpr, // a condition expression which returns a boolean
    StartIf, // the start of an if statement
    ElseIf, // else if()
    Else, // else
    EndIf, // end of if
    CondResult // assignment inside an if branch
};

inline std::ostream& operator<<(std::ostream& os, const CGOpCode& op) {
    switch (op) {
        case CGOpCode::Assign:
            os << "$1 = $2";
            break;
        case CGOpCode::Abs:
            os << "abs($1)";
            break;
        case CGOpCode::Acos:
            os << "acos($1)";
            break;
        case CGOpCode::Acosh:
            os << "acosh($1)";
            break;
        case CGOpCode::Add:
            os << "$1 + $2";
            break;
        case CGOpCode::Alias:
            os << "alias($1)";
            break;
        case CGOpCode::ArrayCreation:
            os << "new array[size]";
            break;
        case CGOpCode::SparseArrayCreation:
            os << "new sparseArray[size]";
            break;
        case CGOpCode::ArrayElement:
            os << "array[i]";
            break;
        case CGOpCode::Asin:
            os << "asin($1)";
            break;
        case CGOpCode::Asinh:
            os << "asinh($1)";
            break;
        case CGOpCode::Atan:
            os << "atan($1)";
            break;
        case CGOpCode::Atanh:
            os << "atanh($1)";
            break;
        case CGOpCode::AtomicForward:
            os << "atomicFunction.forward(q, p, vx, vy, tx, ty)";
            break;
        case CGOpCode::AtomicReverse:
            os << "atomicFunction.reverse(p, tx, ty, px, py)";
            break;
        case CGOpCode::ComLt:
            os << "result = ($1 < $2)? $3 : $4";
            break;
        case CGOpCode::ComLe:
            os << "result = ($1 <= $2)? $3 : $4";
            break;
        case CGOpCode::ComEq:
            os << "result = ($1 == $2)? $3 : $4";
            break;
        case CGOpCode::ComGe:
            os << "result = ($1 >= $2)? $3 : $4";
            break;
        case CGOpCode::ComGt:
            os << "result = ($1 > $2)? $3 : $4";
            break;
        case CGOpCode::ComNe:
            os << "result = ($1 != $2)? $3 : $4";
            break;
        case CGOpCode::Cosh:
            os << "cosh($1)";
            break;
        case CGOpCode::Cos:
            os << "cos($1)";
            break;
        case CGOpCode::Div:
            os << "$1 / $2";
            break;
        case CGOpCode::Erf:
            os << "erf($1)";
            break;
        case CGOpCode::Exp:
            os << "exp($1)";
            break;
        case CGOpCode::Expm1:
            os << "expm1($1)";
            break;
        case CGOpCode::Inv:
            os << "independent()";
            break;
        case CGOpCode::Log:
            os << "log($1)";
            break;
        case CGOpCode::Log1p:
            os << "log1p($1)";
            break;
        case CGOpCode::Mul:
            os << "$1 * $2";
            break;
        case CGOpCode::Pow:
            os << "pow($1, $2)";
            break;
        case CGOpCode::Pri:
            os << "print($1)";
            break;
        case CGOpCode::Sign:
            os << "sign($1)";
            break;
        case CGOpCode::Sinh:
            os << "sinh($1)";
            break;
        case CGOpCode::Sin:
            os << "sin($1)";
            break;
        case CGOpCode::Sqrt:
            os << "sqrt($1)";
            break;
        case CGOpCode::Sub:
            os << "$1 - $2";
            break;
        case CGOpCode::Tanh:
            os << "tanh($1)";
            break;
        case CGOpCode::Tan:
            os << "tan($1)";
            break;
        case CGOpCode::UnMinus:
            os << "-($1)";
            break;
        case CGOpCode::DependentMultiAssign:
            os << "dep($1) = ($2) + ...";
            break;
        case CGOpCode::DependentRefRhs:
            os << "depref($1)";
            break;
        case CGOpCode::IndexDeclaration:
            os << "index declaration";
            break;
        case CGOpCode::Index:
            os << "index";
            break;
        case CGOpCode::IndexAssign:
            os << "index = expression()";
            break;
        case CGOpCode::LoopStart:
            os << "for";
            break;
        case CGOpCode::LoopIndexedIndep:
            os << "loopIndexedIndep";
            break;
        case CGOpCode::LoopIndexedDep:
            os << "loopIndexedDep";
            break;
        case CGOpCode::LoopIndexedTmp:
            os << "loopIndexedTmp";
            break;
        case CGOpCode::LoopEnd:
            os << "endfor";
            break;
        case CGOpCode::TmpDcl:
            os << "declare tempVar";
            break;
        case CGOpCode::Tmp:
            os << "tempVar";
            break;
        case CGOpCode::IndexCondExpr:
            os << "bool(index expression)";
            break;
        case CGOpCode::StartIf:
            os << "if()";
            break;
        case CGOpCode::ElseIf:
            os << "else if()";
            break;
        case CGOpCode::Else:
            os << "else";
            break;
        case CGOpCode::EndIf:
            os << "endif";
            break;
        case CGOpCode::CondResult:
            os << "ifResult =";
            break;
        default:
            os << "\?\?\?\?()";
            break;
    }
    return os;
}

} // END cg namespace
} // END CppAD namespace

#endif