#ifndef CPPAD_CG_ARITHMETIC_INCLUDED
#define CPPAD_CG_ARITHMETIC_INCLUDED
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

template<class Base>
CodeHandler<Base>* getHandler(const CG<Base> &left,
                              const CG<Base> &right) throw (CGException) {

    CPPADCG_ASSERT_UNKNOWN(!left.isParameter() || !right.isParameter());

    CodeHandler<Base>* handler;
    if (left.isParameter()) {
        handler = right.getCodeHandler();
    } else if (right.isParameter()) {
        handler = left.getCodeHandler();
    } else {
        if (left.getCodeHandler() != right.getCodeHandler()) {
            throw CGException("Attempting to use several source code generation handlers in the same source code generation");
        }
        handler = left.getCodeHandler();
    }
    return handler;
}

template<class Base>
inline CG<Base> operator+(const CG<Base> &left, const CG<Base> &right) {
    if (left.isParameter() && right.isParameter()) {
        return CG<Base> (left.getValue() + right.getValue());

    } else {
        if (left.isParameter()) {
            if (left.IdenticalZero()) {
                return right;
            }
        } else if (right.isParameter()) {
            if (right.IdenticalZero()) {
                return left;
            }
        }

        CodeHandler<Base>* handler = getHandler(left, right);

        CG<Base> result(*handler, new OperationNode<Base > (CGAddOp,{left.argument(), right.argument()}));
        if (left.isValueDefined() && right.isValueDefined()) {
            result.setValue(left.getValue() + right.getValue());
        }
        return result;
    }
}

template<class Base>
inline CG<Base> operator-(const CG<Base> &left, const CG<Base> &right) {
    if (left.isParameter() && right.isParameter()) {
        return CG<Base> (left.getValue() - right.getValue());

    } else {
        if (right.isParameter()) {
            if (right.IdenticalZero()) {
                return left;
            }
        }

        CodeHandler<Base>* handler = getHandler(left, right);

        CG<Base> result(*handler, new OperationNode<Base> (CGSubOp,{left.argument(), right.argument()}));
        if (left.isValueDefined() && right.isValueDefined()) {
            result.setValue(left.getValue() - right.getValue());
        }
        return result;
    }
}

template<class Base>
inline CG<Base> operator*(const CG<Base> &left, const CG<Base> &right) {
    if (left.isParameter() && right.isParameter()) {
        return CG<Base> (left.getValue() * right.getValue());

    } else {
        if (left.isParameter()) {
            if (left.IdenticalZero()) {
                return CG<Base> (Base(0.0)); // does not consider the possibility of right being infinity
            } else if (left.IdenticalOne()) {
                return right;
            }
        } else if (right.isParameter()) {
            if (right.IdenticalZero()) {
                return CG<Base> (Base(0.0)); // does not consider the possibility of left being infinity
            } else if (right.IdenticalOne()) {
                return left;
            }
        }

        CodeHandler<Base>* handler = getHandler(left, right);

        CG<Base> result(*handler, new OperationNode<Base> (CGMulOp,{left.argument(), right.argument()}));
        if (left.isValueDefined() && right.isValueDefined()) {
            result.setValue(left.getValue() * right.getValue());
        }
        return result;
    }
}

template<class Base>
inline CG<Base> operator/(const CG<Base> &left, const CG<Base> &right) {
    if (left.isParameter() && right.isParameter()) {
        return CG<Base> (left.getValue() / right.getValue());

    } else {
        if (left.isParameter()) {
            if (left.IdenticalZero()) {
                return CG<Base> (Base(0.0)); // does not consider the possibility of right being infinity or zero
            }
        } else if (right.isParameter()) {
            if (right.IdenticalOne()) {
                return left;
            }
        }

        CodeHandler<Base>* handler = getHandler(left, right);

        CG<Base> result(*handler, new OperationNode<Base> (CGDivOp,{left.argument(), right.argument()}));
        if (left.isValueDefined() && right.isValueDefined()) {
            result.setValue(left.getValue() / right.getValue());
        }
        return result;
    }
}

template<class Base>
inline CG<Base> operator+(const Base &left, const CG<Base> &right) {
    return CG<Base> (left) + right;
}

template<class Base>
inline CG<Base> operator+(const CG<Base> &left, const Base &right) {
    return left + CG<Base > (right);
}

template<class Base>
inline CG<Base> operator-(const Base &left, const CG<Base> &right) {
    return CG<Base> (left) - right;
}

template<class Base>
inline CG<Base> operator-(const CG<Base> &left, const Base &right) {
    return left - CG<Base> (right);
}

template<class Base>
inline CG<Base> operator/(const Base &left, const CG<Base> &right) {
    return CG<Base> (left) / right;
}

template<class Base>
inline CG<Base> operator/(const CG<Base> &left, const Base &right) {
    return left / CG<Base> (right);
}

template<class Base>
inline CG<Base> operator*(const Base &left, const CG<Base> &right) {
    return CG<Base> (left) * right;
}

template<class Base>
inline CG<Base> operator*(const CG<Base> &left, const Base &right) {
    return left * CG<Base> (right);
}

} // END cg namespace
} // END CppAD namespace

#endif

