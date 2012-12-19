#ifndef CPPAD_CG_ARITHMETIC_INCLUDED
#define	CPPAD_CG_ARITHMETIC_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<class Base>
    CodeHandler<Base>* getHandler(const CG<Base> &left,
                                  const CG<Base> &right) {

        assert(!left.isParameter() || !right.isParameter());

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
            return CG<Base > (left.getParameterValue() + right.getParameterValue());

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

            return CG<Base > (*handler, new SourceCodeFragment<Base > (CGAddOp, left.argument(), right.argument()));
        }
    }

    template<class Base>
    inline CG<Base> operator-(const CG<Base> &left, const CG<Base> &right) {
        if (left.isParameter() && right.isParameter()) {
            return CG<Base > (left.getParameterValue() - right.getParameterValue());

        } else {
            if (right.isParameter()) {
                if (right.IdenticalZero()) {
                    return left;
                }
            }

            CodeHandler<Base>* handler = getHandler(left, right);

            return CG<Base > (*handler, new SourceCodeFragment<Base > (CGSubOp, left.argument(), right.argument()));
        }
    }

    template<class Base>
    inline CG<Base> operator*(const CG<Base> &left, const CG<Base> &right) {
        if (left.isParameter() && right.isParameter()) {
            return CG<Base > (left.getParameterValue() * right.getParameterValue());

        } else {
            if (left.isParameter()) {
                if (left.IdenticalZero()) {
                    return CG<Base > (Base(0.0)); // does not consider the possibility of right being infinity
                } else if (left.IdenticalOne()) {
                    return right;
                }
            } else if (right.isParameter()) {
                if (right.IdenticalZero()) {
                    return CG<Base > (Base(0.0)); // does not consider the possibility of left being infinity
                } else if (right.IdenticalOne()) {
                    return left;
                }
            }

            CodeHandler<Base>* handler = getHandler(left, right);

            return CG<Base > (*handler, new SourceCodeFragment<Base > (CGMulOp, left.argument(), right.argument()));
        }
    }

    template<class Base>
    inline CG<Base> operator/(const CG<Base> &left, const CG<Base> &right) {
        if (left.isParameter() && right.isParameter()) {
            return CG<Base > (left.getParameterValue() / right.getParameterValue());

        } else {
            if (left.isParameter()) {
                if (left.IdenticalZero()) {
                    return CG<Base > (Base(0.0)); // does not consider the possibility of right being infinity or zero
                }
            } else if (right.isParameter()) {
                if (right.IdenticalOne()) {
                    return left;
                }
            }

            CodeHandler<Base>* handler = getHandler(left, right);

            return CG<Base > (*handler, new SourceCodeFragment<Base > (CGDivOp, left.argument(), right.argument()));
        }
    }

    template<class Base>
    inline CG<Base> operator+(const Base &left, const CG<Base> &right) {
        return CG<Base > (left) + right;
    }

    template<class Base>
    inline CG<Base> operator+(const CG<Base> &left, const Base &right) {
        return left + CG<Base > (right);
    }

    template<class Base>
    inline CG<Base> operator-(const Base &left, const CG<Base> &right) {
        return CG<Base > (left) - right;
    }

    template<class Base>
    inline CG<Base> operator-(const CG<Base> &left, const Base &right) {
        return left - CG<Base > (right);
    }

    template<class Base>
    inline CG<Base> operator/(const Base &left, const CG<Base> &right) {
        return CG<Base > (left) / right;
    }

    template<class Base>
    inline CG<Base> operator/(const CG<Base> &left, const Base &right) {
        return left / CG<Base > (right);
    }

    template<class Base>
    inline CG<Base> operator*(const Base &left, const CG<Base> &right) {
        return CG<Base > (left) * right;
    }

    template<class Base>
    inline CG<Base> operator*(const CG<Base> &left, const Base &right) {
        return left * CG<Base > (right);
    }

}

#endif

