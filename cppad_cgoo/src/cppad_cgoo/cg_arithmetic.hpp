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
    CodeHandler<Base>* getOperations(const CG<Base> &left,
                                     const CG<Base> &right,
                                     std::string& leftOps,
                                     std::string& rightOps) {

        assert(!left.isParameter() || !right.isParameter());

        CodeHandler<Base>* handler;
        if (left.isParameter()) {
            handler = right.getCodeHandler();
            leftOps = handler->baseToString(left.getParameterValue());
            rightOps = right.operations();
        } else if (right.isParameter()) {
            handler = left.getCodeHandler();
            leftOps = left.operations();
            rightOps = handler->baseToString(right.getParameterValue());
        } else {
            if (left.getCodeHandler() != right.getCodeHandler()) {
                throw CGException("Attempting to use several source code generation handlers in the same source code generation");
            }
            handler = left.getCodeHandler();
            leftOps = left.operations();
            rightOps = right.operations();
        }
        return handler;
    }

    template<class Base>
    inline CG<Base> operator+(const CG<Base> &left, const CG<Base> &right) {

        CPPAD_CG_CHECK_CG(left);
        CPPAD_CG_CHECK_CG(right);

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

            std::string leftOps;
            std::string rightOps;
            CodeHandler<Base>* handler = getOperations(left, right, leftOps, rightOps);

            std::string operations = leftOps + " + " + rightOps;

            CG<Base> result;
            result.makeTemporaryVariable(*handler, operations, PLUS_MINUS_BINARY);
            return result;
        }
    }

    template<class Base>
    inline CG<Base> operator-(const CG<Base> &left, const CG<Base> &right) {

        CPPAD_CG_CHECK_CG(left);
        CPPAD_CG_CHECK_CG(right);

        if (left.isParameter() && right.isParameter()) {
            return CG<Base > (left.getParameterValue() - right.getParameterValue());

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
            std::string leftOps;
            std::string rightOps;
            CodeHandler<Base>* handler = getOperations(left, right, leftOps, rightOps);

            std::string operations = leftOps + " - ";
            if (right.opTypes_ == NONE || right.opTypes_ == MULT_DIV || right.opTypes_ == FUNCTION) {
                operations += rightOps;
            } else {
                operations += "(" + rightOps + ")";
            }

            CG<Base> result;
            result.makeTemporaryVariable(*handler, operations, PLUS_MINUS_BINARY);
            return result;
        }
    }

    template<class Base>
    inline CG<Base> operator*(const CG<Base> &left, const CG<Base> &right) {

        CPPAD_CG_CHECK_CG(left);
        CPPAD_CG_CHECK_CG(right);

        if (left.isParameter() && right.isParameter()) {
            return CG<Base > (left.getParameterValue() * right.getParameterValue());

        } else {
            if (left.isParameter()) {
                if (left.IdenticalZero()) {
                    return CG<Base > (Base(0.0));
                } else if (left.IdenticalOne()) {
                    return right;
                }
            } else if (right.isParameter()) {
                if (right.IdenticalZero()) {
                    return CG<Base > (Base(0.0));
                } else if (right.IdenticalOne()) {
                    return left;
                }
            }
            std::string leftOps;
            std::string rightOps;
            CodeHandler<Base>* handler = getOperations(left, right, leftOps, rightOps);

            std::string operations;
            if (left.opTypes_ == NONE || left.opTypes_ == MULT_DIV || left.opTypes_ == FUNCTION) {
                operations += leftOps;
            } else {
                operations += "(" + leftOps + ")";
            }

            operations += " * ";

            if (right.opTypes_ == NONE || right.opTypes_ == MULT_DIV || right.opTypes_ == FUNCTION) {
                operations += rightOps;
            } else {
                operations += "(" + rightOps + ")";
            }

            CG<Base> result;
            result.makeTemporaryVariable(*handler, operations, MULT_DIV);
            return result;
        }
    }

    template<class Base>
    inline CG<Base> operator/(const CG<Base> &left, const CG<Base> &right) {

        CPPAD_CG_CHECK_CG(left);
        CPPAD_CG_CHECK_CG(right);

        if (left.isParameter() && right.isParameter()) {
            return CG<Base > (left.getParameterValue() / right.getParameterValue());

        } else {
            if (right.isParameter()) {
                if (right.IdenticalOne()) {
                    return left;
                }
            }
            std::string leftOps;
            std::string rightOps;
            CodeHandler<Base>* handler = getOperations(left, right, leftOps, rightOps);

            std::string operations;
            if (left.opTypes_ == NONE || left.opTypes_ == FUNCTION) {
                operations += leftOps;
            } else {
                operations += "(" + leftOps + ")";
            }

            operations += " / ";

            if (right.opTypes_ == NONE || right.opTypes_ == FUNCTION) {
                operations += rightOps;
            } else {
                operations += "(" + rightOps + ")";
            }

            CG<Base> result;
            result.makeTemporaryVariable(*handler, operations, MULT_DIV);
            return result;
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

//    template<class Base, class T>
//    inline CG<Base> operator+(const T &left, const CG<Base> &right) {
//        return CG<Base > (Base(left)) + right;
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator+(const CG<Base> &left, const T &right) {
//        return left + CG<Base > (Base(right));
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator-(const T &left, const CG<Base> &right) {
//        return CG<Base > (Base(left)) - right;
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator-(const CG<Base> &left, const T &right) {
//        return left - CG<Base > (Base(right));
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator/(const Base &left, const CG<Base> &right) {
//        return CG<Base > (Base(left)) / right;
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator/(const CG<Base> &left, const Base &right) {
//        return left / CG<Base > (Base(right));
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator*(const Base &left, const CG<Base> &right) {
//        return CG<Base > (Base(left)) * right;
//    }
//
//    template<class Base, class T>
//    inline CG<Base> operator*(const CG<Base> &left, const Base &right) {
//        return left * CG<Base > (Base(right));
//    }

}

#endif

