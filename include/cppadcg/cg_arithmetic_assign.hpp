#ifndef CPPAD_CG_ARITHMETIC_ASSIGN_INCLUDED
#define	CPPAD_CG_ARITHMETIC_ASSIGN_INCLUDED
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
    inline CG<Base>& CG<Base>::operator+=(const CG<Base> &right) {
        if (isParameter() && right.isParameter()) {
            *value_ += *right.value_;

        } else {
            CodeHandler<Base>* handler;
            if (isParameter()) {
                if (IdenticalZero()) {
                    *this = right;
                    return *this;
                }

                handler = right.getCodeHandler();

            } else if (right.isParameter()) {
                if (right.IdenticalZero()) {
                    return *this; // nothing to do
                }

                handler = getCodeHandler();

            } else {
                assert(getCodeHandler() == right.getCodeHandler());
                handler = getCodeHandler();
            }

            makeVariable(*handler, new SourceCodeFragment<Base>(CGAddOp, argument(), right.argument()));
        }

        return *this;
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator-=(const CG<Base> &right) {
        if (isParameter() && right.isParameter()) {
            *value_ -= *right.value_;

        } else {
            CodeHandler<Base>* handler;
            if (isParameter()) {
                handler = right.getCodeHandler();

            } else if (right.isParameter()) {
                if (right.IdenticalZero()) {
                    return *this; // nothing to do
                }

                handler = getCodeHandler();

            } else {
                assert(getCodeHandler() == right.getCodeHandler());
                handler = getCodeHandler();
            }

            makeVariable(*handler, new SourceCodeFragment<Base>(CGSubOp, argument(), right.argument()));
        }

        return *this;
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator*=(const CG<Base> &right) {
        if (isParameter() && right.isParameter()) {
            *value_ *= *right.value_;

        } else {
            CodeHandler<Base>* handler;
            if (isParameter()) {
                if (IdenticalZero()) {
                    return *this; // nothing to do (does not consider that right might be infinity)
                } else if (IdenticalOne()) {
                    *this = right;
                    return *this;
                }

                handler = right.getCodeHandler();

            } else if (right.isParameter()) {
                if (right.IdenticalZero()) {
                    makeParameter(Base(0.0)); // does not consider that left might be infinity
                    return *this;
                } else if (right.IdenticalOne()) {
                    return *this; // nothing to do
                }

                handler = getCodeHandler();

            } else {
                assert(getCodeHandler() == right.getCodeHandler());
                handler = getCodeHandler();
            }

            makeVariable(*handler, new SourceCodeFragment<Base>(CGMulOp, argument(), right.argument()));
        }

        return *this;
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator/=(const CG<Base> &right) {
        if (isParameter() && right.isParameter()) {
            *value_ /= *right.value_;

        } else {
            CodeHandler<Base>* handler;
            if (isParameter()) {
                if (IdenticalZero()) {
                    return *this; // nothing to do (does not consider that right might be infinity or zero)
                }

                handler = right.getCodeHandler();

            } else if (right.isParameter()) {
                if (right.IdenticalOne()) {
                    return *this; // nothing to do
                }

                handler = getCodeHandler();

            } else {
                assert(getCodeHandler() == right.getCodeHandler());
                handler = getCodeHandler();
            }

            makeVariable(*handler, new SourceCodeFragment<Base>(CGDivOp, argument(), right.argument()));
        }

        return *this;
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator+=(const Base &right) {
        return operator+=(CG<Base > (right));
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator-=(const Base &right) {
        return operator-=(CG<Base > (right));
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator/=(const Base &right) {
        return operator/=(CG<Base > (right));
    }

    template<class Base>
    inline CG<Base>& CG<Base>::operator*=(const Base &right) {
        return operator*=(CG<Base > (right));
    }

    template<class Base>
    template<class T>
    inline CG<Base>& CG<Base>::operator+=(const T &right) {
        return operator+=(CG<Base > (right));
    }

    template<class Base>
    template<class T>
    inline CG<Base>& CG<Base>::operator-=(const T &right) {
        return operator-=(CG<Base > (right));
    }

    template<class Base>
    template<class T>
    inline CG<Base>& CG<Base>::operator/=(const T &right) {
        return operator/=(CG<Base > (right));
    }

    template<class Base>
    template<class T>
    inline CG<Base>& CG<Base>::operator*=(const T &right) {
        return operator*=(CG<Base > (right));
    }

}

#endif

