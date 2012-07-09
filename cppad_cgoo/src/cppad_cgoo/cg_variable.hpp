#ifndef CPPAD_CG_VARIABLE_INCLUDED
#define	CPPAD_CG_VARIABLE_INCLUDED
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
    inline bool CG<Base>::isVariable() const {
        return sourceCode_ != NULL;
    }

    template<class Base>
    inline bool CG<Base>::isParameter() const {
        return value_ != NULL;
    }

    template<class Base>
    inline const Base& CG<Base>::getParameterValue() const throw (CGException) {
        if (!isParameter()) {
            throw CGException("getParameterValue() can only be used for parameters");
        }

        return *value_;
    }

    template<class Base>
    inline bool CG<Base>::IdenticalZero() const throw (CGException) {
        return CppAD::IdenticalZero(getParameterValue());
    }

    template<class Base>
    inline bool CG<Base>::IdenticalOne() const throw (CGException) {
        return CppAD::IdenticalOne(getParameterValue());
    }

    template<class Base>
    inline void CG<Base>::makeParameter(const Base &b) {
        sourceCode_ = NULL;
        handler_ = NULL;
        if (value_ != NULL) {
            *value_ = b;
        } else {
            value_ = new Base(b);
        }
    }

    template<class Base>
    inline void CG<Base>::makeVariable(CodeHandler<Base>& handler, SourceCodeFragment<Base>* operation) {
        assert(operation != NULL);
        sourceCode_ = operation;
        handler_ = &handler;
        delete value_;
        value_ = NULL;
        
        handler.manageSourceCodeBlock(operation);
    }

    template<class Base>
    inline SourceCodeFragment<Base>* CG<Base>::getSourceCodeFragment() const {
        return sourceCode_;
    }

    template<class Base>
    inline Argument<Base> CG<Base>::argument() const {
        if (sourceCode_ != NULL)
            return Argument<Base > (*sourceCode_);
        else
            return Argument<Base > (*value_);
    }
}

#endif

