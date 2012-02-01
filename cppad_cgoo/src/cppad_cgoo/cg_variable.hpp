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

#include <stddef.h>

namespace CppAD {

    template<class Base>
    inline bool CG<Base>::isVariable() const {
        return id_ > 0;
    }

    template<class Base>
    inline bool CG<Base>::isTemporaryVariable() const {
        return !isVariable() && !isParameter();
    }

    template<class Base>
    inline bool CG<Base>::isParameter() const {
        return value_ != NULL;
    }

    template<class Base>
    inline size_t CG<Base>::getVariableID() const {
        return id_;
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
        if (!isParameter()) {
            throw CGException("IdenticalZero() can only be used for parameters");
        }

        return CppAD::IdenticalZero(*value_);
    }

    template<class Base>
    inline bool CG<Base>::IdenticalOne() const throw (CGException) {
        if (!isParameter()) {
            throw CGException("IdenticalOne() can only be used for parameters");
        }

        return CppAD::IdenticalOne(*value_);
    }

    template <class Base>
    inline void CG<Base>::makeParameter(const Base &b) {
        id_ = 0;
        handler_ = NULL;
        operations_.clear();
        opTypes_ = NONE;

        if (value_ == NULL) {
            value_ = new Base(b);
        } else {
            *value_ = b;
        }
    }

    template <class Base>
    inline void CG<Base>::makeVariable(CodeHandler<Base>& handler) {
        operations_.clear();
        opTypes_ = NONE;
        delete value_;
        value_ = NULL;
        handler_ = &handler;
        id_ = handler.createID();
    }

    template <class Base>
    inline void CG<Base>::makeVariable(CodeHandler<Base>& handler, size_t id) {
        operations_.clear();
        opTypes_ = NONE;
        delete value_;
        value_ = NULL;
        handler_ = &handler;
        id_ = id;
    }

    template <class Base>
    inline void CG<Base>::makeTemporaryVariable(CodeHandler<Base>& handler, const std::string& operations, OpContainement op) {
        delete value_;
        value_ = NULL;
        operations_ = operations;
        opTypes_ = op;
        handler_ = &handler;
        id_ = 0;
    }

    template<class Base>
    inline std::string CG<Base>::operations() const throw (CGException) {
        if (isParameter()) {
            throw CGException("Cannot call operations() for a parameter");
        }

        if (isTemporaryVariable()) {
            return operations_;
        } else {
            // named variable
            return createVariableName();
        }
    }

}

#endif

