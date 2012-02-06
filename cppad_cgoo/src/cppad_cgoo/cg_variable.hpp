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
    inline bool CG<Base>::isReference() const {
        return referenceTo_ != NULL;
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
        if (isVariable()) {
            if (isReference()) {
                handler_->removeVariableReference(*this);
            } else {
                handler_->removePureVariable(*this);
            }
        }

        makeParameterNoChecks(b);
    }

    template <class Base>
    inline void CG<Base>::makeParameterNoChecks(const Base &b) {
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
        bool wasReference = isReference();
        bool isNewHandler = &handler != handler_;

        if (!isNewHandler && wasReference) {
            handler_->removeVariableReference(*this);
        }

        operations_.clear();
        opTypes_ = NONE;
        delete value_;
        value_ = NULL;
        handler_ = &handler;
        id_ = handler.createID(); // generate a new variable ID (even if it was already a variable)
        referenceTo_ = NULL;

        if (isNewHandler || wasReference) {
            handler_->addPureVariable(*this);
        }
    }

    template <class Base>
    inline void CG<Base>::makeVariableProxy(const CG<Base>& referenceTo) {
        assert(referenceTo.isVariable());
        bool wasVariable = isVariable();
        bool wasReference = isReference();

        if (handler_ == referenceTo.getCodeHandler()) {
            if (wasReference) {
                handler_->removeVariableReference(*this);
            } else if (wasVariable) {
                handler_->removePureVariable(*this);
            }
        }

        operations_.clear();
        opTypes_ = NONE;
        delete value_;
        value_ = NULL;

        if (referenceTo.isReference()) {
            referenceTo_ = referenceTo.referenceTo_; // there are no 2nd level proxies
        } else {
            referenceTo_ = &referenceTo;
        }
        assert(referenceTo_ != this);

        id_ = referenceTo_->id_;
        handler_ = referenceTo_->handler_;
        handler_->putVariableReference(referenceTo, *this);
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

    template <class Base>
    inline void CG<Base>::variableValueWillChange() {
        assert(isVariable());

        if (isReference()) {
            const CG<Base>* oldRef = referenceTo_;
            makeVariable(*handler_); // assign it a new variable ID

            // variable :: print operation
            std::string name = createVariableName();
            printOperationAssig(name, handler_->operations(*oldRef));

        } else {
            handler_->releaseReferences(*this);
        }
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

