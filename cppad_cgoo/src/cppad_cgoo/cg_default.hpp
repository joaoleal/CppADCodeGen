#ifndef CPPAD_CG_DEFAULT_INCLUDED
#define	CPPAD_CG_DEFAULT_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    /**
     * Creates a parameter with a zero value
     */
    template <class Base>
    inline CG<Base>::CG() :
    handler_(NULL),
    sourceCode_(NULL),
    value_(new Base(0.0)) {
    }

    template <class Base>
    inline CG<Base>::CG(CodeHandler<Base>& handler, SourceCodeFragment<Base>* sourceCode) :
    handler_(&handler),
    sourceCode_(sourceCode),
    value_(NULL) {
        assert(sourceCode != NULL);
        handler.manageSourceCodeBlock(sourceCode);
    }

    template <class Base>
    inline CG<Base>::CG(CodeHandler<Base>& handler, const Argument<Base>& arg) :
    handler_(&handler),
    sourceCode_(arg.operation()),
    value_(arg.parameter() != NULL ? new Base(*arg.parameter()) : NULL) {

    }

    /**
     * Creates a parameter with the given value
     */
    template <class Base>
    inline CG<Base>::CG(const Base &b) :
    handler_(NULL),
    sourceCode_(NULL),
    value_(new Base(b)) {
    }

    /**
     * Copy constructor
     */
    template <class Base>
    inline CG<Base>::CG(const CG<Base>& orig) :
    handler_(orig.handler_),
    sourceCode_(orig.sourceCode_),
    value_(orig.value_ != NULL ? new Base(*orig.value_) : NULL) {
    }

    /**
     * Creates a parameter with the given value
     */
    template <class Base>
    inline CG<Base>& CG<Base>::operator=(const Base &b) {
        handler_ = NULL;
        sourceCode_ = NULL;
        if (value_ != NULL) {
            *value_ = b;
        } else {
            value_ = new Base(b);
        }
        return *this;
    }

    template <class Base>
    inline CG<Base>& CG<Base>::operator=(const CG<Base> &rhs) {
        if (&rhs == this) {
            return *this;
        }
        handler_ = rhs.handler_;
        sourceCode_ = rhs.sourceCode_;
        if (rhs.value_ != NULL) {
            if (value_ != NULL) {
                *value_ = *rhs.value_;
            } else {
                value_ = new Base(*rhs.value_);
            }
        } else {
            delete value_;
            value_ = NULL;
        }

        return *this;
    }

    template <class Base>
    CG<Base>::~CG() {
        delete value_;
    }

}

#endif

