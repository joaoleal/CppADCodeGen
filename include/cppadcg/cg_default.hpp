#ifndef CPPAD_CG_DEFAULT_INCLUDED
#define CPPAD_CG_DEFAULT_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

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

