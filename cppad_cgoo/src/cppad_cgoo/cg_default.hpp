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

#include <stddef.h>

namespace CppAD {

    /**
     * Creates a parameter with a zero value
     */
    template <class Base>
    inline CG<Base>::CG() :
    opTypes_(NONE),
    handler_(NULL),
    value_(new Base(0.0)),
    id_(0),
    referenceTo_(NULL) {
    }

    /**
     * Creates a temporary variable
     */
    template <class Base>
    inline CG<Base>::CG(CodeHandler<Base>& handler, const std::string& ops, OpContainement contain) :
    opTypes_(contain),
    handler_(&handler),
    operations_(ops),
    value_(NULL),
    id_(0),
    referenceTo_(NULL) {
    }

    /**
     * Default copy constructor
     */
    template <class Base>
    inline CG<Base>::CG(const CG<Base>& orig) :
    opTypes_(CppAD::NONE),
    handler_(NULL),
    id_(0),
    referenceTo_(NULL),
    value_(NULL) {
        if (orig.isParameter()) {
            // parameter
            makeParameter(*orig.value_);
        } else if (orig.isVariable()) {
            // make this variable a reference to orig
            makeVariableProxy(orig);
        } else {
            // the other is a temporary variable
            makeTemporaryVariable(*orig.handler_, orig.operations_, orig.opTypes_);
        }
    }

    /**
     * Creates a parameter with the given value
     */
    template <class Base>
    inline CG<Base>::CG(const Base &b) :
    opTypes_(NONE),
    handler_(NULL),
    value_(NULL),
    id_(0),
    referenceTo_(NULL) {
        // make it a parameter
        makeParameter(b);
    }

    template <class Base>
    inline CG<Base>& CG<Base>::operator =(const CG<Base> &other) {
        if (this == &other) {
            return *this;
        }

        if (isVariable() && other.getCodeHandler() == handler_) {
            variableValueWillChange();
        }

        if (other.isParameter()) {
            // parameter
            makeParameter(*other.value_);
        } else if (other.isVariable()) {
            // make this variable a reference to other
            makeVariableProxy(other);

        } else {
            // the other is a temporary variable
            if (!isVariable() || other.getCodeHandler() != getCodeHandler()) {
                // this can be used to switch between source code handlers
                makeVariable(*other.getCodeHandler());
            }

            // variable :: print operation
            std::string name = createVariableName();
            printOperationAssig(name, handler_->operations(other));
        }

        return *this;
    }

    template <class Base>
    inline CG<Base>& CG<Base>::operator=(const Base &b) {
        if (isVariable()) {
            variableValueWillChange();
        }

        // make it a parameter
        makeParameter(b);

        return *this;
    }

    template <class Base>
    CG<Base>::~CG() {
        if (isVariable()) {
            if (isReference()) {
                handler_->removeVariableReference(*this);
            } else {
                handler_->removePureVariable(*this);
            }
        }

        delete value_;
        value_ = NULL; // not really required
    }

}

#endif

