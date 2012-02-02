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

#define CPPAD_CG_CONSTRUCTION_CHECK() \
           deleted_ = false;\
           CPPAD_CG_CHECK_CG_THIS();

namespace CppAD {

    template <class Base>
    inline CG<Base>::CG() :
    opTypes_(NONE),
    handler_(NULL),
    value_(new Base(0.0)),
    id_(0) {
        CPPAD_CG_CONSTRUCTION_CHECK();
    }

    template <class Base>
    inline CG<Base>::CG(CodeHandler<Base>& handler, const std::string& ops, OpContainement contain) :
    opTypes_(contain),
    handler_(&handler),
    operations_(ops),
    value_(NULL),
    id_(0) {
        CPPAD_CG_CONSTRUCTION_CHECK();
    }

    template <class Base>
    inline CG<Base>::CG(const CG<Base>& orig) :
    opTypes_(orig.opTypes_),
    handler_(orig.handler_),
    operations_(orig.operations_),
    id_(orig.id_) {
        if (orig.value_ != NULL) {
            value_ = new Base(*orig.value_);
        } else {
            value_ = NULL;
        }

        CPPAD_CG_CONSTRUCTION_CHECK();
    }

    template <class Base>
    inline CG<Base>::CG(const Base &b) :
    opTypes_(NONE),
    handler_(NULL),
    value_(NULL),
    id_(0) {
        // make it a parameter
        makeParameter(b);

        CPPAD_CG_CONSTRUCTION_CHECK();
    }

    template <class Base>
    inline CG<Base>& CG<Base>::operator =(const CG<Base> &other) {
        if (this == &other) {
            return *this;
        }

        if (other.isParameter()) {
            // parameter
            makeParameter(*other.value_);
        } else {
            if (other.isVariable()) {
                // variable
                makeVariable(*other.getCodeHandler(), other.id_);

            } else {
                // variable :: print operation
                //get id and name
                makeVariable(*other.getCodeHandler());

                std::string name = createVariableName();
                printOperationAssig(name, other.operations());
            }
        }

        CPPAD_CG_CONSTRUCTION_CHECK();
        return *this;
    }

    template <class Base>
    inline CG<Base>& CG<Base>::operator=(const Base &b) {
        // make it a parameter
        makeParameter(b);

        CPPAD_CG_CONSTRUCTION_CHECK();
        return *this;
    }

    template <class Base>
    CG<Base>::~CG() {
        delete value_;
        value_ = NULL; // not really required
        deleted_ = true;
    }

}

#endif

