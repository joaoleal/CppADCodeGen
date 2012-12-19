#ifndef CPPAD_CG_ARGUMENT_INCLUDED
#define	CPPAD_CG_ARGUMENT_INCLUDED
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
    class SourceCodeFragment;

    /**
     * An argument used by an operation wich can be either a constant value
     * or the result of another operation
     * 
     * \author Joao Leal
     */
    template<class Base>
    class Argument {
    private:
        SourceCodeFragment<Base>* operation_;
        Base* parameter_;
    public:

        inline Argument() :
            operation_(NULL),
            parameter_(NULL) {
        }

        inline Argument(SourceCodeFragment<Base>& operation) :
            operation_(&operation),
            parameter_(NULL) {
        }

        inline Argument(const Base& parameter) :
            operation_(NULL),
            parameter_(new Base(parameter)) {
        }

        inline Argument(const Argument& orig) :
            operation_(orig.operation_),
            parameter_(orig.parameter_ != NULL ? new Base(*orig.parameter_) : NULL) {
        }

        inline Argument& operator=(const Argument& rhs) {
            if (&rhs == this) {
                return *this;
            }
            if (rhs.operation_ != NULL) {
                operation_ = rhs.operation_;
                delete parameter_;
                parameter_ = NULL;
            } else {
                operation_ = NULL;
                if (parameter_ != NULL) {
                    *parameter_ = *rhs.parameter_;
                } else {
                    parameter_ = new Base(*rhs.parameter_);
                }
            }
            return *this;
        }

        inline SourceCodeFragment<Base>* operation() const {
            return operation_;
        }

        inline Base* parameter() const {
            return parameter_;
        }

        virtual ~Argument() {
            delete parameter_;
        }
    };
}

#endif
