#ifndef CPPAD_CG_ARGUMENT_INCLUDED
#define CPPAD_CG_ARGUMENT_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    template<class Base>
    class OperationNode;

    /**
     * An argument used by an operation which can be either a constant value
     * or the result of another operation
     * 
     * @author Joao Leal
     */
    template<class Base>
    class Argument {
    private:
        OperationNode<Base>* operation_;
        Base* parameter_;
    public:

        inline Argument() :
            operation_(NULL),
            parameter_(NULL) {
        }

        inline Argument(OperationNode<Base>& operation) :
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

        inline OperationNode<Base>* getOperation() const {
            return operation_;
        }

        inline Base* getParameter() const {
            return parameter_;
        }

        virtual ~Argument() {
            delete parameter_;
        }
    };
}

#endif
