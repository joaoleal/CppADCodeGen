#ifndef CPPAD_CG_SOURCE_CODE_FRAGMENT_INCLUDED
#define	CPPAD_CG_SOURCE_CODE_FRAGMENT_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <vector>

namespace CppAD {

    /**
     * An operation
     * 
     * \author Joao Leal
     */
    template<class Base>
    class SourceCodeFragment {
    private:
        // the operations used to create this variable (temporary variables only)
        CGOpCode operation_;
        // the code blocks this block depends upon (empty for independent 
        // variables and possibly for the 1st assignment of a dependent variable)
        std::vector<Argument<Base> > arguments_;
        // variable ID that was altered/assigned in this source code
        // (zero means that no variable is assigned)
        size_t var_id_;
        // the number of times the result of this operation is used
        size_t use_count_;

    public:

        inline SourceCodeFragment(CGOpCode op) :
            operation_(op),
            var_id_(0),
            use_count_(0) {
        }

        inline SourceCodeFragment(CGOpCode op,
                                  const Argument<Base>& arg) :
            operation_(op),
            arguments_(1),
            var_id_(0),
            use_count_(0) {
            assert(arg.operation() != NULL);
            arguments_[0] = arg;
        }

        inline SourceCodeFragment(CGOpCode op,
                                  const Argument<Base>& arg1,
                                  const Argument<Base>& arg2) :
            operation_(op),
            arguments_(2),
            var_id_(0),
            use_count_(0) {
            assert(arg1.operation() != NULL || arg2.operation() != NULL);
            arguments_[0] = arg1;
            arguments_[1] = arg2;
        }

        inline SourceCodeFragment(CGOpCode op,
                                  const Argument<Base>& arg1,
                                  const Argument<Base>& arg2,
                                  const Argument<Base>& arg3,
                                  const Argument<Base>& arg4) :
            operation_(op),
            arguments_(4),
            var_id_(0),
            use_count_(0) {
            assert(arg1.operation() != NULL || arg2.operation() != NULL ||
                   arg3.operation() != NULL || arg4.operation() != NULL);
            arguments_[0] = arg1;
            arguments_[1] = arg2;
            arguments_[2] = arg3;
            arguments_[3] = arg4;
        }

        SourceCodeFragment(const SourceCodeFragment& orig) :
            operation_(orig.operation_),
            arguments_(orig.arguments_),
            var_id_(0),
            use_count_(0) {
        }

        inline CGOpCode operation() const {
            return operation_;
        }

        /**
         * Provides the arguments used in the operation represnted by this
         * code fragment.
         * \return the arguments for the operation in this code fragment
         */
        inline const std::vector<Argument<Base> >& arguments() const {
            return arguments_;
        }

        /**
         * Provides the variable ID that was altered/assigned in this source 
         * code (zero means that no variable is assigned).
         * \return the variable ID
         */
        inline size_t variableID() const {
            return var_id_;
        }

        /**
         * Specifies a variable ID to the result of this source code
         * (zero means that no variable is created).
         */
        inline void setVariableID(size_t var_id) {
            var_id_ = var_id;
        }

        /**
         * Provides the number of times this code fragement is marked as being 
         * used as an argument of another code fragement.
         * \return the usage count
         */
        inline size_t usageCount() const {
            return use_count_;
        }

        virtual ~SourceCodeFragment() {
        }


        friend class CodeHandler<Base>;

    };

    template<class Base>
    inline std::ostream& operator <<(
    std::ostream& os, //< stream to write to
    const CppAD::SourceCodeFragment<Base>& c) {
        switch (c.operation()) {
            case CGAbsOp:
                os << "abs( $1 )";
                break;
            case CGAcosOp:
                os << "acos( $1 )";
                break;
            case CGAddOp:
                os << "$1 + $2";
                break;
            case CGAsinOp:
                os << "asin( $1 )";
                break;
            case CGAtanOp:
                os << "atan( $1 )";
                break;
            case CGComOpLt:
                os << "($1 < $2)? $3 : $4";
                break;
            case CGComOpLe:
                os << "($1 <= $2)? $3 : $4";
                break;
            case CGComOpEq:
                os << "($1 == $2)? $3 : $4";
                break;
            case CGComOpGe:
                os << "($1 > $2)? $3 : $4";
                break;
            case CGComOpGt:
                os << "($1 >= $2)? $3 : $4";
                break;
            case CGComOpNe:
                os << "($1 != $2)? $3 : $4";
                break;
            case CGCoshOp:
                os << "cosh( $1 )";
                break;
            case CGCosOp:
                os << "cosh( $1 )";
                break;
            case CGDivOp:
                os << "$1 / $2";
                break;
            case CGExpOp:
                os << "e^$1";
                break;
            case CGInvOp:
                os << "independent( $1 )";
                break;
            case CGLogOp:
                os << "log( $1 )";
                break;
            case CGMulOp:
                os << "$1 * $2";
                break;
            case CGPowOp:
                os << "$1^$2";
                break;
            case CGSignOp:
                os << "if($1 > 0) { 1 } else if($1 == 0) { 0 } else { -1 }";
                break;
            case CGSinhOp:
                os << "sinh( $1 )";
                break;
            case CGSinOp:
                os << "sin( $1 )";
                break;
            case CGSqrtOp:
                os << "sqrt( $1 )";
                break;
            case CGSubOp:
                os << "$1 - $2";
                break;
            case CGTanhOp:
                os << "tanh( $1 )";
                break;
            case CGTanOp:
                os << "tan( $1 )";
                break;
            case CGUnMinusOp:
                os << "-$1";
                break;
            default:
                os << "???";
        }

        return os;
    }

}

#endif
