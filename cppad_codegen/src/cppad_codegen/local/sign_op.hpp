#ifndef CPPAD_CODEGEN_SIGN_OP_INCLUDED
#define CPPAD_CODEGEN_SIGN_OP_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <cppad_codegen/local/ad_code_gen_name_provider.hpp>
#include <cppad_codegen/local/sign.hpp>

CPPAD_BEGIN_NAMESPACE

template <class Base>
inline void forward_code_gen_sign_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SignOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SignOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_d = n.generateVarName(d, i_x);
    std::string sz_d = n.generateVarName(d, i_z);


    if (d == 0)
        s_out << sz_d << " = " << code_gen_sign(n, sx_d) << n.endl();
    else
        s_out << sz_d << " = " << n.zero() << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = SignOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sign(x)
\endverbatim

\copydetails forward_unary1_op_0
 */
template <class Base>
inline void forward_code_gen_sign_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SignOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SignOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sz = n.generateVarName(0, i_z);

    s_out << sz << " = " << code_gen_sign(n, sx_0) << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = SignOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sign(x)
\endverbatim

\copydetails reverse_unary1_op
 */

template <class Base>
inline void reverse_code_gen_sign_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SignOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SignOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // nothing to do because partials of sign are zero
    return;
}

CPPAD_END_NAMESPACE

#endif	/* SIGN_OP_HPP */

