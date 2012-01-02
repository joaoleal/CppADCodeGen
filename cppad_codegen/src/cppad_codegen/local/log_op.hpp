#ifndef CPPAD_CODEGEN_LOG_OP_INCLUDED
#define CPPAD_CODEGEN_LOG_OP_INCLUDED

#include "ad_code_gen_name_provider.hpp"

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

CPPAD_BEGIN_NAMESPACE
/*!
\file log_op.hpp
Forward and reverse mode calculations for z = log(x).
 */

/*!
Compute forward mode Taylor coefficient for result of op = LogOp.

The C++ source code corresponding to this operation is
\verbatim
        z = log(x)
\endverbatim

\copydetails forward_unary1_op
 */
template <class Base>
inline void forward_code_gen_log_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    size_t k;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(LogOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(LogOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sz_d = n.generateVarName(d, i_z);
    std::string sx_0 = n.generateVarName(0, i_x);

    if (d == 0) {
        s_out << sz_d << " = log(" << sx_0 << ")" << n.endl();
    } else if (d == 1) {
        std::string sx_1 = n.generateVarName(1, i_x);
        s_out << sz_d << " = " << sx_1 << " / " << sx_0 << n.endl();
    } else {
        std::string sz_1 = n.generateVarName(1, i_z);
        std::string sx_d1 = n.generateVarName(d - 1, i_x);
        std::string sx_d = n.generateVarName(d, i_x);

        s_out << sz_d << " = ( ( -" << sz_1 << " * " << sx_d1;

        for (k = 2; k < d; k++) {
            std::string sz_k = n.generateVarName(k, i_z);
            std::string sx_dk = n.generateVarName(d - k, i_x);
            s_out << " -" << k << " * " << sz_k << " * " << sx_dk;
        }
        s_out << ") / " << n.CastToBase(d) << " + " << sx_d << ") / " << sx_0 << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = LogOp.

The C++ source code corresponding to this operation is
\verbatim
        z = log(x)
\endverbatim

\copydetails forward_unary1_op_0
 */
template <class Base>
inline void forward_code_gen_log_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(LogOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(LogOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sz_0 = n.generateVarName(0, i_z);
    std::string sx_0 = n.generateVarName(0, i_x);

    s_out << sz_0 << " = log(" << sx_0 << ")" << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = LogOp.

The C++ source code corresponding to this operation is
\verbatim
        z = log(x)
\endverbatim

\copydetails reverse_unary1_op
 */

template <class Base>
inline void reverse_code_gen_log_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(LogOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(LogOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    std::string px_0 = n.generatePartialName(0, i_x);
    std::string pz_0 = n.generatePartialName(0, i_z);
    std::string sx_0 = n.generateVarName(0, i_x);

    size_t j = d;
    while (j) { // scale partial w.r.t z[j]
        std::string px_j = n.generatePartialName(j, i_x);
        std::string pz_j = n.generatePartialName(j, i_z);
        std::string sz_j = n.generateVarName(j, i_z);

        s_out << pz_j << " /= " << sx_0 << n.endl();
        s_out << px_0 << " -= " << pz_j << " * " << sz_j << n.endl();
        s_out << px_j << " += " << pz_j << n.endl();

        // further scale partial w.r.t. z[j]
        s_out << pz_j << " /= " << n.CastToBase(j) << n.endl();

        for (size_t k = 1; k < j; k++) {
            std::string px_jk = n.generatePartialName(j - k, i_x);
            std::string pz_k = n.generatePartialName(k, i_z);
            std::string sz_k = n.generateVarName(k, i_z);
            std::string sx_jk = n.generateVarName(j - k, i_x);

            s_out << pz_k << " -= " << pz_j << " * " << n.CastToBase(k) << " * " << sx_jk << n.endl();
            s_out << px_jk << " -= " << pz_j << " * " << n.CastToBase(k) << " * " << sz_k << n.endl();
        }
        --j;
    }

    s_out << px_0 << " += " << pz_0 << " / " << sx_0 << n.endl();
}

CPPAD_END_NAMESPACE
#endif
