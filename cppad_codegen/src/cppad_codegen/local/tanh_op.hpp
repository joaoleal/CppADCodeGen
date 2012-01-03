#ifndef CPPAD_CODEGEN_TANH_OP_INCLUDED
#define CPPAD_CODEGEN_TANH_OP_INCLUDED

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
\file tanh_op.hpp
Forward and reverse mode calculations for z = tanh(x).
 */


/*!
Compute forward mode Taylor coefficient for result of op = TanOp.

The C++ source code corresponding to this operation is
\verbatim
        z = tanh(x)
\endverbatim
The auxillary result is
\verbatim
        y = tanh(x)^2
\endverbatim
The value of y, and its derivatives, are computed along with the value
and derivatives of z.

\copydetails forward_unary2_op
 */
template <class Base>
inline void forward_code_gen_tanh_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(TanOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(TanOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1;

    std::string sx_d = n.generateVarName(d, i_x);
    std::string sy_d = n.generateVarName(d, i_y);
    std::string sz_d = n.generateVarName(d, i_z);

    if (d == 0) {
        s_out << sz_d << " = tanh(" << sx_d << ")" << n.endl();
        s_out << sy_d << " = " << sz_d << " * " << sz_d << n.endl();
    } else {
        s_out << sz_d << " = " << sx_d;
        for (size_t k = 1; k <= d; k++) {
            std::string sx_k = n.generateVarName(k, i_x);
            std::string sy_dk = n.generateVarName(d - k, i_y);

            s_out << " - " << n.CastToBase(k) << " * " << sx_k << " * " << sy_dk << " / " << n.CastToBase(d);
        }
        s_out << n.endl();

        std::string sz_0 = n.generateVarName(0, i_z);

        s_out << sy_d << " = " << sz_0 << " * " << sz_d;
        for (size_t k = 1; k <= d; k++) {
            std::string sz_k = n.generateVarName(k, i_z);
            std::string sz_dk = n.generateVarName(d - k, i_z);

            s_out << " + " << sz_k << " * " << sz_dk;
        }
        s_out << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = TanOp.

The C++ source code corresponding to this operation is
\verbatim
        z = tanh(x)
\endverbatim
The auxillary result is
\verbatim
        y = cos(x)
\endverbatim
The value of y is computed along with the value of z.

\copydetails forward_unary2_op_0
 */
template <class Base>
inline void forward_code_gen_tanh_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(TanOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(TanOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1;

    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = tanh(" << sx_0 << ")" << n.endl();
    s_out << sy_0 << " = " << sz_0 << " * " << sz_0 << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = TanOp.

The C++ source code corresponding to this operation is
\verbatim
        z = tanh(x)
\endverbatim
The auxillary result is
\verbatim
        y = cos(x)
\endverbatim
The value of y is computed along with the value of z.

\copydetails reverse_unary2_op
 */

template <class Base>
inline void reverse_code_gen_tanh_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(TanOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(TanOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    size_t i_y = i_z - 1;

    size_t j = d;
    size_t k;
    std::string base_two = n.CastToBase(2);
    while (j) {
        std::string px_j = n.generatePartialName(j, i_x);
        std::string sy_jk = n.generateVarName(j - k, i_y);
        std::string pz_j = n.generatePartialName(j, i_z);

        s_out << px_j << " += " << pz_j << n.endl();
        s_out << pz_j << " /= " << n.CastToBase(j) << n.endl();

        for (k = 1; k <= j; k++) {
            std::string sx_k = n.generateVarName(k, i_x);
            std::string px_k = n.generatePartialName(k, i_x);
            std::string py_jk = n.generatePartialName(j - k, i_y);

            s_out << px_k << " -= " << pz_j << " * " << sy_jk << " * " << k << n.endl();
            s_out << py_jk << " -= " << pz_j << " * " << sx_k << " * " << k << n.endl();
        }

        for (k = 0; k < j; k++) {
            std::string py_j1 = n.generatePartialName(j - 1, i_y);
            std::string z_jk1 = n.generateVarName(j - k - 1, i_z);
            std::string pz_k = n.generatePartialName(k, i_z);
            s_out << pz_k << " += " << py_j1 << " * " << z_jk1 << " * " << base_two << n.endl();
        }

        --j;
    }
    std::string px_0 = n.generatePartialName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string pz_0 = n.generatePartialName(0, i_z);

    s_out << px_0 << " += " << pz_0 << " * (" << n.CastToBase(1) << " - " << sy_0 << ")" << n.endl();
}

CPPAD_END_NAMESPACE
#endif
