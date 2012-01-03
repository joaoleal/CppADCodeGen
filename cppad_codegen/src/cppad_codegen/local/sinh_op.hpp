#ifndef CPPAD_CODEGEN_SINH_OP_INCLUDED
#define CPPAD_CODEGEN_SINH_OP_INCLUDED

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
\file sinh_op.hpp
Forward and reverse mode calculations for z = sinh(x).
 */


/*!
Compute forward mode Taylor coefficient for result of op = SinhOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sinh(x)
\endverbatim
The auxillary result is
\verbatim
        y = cosh(x)
\endverbatim
The value of y, and its derivatives, are computed along with the value
and derivatives of z.

\copydetails forward_unary2_op
 */
template <class Base>
inline void forward_code_gen_sinh_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SinhOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SinhOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1;

    std::string sx_d = n.generateVarName(d, i_x);
    std::string sy_d = n.generateVarName(d, i_y);
    std::string sz_d = n.generateVarName(d, i_z);

    // rest of this routine is identical for the following cases:
    // forward_sin_op, forward_cos_op, forward_sinh_op, forward_cosh_op.
    if (d == 0) {
        s_out << sz_d << " = sinh(" << sx_d << ")" << n.endl();
        s_out << sy_d << " = cosh(" << sx_d << ")" << n.endl();
    } else {
        s_out << sz_d << " = " << n.zero() << n.endl();
        s_out << sy_d << " = " << n.zero() << n.endl();
        for (size_t k = 1; k <= d; k++) {
            std::string sx_k = n.generateVarName(k, i_x);
            std::string sy_dk = n.generateVarName(d - k, i_y);
            std::string sz_dk = n.generateVarName(d - k, i_z);

            s_out << sz_d << " += " << k << " * " << sx_k << " * " << sy_dk << n.endl();
            s_out << sy_d << " += " << k << " * " << sx_k << " * " << sz_dk << n.endl();
        }
        s_out << sz_d << " /= " << n.CastToBase(d) << n.endl();
        s_out << sy_d << " /= " << n.CastToBase(d) << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = SinhOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sinh(x)
\endverbatim
The auxillary result is
\verbatim
        y = cosh(x)
\endverbatim
The value of y is computed along with the value of z.

\copydetails forward_unary2_op_0
 */
template <class Base>
inline void forward_code_gen_sinh_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SinhOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SinhOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1;

    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = sinh(" << sx_0 << ")" << n.endl();
    s_out << sy_0 << " = cosh(" << sx_0 << ")" << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = SinhOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sinh(x)
\endverbatim
The auxillary result is
\verbatim
        y = cosh(x)
\endverbatim
The value of y is computed along with the value of z.

\copydetails reverse_unary2_op
 */

template <class Base>
inline void reverse_code_gen_sinh_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SinhOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SinhOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    size_t i_y = i_z - 1;

    // rest of this routine is identical for the following cases:
    // reverse_sin_op, reverse_cos_op, reverse_sinh_op, reverse_cosh_op.
    size_t j = d;
    size_t k;
    while (j) {
        std::string py_j = n.generatePartialName(j, i_y);
        std::string sy_jk = n.generateVarName(j - k, i_y);
        std::string pz_j = n.generatePartialName(j, i_z);
        std::string sz_jk = n.generateVarName(j - k, i_z);

        s_out << pz_j << " /= " << n.CastToBase(j) << n.endl();
        s_out << py_j << " /= " << n.CastToBase(j) << n.endl();

        for (k = 1; k <= j; k++) {
            std::string sx_k = n.generateVarName(k, i_x);
            std::string px_k = n.generatePartialName(k, i_x);
            std::string py_jk = n.generatePartialName(j - k, i_y);
            std::string pz_jk = n.generatePartialName(j - k, i_z);

            s_out << px_k << " += " << pz_j << " * " << k << " * " << sy_jk
                    << " + " << py_j << " * " << k << " * " << sz_jk << n.endl();

            s_out << pz_jk << " += " << py_j << " * " << k << " * " << sx_k << n.endl();
            s_out << py_jk << " += " << pz_j << " * " << k << " * " << sx_k << n.endl();
        }
        --j;
    }

    std::string px_0 = n.generatePartialName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string py_0 = n.generatePartialName(0, i_y);
    std::string sz_0 = n.generateVarName(0, i_z);
    std::string pz_0 = n.generatePartialName(0, i_z);

    s_out << px_0 << " += " << py_0 << " * " << sz_0
            << " + " << pz_0 << " * " << sy_0 << n.endl();
}

CPPAD_END_NAMESPACE
#endif
