#ifndef CPPAD_CODEGEN_ATAN_OP_INCLUDED
#define CPPAD_CODEGEN_ATAN_OP_INCLUDED

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
\file atan_op.hpp
Forward and reverse mode calculations for z = atan(x).
 */


/*!
Forward mode Taylor coefficient for result of op = AtanOp.

The C++ source code corresponding to this operation is
\verbatim
        z = atan(x)
\endverbatim
The auxillary result is
\verbatim
        y = 1 + x * x
\endverbatim
The value of y, and its derivatives, are computed along with the value
and derivatives of z.

\copydetails forward_unary2_op
 */
template <class Base>
inline void forward_code_gen_atan_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AtanOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AtanOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1; // ????

    std::string sx_d = n.generateVarName(d, i_x);
    std::string sy_d = n.generateVarName(d, i_y);
    std::string sz_d = n.generateVarName(d, i_z);

    if (d == 0) {
        s_out << sz_d << " = atan(" << sx_d << ")" << n.endl();
        s_out << sy_d << " = " << n.one() << " + " << sx_d << " * " << sx_d << n.endl();
    } else {
        std::string sx_0 = n.generateVarName(0, i_x);

        s_out << sy_d << " = 2 * " << sx_0 << " * " << sx_d << n.endl();
        s_out << sz_d << " = " << n.zero() << n.endl();

        for (size_t k = 1; k < d; k++) {
            std::string sy_dk = n.generateVarName(d - k, i_y);
            std::string sx_k = n.generateVarName(k, i_x);
            std::string sx_dk = n.generateVarName(d - k, i_x);
            std::string sz_k = n.generateVarName(k, i_z);

            s_out << sy_d << " += " << sx_k << " * " << sx_dk << n.endl();
            s_out << sz_d << " -= " << k << " * " << sz_k << " * " << sy_dk << n.endl();
        }

        std::string sy_0 = n.generateVarName(0, i_y);

        s_out << sz_d << " /= " << n.CastToBase(d) << n.endl();
        s_out << sz_d << " += " << sx_d << n.endl();
        s_out << sz_d << " /= " << sy_0 << n.endl();
    }
}

/*!
Zero order forward mode Taylor coefficient for result of op = AtanOp.

The C++ source code corresponding to this operation is
\verbatim
        z = atan(x)
\endverbatim
The auxillary result is
\verbatim
        y = 1 + x * x
\endverbatim
The value of y is computed along with the value of z.

\copydetails forward_unary2_op_0
 */
template <class Base>
inline void forward_code_gen_atan_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AtanOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AtanOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1; // ????

    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = atan(" << sx_0 << ")" << n.endl();
    s_out << sy_0 << " = " << n.one() << " + " << sx_0 << " * " << sx_0 << n.endl();
}

/*!
Reverse mode partial derivatives for result of op = AtanOp.

The C++ source code corresponding to this operation is
\verbatim
        z = atan(x)
\endverbatim
The auxillary result is
\verbatim
        y = 1 + x * x
\endverbatim
The value of y is computed along with the value of z.

\copydetails reverse_unary2_op
 */

template <class Base>
inline void reverse_code_gen_atan_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AtanOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AtanOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    size_t i_y = i_z - 1;

    std::string px_0 = n.generatePartialName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string py_0 = n.generatePartialName(0, i_y);
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string pz_0 = n.generatePartialName(0, i_z);

    // number of indices to access
    size_t j = d;
    size_t k;
    while (j) {
        std::string sx_j = n.generateVarName(j, i_x);
        std::string px_j = n.generatePartialName(j, i_x);
        std::string py_j = n.generatePartialName(j, i_y);
        std::string pz_j = n.generatePartialName(j, i_z);
        std::string sz_j = n.generateVarName(j, i_z);

        // scale partials w.r.t z[j] and b[j]
        s_out << pz_j << " /= " << sy_0 << n.endl();

        s_out << py_j << " *= " << n.CastToBase(2) << n.endl();

        s_out << py_0 << " -= " << pz_j << " * " << sz_j << n.endl();
        s_out << px_j << " += " << pz_j << " + " << py_j << " * " << sx_0 << n.endl();
        s_out << px_0 << " += " << py_j << " * " << sx_j << n.endl();

        // more scaling of partials w.r.t z[j]
        s_out << pz_j << " /= " << n.CastToBase(j) << n.endl();

        for (k = 1; k < j; k++) { // update partials w.r.t b^(j-k)
            std::string px_k = n.generatePartialName(k, i_x);
            std::string sx_jk = n.generateVarName(j - k, i_x);
            std::string py_jk = n.generatePartialName(j - k, i_y);
            std::string sy_jk = n.generateVarName(j - k, i_y);
            std::string sz_k = n.generateVarName(k, i_z);
            std::string pz_k = n.generatePartialName(k, i_z);

            s_out << py_jk << " -= " << pz_j << " * " << n.CastToBase(k) << " * " << sz_k << n.endl();
            s_out << pz_k << " -= " << pz_j << " * " << n.CastToBase(k) << " * " << sy_jk << n.endl();
            s_out << px_k << " += " << py_j << " * " << sx_jk << n.endl();
        }
        --j;
    }

    // j == 0 case
    s_out << px_0 << " += " << pz_0 << " / " << sy_0 << " + " << py_0 << " * " << n.CastToBase(2) << " * " << sx_0 << n.endl();
}

CPPAD_END_NAMESPACE
#endif
