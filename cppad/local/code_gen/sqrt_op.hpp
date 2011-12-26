#ifndef CPPAD_CODE_GEN_SQRT_OP_INCLUDED
#define CPPAD_CODE_GEN_SQRT_OP_INCLUDED

#include "ad_code_gen_name_provider.hpp"

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-10 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */


CPPAD_BEGIN_NAMESPACE
/*!
\file sqrt_op.hpp
Forward and reverse mode calculations for z = sqrt(x).
 */


/*!
Compute forward mode Taylor coefficient for result of op = SqrtOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sqrt(x)
\endverbatim

\copydetails forward_unary1_op
 */
template <class Base>
inline void forward_code_gen_sqrt_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SqrtOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SqrtOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sz_d = n.generateVarName(d, i_z);

    if (d == 0) {
        std::string sx_0 = n.generateVarName(0, i_x);
        s_out << sz_d << " = sqrt(" << sx_0 << ")" << n.endl();
    } else {
        std::string sz_0 = n.generateVarName(0, i_z);
        std::string sx_d = n.generateVarName(d, i_x);

        s_out << sz_d << " = (";
        if (d > 1) {
            s_out << "(";
            for (size_t k = 1; k < d; k++) {
                std::string sz_k = n.generateVarName(k, i_z);
                std::string sz_dk = n.generateVarName(d - k, i_z);

                s_out << sz_d << " - " << k << " * " << sz_k << " * " << sz_dk;
            }
            s_out << ") / " << n.CastToBase(d) << " + ";
        }

        s_out << sx_d << " / " << n.CastToBase(2);
        s_out << ") / " << sz_0 << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = SqrtOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sqrt(x)
\endverbatim

\copydetails forward_unary1_op_0
 */
template <class Base>
inline void forward_code_gen_sqrt_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SqrtOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(SqrtOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = sqrt(" << sx_0 << ")" << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = SqrtOp.

The C++ source code corresponding to this operation is
\verbatim
        z = sqrt(x)
\endverbatim

\copydetails reverse_unary1_op
 */

template <class Base>
inline void reverse_code_gen_sqrt_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x,
size_t nc_taylor,
const Base* taylor,
size_t nc_partial,
Base* partial) {
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(SqrtOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(SqrtOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(i_x < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Taylor coefficients and partials corresponding to argument
    //    Base* px = partial + i_x * nc_partial;
    //
    //    // Taylor coefficients and partials corresponding to result
    //    const Base* z = taylor + i_z * nc_taylor;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // number of indices to access
    //    size_t j = d;
    //    size_t k;
    //    while (j) { // scale partial w.r.t. z[j]
    //        pz[j] /= z[0];
    //
    //        pz[0] -= pz[j] * z[j];
    //        px[j] += pz[j] / Base(2);
    //        for (k = 1; k < j; k++)
    //            pz[k] -= pz[j] * z[j - k];
    //        --j;
    //    }
    //    px[0] += pz[0] / (Base(2) * z[0]);
    throw "not implemented yet";
}

CPPAD_END_NAMESPACE
#endif
