#ifndef CPPAD_CODEGEN_ACOS_OP_INCLUDED
#define CPPAD_CODEGEN_ACOS_OP_INCLUDED

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
\file acos_op.hpp
Forward and reverse mode calculations for z = acos(x).
 */


/*!
Compute forward mode Taylor coefficient for result of op = AcosOp.

The C++ source code corresponding to this operation is
\verbatim
        z = acos(x)
\endverbatim
The auxillary result is
\verbatim
        y = sqrt(1 - x * x)
\endverbatim
The value of y, and its derivatives, are computed along with the value
and derivatives of z.

\copydetails forward_unary2_op
 */
template<class Base>
inline void forward_code_gen_acos_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AcosOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AcosOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1; // ????

    std::string sx_d = n.generateVarName(d, i_x);
    std::string sy_d = n.generateVarName(d, i_y);
    std::string sz_d = n.generateVarName(d, i_z);

    size_t k;
    if (d == 0) {
        s_out << sz_d << " = acos(" << sx_d << ")" << n.endl();
        s_out << sy_d << " = sqrt(" << n.one() << " - " << sx_d << " * " << sx_d << ")" << n.endl();
    } else {
        s_out << sy_d << " = " << n.zero() << n.endl();
        s_out << sz_d << " = " << n.zero() << n.endl();

        for (k = 1; k < d; k++) {
            std::string sy_k = n.generateVarName(k, i_y);
            std::string sy_dk = n.generateVarName(d - k, i_y);
            std::string sz_k = n.generateVarName(k, i_z);

            s_out << sy_d << " -= " << k << " * " << sy_k << " * " << sy_dk << n.endl();
            s_out << sz_d << " -= " << k << " * " << sz_k << " * " << sy_dk << n.endl();
        }
        s_out << sy_d << " /= " << n.CastToBase(d) << n.endl();
        s_out << sz_d << " /= " << n.CastToBase(d) << n.endl();

        std::string qj;
        for (k = 0; k <= d; k++) {
            std::string sx_k = n.generateVarName(k, i_x);
            std::string sx_dk = n.generateVarName(d - k, i_x);
            qj += std::string("- ") + sx_k + " * " + sx_dk;
        }

        //
        s_out << sy_d << " += (" << qj << ") / " << n.CastToBase(2) << n.endl();
        s_out << sz_d << " -= " << sx_d << n.endl();
        //
        std::string sy_0 = n.generateVarName(0, i_y);
        s_out << sy_d << " /= " << sy_0 << n.endl();
        s_out << sz_d << " /= " << sy_0 << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = AcosOp.

The C++ source code corresponding to this operation is
\verbatim
        z = acos(x)
\endverbatim
The auxillary result is
\verbatim
        y = sqrt( 1 - x * x )
\endverbatim
The value of y is computed along with the value of z.

\copydetails forward_unary2_op_0
 */
template<class Base>
inline void forward_code_gen_acos_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AcosOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AcosOp) == 2);
    CPPAD_ASSERT_UNKNOWN(i_x + 1 < i_z);

    // Taylor coefficients corresponding to argument and result
    size_t i_y = i_z - 1; // ????
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = acos(" << sx_0 << ")" << n.endl();
    s_out << sy_0 << " = sqrt(" << n.one() << " - " << sx_0 << " * " << sx_0 << ")" << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = AcosOp.

The C++ source code corresponding to this operation is
\verbatim
        z = acos(x)
\endverbatim
The auxillary result is
\verbatim
        y = sqrt( 1 - x * x )
\endverbatim
The value of y is computed along with the value of z.

\copydetails reverse_unary2_op
 */

template <class Base>
inline void reverse_code_gen_acos_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AcosOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AcosOp) == 2);
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
        std::string sy_j = n.generateVarName(j, i_y);
        std::string sz_j = n.generateVarName(j, i_z);

        // scale partials w.r.t b[j] by 1 / b[0]
        s_out << py_j << " /= " << sy_0 << n.endl();

        // scale partials w.r.t z[j] by 1 / b[0]
        s_out << pz_j << " /= " << sy_0 << n.endl();

        // update partials w.r.t b^0 
        s_out << py_0 << " -= " << pz_j << " * " << sz_j << " + " << py_j << " * " << sy_j << n.endl();

        // update partial w.r.t. x^0
        s_out << px_0 << " -= " << py_j << " * " << sx_j << n.endl();

        // update partial w.r.t. x^j
        s_out << px_j << " -= " << pz_j << " + " << py_j << " * " << sx_0 << n.endl();

        // further scale partial w.r.t. z[j] by 1 / j
        s_out << pz_j << " /= " << n.CastToBase(j) << n.endl();

        for (k = 1; k < j; k++) { // update partials w.r.t b^(j-k)
            std::string px_k = n.generatePartialName(k, i_x);
            std::string sx_jk = n.generateVarName(j - k, i_x);
            std::string py_jk = n.generatePartialName(j - k, i_y);
            std::string sy_k = n.generateVarName(k, i_y);
            std::string sy_jk = n.generateVarName(j - k, i_y);
            std::string sz_k = n.generateVarName(k, i_z);
            std::string pz_k = n.generatePartialName(k, i_z);

            s_out << py_jk << " -= " << n.CastToBase(k) << " * " << pz_j << " * " << sz_k << " + " << py_j << " * " << sy_k << n.endl();

            // update partials w.r.t. x^k 
            s_out << px_k << " -= " << py_j << " * " << sx_jk << n.endl();

            // update partials w.r.t. z^k
            s_out << pz_k << " -= " << pz_j << " * " << n.CastToBase(k) << " * " << sy_jk << n.endl();
        }
        --j;
    }

    // j == 0 case
    s_out << px_0 << " -= (" << pz_0 << " + " << py_0 << " * " << sx_0 << ") / " << sy_0 << n.endl();
}

CPPAD_END_NAMESPACE
#endif
