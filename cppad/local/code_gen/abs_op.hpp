#ifndef CPPAD_CODE_GEN_ABS_OP_INCLUDED
#define CPPAD_CODE_GEN_ABS_OP_INCLUDED

#include <cppad/local/code_gen/ad_code_gen_name_provider.hpp>

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
\file abs_op.hpp
Forward and reverse mode calculations for z = abs(x).
 */

/*!
Compute forward mode Taylor coefficient for result of op = AbsOp.

The C++ source code corresponding to this operation is
\verbatim
        z = abs(x)
\endverbatim

\copydetails forward_unary1_op
 */
template<class Base>
inline void forward_code_gen_abs_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x) {
    size_t k;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AbsOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AbsOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_d = n.generateVarName(d, i_x);
    std::string sz_d = n.generateVarName(d, i_z);

    std::string aux;

    if (d > 1) {
        aux = n.tempBaseVarName();
        // order that decides positive, negative or zero
        for (size_t j = 0; j < d - 1; j++) {
            std::string sx_j = n.generateVarName(j, i_x);
            if (j > 0) s_out << "else ";
            s_out << "if (" << sx_j << " == " << n.zero() << ") "
                    << aux << " = " << sx_j << n.endl();
        }

        s_out << "else "
                << aux << " = " << sx_d << n.endl();
    } else {
        aux = sx_d;
    }

    s_out << "if(" << aux << " < " << n.zero() << ") "
            << sz_d << " = -" << sx_d << n.endl() <<
            "else "
            << sz_d << " = " << sx_d << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = AbsOp.

The C++ source code corresponding to this operation is
\verbatim
        z = abs(x)
\endverbatim

\copydetails forward_unary1_op_0
 */
template<class Base>
inline void forward_code_gen_abs_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
size_t i_x) {

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(AbsOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(AbsOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx = n.generateVarName(0, i_x);
    std::string sz = n.generateVarName(0, i_z);

    s_out << "if(" << sx << " < " << n.zero() << ") "
            << sz << " = -" << sx << n.endl() <<
            "else "
            << sz << " = " << sx << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = AbsOp.

The C++ source code corresponding to this operation is
\verbatim
        z = abs(x)
\endverbatim

\copydetails reverse_unary1_op
 */
template<typename Base>
inline void reverse_code_gen_abs_op(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t d,
        size_t i_z,
        size_t i_x,
        size_t nc_partial,
        Base* partial) {
    //    size_t j, k;
    //    Base zero(0.);
    //
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(AbsOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(AbsOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(i_x < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Taylor coefficients and partials corresponding to argument
    //    const Base* x = taylor + i_x * nc_taylor;
    //    Base* px = partial + i_x * nc_partial;
    //
    //    // Taylor coefficients and partials corresponding to result
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // order that decides positive, negative or zero
    //    k = 0;
    //    while ((k < d) & (x[k] == zero))
    //        k++;
    //
    //    if (GreaterThanZero(x[k])) { // partial of z w.r.t y is +1
    //        for (j = k; j <= d; j++)
    //            px[j] += pz[j];
    //    } else if (LessThanZero(x[k])) { // partial of z w.r.t y is -1
    //        for (j = k; j <= d; j++)
    //            px[j] -= pz[j];
    //    }

    throw "not implemented yet";
}

CPPAD_END_NAMESPACE
#endif
