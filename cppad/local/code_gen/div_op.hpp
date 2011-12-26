#ifndef CPPAD_CODE_GEN_DIV_OP_INCLUDED
#define CPPAD_CODE_GEN_DIV_OP_INCLUDED

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
\file div_op.hpp
Forward and reverse mode calculations for z = x / y.
 */

// --------------------------- Divvv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = DivvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op
 */
template <class Base>
inline void forward_code_gen_divvv_op(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t d,
        size_t i_z,
        const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    size_t i_x = arg[0];
    size_t i_y = arg[1];

    // Taylor coefficients corresponding to arguments and result
    std::string sx_d = n.generateVarName(d, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_d = n.generateVarName(d, i_z);

    // Using CondExp, it can make sense to divide by zero,
    // so do not make it an error.
    s_out << sz_d << " = (" << sx_d;
    for (size_t k = 1; k <= d; k++) {
        std::string sz_dk = n.generateVarName(d - k, i_z);
        std::string sy_k = n.generateVarName(k, i_y);

        s_out << " - " << sz_dk << " * " << sy_k;
    }
    s_out << ") / " << sy_0 << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = DivvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_divvv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    size_t i_x = arg[0];
    size_t i_y = arg[1];

    // Taylor coefficients corresponding to arguments and result
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = " << sx_0 << " / " << sy_0 << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = DivvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_divvv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter,
size_t nc_taylor,
const Base* taylor,
size_t nc_partial,
Base* partial) {
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(DivvvOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(DivvvOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Arguments
    //    const Base* y = taylor + arg[1] * nc_taylor;
    //    const Base* z = taylor + i_z * nc_taylor;
    //
    //    // Partial derivatives corresponding to arguments and result
    //    Base* px = partial + arg[0] * nc_partial;
    //    Base* py = partial + arg[1] * nc_partial;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // Using CondExp, it can make sense to divide by zero
    //    // so do not make it an error.
    //
    //    size_t k;
    //    // number of indices to access
    //    size_t j = d + 1;
    //    while (j) {
    //        --j;
    //        // scale partial w.r.t. z[j]
    //        pz[j] /= y[0];
    //
    //        px[j] += pz[j];
    //        for (k = 1; k <= j; k++) {
    //            pz[j - k] -= pz[j] * y[k];
    //            py[k] -= pz[j] * z[j - k];
    //        }
    //        py[0] -= pz[j] * z[j];
    //    }
    throw "not implemented yet";
}

// --------------------------- Divpv -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = DivpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op
 */

template <class Base>
inline void forward_code_gen_divpv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivpvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    size_t i_y = arg[1];

    // Taylor coefficients corresponding to arguments and result
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_d = n.generateVarName(d, i_z);

    // Using CondExp, it can make sense to divide by zero,
    // so do not make it an error.
    s_out << sz_d << " = (";
    if (d == 0) {
        // Parameter value
        Base x = parameter[ arg[0] ];
        
        s_out << n.PrintBase(x);
    } else {
        for (size_t k = 1; k <= d; k++) {
            std::string sz_dk = n.generateVarName(d - k, i_z);
            std::string sy_k = n.generateVarName(k, i_y);

            s_out << " - " << sz_dk << " * " << sy_k;
        }
    }

    s_out << ") / " << sy_0 << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = DivpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_divpv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivpvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    size_t i_y = arg[1];

    // Parameter value
    Base x = parameter[ arg[0] ];

    // Taylor coefficients corresponding to arguments and result
    std::string sy_0 = n.generateVarName(0, i_y);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = " << n.PrintBase(x) << " / " << sy_0 << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = DivpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_divpv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter,
size_t nc_taylor,
const Base* taylor,
size_t nc_partial,
Base* partial) {
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(DivvvOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(DivvvOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Arguments
    //    const Base* y = taylor + arg[1] * nc_taylor;
    //    const Base* z = taylor + i_z * nc_taylor;
    //
    //    // Partial derivatives corresponding to arguments and result
    //    Base* py = partial + arg[1] * nc_partial;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // Using CondExp, it can make sense to divide by zero so do not
    //    // make it an error.
    //
    //    size_t k;
    //    // number of indices to access
    //    size_t j = d + 1;
    //    while (j) {
    //        --j;
    //        // scale partial w.r.t z[j]
    //        pz[j] /= y[0];
    //
    //        for (k = 1; k <= j; k++) {
    //            pz[j - k] -= pz[j] * y[k];
    //            py[k] -= pz[j] * z[j - k];
    //        }
    //        py[0] -= pz[j] * z[j];
    //    }
    throw "not implemented yet";
}


// --------------------------- Divvp -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = DivvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_binary_op
 */

template <class Base>
inline void forward_code_gen_divvp_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    size_t i_x = arg[0];

    // Parameter value
    Base y = parameter[ arg[1] ];

    // Taylor coefficients corresponding to arguments and result
    std::string sx_d = n.generateVarName(d, i_x);
    std::string sz_d = n.generateVarName(d, i_z);

    // Using CondExp and multiple levels of AD, it can make sense 
    // to divide by zero so do not make it an error.
    s_out << sz_d << " = " << sx_d << " / " << n.PrintBase(y) << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = DivvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_divvp_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    size_t i_x = arg[0];

    // Parameter value
    Base y = parameter[ arg[1] ];

    // Taylor coefficients corresponding to arguments and result
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = " << sx_0 << " / " << n.PrintBase(y) << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = DivvpOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x / y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_divvp_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter,
size_t nc_taylor,
const Base* taylor,
size_t nc_partial,
Base* partial) {
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(DivvpOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(DivvpOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Argument values
    //    Base y = parameter[ arg[1] ];
    //
    //    // Partial derivatives corresponding to arguments and result
    //    Base* px = partial + arg[0] * nc_partial;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // Using CondExp, it can make sense to divide by zero
    //    // so do not make it an error.
    //
    //    // number of indices to access
    //    size_t j = d + 1;
    //    while (j) {
    //        --j;
    //        px[j] += pz[j] / y;
    //    }
    throw "not implemented yet";
}

CPPAD_END_NAMESPACE
#endif
