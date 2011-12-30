#ifndef CPPAD_CODEGEN_SUB_OP_INCLUDED
#define CPPAD_CODEGEN_SUB_OP_INCLUDED

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
\file sub_op.hpp
Forward and reverse mode calculations for z = x - y.
 */

// --------------------------- Subvv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op
 */

template <class Base>
inline void forward_code_gen_subvv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SubvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(SubvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sx = n.generateVarName(d, arg[0]);
    std::string sy = n.generateVarName(d, arg[1]);
    std::string sz = n.generateVarName(d, i_z);

    s_out << sz << " = " << sx << " - " << sy << ";";
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_subvv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SubvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(SubvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sx = n.generateVarName(0, arg[0]);
    std::string sy = n.generateVarName(0, arg[1]);
    std::string sz = n.generateVarName(0, i_z);

    s_out << sz << " = " << sx << " - " << sy << ";";
}

/*!
Compute reverse mode partial derivatives for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_subvv_op(
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
    //    CPPAD_ASSERT_UNKNOWN(NumArg(SubvvOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(SubvvOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Partial derivatives corresponding to arguments and result
    //    Base* px = partial + arg[0] * nc_partial;
    //    Base* py = partial + arg[1] * nc_partial;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // number of indices to access
    //    size_t i = d + 1;
    //    while (i) {
    //        --i;
    //        px[i] += pz[i];
    //        py[i] -= pz[i];
    //    }
    throw "not implemented yet";
}

// --------------------------- Subpv -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = SubpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op
 */

template <class Base>
inline void forward_code_gen_subpv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SubpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(SubpvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sy = n.generateVarName(d, arg[1]);
    std::string sz = n.generateVarName(d, i_z);

    // Parameter value
    if (d == 0) {
        Base x = parameter[ arg[0] ];
        s_out << sz << " = " << n.PrintBase(x) << " - " << sy << n.endl();
    } else {
        s_out << sz << " = -" << sy << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = SubpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_subpv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SubpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(SubpvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sy = n.generateVarName(0, arg[1]);
    std::string sz = n.generateVarName(0, i_z);

    // Parameter value
    Base x = parameter[ arg[0] ];
    s_out << sz << " = " << n.PrintBase(x) << " - " << sy << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = SubpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_subpv_op(
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
    //    CPPAD_ASSERT_UNKNOWN(NumArg(SubvvOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(SubvvOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Partial derivatives corresponding to arguments and result
    //    Base* py = partial + arg[1] * nc_partial;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // number of indices to access
    //    size_t i = d + 1;
    //    while (i) {
    //        --i;
    //        py[i] -= pz[i];
    //    }
    throw "not implemented yet";
}

// --------------------------- Subvp -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_binary_op
 */

template <class Base>
inline void forward_code_gen_subvp_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SubvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(SubvpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sx = n.generateVarName(d, arg[0]);
    std::string sz = n.generateVarName(d, i_z);

    // Parameter value
    if (d == 0) {
        Base y = parameter[ arg[1] ];
        s_out << sz << " = " << sx << " - " << n.PrintBase(y) << n.endl();
    } else {
        s_out << sz << " = " << sx << n.endl();
    }
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = SubvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_subvp_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(SubvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(SubvpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sx = n.generateVarName(0, arg[0]);
    std::string sz = n.generateVarName(0, i_z);

    Base y = parameter[ arg[1] ];
    s_out << sz << " = " << sx << " - " << n.PrintBase(y) << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = SubvpOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x - y
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_subvp_op(
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
    //    CPPAD_ASSERT_UNKNOWN(NumArg(SubvpOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(SubvpOp) == 1);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // Partial derivatives corresponding to arguments and result
    //    Base* px = partial + arg[0] * nc_partial;
    //    Base* pz = partial + i_z * nc_partial;
    //
    //    // number of indices to access
    //    size_t i = d + 1;
    //    while (i) {
    //        --i;
    //        px[i] += pz[i];
    //    }

    throw "not implemented yet";
}

CPPAD_END_NAMESPACE
#endif
