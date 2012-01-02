#ifndef CPPAD_CODEGEN_DIV_OP_INCLUDED
#define CPPAD_CODEGEN_DIV_OP_INCLUDED

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
const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Using CondExp, it can make sense to divide by zero
    // so do not make it an error.
    size_t i_x = arg[0];
    size_t i_y = arg[1];

    std::string py_0 = n.generatePartialName(0, i_y);
    std::string sy_0 = n.generateVarName(0, i_y);

    size_t j = d + 1;
    while (j) {
        --j;
        std::string pz_j = n.generatePartialName(j, i_z);
        std::string sz_j = n.generateVarName(j, i_z);
        std::string px_j = n.generatePartialName(j, i_x);

        // scale partial w.r.t. z[j]
        s_out << pz_j << " /= " << sy_0 << n.endl();

        s_out << px_j << " += " << pz_j << n.endl();
        for (size_t k = 1; k <= j; k++) {
            std::string pz_jk = n.generatePartialName(j - k, i_z);
            std::string sy_k = n.generateVarName(k, i_y);
            std::string sz_jk = n.generateVarName(j - k, i_z);
            std::string py_k = n.generatePartialName(k, i_y);

            s_out << pz_jk << " -= " << pz_j << " * " << sy_k << n.endl();
            s_out << py_k << " -= " << pz_j << " * " << sz_jk << n.endl();
        }

        s_out << py_0 << " -= " << pz_j << " * " << sz_j << n.endl();
    }
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
const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Using CondExp, it can make sense to divide by zero
    // so do not make it an error.
    size_t i_y = arg[1];

    std::string py_0 = n.generatePartialName(0, i_y);
    std::string sy_0 = n.generateVarName(0, i_y);

    size_t j = d + 1;
    while (j) {
        --j;
        std::string pz_j = n.generatePartialName(j, i_z);
        std::string sz_j = n.generateVarName(j, i_z);

        // scale partial w.r.t. z[j]
        s_out << pz_j << " /= " << sy_0 << n.endl();

        for (size_t k = 1; k <= j; k++) {
            std::string pz_jk = n.generatePartialName(j - k, i_z);
            std::string sy_k = n.generateVarName(k, i_y);
            std::string sz_jk = n.generateVarName(j - k, i_z);
            std::string py_k = n.generatePartialName(k, i_y);

            s_out << pz_jk << " -= " << pz_j << " * " << sy_k << n.endl();
            s_out << py_k << " -= " << pz_j << " * " << sz_jk << n.endl();
        }

        s_out << py_0 << " -= " << pz_j << " * " << sz_j << n.endl();
    }
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
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DivvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DivvpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    // Using CondExp, it can make sense to divide by zero
    // so do not make it an error.

    size_t i_x = arg[0];
    
    // Argument values
    std::string sy = n.PrintBase(parameter[ arg[1] ]);
    
    size_t j = d + 1;
    while (j) {
        --j;
        std::string px_j = n.generatePartialName(j, i_x);
        std::string pz_j = n.generatePartialName(j, i_z);

        s_out << px_j << " += " << pz_j << " / " << sy << n.endl();
    }
}

CPPAD_END_NAMESPACE
#endif
