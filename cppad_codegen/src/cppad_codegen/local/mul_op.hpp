#ifndef CPPAD_CODEGEN_MUL_OP_INCLUDED
#define CPPAD_CODEGEN_MUL_OP_INCLUDED

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
\file mul_op.hpp
Forward and reverse mode calculations for z = x * y.
 */

// --------------------------- Mulvv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = MulvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x * y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op
 */
template <class Base>
inline void forward_code_gen_mulvv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(MulvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(MulvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    size_t i_x = arg[0];
    size_t i_y = arg[1];

    std::string sz_d = n.generateVarName(d, i_z);

    s_out << sz_d << " = ";
    for (size_t k = 0; k <= d; k++) {
        std::string sx_dk = n.generateVarName(d - k, i_x);
        std::string sy_k = n.generateVarName(k, i_y);

        s_out << " + " << sx_dk << " * " << sy_k;
    }
    s_out << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = MulvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x * y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_binary_op_0
 */
template <class Base>
inline void forward_code_gen_mulvv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(MulvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(MulvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    size_t i_x = arg[0];
    size_t i_y = arg[1];

    std::string sz_0 = n.generateVarName(0, i_z);
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sy_0 = n.generateVarName(0, i_y);

    s_out << sz_0 << " = " << sx_0 << " * " << sy_0 << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = MulvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x * y
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_mulvv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(MulvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(MulvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Arguments
    size_t i_x = arg[0];
    size_t i_y = arg[1];

    // number of indices to access
    size_t j = d + 1;
    while (j) {
        --j;
        for (size_t k = 0; k <= j; k++) {
            std::string px_jk = n.generatePartialName(j - k, i_x);
            std::string pz_j = n.generatePartialName(j, i_z);
            std::string sy_k = n.generateVarName(k, i_y);
            std::string sx_jk = n.generateVarName(j - k, i_x);
            std::string py_k = n.generatePartialName(k, i_y);

            s_out << px_jk << " += " << pz_j << " * " << sy_k << n.endl();
            s_out << py_k << " += " << pz_j << " * " << sx_jk << n.endl();
        }
    }
}
// --------------------------- Mulpv -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = MulpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x * y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op
 */

template <class Base>
inline void forward_code_gen_mulpv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(MulpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(MulpvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Parameter value
    Base x = parameter[ arg[0] ];

    // Taylor coefficients corresponding to arguments and result
    size_t i_y = arg[1];

    std::string sz_d = n.generateVarName(d, i_z);
    std::string sy_d = n.generateVarName(d, i_y);

    s_out << sz_d << " = " << n.PrintBase(x) << " * " << sy_d << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = MulpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x * y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_binary_op_0
 */

template <class Base>
inline void forward_code_gen_mulpv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(MulpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(MulpvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Parameter value
    Base x = parameter[ arg[0] ];

    // Taylor coefficients corresponding to arguments and result
    size_t i_y = arg[1];

    std::string sz_0 = n.generateVarName(0, i_z);
    std::string sy_0 = n.generateVarName(0, i_y);

    s_out << sz_0 << " = " << n.PrintBase(x) << " * " << sy_0 << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = MulpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = x * y
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails reverse_binary_op
 */

template <class Base>
inline void reverse_code_gen_mulpv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(MulvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(MulvvOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    size_t i_y = arg[1];

    // Arguments
    std::string x = n.PrintBase(parameter[ arg[0] ]);

    // number of indices to access
    size_t j = d + 1;
    while (j) {
        --j;
        std::string pz_j = n.generatePartialName(j, i_z);
        std::string py_j = n.generatePartialName(j, i_y);

        s_out << py_j << " += " << pz_j << " * " << x << n.endl();

    }
}


CPPAD_END_NAMESPACE
#endif

