#ifndef CPPAD_CODEGEN_POW_OP_INCLUDED
#define CPPAD_CODEGEN_POW_OP_INCLUDED

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
\file pow_op.hpp
Forward and reverse mode calculations for z = pow(x, y).
 */

// --------------------------- Powvv -----------------------------------------
/*!
Compute forward mode Taylor coefficients for result of op = PowvvOp.

In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_pow_op
 */
template <class Base>
inline void forward_code_gen_powvv_op(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t d,
        size_t i_z,
        const addr_t* arg) {
    // convert from final result to first result
    i_z -= 2; // NumRes(PowvvOp) - 1;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(PowvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(PowvvOp) == 3);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // z_0 = log(x)
    forward_code_gen_log_op(s_out, n, d, i_z, arg[0]);

    // z_1 = z_0 * y
    addr_t adr[2];
    adr[0] = i_z;
    adr[1] = arg[1];
    forward_code_gen_mulvv_op(s_out, n, d, i_z + 1, adr);

    // z_2 = exp(z_1)
#if CPPAD_USE_FORWARD0SWEEP
    CPPAD_ASSERT_UNKNOWN(d > 0);
    forward_code_gen_exp_op(s_out, n, d, i_z + 2, i_z + 1);
#else
    // final result for zero order case is exactly the same as for Base
    if (d == 0) {
        // Taylor coefficients corresponding to arguments and result
        std::string sx_0 = n.generateVarName(0, arg[0]);
        std::string sy_0 = n.generateVarName(0, arg[1]);
        std::string sz2_0 = n.generateVarName(0, i_z + 2);

        s_out << sz2_0 << " = pow(" << sx_0 << ", " << sy_0 << ")" << n.endl();
    } else {
        forward_code_gen_exp_op(s_out, n, d, i_z + 2, i_z + 1);
    }
#endif
}

/*!
Compute zero order forward mode Taylor coefficients for result of op = PowvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails forward_pow_op_0
 */
template <class Base>
inline void forward_code_gen_powvv_op_0(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t i_z,
        const addr_t* arg) {
    // convert from final result to first result
    i_z -= 2; // NumRes(PowvvOp) - 1;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(PowvvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(PowvvOp) == 3);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sx_0 = n.generateVarName(0, arg[0]);
    std::string sy_0 = n.generateVarName(0, arg[1]);
    std::string sz0_0 = n.generateVarName(0, i_z);
    std::string sz1_0 = n.generateVarName(0, i_z + 1);
    std::string sz2_0 = n.generateVarName(0, i_z + 2);

    s_out << sz0_0 << " = log(" << sx_0 << ")" << n.endl();
    s_out << sz1_0 << " = " << sz0_0 << " * " << sy_0 << n.endl();
    s_out << sz2_0 << " = pow(" << sx_0 << ", " << sy_0 << ")" << n.endl();
}

/*!
Compute reverse mode partial derivatives for result of op = PowvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where both x and y are variables
and the argument \a parameter is not used.

\copydetails reverse_pow_op
 */

template <class Base>
inline void reverse_code_gen_powvv_op(
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
    //	// convert from final result to first result
    //	i_z -= 2; // NumRes(PowvvOp) - 1;
    //
    //	// check assumptions
    //	CPPAD_ASSERT_UNKNOWN( NumArg(PowvvOp) == 2 );
    //	CPPAD_ASSERT_UNKNOWN( NumRes(PowvvOp) == 3 );
    //	CPPAD_ASSERT_UNKNOWN( size_t(arg[0]) < i_z );
    //	CPPAD_ASSERT_UNKNOWN( size_t(arg[1]) < i_z );
    //	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
    //	CPPAD_ASSERT_UNKNOWN( d < nc_partial );
    //
    //	// z_2 = exp(z_1)
    //	reverse_exp_op(
    //		d, i_z+2, i_z+1, nc_taylor, taylor, nc_partial, partial
    //	);
    //
    //	// z_1 = z_0 * y
    //	addr_t adr[2];
    //	adr[0] = i_z;
    //	adr[1] = arg[1];
    //	reverse_mulvv_op(
    //	d, i_z+1, adr, parameter, nc_taylor, taylor, nc_partial, partial
    //	);
    //
    //	// z_0 = log(x)
    //	reverse_log_op(
    //		d, i_z, arg[0], nc_taylor, taylor, nc_partial, partial
    //	);
    throw "not implemented yet";
}

// --------------------------- Powpv -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = PowpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_pow_op
 */

template <class Base>
inline void forward_code_gen_powpv_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // convert from final result to first result
    i_z -= 2; // NumRes(PowpvOp) - 1;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(PowpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(PowpvOp) == 3);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to arguments and result
    std::string sz0_0 = n.generateVarName(0, i_z);
    std::string sz0_d = n.generateVarName(d, i_z);
    std::string sz1_d = n.generateVarName(d, i_z + 1);
    std::string sy_d = n.generateVarName(d, arg[1]);

    // z_0 = log(x)
    Base x = parameter[ arg[0] ];
    if (d == 0) {
        s_out << sz0_0 << " = log(" << n.PrintBase(x) << ")" << n.endl();
    } else {
        s_out << sz0_d << " = " << n.zero() << n.endl();
    }

    // z_1 = z_0 * y
    s_out << sz1_d << " = " << sz0_0 << " * " << sy_d << n.endl();
    /**
     * todo: CHECK THIS!!!!!!!!!!!!!
     */

    // z_2 = exp(z_1)
#if CPPAD_USE_FORWARD0SWEEP
    forward_code_gen_exp_op(s_out, n, d, i_z + 2, i_z + 1);
#else
    // zero order case exactly same as Base type operation
    if (d == 0) {
        std::string sz2_0 = n.generateVarName(0, i_z + 2);
        std::string sy_0 = n.generateVarName(0, arg[1]);

        s_out << sz2_0 << " = pow(" << x << ", " << sy_0 << ")" << n.endl();
    } else {
        forward_code_gen_exp_op(s_out, n, d, i_z + 2, i_z + 1);
    }
#endif
}

/*!
Compute zero order forward mode Taylor coefficient for result of op = PowpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails forward_pow_op_0
 */

template <class Base>
inline void forward_code_gen_powpv_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // convert from final result to first result
    i_z -= 2; // NumRes(PowpvOp) - 1;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(PowpvOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(PowpvOp) == 3);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Paraemter value
    Base x = parameter[ arg[0] ];

    // Taylor coefficients corresponding to arguments and result
    std::string sy_0 = n.generateVarName(0, arg[1]);
    std::string sz0_0 = n.generateVarName(0, i_z);
    std::string sz1_0 = n.generateVarName(0, i_z + 1);
    std::string sz2_0 = n.generateVarName(0, i_z + 2);

    // z_0 = log(x)
    s_out << sz0_0 << " = log(" << n.PrintBase(x) << ")" << n.endl();

    // z_1 = z_0 * y
    s_out << sz1_0 << " = " << sz0_0 << " * " << sy_0 << n.endl();

    // z_2 = exp(z_1)
    // zero order case exactly same as Base type operation
    s_out << sz2_0 << " = pow(" << n.PrintBase(x) << ", " << sy_0 << ")" << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = PowpvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where x is a parameter and y is a variable.

\copydetails reverse_pow_op
 */

template <class Base>
inline void reverse_code_gen_powpv_op(
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
    //    // convert from final result to first result
    //    i_z -= 2; // NumRes(PowpvOp) - 1;
    //
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(PowvvOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(PowvvOp) == 3);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // z_2 = exp(z_1)
    //    reverse_exp_op(
    //            d, i_z + 2, i_z + 1, nc_taylor, taylor, nc_partial, partial
    //            );
    //
    //    // z_1 = z_0 * y
    //    addr_t adr[2];
    //    adr[0] = i_z * nc_taylor; // offset of z_0[0] in taylor 
    //    adr[1] = arg[1]; // index of y in taylor and partial
    //    // use taylor both for parameter and variable values
    //    reverse_mulpv_op(
    //            d, i_z + 1, adr, taylor, nc_taylor, taylor, nc_partial, partial
    //            );
    //
    //    // z_0 = log(x)
    //    // x is a parameter
    throw "not implemented yet";
}

// --------------------------- Powvp -----------------------------------------

/*!
Compute forward mode Taylor coefficients for result of op = PowvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_pow_op
 */

template <class Base>
inline void forward_code_gen_powvp_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // convert from final result to first result
    i_z -= 2; // NumRes(PowvpOp) - 1;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(PowvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(PowvpOp) == 3);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    // z_0 = log(x)
    forward_code_gen_log_op(s_out, n, d, i_z, arg[0]);

    // z_1 = y * z_0
    addr_t adr[2];
    adr[0] = arg[1];
    adr[1] = i_z;
    forward_code_gen_mulpv_op(s_out, n, d, i_z + 1, adr, parameter);

    // z_2 = exp(z_1)
    // zero order case exactly same as Base type operation
#if CPPAD_USE_FORWARD0SWEEP
    CPPAD_ASSERT_UNKNOWN(d > 0);
    forward_code_gen_exp_op(s_out, n, d, i_z + 2, i_z + 1);
#else
    if (d == 0) {
        Base y = parameter[ arg[1] ];
        std::string sz2_0 = n.generateVarName(0, i_z + 2);
        std::string sx_0 = n.generateVarName(0, arg[0]);

        s_out << sz2_0 << " = pow(" << sx_0 << ", " << n.PrintBase(y) << ")" << n.endl();

    } else {
        forward_code_gen_exp_op(s_out, n, d, i_z + 2, i_z + 1);
    }
#endif

}

/*!
Compute zero order forward mode Taylor coefficients for result of op = PowvvOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails forward_pow_op_0
 */

template <class Base>
inline void forward_code_gen_powvp_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
const Base* parameter) {
    // convert from final result to first result
    i_z -= 2; // NumRes(PowvpOp) - 1;

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(PowvpOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(PowvpOp) == 3);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);

    // Parameter value
    Base y = parameter[ arg[1] ];

    // Taylor coefficients corresponding to arguments and result
    std::string sx_0 = n.generateVarName(0, arg[0]);
    std::string sz0_0 = n.generateVarName(0, i_z);
    std::string sz1_0 = n.generateVarName(0, i_z + 1);
    std::string sz2_0 = n.generateVarName(0, i_z + 2);

    // z_0 = log(x)
    s_out << sz0_0 << " = log(" << sx_0 << ")" << n.endl();

    // z_1 = z_0 * y
    s_out << sz1_0 << " = " << sz0_0 << " * " << n.PrintBase(y) << n.endl();

    // z_2 = exp(z_1)
    // zero order case exactly same as Base type operation
    s_out << sz2_0 << " = pow(" << sx_0 << ", " << n.PrintBase(y) << ")" << n.endl();
}

/*!
Compute reverse mode partial derivative for result of op = PowvpOp.

The C++ source code corresponding to this operation is
\verbatim
        z = pow(x, y)
\endverbatim
In the documentation below,
this operations is for the case where x is a variable and y is a parameter.

\copydetails reverse_pow_op
 */

template <class Base>
inline void reverse_code_gen_powvp_op(
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
    //    // convert from final result to first result
    //    i_z -= 2; // NumRes(PowvpOp) - 1;
    //
    //    // check assumptions
    //    CPPAD_ASSERT_UNKNOWN(NumArg(PowvpOp) == 2);
    //    CPPAD_ASSERT_UNKNOWN(NumRes(PowvpOp) == 3);
    //    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < i_z);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_taylor);
    //    CPPAD_ASSERT_UNKNOWN(d < nc_partial);
    //
    //    // z_2 = exp(z_1)
    //    reverse_exp_op(
    //            d, i_z + 2, i_z + 1, nc_taylor, taylor, nc_partial, partial
    //            );
    //
    //    // z_1 = y * z_0
    //    addr_t adr[2];
    //    adr[0] = arg[1];
    //    adr[1] = i_z;
    //    reverse_mulpv_op(
    //            d, i_z + 1, adr, parameter, nc_taylor, taylor, nc_partial, partial
    //            );
    //
    //    // z_0 = log(x)
    //    reverse_log_op(
    //            d, i_z, arg[0], nc_taylor, taylor, nc_partial, partial
    //            );
    throw "not implemented yet";
}

CPPAD_END_NAMESPACE
#endif
