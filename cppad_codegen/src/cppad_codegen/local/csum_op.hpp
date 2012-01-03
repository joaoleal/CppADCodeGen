#ifndef CPPAD_CODEGEN_CSUM_OP_INCLUDED
#define CPPAD_CODEGEN_CSUM_OP_INCLUDED

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
\file csum_op.hpp
Forward, reverse and sparsity calculations for cummulative summation.
 */

/*!
Compute forward mode Taylor coefficients for result of op = CsumOp.

This operation is 
\verbatim
        z = p + x(1) + ... + x(m) - y(1) - ... - y(n).
\endverbatim

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type
\a Base.

\param d
order of the Taylor coefficient that we are computing.

\param i_z
variable index corresponding to the result for this operation;
i.e. the row index in \a taylor corresponding to z.

\param arg
\a arg[0] 
is the number of addition variables in this cummulative summation; i.e.,
<tt>m</tt>.
\n
\a arg[1] 
is the number of subtraction variables in this cummulative summation; i.e.,
\c m.
\n
<tt>parameter[ arg[2] ]</tt>
is the parameter value \c p in this cummunative summation.
\n
<tt>arg[2+i]</tt>
for <tt>i = 1 , ... , m</tt> is the value <tt>x(i)</tt>. 
\n
<tt>arg[2+arg[0]+i]</tt>
for <tt>i = 1 , ... , n</tt> is the value <tt>y(i)</tt>. 

\param num_par
is the number of parameters in \a parameter.

\param parameter
is the parameter vector for this operation sequence.

\param nc_taylor
number of colums in the matrix containing all the Taylor coefficients.

\param taylor
\b Input: <tt>taylor [ arg[2+i] * nc_taylor + k ]</tt>
for <tt>i = 1 , ... , m</tt> 
and <tt>k = 0 , ... , d</tt>
is the k-th order Taylor coefficient corresponding to <tt>x(i)</tt>
\n
\b Input: <tt>taylor [ arg[2+m+i] * nc_taylor + k ]</tt>
for <tt>i = 1 , ... , n</tt> 
and <tt>k = 0 , ... , d</tt>
is the k-th order Taylor coefficient corresponding to <tt>y(i)</tt>
\n
\b Input: <tt>taylor [ i_z * nc_taylor + k ]</tt>
for k = 0 , ... , \a d - 1
is the k-th order Taylor coefficient corresponding to z.
\n
\b Output: <tt>taylor [ i_z * nc_taylor + d ]</tt>
is the \a d-th order Taylor coefficient corresponding to z.
 */
template <class Base>
inline void forward_code_gen_csum_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
size_t num_par,
const Base* parameter) {

    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumRes(CSumOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < num_par);
    CPPAD_ASSERT_UNKNOWN(arg[0] + arg[1] == arg[ arg[0] + arg[1] + 3 ]);

    // Taylor coefficients corresponding to result
    std::string sz_d = n.generateVarName(d, i_z);
    if (d == 0) {
        s_out << sz_d << " = " << n.PrintBase(parameter[ arg[2] ]);
    } else {
        s_out << sz_d << " = " << n.zero();
    }

    size_t i = arg[0];
    size_t j = 2;
    while (i--) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[j + 1]) < i_z);
        std::string sx_d = n.generateVarName(d, arg[++j]);
        s_out << " + " << sx_d;
    }
    i = arg[1];
    while (i--) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[j + 1]) < i_z);
        std::string sx_d = n.generateVarName(d, arg[++j]);
        s_out << " - " << sx_d;
    }
    s_out << n.endl();
}

/*!
Compute reverse mode Taylor coefficients for result of op = CsumOp.

This operation is 
\verbatim
        z = p + x(1) + ... + x(m) - y(1) - ... - y(n).
        H(y, x, w, ...) = G[ z(x, y), y, x, w, ... ] 
\endverbatim

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type
\a Base.

\param d
order the highest order Taylor coefficient that we are computing
the partial derivatives with respect to.

\param i_z
variable index corresponding to the result for this operation;
i.e. the row index in \a taylor corresponding to z.

\param arg
\a arg[0] 
is the number of addition variables in this cummulative summation; i.e.,
<tt>m</tt>.
\n
\a arg[1] 
is the number of subtraction variables in this cummulative summation; i.e.,
\c m.
\n
<tt>parameter[ arg[2] ]</tt>
is the parameter value \c p in this cummunative summation.
\n
<tt>arg[2+i]</tt>
for <tt>i = 1 , ... , m</tt> is the value <tt>x(i)</tt>. 
\n
<tt>arg[2+arg[0]+i]</tt>
for <tt>i = 1 , ... , n</tt> is the value <tt>y(i)</tt>. 

\param nc_partial
number of colums in the matrix containing all the partial derivatives.

\param partial
\b Input: <tt>partial [ arg[2+i] * nc_partial + k ]</tt>
for <tt>i = 1 , ... , m</tt> 
and <tt>k = 0 , ... , d</tt>
is the partial derivative of G(z, y, x, w, ...) with respect to the
k-th order Taylor coefficient corresponding to <tt>x(i)</tt>
\n
\b Input: <tt>partial [ arg[2+m+i] * nc_partial + k ]</tt>
for <tt>i = 1 , ... , n</tt> 
and <tt>k = 0 , ... , d</tt>
is the partial derivative of G(z, y, x, w, ...) with respect to the
k-th order Taylor coefficient corresponding to <tt>y(i)</tt>
\n
\b Input: <tt>partial [ i_z * nc_partial + k ]</tt>
for <tt>i = 1 , ... , n</tt> 
and <tt>k = 0 , ... , d</tt>
is the partial derivative of G(z, y, x, w, ...) with respect to the
k-th order Taylor coefficient corresponding to \c z.
\n
\b Output: <tt>partial [ arg[2+i] * nc_partial + k ]</tt>
for <tt>i = 1 , ... , m</tt> 
and <tt>k = 0 , ... , d</tt>
is the partial derivative of H(y, x, w, ...) with respect to the
k-th order Taylor coefficient corresponding to <tt>x(i)</tt>
\n
\b Output: <tt>partial [ arg[2+m+i] * nc_partial + k ]</tt>
for <tt>i = 1 , ... , n</tt> 
and <tt>k = 0 , ... , d</tt>
is the partial derivative of H(y, x, w, ...) with respect to the
k-th order Taylor coefficient corresponding to <tt>y(i)</tt>
 */

template <class Base>
inline void reverse_code_gen_csum_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumRes(CSumOp) == 1);

    size_t i, j, k;
    size_t d1 = d + 1;

    i = arg[0];
    j = 2;
    while (i--) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[j + 1]) < i_z);
        size_t i_x = arg[++j];
        k = d1;
        while (k--) {
            std::string px_k = n.generatePartialName(k, i_x);
            std::string pz_k = n.generatePartialName(k, i_z);

            s_out << px_k << " += " << pz_k << n.endl();
        }
    }

    i = arg[1];
    while (i--) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[j + 1]) < i_z);
        size_t i_x = arg[++j];
        k = d1;
        while (k--) {
            std::string px_k = n.generatePartialName(k, i_x);
            std::string pz_k = n.generatePartialName(k, i_z);

            s_out << px_k << " -= " << pz_k << n.endl();
        }
    }
}

CPPAD_END_NAMESPACE
#endif
