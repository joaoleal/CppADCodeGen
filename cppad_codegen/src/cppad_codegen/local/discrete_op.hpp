#ifndef CPPAD_CODEGEN_DISCRETE_OP_INCLUDED
#define CPPAD_CODEGEN_DISCRETE_OP_INCLUDED

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
\file discrete_op.hpp
Zero order forward mode for z = f(x) where f is piecewise constant.
 */


/*!
Compute zero order forward mode Taylor coefficient for result of op = DisOp.

The C++ source code corresponding to this operation is
\verbatim
        z = f(x)
\endverbatim
where f is a piecewise constant function (and it's derivative is always
calculated as zero).

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base .

\param i_z
variable index corresponding to the result for this operation; 
i.e. the row index in \a taylor corresponding to z. 

\param arg
\a arg[0]
\n
is the index, in the order of the discrete functions defined by the user, 
for this discrete function.
\n
\n
\a arg[1]
variable index corresponding to the argument for this operator;
i.e. the row index in \a taylor corresponding to x.

\param nc_taylor
number of colums in the matrix containing all the Taylor coefficients.

\param taylor
\b Input: \a taylor [ \a arg[1] * \a nc_taylor + 0 ] 
is the zero order Taylor coefficient corresponding to x. 
\n
\b Output: \a taylor [ \a i_z * \a nc_taylor + 0 ] 
is the zero order Taylor coefficient corresponding to z. 

\par Checked Assertions where op is the unary operator with one result:
\li NumArg(op) == 2
\li NumRes(op) == 1
\li \a arg[1] < \a i_z 
\li \a 0 < \a nc_taylor
 */
template<class Base>
inline void forward_code_gen_dis_op(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t d,
        size_t i_z,
        const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DisOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DisOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_d = n.generateVarName(d, arg[1]);
    std::string sz_d = n.generateVarName(d, i_z);

    if (d == 0) {
        s_out << sz_d << " = " << discrete<Base>::name(arg[0]) << "(" << sx_d << ")" << n.endl();
    } else {
        s_out << sz_d << " = " << n.zero() << n.endl();
    }
}

template<class Base>
inline void forward_code_gen_dis_op_0(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t i_z,
        const addr_t* arg) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(DisOp) == 2);
    CPPAD_ASSERT_UNKNOWN(NumRes(DisOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[1]) < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_0 = n.generateVarName(0, arg[1]);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = " << discrete<Base>::name(arg[0]) << "(" << sx_0 << ")" << n.endl();
}


CPPAD_END_NAMESPACE
#endif
