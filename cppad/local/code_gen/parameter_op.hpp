#ifndef CPPAD_CODE_GEN_PARAMETER_OP_INCLUDED
#define CPPAD_CODE_GEN_PARAMETER_OP_INCLUDED

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
\file par_op.hpp
Forward and reverse mode calculations for parameters.
 */


/*!
Compute zero order forward mode Taylor coefficient for result of op = ParOp.

The C++ source code corresponding to this operation is one of the following
\verbatim
        ADFun<Base> f(x, y)
        f.Dependent(x, y)
\endverbatim
where some of the components of the vector y are parameters.

\tparam Base
base type for the operator; i.e., this operation was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base .

\param i_z
variable index corresponding to the result for this operation; 
i.e. the row index in \a taylor corresponding to the component of y
that is a parameter. 

\param arg
\a arg[0]
\n
index corresponding to the parameter value for this operator.

\param num_par
is the number of parameters in \a parameter.

\param parameter
\b Input: \a parameter[ \a arg[0] ] is the value of a component
of y that is a parameter. 

\param nc_taylor
number of colums in the matrix containing all the Taylor coefficients.

\param taylor
\b Output: \a taylor [ \a i_z * \a nc_taylor + 0 ] 
is the zero order Taylor coefficient corresponding to z. 

\par Checked Assertions where op is the unary operator with one result:
\li NumArg(op) == 1
\li NumRes(op) == 1
\li \a size_t(arg[0]) < num_par 
\li \a 0 < \a nc_taylor
 */
template <class Base>
inline void forward_code_gen_par_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
size_t num_par,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(ParOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(ParOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);

    std::string sz_d = n.generateVarName(d, i_z);

    if (d == 0) {
        s_out << sz_d << " = " << n.PrintBase(parameter[ arg[0] ]) << n.endl();
    } else {
        s_out << sz_d << " = " << n.zero() << n.endl();
    }
}

template <class Base>
inline void forward_code_gen_par_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
size_t num_par,
const Base* parameter) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(ParOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(ParOp) == 1);
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < num_par);

    std::string sz_0 = n.generateVarName(0, i_z);
    s_out << sz_0 << " = " << n.PrintBase(parameter[ arg[0] ]) << n.endl();
}

CPPAD_END_NAMESPACE
#endif
