#ifndef CPPAD_CODEGEN_EXP_OP_INCLUDED
#define CPPAD_CODEGEN_EXP_OP_INCLUDED

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
\file exp_op.hpp
Forward and reverse mode calculations for z = exp(x).
 */


/*!
Forward mode Taylor coefficient for result of op = ExpOp.

The C++ source code corresponding to this operation is
\verbatim
        z = exp(x)
\endverbatim

\copydetails forward_unary1_op
 */
template <class Base>
inline void forward_code_gen_exp_op(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t d,
        size_t i_z,
        size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(ExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(ExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sz_d = n.generateVarName(d, i_z);

    if (d == 0) {
        std::string sx_0 = n.generateVarName(0, i_x);
        s_out << sz_d << " = exp(" << sx_0 << ")" << n.endl();
    } else {
        std::string sx_1 = n.generateVarName(1, i_x);
        std::string sz_d1 = n.generateVarName(d - 1, i_z);

        s_out << sz_d << " = (" << sx_1 << " * " << sz_d1;
        for (size_t k = 2; k <= d; k++) {
            std::string sx_k = n.generateVarName(k, i_x);
            std::string sz_dk = n.generateVarName(d - k, i_z);

            s_out << " + " << k << " * " << sx_k << " * " << sz_dk;
        }
        s_out << ") / " << n.CastToBase(d) << n.endl();
    }
}

/*!
Zero order forward mode Taylor coefficient for result of op = ExpOp.

The C++ source code corresponding to this operation is
\verbatim
        z = exp(x)
\endverbatim

\copydetails forward_unary1_op_0
 */
template <class Base>
inline void forward_code_gen_exp_op_0(
        std::ostream& s_out,
        CodeGenNameProvider<Base>& n,
        size_t i_z,
        size_t i_x) {
    // check assumptions
    CPPAD_ASSERT_UNKNOWN(NumArg(ExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(NumRes(ExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(i_x < i_z);

    // Taylor coefficients corresponding to argument and result
    std::string sx_0 = n.generateVarName(0, i_x);
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << sz_0 << " = exp(" << sx_0 << ")" << n.endl();
}

/*!
Reverse mode partial derivatives for result of op = ExpOp.

The C++ source code corresponding to this operation is
\verbatim
        z = exp(x)
\endverbatim

\copydetails reverse_unary1_op
 */

template <class Base>
inline void reverse_code_gen_exp_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
size_t i_x,
size_t nc_taylor,
const Base* taylor,
size_t nc_partial,
Base* partial) {
    //	// check assumptions
    //	CPPAD_ASSERT_UNKNOWN( NumArg(ExpOp) == 1 );
    //	CPPAD_ASSERT_UNKNOWN( NumRes(ExpOp) == 1 );
    //	CPPAD_ASSERT_UNKNOWN( i_x < i_z );
    //	CPPAD_ASSERT_UNKNOWN( d < nc_taylor );
    //	CPPAD_ASSERT_UNKNOWN( d < nc_partial );
    //
    //	// Taylor coefficients and partials corresponding to argument
    //	const Base* x  = taylor  + i_x * nc_taylor;
    //	Base* px       = partial + i_x * nc_partial;
    //
    //	// Taylor coefficients and partials corresponding to result
    //	const Base* z  = taylor  + i_z * nc_taylor;
    //	Base* pz       = partial + i_z * nc_partial;
    //
    //	// lopp through orders in reverse
    //	size_t j, k;
    //	j = d;
    //	while(j)
    //	{	// scale partial w.r.t z[j]
    //		pz[j] /= Base(j);
    //
    //		for(k = 1; k <= j; k++)
    //		{	px[k]   += pz[j] * Base(k) * z[j-k]; 	
    //			pz[j-k] += pz[j] * Base(k) * x[k];
    //		}
    //		--j;
    //	}
    //	px[0] += pz[0] * z[0];
    //        
    throw "not implemented yet";
}

CPPAD_END_NAMESPACE
#endif
