#ifndef CPPAD_CODEGEN_COND_OP_INCLUDED
#define CPPAD_CODEGEN_COND_OP_INCLUDED

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
\file cond_op.hpp
Forward, reverse, and sparse operations for conditional expressions.
 */

inline std::string getComparisonString(CompareOp op) {
    switch (op) {
        case CompareLt:
            return "<";

        case CompareLe:
            return "<=";

        case CompareEq:
            return "==";

        case CompareGe:
            return ">=";

        case CompareGt:
            return ">";

        case CompareNe:
            return "!=";

        default:
            CPPAD_ASSERT_UNKNOWN(0);
    }
}

/*!
Compute forward mode Taylor coefficients for op = CExpOp.

\copydetails conditional_exp_op

\param d
is the order of the Taylor coefficient of z that we are computing.

\param taylor
\b Input:
For j = 0, 1, 2, 3 and k = 0 , ... , \a d,
if y_j is a variable then
\a taylor [ \a arg[2+j] * nc_taylor + k ]
is the k-th order Taylor coefficient corresponding to y_j.
\n
\b Input: \a taylor [ \a i_z * \a nc_taylor + k ] 
for k = 0 , ... , \a d - 1
is the k-th order Taylor coefficient corresponding to z.
\n
\b Output: \a taylor [ \a i_z * \a nc_taylor + \a d ] 
is the \a d-th order Taylor coefficient corresponding to z. 

 */
template <class Base>
inline void forward_code_gen_cond_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
size_t num_par,
const Base* parameter) {
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < static_cast<size_t> (CompareNe));
    CPPAD_ASSERT_UNKNOWN(NumArg(CExpOp) == 6);
    CPPAD_ASSERT_UNKNOWN(NumRes(CExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(arg[1] != 0);

    std::string y_0, y_1, y_2, y_3;

    if (arg[1] & 1) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < i_z);
        y_0 = n.generateVarName(0, arg[2]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < num_par);
        y_0 = n.PrintBase(parameter[ arg[2] ]);
    }
    if (arg[1] & 2) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[3]) < i_z);
        y_1 = n.generateVarName(0, arg[3]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[3]) < num_par);
        y_1 = n.PrintBase(parameter[ arg[3] ]);
    }
#if CPPAD_USE_FORWARD0SWEEP
    CPPAD_ASSERT_UNKNOWN(d > 0);
#else
    if (d == 0) {
        if (arg[1] & 4) {
            CPPAD_ASSERT_UNKNOWN(size_t(arg[4]) < i_z);
            y_2 = n.generateVarName(0, arg[4]);
        } else {
            CPPAD_ASSERT_UNKNOWN(size_t(arg[4]) < num_par);
            y_2 = n.PrintBase(parameter[ arg[4] ]);
        }
        if (arg[1] & 8) {
            CPPAD_ASSERT_UNKNOWN(size_t(arg[5]) < i_z);
            y_3 = n.generateVarName(0, arg[5]);
        } else {
            CPPAD_ASSERT_UNKNOWN(size_t(arg[5]) < num_par);
            y_3 = n.PrintBase(parameter[ arg[5] ]);
        }
    } else
#endif
    {
        if (arg[1] & 4) {
            CPPAD_ASSERT_UNKNOWN(size_t(arg[4]) < i_z);
            y_2 = n.generateVarName(d, arg[4]);
        } else {
            y_2 = n.zero();
        }
        if (arg[1] & 8) {
            CPPAD_ASSERT_UNKNOWN(size_t(arg[5]) < i_z);
            y_3 = n.generateVarName(d, arg[5]);
        } else {
            y_3 = n.zero();
        }
    }

    std::string opCode = getComparisonString(CompareOp(arg[0]));
    std::string sz_d = n.generateVarName(d, i_z);

    s_out << "if(" << y_0 << " " << opCode << " " << y_1 << ") "
            << sz_d << " = " << y_2 << n.endl()
            << " else "
            << sz_d << " = " << y_3 << n.endl();
}

/*!
Compute zero order forward mode Taylor coefficients for op = CExpOp.

\copydetails conditional_exp_op

\param taylor
\b Input:
For j = 0, 1, 2, 3,
if y_j is a variable then
\a taylor [ \a arg[2+j] * nc_taylor + 0 ]
is the zero order Taylor coefficient corresponding to y_j.
\n
\b Output: \a taylor [ \a i_z * \a nc_taylor + 0 ] 
is the zero order Taylor coefficient corresponding to z. 
 */
template <class Base>
inline void forward_code_gen_cond_op_0(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t i_z,
const addr_t* arg,
size_t num_par,
const Base* parameter) {
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < static_cast<size_t> (CompareNe));
    CPPAD_ASSERT_UNKNOWN(NumArg(CExpOp) == 6);
    CPPAD_ASSERT_UNKNOWN(NumRes(CExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(arg[1] != 0);

    std::string y_0, y_1, y_2, y_3;

    if (arg[1] & 1) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < i_z);
        y_0 = n.generateVarName(0, arg[2]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < num_par);
        y_0 = n.PrintBase(parameter[ arg[2] ]);
    }
    if (arg[1] & 2) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[3]) < i_z);
        y_1 = n.generateVarName(0, arg[3]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[3]) < num_par);
        y_1 = n.PrintBase(parameter[ arg[3] ]);
    }
    if (arg[1] & 4) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[4]) < i_z);
        y_2 = n.generateVarName(0, arg[4]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[4]) < num_par);
        y_2 = n.PrintBase(parameter[ arg[4] ]);
    }
    if (arg[1] & 8) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[5]) < i_z);
        y_3 = n.generateVarName(0, arg[5]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[5]) < num_par);
        y_3 = n.PrintBase(parameter[ arg[5] ]);
    }

    std::string opCode = getComparisonString(CompareOp(arg[0]));
    std::string sz_0 = n.generateVarName(0, i_z);

    s_out << "if(" << y_0 << " " << opCode << " " << y_1 << ") "
            << sz_0 << " = " << y_2 << n.endl()
            << " else "
            << sz_0 << " = " << y_3 << n.endl();
}

/*!
Compute reverse mode Taylor coefficients for op = CExpOp.

\copydetails conditional_exp_op

This routine is given the partial derivatives of a function 
G( z , y , x , w , ... )
and it uses them to compute the partial derivatives of 
\verbatim
        H( y , x , w , u , ... ) = G[ z(y) , y , x , w , u , ... ]
\endverbatim
where y above represents y_0, y_1, y_2, y_3.

\param d
is the order of the Taylor coefficient of z that we are  computing.

\param taylor
\b Input:
For j = 0, 1, 2, 3 and k = 0 , ... , \a d,
if y_j is a variable then
\a taylor [ \a arg[2+j] * nc_taylor + k ]
is the k-th order Taylor coefficient corresponding to y_j.
\n
\a taylor [ \a i_z * \a nc_taylor + k ] 
for k = 0 , ... , \a d
is the k-th order Taylor coefficient corresponding to z.

\param nc_partial
number of columns in the matrix containing the Taylor coefficients.

\param partial
\b Input:
For j = 0, 1, 2, 3 and k = 0 , ... , \a d,
if y_j is a variable then
\a partial [ \a arg[2+j] * nc_partial + k ]
is the partial derivative of G( z , y , x , w , u , ... )
with respect to the k-th order Taylor coefficient corresponding to y_j.
\n
\b Input: \a partial [ \a i_z * \a nc_taylor + k ] 
for k = 0 , ... , \a d
is the partial derivative of G( z , y , x , w , u , ... )
with respect to the k-th order Taylor coefficient corresponding to z.
\n
\b Output:
For j = 0, 1, 2, 3 and k = 0 , ... , \a d,
if y_j is a variable then
\a partial [ \a arg[2+j] * nc_partial + k ]
is the partial derivative of H( y , x , w , u , ... )
with respect to the k-th order Taylor coefficient corresponding to y_j.

 */
template <class Base>
inline void reverse_code_gen_cond_op(
std::ostream& s_out,
CodeGenNameProvider<Base>& n,
size_t d,
size_t i_z,
const addr_t* arg,
size_t num_par,
const Base* parameter) {
    CPPAD_ASSERT_UNKNOWN(size_t(arg[0]) < static_cast<size_t> (CompareNe));
    CPPAD_ASSERT_UNKNOWN(NumArg(CExpOp) == 6);
    CPPAD_ASSERT_UNKNOWN(NumRes(CExpOp) == 1);
    CPPAD_ASSERT_UNKNOWN(arg[1] != 0);


    std::string y_0, y_1;

    std::string opCode = getComparisonString(CompareOp(arg[0]));
    if (arg[1] & 1) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < i_z);
        y_0 = n.generateVarName(0, arg[2]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[2]) < num_par);
        y_0 = n.PrintBase(parameter[ arg[2] ]);
    }
    if (arg[1] & 2) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[3]) < i_z);
        y_1 = n.generateVarName(0, arg[3]);
    } else {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[3]) < num_par);
        y_1 = n.PrintBase(parameter[ arg[3] ]);
    }
    if (arg[1] & 4) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[4]) < i_z);
        size_t i_y2 = arg[4];
        s_out << "if(" << y_0 << " " << opCode << " " << y_1 << ") {\n";
        size_t j = d + 1;
        while (j--) {
            std::string py2_j = n.generatePartialName(j, i_y2);
            std::string pz_j = n.generatePartialName(j, i_z);
            s_out << py2_j << " += " << pz_j << n.endl();
        }
        s_out << "} else {\n";
        j = d + 1;
        while (j--) {
            std::string py2_j = n.generatePartialName(j, i_y2);
            s_out << py2_j << " += " << n.zero() << n.endl();
        }
        s_out << "}\n";
    }
    if (arg[1] & 8) {
        CPPAD_ASSERT_UNKNOWN(size_t(arg[5]) < i_z);
        size_t i_y3 = arg[5];

        s_out << "if(" << y_0 << " " << opCode << " " << y_1 << ") {\n";
        size_t j = d + 1;
        while (j--) {
            std::string py3_j = n.generatePartialName(j, i_y3);
            s_out << py3_j << " += " << n.zero() << n.endl();
        }
        s_out << "} else {\n";
        j = d + 1;
        while (j--) {
            std::string py3_j = n.generatePartialName(j, i_y3);
            std::string pz_j = n.generatePartialName(j, i_z);
            s_out << py3_j << " += " << pz_j << n.endl();
        }
        s_out << "}\n";
    }
}

CPPAD_END_NAMESPACE
#endif
