#ifndef CPPAD_CODEGEN_REVERSE_CODEGEN_INCLUDED
#define CPPAD_CODEGEN_REVERSE_CODEGEN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

/*
$begin Reverse$$
$spell
        typename
        xk
        xp
        dw
        Ind
        uj
        std
        arg
        const
        taylor_
$$

$section Reverse Mode$$ 

$childtable%
        omh/reverse.omh
%$$

$end
-----------------------------------------------------------------------------
 */
#include <algorithm>
#include <cppad/local/pod_vector.hpp>

CPPAD_BEGIN_NAMESPACE
/*!
\file reverse.hpp
Compute derivatives using reverse mode.
 */


/*!
Use reverse mode to compute derivative of forward mode Taylor coefficients.

The function 
\f$ X : {\rm R} \times {\rm R}^{n \times p} \rightarrow {\rm R} \f$ 
is defined by
\f[
X(t , u) = \sum_{k=0}^{p-1} u^{(k)} t^k
\f]
The function 
\f$ Y : {\rm R} \times {\rm R}^{n \times p} \rightarrow {\rm R} \f$ 
is defined by
\f[
Y(t , u) = F[ X(t, u) ]
\f]
The function 
\f$ W : {\rm R}^{n \times p} \rightarrow {\rm R} \f$ is defined by
\f[
W(u) = \sum_{k=0}^{p-1} ( w^{(k)} )^{\rm T} 
\frac{1}{k !} \frac{ \partial^k } { t^k } Y(0, u)
\f]

\tparam Base
base type for the operator; i.e., this operation sequence was recorded
using AD< \a Base > and computations by this routine are done using type 
\a Base.

\tparam VectorBase
is a Simple Vector class with elements of type \a Base.

\param p
is the number of the number of Taylor coefficients that are being
differentiated (per variable).

\param w
is the weighting for each of the Taylor coefficients corresponding
to dependent variables.
If the argument \a w has size <tt>m * p </tt>,
for \f$ k = 0 , \ldots , p-1 \f$ and \f$ i = 0, \ldots , m-1 \f$,
\f[
        w_i^{(k)} = w [ i * p + k ]
\f]
If the argument \a w has size \c m ,
for \f$ k = 0 , \ldots , p-1 \f$ and \f$ i = 0, \ldots , m-1 \f$,
\f[
w_i^{(k)} = \left\{ \begin{array}{ll}
        w [ i ] & {\rm if} \; k = p-1
        \\
        0       & {\rm otherwise}
\end{array} \right.
\f]

\return
Is a vector \f$ dw \f$ such that
for \f$ j = 0 , \ldots , n-1 \f$ and
\f$ k = 0 , \ldots , p-1 \f$
\f[ 
        dw[ j * p + k ] = W^{(1)} ( x )_{j,k}
\f]
where the matrix \f$ x \f$ is the value for \f$ u \f$
that corresponding to the forward mode Taylor coefficients
for the independent variables as specified by previous calls to Forward.

 */
template <typename Base>
template <typename VectorBase>
void ADFunCodeGen<Base>::ReverseCodeGen(size_t p,
        const VectorBase &w,
        std::ostream& s_out) {
    // temporary indices
    size_t i, j, k;

    // number of independent variables
    size_t n = ADFun<Base>::Domain();

    // number of dependent variables
    size_t m = ADFun<Base>::Range();

    player<Base>& play = ADFun<Base>::Play();

    size_t total_num_var = ADFun<Base>::size_var();

    // update maximum memory requirement
    // memoryMax = std::max( memoryMax, 
    // 	Memory() + total_num_var_ * p * sizeof(Base)
    // );

    CPPAD_ASSERT_KNOWN(
            w.size() == m || w.size() == (m * p),
            "Argument w to Reverse does not have length equal to\n"
            "the dimension of the range for the corresponding ADFun."
            );
    CPPAD_ASSERT_KNOWN(
            p > 0,
            "The first argument to Reverse must be greater than zero."
            );
    CPPAD_ASSERT_KNOWN(
            taylor_per_var_ >= p,
            "Less that p taylor_ coefficients are currently stored"
            " in this ADFun object."
            );

    // initialize entire Partial matrix to zero
    for (i = 0; i < total_num_var; i++) {
        for (j = 0; j < p; j++) {
            std::string partial = nameGen_->generatePartialName(j, i);
            s_out << partial << " = " << nameGen_->zero() << nameGen_->endl();
        }
    }

    const CppAD::vector<size_t>& dep_taddr = ADFun<Base>::DependentTapeAddr();
    const CppAD::vector<size_t>& ind_taddr = ADFun<Base>::IndependentTapeAddr();

    // set the dependent variable direction
    // (use += because two dependent variables can point to same location)
    for (i = 0; i < m; i++) {
        CPPAD_ASSERT_UNKNOWN(dep_taddr[i] < total_num_var);
        if (w.size() == m) {
            std::string partial = nameGen_->generatePartialName(p - 1, dep_taddr[i]);
            s_out << partial << " = " << nameGen_->PrintBase(w[i]) << nameGen_->endl();
            //Partial[dep_taddr[i] * p + p - 1] += w[i];
        } else {
            for (k = 0; k < p; k++) {
                // ? should use += here, first make test to demonstrate bug
                std::string partial = nameGen_->generatePartialName(k, dep_taddr[i]);
                s_out << partial << " = " << nameGen_->PrintBase(w[i * p + k]) << nameGen_->endl();
                //Partial[ dep_taddr[i] * p + k ] = w[i * p + k ];
            }
        }
    }

    // evaluate the derivatives
    reverse_code_gen_sweep(s_out, *nameGen_, p - 1, n, total_num_var, &play);

    for (j = 0; j < n; j++) {
        CPPAD_ASSERT_UNKNOWN(ind_taddr[j] < total_num_var);

        // independent variable taddr equals its operator taddr 
        CPPAD_ASSERT_UNKNOWN(play.GetOp(ind_taddr[j]) == InvOp);

        // by the Reverse Identity Theorem 
        // partial of y^{(k)} w.r.t. u^{(0)} is equal to
        // partial of y^{(p-1)} w.r.t. u^{(p - 1 - k)}

        //        if (w.size() == m) {
        //            for (k = 0; k < p; k++) {
        //                std::string partial = nameGen_->generatePartialName(p - 1 - k, ind_taddr[j]);
        //                //value[j * p + k ] = Partial[ind_taddr[j] * p + p - 1 - k];
        //            }
        //        } else {
        //            for (k = 0; k < p; k++) {
        //                std::string partial = nameGen_->generatePartialName(k, ind_taddr[j]);
        //                //value[j * p + k ] = Partial[ind_taddr[j] * p + k];
        //            }
        //        }
    }
}


} // END CppAD namespace


#endif
