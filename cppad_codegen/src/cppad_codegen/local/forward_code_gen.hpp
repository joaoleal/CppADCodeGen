/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#ifndef CPPAD_CODEGEN_FORWARD_CODE_GEN_INCLUDED
#define	CPPAD_CODEGEN_FORWARD_CODE_GEN_INCLUDED

#include "ad_fun_code_gen.hpp"

CPPAD_BEGIN_NAMESPACE

template <typename Base>
void ADFunCodeGen<Base>::ForwardCodeGen(
        size_t d,
        std::ostream& s_out) {

    // number of independent variables
    size_t n = ADFun<Base>::Domain();

    // number of dependent variables is not used
    // size_t m = ADFun<Base>::Range();

    player<Base>& play = ADFun<Base>::Play();

    size_t total_num_var = ADFun<Base>::size_var();

    const CppAD::vector<size_t>& ind_taddr = ADFun<Base>::IndependentTapeAddr();

    //    if (zeroTaylor_.capacity() == 0) {
    //        zeroTaylor_ = Matrix(d, total_num_var);
    //    } else {
    //        zeroTaylor_.resizeRows(d);
    //        CPPAD_ASSERT_UNKNOWN(zeroTaylor_.columns() == total_num_var);
    //    }

    CPPAD_ASSERT_KNOWN(
            d <= taylor_per_var_,
            "The number of taylor_ coefficient currently generated\n"
            "is less than d."
            );

    // set the p-th order taylor_ coefficients for independent variables
    for (size_t j = 0; j < n; j++) {
        CPPAD_ASSERT_UNKNOWN(ind_taddr[j] < total_num_var);

        // ind_taddr_[j] is operator taddr for j-th independent variable
        CPPAD_ASSERT_UNKNOWN(play.GetOp(ind_taddr[j]) == InvOp);

        // It is also variable taddr for j-th independent variable
        //taylor_[ind_taddr_[j] * taylor_col_dim_ + p] = x_p[j];
    }

    s_out << nameGen_->Comment("forward sweep (order: " + nameGen_->toString(d) + ")");

    // generate the code
    forward_code_gen_sweep(s_out, *nameGen_, false, d, n, total_num_var, &play);

    // now we have p + 1  taylor_ coefficients per variable
    taylor_per_var_ = d + 1;
}

CPPAD_END_NAMESPACE
#endif	/* FORWARDCODEGEN_H */

