#ifndef CPPAD_CODEGEN_SIGN_INCLUDED
#define	CPPAD_CODEGEN_SIGN_INCLUDED

/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template <class Base>
    inline std::string code_gen_sign(CodeGenNameProvider<Base>& n,
    const std::string& sx) {
        return std::string("(") +
                sx + " > " + n.zero() + "? " + n.CastToBase(1) + " : ("
                + sx + " == " + n.zero() + "?" + n.zero() + " : " + n.CastToBase(-1) +
                ")"
                ")";
    }

    template <class Base>
    inline std::string code_gen_sign(CodeGenNameProvider<Base>& n,
    const std::string& sx_0,
    const std::string& sx_d) {
        return std::string("(")
                + sx_0 + " > " + n.zero() + "? " + sx_d + " : ("
                + sx_0 + " == " + n.zero() + "?" + n.zero() + " : "
                "-" + sx_d +
                ")"
                ")";
    }
}

#endif

