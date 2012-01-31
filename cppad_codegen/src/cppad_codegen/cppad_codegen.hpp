#ifndef CPPAD_CODEGEN_CPPAD_CODEGEN_INCLUDED
#define CPPAD_CODEGEN_CPPAD_CODEGEN_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
/*!
\file cppad_codegen.hpp
\brief includes the entire CppAD package in the necessary order for source code generation.

\namespace CppAD
\brief contains all the variables and functions defined by the CppAD package.
 */

#include <cppad/cppad.hpp>

//#include <cppad/local/ad_fun.hpp>

// non-user interfaces
#include <cppad_codegen/local/ad_code_gen_name_provider.hpp>
#include <cppad_codegen/local/ad_fun_code_gen.hpp> 
#include <cppad_codegen/local/op.hpp> // executes taped operations
#include <cppad_codegen/local/forward_code_gen_sweep.hpp>
#include <cppad_codegen/local/forward_code_gen.hpp>
#include <cppad_codegen/local/reverse_code_gen_sweep.hpp>
#include <cppad_codegen/local/reverse_code_gen.hpp>

#include <cppad_codegen/local/sparse_jacobian.hpp>

#endif
