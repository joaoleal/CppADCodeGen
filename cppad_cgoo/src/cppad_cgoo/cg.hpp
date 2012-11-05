#ifndef CPPAD_CG_INCLUDED
#define	CPPAD_CG_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
// all base type requirements
#include <cppad/base_require.hpp> 

// --------------------------------------------------------------------------
// System routines that can be used by rest of CppAD with out including 

#include <algorithm>
#include <assert.h>
#include <cstddef>
#include <errno.h>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <set>
#include <stddef.h>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <string.h>

// ---------------------------------------------------------------------------
// definitions needed by rest of includes

// definitions that come from the installation
#include <cppad/configure.hpp>

// definitions that are local to the CppAD include files
#include <cppad/local/define.hpp>

// vectors used with CppAD
#include <cppad/local/testvector.hpp>

// deprecated vectors used with CppAD
#include <cppad/local/test_vector.hpp>

// Declare classes and fucntions that are used before defined
#include <cppad/local/declare_ad.hpp>

// ---------------------------------------------------------------------------
#include <cppad_cgoo/cg_cppadcg_assert.hpp>
#include <cppad_cgoo/cg_exception.hpp>
#include <cppad_cgoo/cg_declare_cg.hpp>

// ---------------------------------------------------------------------------
// CppAD
#include <cppad/cppad.hpp>

// ---------------------------------------------------------------------------
// system dependent files
#include <cppad_cgoo/dynamic_lib/cg_system.hpp>
#include <cppad_cgoo/dynamic_lib/linux/cg_linux_system.hpp>

// ---------------------------------------------------------------------------
// core files
#include <cppad_cgoo/cg_debug.hpp>
#include <cppad_cgoo/cg_operation.hpp>
#include <cppad_cgoo/cg_argument.hpp>
#include <cppad_cgoo/cg_source_code_fragment.hpp>
#include <cppad_cgoo/cg_cg.hpp>
#include <cppad_cgoo/cg_default.hpp>
#include <cppad_cgoo/cg_variable.hpp>
#include <cppad_cgoo/cg_identical.hpp>
#include <cppad_cgoo/cg_variable_name_generator.hpp>
#include <cppad_cgoo/cg_language.hpp>
#include <cppad_cgoo/cg_code_handler.hpp>
#include <cppad_cgoo/cg_arithmetic.hpp>
#include <cppad_cgoo/cg_arithmetic_assign.hpp>
#include <cppad_cgoo/cg_math.hpp>
#include <cppad_cgoo/cg_math_other.hpp>
#include <cppad_cgoo/cg_cond_exp_op.hpp>
#include <cppad_cgoo/cg_compare.hpp>
#include <cppad_cgoo/cg_ordered.hpp>
#include <cppad_cgoo/cg_unary.hpp>
#include <cppad_cgoo/cg_base_double.hpp>

// ---------------------------------------------------------------------------
// additional utilities
#include <cppad_cgoo/cg_util.hpp>
#include <cppad_cgoo/cg_evaluator.hpp>
#include <cppad_cgoo/cg_source_code_path.hpp>
#include <cppad_cgoo/cg_solver.hpp>
#include <cppad_cgoo/cg_graph_mod.hpp>

// ---------------------------------------------------------------------------
// C source code generation
#include <cppad_cgoo/c/cg_c_language.hpp>
#include <cppad_cgoo/c/cg_c_language_double.hpp>
#include <cppad_cgoo/c/cg_c_lang_default_var_name_gen.hpp>
#include <cppad_cgoo/c/cg_c_lang_default_hessian_var_name_gen.hpp>
#include <cppad_cgoo/c/cg_c_lang_custom_var_name_gen.hpp>

// automated dynamic library creation
#include <cppad_cgoo/dynamic_lib/cg_dynamiclib.hpp>
#include <cppad_cgoo/dynamic_lib/cg_dynamiclib_model.hpp>
#include <cppad_cgoo/dynamic_lib/cg_c_lang_compiler.hpp>
#include <cppad_cgoo/dynamic_lib/cg_c_lang_compile_model_helper.hpp>
#include <cppad_cgoo/dynamic_lib/cg_c_lang_compile_dynamic_helper.hpp>
#include <cppad_cgoo/dynamic_lib/cg_c_lang_compile_model_helper_impl.hpp>
#include <cppad_cgoo/dynamic_lib/cg_c_lang_compile_dynamic_helper_impl.hpp>
#include <cppad_cgoo/dynamic_lib/cg_gcc_compiler.hpp>

// ---------------------------------------------------------------------------
// automated dynamic library creation for Linux
#include <cppad_cgoo/dynamic_lib/linux/cg_linux_dynamiclib_model.hpp>
#include <cppad_cgoo/dynamic_lib/linux/cg_linux_dynamiclib.hpp>
#include <cppad_cgoo/dynamic_lib/linux/cg_linux_c_lang_compile_dynamic_helper.hpp>

// ---------------------------------------------------------------------------
// automated DAE differential index reduction
#define CPPAD_CG_DAE_VERBOSE

#endif

