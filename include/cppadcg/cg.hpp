#ifndef CPPAD_CG_INCLUDED
#define CPPAD_CG_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Eclipse Public License Version 1.0 (EPL1), and
 *   - GNU General Public License Version 3 (GPL3).
 *
 *  EPL1 terms and conditions can be found in the file "epl-v10.txt", while
 *  terms and conditions for the GPL3 can be found in the file "gpl3.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
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
#include <list>
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
#include <cppadcg/cg_cppadcg_assert.hpp>
#include <cppadcg/cg_exception.hpp>
#include <cppadcg/cg_declare_cg.hpp>

// ---------------------------------------------------------------------------
// CppAD
#include <cppad/cppad.hpp>
// addons
#include <cppad/extra/extra.hpp>

// ---------------------------------------------------------------------------
// system dependent files
#include <cppadcg/dynamic_lib/cg_system.hpp>
#include <cppadcg/dynamic_lib/linux/cg_linux_system.hpp>

// ---------------------------------------------------------------------------
// core files
#include <cppadcg/cg_debug.hpp>
#include <cppadcg/cg_operation.hpp>
#include <cppadcg/cg_argument.hpp>
#include <cppadcg/cg_operation_node.hpp>
#include <cppadcg/nodes/cg_index_operation_node.hpp>
#include <cppadcg/nodes/cg_index_dclr_operation_node.hpp>
#include <cppadcg/nodes/cg_index_assign_operation_node.hpp>
#include <cppadcg/nodes/cg_loop_start_operation_node.hpp>
#include <cppadcg/nodes/cg_loop_end_operation_node.hpp>
#include <cppadcg/cg_cg.hpp>
#include <cppadcg/cg_default.hpp>
#include <cppadcg/cg_variable.hpp>
#include <cppadcg/cg_identical.hpp>
#include <cppadcg/cg_variable_name_generator.hpp>
#include <cppadcg/cg_job_time.hpp>
#include <cppadcg/cg_language.hpp>
#include <cppadcg/cg_scope_path_element.hpp>
#include <cppadcg/cg_code_handler.hpp>
#include <cppadcg/patterns/cg_loop_position.hpp>
#include <cppadcg/cg_code_handler_loops.hpp>
#include <cppadcg/cg_arithmetic.hpp>
#include <cppadcg/cg_arithmetic_assign.hpp>
#include <cppadcg/cg_math.hpp>
#include <cppadcg/cg_math_other.hpp>
#include <cppadcg/cg_nan.hpp>
#include <cppadcg/cg_cond_exp_op.hpp>
#include <cppadcg/cg_compare.hpp>
#include <cppadcg/cg_ordered.hpp>
#include <cppadcg/cg_unary.hpp>
#include <cppadcg/cg_base_double.hpp>

// ---------------------------------------------------------------------------
// additional utilities
#include <cppadcg/cg_util.hpp>
#include <cppadcg/cg_evaluator.hpp>
#include <cppadcg/cg_operation_path.hpp>
#include <cppadcg/cg_solver.hpp>
#include <cppadcg/cg_graph_mod.hpp>

// ---------------------------------------------------------------------------
// atomic function utilities
#include <cppadcg/cg_custom_position.hpp>
#include <cppadcg/cg_base_abstract_atomic_fun.hpp>
#include <cppadcg/cg_abstract_atomic_fun.hpp>
#include <cppadcg/cg_atomic_fun.hpp>
#include <cppadcg/cg_atomic_fun_bridge.hpp>
#include <cppadcg/dynamic_lib/cg_atomic_lib_model.hpp>

// ---------------------------------------------------------------------------
// loop/pattern detection
#include <cppadcg/patterns/index/cg_index_pattern.hpp>
#include <cppadcg/patterns/index/cg_linear_index_pattern.hpp>
#include <cppadcg/patterns/index/cg_sectioned_index_pattern.hpp>
#include <cppadcg/patterns/index/cg_plane_2d_index_pattern.hpp>
#include <cppadcg/patterns/index/cg_random_index_pattern.hpp>
#include <cppadcg/patterns/index/cg_index_pattern_impl.hpp>
#include <cppadcg/patterns/cg_loop_model.hpp>
#include <cppadcg/patterns/cg_loop_free_model.hpp>
#include <cppadcg/patterns/cg_equation_pattern.hpp>
#include <cppadcg/patterns/cg_loop.hpp>
#include <cppadcg/patterns/cg_dependent_pattern_matcher.hpp>

// ---------------------------------------------------------------------------
// C source code generation
#include <cppadcg/c/cg_c_lang_atomic_fun.hpp>
#include <cppadcg/c/cg_c_language.hpp>
#include <cppadcg/c/cg_c_language_index_patterns.hpp>
#include <cppadcg/c/cg_c_language_double.hpp>
#include <cppadcg/c/cg_c_lang_default_var_name_gen.hpp>
#include <cppadcg/c/cg_c_lang_default_hessian_var_name_gen.hpp>
#include <cppadcg/c/cg_c_lang_default_reverse2_var_name_gen.hpp>
#include <cppadcg/c/cg_c_lang_custom_var_name_gen.hpp>

// automated static library creation
#include <cppadcg/dynamic_lib/cg_archiver.hpp>
#include <cppadcg/dynamic_lib/cg_ar_archiver.hpp>

// automated dynamic library creation
#include <cppadcg/dynamic_lib/cg_dynamiclib.hpp>
#include <cppadcg/dynamic_lib/cg_dynamiclib_model.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compiler.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_dynamic_helper.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_impl.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_for0.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_for1.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_rev1.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_rev2.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_jac.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_model_helper_hes.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_for0.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_for1.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_jac.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_jac_fr1.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_hess.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_rev1.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_rev2.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_c_lang_compile_model_helper_loops_hess_r2.hpp>
#include <cppadcg/dynamic_lib/patterns/cg_hessian_with_loops_info.hpp>
#include <cppadcg/dynamic_lib/cg_c_lang_compile_dynamic_helper_impl.hpp>
#include <cppadcg/dynamic_lib/cg_gcc_compiler.hpp>

// ---------------------------------------------------------------------------
// automated dynamic library creation for Linux
#include <cppadcg/dynamic_lib/linux/cg_linux_dynamiclib_model.hpp>
#include <cppadcg/dynamic_lib/linux/cg_linux_dynamiclib.hpp>
#include <cppadcg/dynamic_lib/linux/cg_linux_c_lang_compile_dynamic_helper.hpp>

// ---------------------------------------------------------------------------
// automated DAE differential index reduction
//#define CPPAD_CG_DAE_VERBOSE

#endif

