#ifndef CPPAD_CG_LINUX_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
#define CPPAD_CG_LINUX_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
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
#ifdef __linux__

namespace CppAD {
namespace cg {

template<class Base>
DynamicLib<Base>* DynamicModelLibraryProcessor<Base>::loadDynamicLibrary() {
    return new LinuxDynamicLib<Base> (_libraryName + system::SystemInfo<>::DYNAMIC_LIB_EXTENSION);
}

} // END cg namespace
} // END CppAD namespace

#endif
#endif