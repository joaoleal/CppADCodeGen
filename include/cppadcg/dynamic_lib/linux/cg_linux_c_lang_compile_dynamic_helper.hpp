#ifndef CPPAD_CG_LINUX_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
#define CPPAD_CG_LINUX_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */
#ifdef __linux__

namespace CppAD {

    template<class Base>
    DynamicLib<Base>* CLangCompileDynamicHelper<Base>::loadDynamicLibrary() {
        return new LinuxDynamicLib<Base> (_libraryName + system::SystemInfo<>::DYNAMIC_LIB_EXTENSION);
    }

}
#endif

#endif
