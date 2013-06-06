#ifndef CPPAD_CG_C_LANG_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_C_LANG_ATOMIC_FUN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

extern "C" {

    struct CLangAtomicFun {
        void* libModel;
        int (*forward)(void* libModel,
                int atomicIndex,
                int q,
                int p,
                const void* tx,
                unsigned long int txSize,
                void* ty,
                unsigned long int tySize);

        int (*reverse)(void* libModel,
                int atomicIndex,
                int p,
                const void* tx,
                const void* ty,
                void* px,
                const void* py,
                unsigned long int xSize,
                unsigned long int ySize);
    };

}

#endif