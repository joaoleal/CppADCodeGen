#ifndef CPPAD_CG_C_LANG_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_C_LANG_ATOMIC_FUN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

extern "C" {

    /**
     * Holds function pointers that the compiled code uses to call atomic functions.
     */
    struct CLangAtomicFun {
        /**
         * A pointer to the compiled model object (e.g. LinuxDynamicLibModel)
         */
        void* libModel;
        
        int (*forward)(void* libModel,
                int atomicIndex,
                int q,
                int p,
                const void* tx,
                unsigned long txSize,
                void* ty,
                unsigned long tySize);

        int (*reverse)(void* libModel,
                int atomicIndex,
                int p,
                const void* tx,
                const void* ty,
                void* px,
                const void* py,
                unsigned long xSize,
                unsigned long ySize);
    };

}

#endif