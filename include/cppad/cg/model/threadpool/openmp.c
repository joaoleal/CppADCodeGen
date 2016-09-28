/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2016 Ciengis
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

static volatile int cppadcg_openmp_enabled = 1; // false

void cppadcg_openmp_set_disabled(int disabled) {
    cppadcg_openmp_enabled = !disabled;
}

int cppadcg_openmp_is_disabled() {
    return !cppadcg_openmp_enabled;
}