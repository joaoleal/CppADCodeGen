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

#include <omp.h>
#include <stdio.h>

static volatile int cppadcg_openmp_enabled = 1; // false
static volatile int cppadcg_openmp_verbose = 1; // false
static volatile unsigned int cppadcg_openmp_n_threads = 2;

void cppadcg_openmp_set_disabled(int disabled) {
    cppadcg_openmp_enabled = !disabled;
}

int cppadcg_openmp_is_disabled() {
    return !cppadcg_openmp_enabled;
}

void cppadcg_openmp_set_verbose(int v) {
    cppadcg_openmp_verbose = v;
}

int cppadcg_openmp_is_verbose() {
    return cppadcg_openmp_verbose;
}

void cppadcg_openmp_set_threads(unsigned int n) {
    if(n <= 0)
        n = 1;
    cppadcg_openmp_n_threads = n;
}

unsigned int cppadcg_openmp_get_threads() {
    return cppadcg_openmp_n_threads;
}