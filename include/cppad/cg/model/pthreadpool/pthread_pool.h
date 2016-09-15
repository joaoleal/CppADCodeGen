#ifndef CPPADCG_PTHREAD_POOL_H
#define CPPADCG_PTHREAD_POOL_H
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

#ifdef __cplusplus
extern "C" {
#endif

void cppadcg_thpool_set_threads(int n);

int cppadcg_thpool_get_threads();

void cppadcg_thpool_prepare();

void cppadcg_thpool_add_job(void (*function_p)(void*), void* arg_p);

void cppadcg_thpool_wait();

void cppadcg_thpool_shutdown();

#ifdef __cplusplus
}
#endif

#endif
