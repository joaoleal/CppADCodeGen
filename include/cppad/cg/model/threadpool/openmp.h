#ifndef CPPADCG_OPENMP_H
#define CPPADCG_OPENMP_H
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

enum ScheduleStrategy {SCHED_SINGLE_JOB, // omp_sched_dynamic with chunk size 1
                       SCHED_MULTI_JOB, // omp_sched_guided
                       SCHED_STATIC}; // omp_sched_static

void cppadcg_openmp_set_disabled(int disabled);

int cppadcg_openmp_is_disabled();


void cppadcg_openmp_set_threads(unsigned int n);

unsigned int cppadcg_openmp_get_threads();


#ifdef __cplusplus
}
#endif

#endif