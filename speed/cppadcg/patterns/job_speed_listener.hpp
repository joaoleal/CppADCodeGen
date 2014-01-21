#ifndef CPPAD_CG_JOBSPEEDLISTENER_INCLUDED
#define	CPPAD_CG_JOBSPEEDLISTENER_INCLUDED
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

#include <cppadcg/cg.hpp>

namespace CppAD {

    class JobSpeedListener : public JobListener {
    public:
        // pattern detection
        double patternDection;
        // pattern detection
        double graphGen;
        // total time used for source code generation
        double srcCodeGen;
        // source code compilation
        double srcCodeComp;
        // compilation of the dynamic library
        double dynLibComp;
        // JIT preparation time
        double jit;
        // total time used to compile the sources and generate the library
        double totalLibrary;
    public:
        JobSpeedListener();

        inline void reset() {
            patternDection = 0;
            graphGen = 0;
            srcCodeGen = 0;
            srcCodeComp = 0;
            dynLibComp = 0;
            jit = 0;
            totalLibrary = 0;
        }

        virtual void jobStarted(const std::vector<Job>& job);

        virtual void jobEndended(const std::vector<Job>& job,
                                 double elapsed);
    };

}

#endif