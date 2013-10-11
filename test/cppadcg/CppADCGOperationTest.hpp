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

#include <dlfcn.h>
#include <cppadcg/cg.hpp>

#include "CppADCGTest.hpp"
#include "TestException.hpp"

#ifndef CPPAD_CG_CPPADCGOPERATIONTEST_HPP
#define	CPPAD_CG_CPPADCGOPERATIONTEST_HPP

namespace CppAD {

    class CppADCGOperationTest : public CppADCGTest {
    public:

        inline CppADCGOperationTest(bool verbose = false, bool printValues = false) :
            CppADCGTest(verbose, printValues) {
        }

        virtual void TearDown() {
            ASSERT_FALSE(CppAD::memory_leak());
        }

    protected:
        void* loadLibrary(const std::string& library) throw (CppAD::TestException);

        void* getFunction(void * libHandle, const std::string& functionName) throw (CppAD::TestException);

        void closeLibrary(void* libHandle);

        void compile(const std::string& source, const std::string& library) throw (CppAD::TestException);

        std::vector<std::vector<double> > runDefault0(CppAD::ADFun<double>& f,
                                                      const std::vector<std::vector<double> >& ind);

        std::vector<std::vector<double> > runDefault0(CppAD::ADFun<double>& f,
                                                      const std::vector<std::vector<double> >& indV,
                                                      int& comparisons);

        std::vector<std::vector<double> > run0(CppAD::ADFun<CppAD::CG<double> >& f,
                                               const std::string& library, const std::string& function,
                                               const std::vector<std::vector<double> >& ind);

        std::vector<std::vector<double> > run0(CppAD::ADFun<CppAD::CG<double> >& f,
                                               const std::string& library, const std::string& function,
                                               const std::vector<std::vector<double> >& indV,
                                               int& comparisons);

        std::vector<std::vector<double> > run0TapeWithValues(CppAD::ADFun<CppAD::CG<double> >& f,
                                                             const std::vector<std::vector<double> >& indV);

        std::vector<std::vector<double> > runSparseJacDefault(CppAD::ADFun<double>& f,
                                                              const std::vector<std::vector<double> >& ind);

        std::vector<std::vector<double> > runSparseJac(CppAD::ADFun<CppAD::CG<double> >& f,
                                                       const std::string& library,
                                                       const std::string& functionJac,
                                                       const std::vector<std::vector<double> >& indV);

        void test0(const std::string& test,
                   CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
                   CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
                   const std::vector<std::vector<double> >& indV,
                   int& comparisons,
                   double epsilonR = 1e-14, double epsilonA = 1e-14);

        void test0nJac(const std::string& test,
                       CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
                       CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
                       const std::vector<double>& ind,
                       double epsilonR = 1e-14, double epsilonA = 1e-14);

        void test0nJac(const std::string& test,
                       CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
                       CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
                       const std::vector<std::vector<double> >& indV,
                       double epsilonR = 1e-14, double epsilonA = 1e-14);

        void prepareADFun(const std::vector<double>& indep,
                          CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
                          CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
                          CppAD::ADFun<double>*& f1,
                          CppAD::ADFun<CppAD::CG<double> >*& f2);
    };

}

#include "CppADCGOperationTestImpl.hpp"

#endif
