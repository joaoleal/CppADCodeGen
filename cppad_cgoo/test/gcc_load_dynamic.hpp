#ifndef GCC_LOAD_DYNAMIC_HPP
#define	GCC_LOAD_DYNAMIC_HPP

#include <string>
#include <dlfcn.h>
#include <cppad_cgoo/cg.hpp>

#include "TestException.hpp"


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

std::vector<std::vector<double> > runSparseJacDefault(CppAD::ADFun<double>& f,
                                                      const std::vector<std::vector<double> >& ind);

std::vector<std::vector<double> > runSparseJac(CppAD::ADFun<CppAD::CG<double> >& f,
                                               const std::string& library,
                                               const std::string& functionJac,
                                               const std::vector<std::vector<double> >& indV);

bool test0(const std::string& test,
           CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
           CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
           const std::vector<std::vector<double> >& indV,
           int& comparisons,
           double epsilonR = 1e-14, double epsilonA = 1e-14);

bool test0nJac(const std::string& test,
               CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
               CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
               const std::vector<double>& ind,
               double epsilonR = 1e-14, double epsilonA = 1e-14);

bool test0nJac(const std::string& test,
               CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
               CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
               const std::vector<std::vector<double> >& indV,
               double epsilonR = 1e-14, double epsilonA = 1e-14);

void prepareADFun(const std::vector<double>& indep,
                  CppAD::ADFun<double>* (*func1)(const std::vector<CppAD::AD<double> >&),
                  CppAD::ADFun<CppAD::CG<double> >* (*func2)(const std::vector<CppAD::AD<CppAD::CG<double> > >&),
                  CppAD::ADFun<double>*& f1,
                  CppAD::ADFun<CppAD::CG<double> >*& f2);

bool compareValues(const std::string& testType,
                   const std::vector<std::vector<double> >& depCGen,
                   const std::vector<std::vector<double> >& dep,
                   double epsilonR = 1e-14, double epsilonA = 1e-14);

template<class T>
bool compareValues(std::vector<double> cgen, T orig,
                   const std::string& nameCgen, const std::string& nameOrig,
                   double epsilonR = 1e-14, double epsilonA = 1e-14) {
    for (size_t i = 0; i < cgen.size(); i++) {
        if (!NearEqual(cgen[i], orig[i], epsilonR, epsilonA)) {
            std::cerr << nameCgen << "[" << i << "] = " << std::setprecision(8) << cgen[i] <<
                    " != "
                    << nameOrig << "[" << i << "] = " << std::setprecision(8) << orig[i] << std::endl;
            return false;
        }
    }

    return true;
}

#endif	/* GCC_LOAD_DYNAMIC_HPP */

