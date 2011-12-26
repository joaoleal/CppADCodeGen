#ifndef GCC_LOAD_DYNAMIC_HPP
#define	GCC_LOAD_DYNAMIC_HPP

#include <string>
#include <dlfcn.h>

#include "TestException.hpp"


void* loadLibrary(const std::string& library) throw (CppAD::TestException);

void* getFunction(void * libHandle, const std::string& functionName) throw (CppAD::TestException);

void closeLibrary(void* libHandle);

void compile(const std::string& source, const std::string& library) throw (CppAD::TestException);

bool run0(CppAD::ADFunCodeGen<double>& f, const std::string& library, const std::string& function,
        const std::vector<double>& ind,
        int& comparisons, std::vector<double>& depCGen);

bool runTest0(CppAD::ADFunCodeGen<double>& f, const std::string& library, const std::string& function, const std::vector<double>& ind, const CPPAD_TEST_VECTOR< CppAD::AD<double> >& U, double epsilonR = 1e-14, double epsilonA = 1e-14);

bool runTestSparseJac(CppAD::ADFunCodeGen<double>& f, const std::string& library, const std::string& function, const std::vector<double>& ind, const std::vector<double>& jac, double epsilonR = 1e-14, double epsilonA = 1e-14);

bool test0nJac(const std::string& test, CppAD::ADFunCodeGen<double>& f, const std::vector<double>& ind, const CPPAD_TEST_VECTOR< CppAD::AD<double> >& dep, double epsilonR = 1e-14, double epsilonA = 1e-14);

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

