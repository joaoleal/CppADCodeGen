/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-11 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

// system include files used for I/O
#include <iostream>
#include <vector>

// memory leak checker
#include <cppad/memory_leak.hpp>
#include <cppad_cgoo/cg.hpp>

using namespace CppAD;

void prepareADFun(const std::vector<double>& indep,
                  ADFun<double> (*func1)(const std::vector<AD<double> >&),
                  ADFun<CG<double> > (*func2)(const std::vector<AD<CG<double> > >&),
                  ADFun<double>*& f1,
                  ADFun<CG<double> >*& f2);

// various examples / tests
//extern bool Abs();
//extern bool Acos();
extern bool Add();
//extern bool Asin();
//extern bool Atan();
//extern bool Atan2();
//extern bool CompareChange();
//extern bool CondExp();
//extern bool Cos();
//extern bool Cosh();
//extern bool Div();
//extern bool Exp();
//extern bool Log();
//extern bool Log10();
//extern bool Mul();
//extern bool Parameter();
//extern bool Pow();
//extern bool Sin();
//extern bool Sub();
//extern bool Tan();

namespace {
    // function that runs one test
    static size_t Run_ok_count = 0;
    static size_t Run_error_count = 0;
    static size_t Run_warn_count = 0;

    bool Run(bool TestOk(void), std::string name) {
        bool ok = true;
        std::streamsize width = 20;
        //
        ok &= name.size() < size_t(width);
        ok &= TestOk();

        std::cout.width(width);
        std::cout.setf(std::ios_base::left);
        std::cout << name;
        if (ok) {
            std::cout << "OK" << std::endl;
            Run_ok_count++;
        } else {
            std::cout << "Error" << std::endl;
            Run_error_count++;
        }
        return ok;
    }
}

// main program that runs all the tests

int main(void) {
    bool ok = true;
    using namespace std;

    // This line is used by test_one.sh 
//    ok &= Run(Abs, "Abs");
//    ok &= Run(Acos, "Acos");
    ok &= Run(Add, "Add");
//    ok &= Run(Asin, "Asin");
//    ok &= Run(Atan, "Atan");
//    ok &= Run(Atan2, "Atan2");
//    ok &= Run(CompareChange, "CompareChange");
//    ok &= Run(CondExp, "CondExp");
//    ok &= Run(Cos, "Cos");
//    ok &= Run(Cosh, "Cosh");
//    ok &= Run(Div, "Div");
//    ok &= Run(Exp, "Exp");
//    ok &= Run(Log, "Log");
//    ok &= Run(Log10, "Log10");
//    ok &= Run(Mul, "Mul");
//    ok &= Run(Parameter, "Parameter");
//    ok &= Run(Pow, "Pow");
//    ok &= Run(Sin, "Sin");
//    ok &= Run(Sub, "Sub");
//    ok &= Run(Tan, "Tan");


    // check for errors
    using std::cout;
    using std::endl;
    assert(ok || (Run_error_count > 0));
    if (CppAD::memory_leak()) {
        ok = false;
        Run_error_count++;
        cout << "Error: " << "memory leak detected" << endl;
    } else {
        Run_ok_count++;
        cout << "OK:    " << "No memory leak detected" << endl;
    }
    // convert int(size_t) to avoid warning on _MSC_VER systems
    ok &= Run_warn_count == 0;
    if (ok) {
        cout << int(Run_ok_count) << " tests passed";
    } else {
        cout << int(Run_error_count) << " tests failed." << endl;
    }

    return static_cast<int> (!ok);
}

// END PROGRAM
