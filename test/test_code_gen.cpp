/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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

// system include files used for I/O
#include <iostream>
#include <vector>

// memory leak checker
#include <cppad/memory_leak.hpp>
#include <cppadcg/cg.hpp>

using namespace CppAD;
using namespace std;

// various examples / tests
extern bool Abs();
extern bool Acos();
extern bool Add();
extern bool Asin();
extern bool Assign();
extern bool Atan();
extern bool Atan2();
extern bool CompareChange();
extern bool CondExp();
extern bool Cos();
extern bool Cosh();
extern bool Div();
extern bool Exp();
extern bool Log();
extern bool Log10();
extern bool Mul();
extern bool Parameter();
extern bool Pow();
extern bool Sin();
extern bool Sub();
extern bool Tan();
extern bool Unary();
extern bool Dynamic();
//extern bool HandlerReset();

// solve
extern bool SolveAdd();
extern bool SolveCosh();
extern bool SolveDiv();
extern bool SolveExp();
extern bool SolveLog();
extern bool SolveLog10();
extern bool SolveMul();
extern bool SolvePow();
extern bool SolveSinh();
extern bool SolveSqrt();
extern bool SolveSub();
extern bool SolveTanh();
extern bool SolveUnaryMinus();

// DAE index reduction
extern bool Pantelides();
extern bool DummyDeriv();

bool test_verbose = false;
bool test_printvalues = false;

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
        try {
            ok &= TestOk();
        } catch (const std::exception& ex) {
            std::cerr << "   Error: " << ex.what() << std::endl;
            ok = false;
        } catch (...) {
            ok = false;
        }

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
    using std::cout;
    using std::endl;

    bool ok = true;
    ok &= Run(Assign, "Assign");
    ok &= Run(Abs, "Abs");
    ok &= Run(Acos, "Acos");
    ok &= Run(Add, "Add");
    ok &= Run(Asin, "Asin");
    ok &= Run(Atan, "Atan");
    ok &= Run(Atan2, "Atan2");
    ok &= Run(CompareChange, "CompareChange");
    ok &= Run(CondExp, "CondExp");
    ok &= Run(Cos, "Cos");
    ok &= Run(Cosh, "Cosh");
    ok &= Run(Div, "Div");
    ok &= Run(Exp, "Exp");
    ok &= Run(Log, "Log");
    ok &= Run(Log10, "Log10");
    ok &= Run(Mul, "Mul");
    ok &= Run(Parameter, "Parameter");
    ok &= Run(Pow, "Pow");
    ok &= Run(Sin, "Sin");
    ok &= Run(Sub, "Sub");
    ok &= Run(Tan, "Tan");
    ok &= Run(Unary, "Unary");
    //ok &= Run(HandlerReset, "HandlerReset");

    ok &= Run(Dynamic, "Dynamic");

    // Solve
    ok &= Run(SolveAdd, "SolveAdd");
    ok &= Run(SolveTanh, "SolveCosh");
    ok &= Run(SolveDiv, "SolveDiv");
    ok &= Run(SolveExp, "SolveExp");
    ok &= Run(SolveLog, "SolveLog");
    ok &= Run(SolveLog10, "SolveLog10");
    ok &= Run(SolveMul, "SolveMul");
    ok &= Run(SolvePow, "SolvePow");
    ok &= Run(SolveTanh, "SolveSinh");
    ok &= Run(SolveSqrt, "SolveSqrt");
    ok &= Run(SolveSub, "SolveSub");
    ok &= Run(SolveTanh, "SolveTanh");
    ok &= Run(SolveUnaryMinus, "SolveUnaryMinus");

    // check for errors
    assert(ok || (Run_error_count > 0));
    if (CppAD::memory_leak()) {
        ok = false;
        Run_error_count++;
        cout << "Error: " << "memory leak detected" << endl;
    } else {
        Run_ok_count++;
        cout << "OK:    " << "No memory leak detected" << endl;
    }

    // DAE index reduction (must be called after memory check because the atomic functions use static data)
    ok &= Run(Pantelides, "Pantelides");
    ok &= Run(DummyDeriv, "DummyDerivatives");

    // convert int(size_t) to avoid warning on _MSC_VER systems
    ok &= Run_warn_count == 0;
    if (ok) {
        cout << int(Run_ok_count) << " tests passed";
    } else {
        cout << int(Run_error_count) << " tests failed.";
    }

    cout << endl;

    return static_cast<int> (!ok);
}

// END PROGRAM
