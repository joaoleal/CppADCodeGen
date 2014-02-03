#ifndef CPPAD_CG_CPPADCGOPERATIONTESTIMPL_HPP
#define	CPPAD_CG_CPPADCGOPERATIONTESTIMPL_HPP
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
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
#include <sys/types.h>
#include <sys/wait.h>
#include <memory>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>

namespace CppAD {
namespace cg {

using namespace CppAD;

void* CppADCGOperationTest::loadLibrary(const std::string& library) throw (TestException) {
    void * libHandle = dlopen(library.c_str(), RTLD_NOW);
    if (!libHandle) {
        throw TestException("Failed to dynamically load library");
    }
    return libHandle;
}

void CppADCGOperationTest::closeLibrary(void* libHandle) {
    dlclose(libHandle);
}

void* CppADCGOperationTest::getFunction(void * libHandle, const std::string& functionName) throw (TestException) {
    void* functor = dlsym(libHandle, functionName.c_str());
    char *error;
    if ((error = dlerror()) != nullptr) {
        throw TestException(error);
    }
    return functor;
}

void CppADCGOperationTest::compile(const std::string& source, const std::string& library) throw (TestException) {
    int fd[2];
    //Create pipe for piping source to gcc
    if (pipe(fd) < 0) {
        throw TestException("Failed to create pipe");
    }

    //Fork a gcc, pipe source to it, wait for gcc to exit
    pid_t pid = fork();
    if (pid < 0) {
        throw TestException("Failed to fork program");
    }

    if (pid == 0) {
        //  Child process
        // close write end of pipe
        close(fd[1]);
        // Send pipe input to stdin
        close(STDIN_FILENO);
        dup2(fd[0], STDIN_FILENO);
        /**
         * Call gcc
         * 
         * Arguments:
         *   -O0                   Optimization level
         *   -x c                  C source
         *   -pipe                 Use pipes between gcc stages
         *   -fPIC -shared         Make shared object
         *   -Wl,-soname, library  Pass suitable options to linker
         * 
         */
        std::string linker = "-Wl,-soname," + library;
        execl("/usr/bin/gcc", "gcc", "-x", "c", "-O0", "-pipe", "-", "-fPIC", "-shared",
              linker.c_str(), "-o", library.c_str(), (char *) nullptr);

        exit(0);
    }

    // Parent process
    // close read end of pipe
    close(fd[0]);
    //Pipe source to gcc
    write(fd[1], source.c_str(), source.size());
    close(fd[1]);

    //Wait for gcc to exit
    int status;
    if (wait(&status) < 0) {
        throw TestException("Failed while waiting for gcc");
    }
}

std::vector<std::vector<double> > CppADCGOperationTest::runDefault0(ADFun<double>& f,
                                                                    const std::vector<std::vector<double> >& ind) {
    int comparisons;

    return runDefault0(f, ind, comparisons);
}

std::vector<std::vector<double> > CppADCGOperationTest::runDefault0(ADFun<double>& f,
                                                                    const std::vector<std::vector<double> >& indV,
                                                                    int& comparisons) {

    using std::vector;

    vector<vector<double> >dep(indV.size());
    for (size_t i = 0; i < indV.size(); i++) {

        // the regular CppAD
        dep[i] = f.Forward(0, indV[i]);

        comparisons = f.CompareChange();
    }

    return dep;
}

std::vector<std::vector<double> > CppADCGOperationTest::run0(ADFun<CG<double> >& f,
                                                             const std::string& library, const std::string& function,
                                                             const std::vector<std::vector<double> >& ind) {

    int comparisons;

    return run0(f, library, function, ind, comparisons);
}

std::vector<std::vector<double> > CppADCGOperationTest::run0(ADFun<CG<double> >& f,
                                                             const std::string& library, const std::string& function,
                                                             const std::vector<std::vector<double> >& indV,
                                                             int& comparisons) {
    using CppAD::vector;

    size_t n = indV.begin()->size();

    CodeHandler<double> handler(10 + n * n);

    vector<CG<double> > indVars(n);
    handler.makeVariables(indVars);

    vector<CG<double> > dep = f.Forward(0, indVars);

    CLanguage<double> langC("double");
    CLangDefaultVariableNameGenerator<double> nameGen;

    std::ostringstream code;
    handler.generateCode(code, langC, dep, nameGen);

    std::string spaces = "   ";
    std::string source = "#include <math.h>\n\n"
            "int " + function + "(const double* x, double* y) {\n";

    // declare variables
    source += langC.generateTemporaryVariableDeclaration();

    source += code.str();

    //source += "return " + n->compareChangeCounter() + n->endl();
    source += spaces + "return 0;\n";
    source += "}";

    if (verbose_) {
        std::cout << std::endl << source << std::endl;
    }

    compile(source, library);

    void* libHandle = loadLibrary(library);

    int (*fn)(const double*, double*) = nullptr;

    try {
        *(void **) (&fn) = getFunction(libHandle, function);
    } catch (const std::exception& ex) {
        closeLibrary(libHandle);
        throw;
    }

    std::vector<std::vector<double> >depCGen(indV.size());
    for (size_t i = 0; i < indV.size(); i++) {

        std::vector<double>& depi = depCGen[i];
        depi.resize(dep.size());

        const std::vector<double>& ind = indV[i];

        // the compiled version
        comparisons = (*fn)(&ind[0], &depi[0]);

        if (verbose_ && printValues_) {
            for (size_t j = 0; j < ind.size(); j++) {
                std::cout << " x[" << j << "] = " << ind[j] << "\n";
            }

            for (size_t j = 0; j < depi.size(); j++) {
                std::cout << " y[" << j << "] = " << depi[j] << "\n";
            }
        }
    }

    closeLibrary(libHandle);

    return depCGen;
}

std::vector<std::vector<double> > CppADCGOperationTest::run0TapeWithValues(ADFun<CG<double> >& f,
                                                                           const std::vector<std::vector<double> >& indV) {
    using std::vector;

    vector<vector<double> > results(indV.size());

    size_t n = indV.begin()->size();

    for (size_t i = 0; i < indV.size(); i++) {
        const vector<double>& ind = indV[i];
        CodeHandler<double> handler(10 + n * n);

        vector<CG<double> > indVars(n);
        handler.makeVariables(indVars);
        for (size_t j = 0; j < n; j++) {
            indVars[j].setValue(ind[j]);
        }

        vector<CG<double> > dep = f.Forward(0, indVars);
        results[i].resize(dep.size());
        for (size_t j = 0; j < dep.size(); j++) {
            results[i][j] = dep[j].getValue();
        }
    }

    return results;
}

std::vector<std::vector<double> > CppADCGOperationTest::runSparseJacDefault(CppAD::ADFun<double>& f,
                                                                            const std::vector<std::vector<double> >& ind) {
    using std::vector;

    vector<vector<double> >jac(ind.size());

    for (size_t i = 0; i < ind.size(); i++) {
        /**
         * Jacobian (regular CppAD)
         */
        jac[i] = f.SparseJacobian(ind[i]);
    }

    return jac;
}

std::vector<std::vector<double> > CppADCGOperationTest::runSparseJac(ADFun<CG<double> >& f,
                                                                     const std::string& library,
                                                                     const std::string& functionJac,
                                                                     const std::vector<std::vector<double> >& indV) {
    using CppAD::vector;

    assert(!indV.empty());

    CodeHandler<double> handler(50 + indV.size() * indV.size());

    vector<CG<double> > indVars(indV[0].size());
    handler.makeVariables(indVars);

    vector<CG<double> > jacCG = f.SparseJacobian(indVars);

    CLanguage<double> langC("double");
    CLangDefaultVariableNameGenerator<double> nameGen;

    std::ostringstream code;
    handler.generateCode(code, langC, jacCG, nameGen);

    std::string spaces = "   ";
    std::string source = "#include <math.h>\n\n"
            "int " + functionJac + "(const double* x, double* y) {\n";

    // declare variables
    source += langC.generateTemporaryVariableDeclaration();

    source += code.str();

    //source += "return " + n->compareChangeCounter() + n->endl();
    source += spaces + "return 0;\n";
    source += "}\n\n";

    if (verbose_) {
        std::cout << std::endl << source << std::endl;
    }

    /**
     * Compile
     */
    compile(source, library);

    void* libHandle = loadLibrary(library);

    int (*fn)(const double*, double*) = nullptr;

    try {
        *(void **) (&fn) = getFunction(libHandle, functionJac);
    } catch (const std::exception& ex) {
        closeLibrary(libHandle);
        throw;
    }

    std::vector<std::vector<double> > jac(indV.size());

    for (size_t i = 0; i < indV.size(); i++) {
        jac[i].resize(f.Range() * f.Domain());

        /**
         * test jacobian (compiled)
         */
        (*fn)(&indV[i][0], &jac[i][0]);
    }

    closeLibrary(libHandle);

    return jac;
}

void CppADCGOperationTest::test0nJac(const std::string& test,
                                     ADFun<double>* (*func1)(const std::vector<AD<double> >&),
                                     ADFun<CG<double> >* (*func2)(const std::vector<AD<CG<double> > >&),
                                     const std::vector<double>& ind,
                                     double epsilonR, double epsilonA) {

    std::vector<std::vector<double> > indV;
    indV.push_back(ind);

    test0nJac(test, func1, func2, indV, epsilonR, epsilonA);
}

void CppADCGOperationTest::test0nJac(const std::string& test,
                                     ADFun<double>* (*func1)(const std::vector<AD<double> >&),
                                     ADFun<CG<double> >* (*func2)(const std::vector<AD<CG<double> > >&),
                                     const std::vector<std::vector<double> >& indV,
                                     double epsilonR, double epsilonA) {

    using std::vector;

    assert(!indV.empty());

    /**
     * Determine the values using the default CppAD
     */
    // independent variable vector, indices, values, and declaration
    std::vector< AD<double> > U1(indV[0].size());
    for (size_t i = 0; i < U1.size(); i++) {
        U1[i] = indV[0][i];
    }
    Independent(U1);

    // create f: U -> Z and vectors used for derivative calculations
    CppAD::ADFun<double>* f1 = (*func1)(U1);

    vector<vector<double> > depsDef = runDefault0(*f1, indV);
    vector<vector<double> > jacDef = runSparseJacDefault(*f1, indV);

    delete f1;

    /**
     * Determine the values using the compiled version
     */
    vector<AD<CG<double> > > u2(indV[0].size());
    // values must be given during tapping in order to avoid NaN
    for (size_t i = 0; i < u2.size(); i++) {
        u2[i] = indV[0][i];
    }
    Independent(u2);

    std::unique_ptr<CppAD::ADFun<CG<double> > > f2((*func2)(u2));

    std::string library = "./tmp/test_" + test + ".so";
    std::string function = "test_" + test;
    vector<vector<double> > depsCG = run0(*f2.get(), library, function, indV);

    vector<vector<double> > depCGTape = run0TapeWithValues(*f2.get(), indV);

    library = "./tmp/test_" + test + "_jac.so";
    std::string functionJac = "test_" + test + "_jac_for";
    vector<vector<double> > jacCG = runSparseJac(*f2.get(), library, functionJac, indV);

    /**
     * compare results
     */
    // Forward 0
    ASSERT_TRUE(compareValues(depsCG, depsDef, epsilonR, epsilonA));
    // Forward 0 (from taped values)
    ASSERT_TRUE(compareValues(depCGTape, depsDef, epsilonR, epsilonA));
    // Jacobian
    ASSERT_TRUE(compareValues(jacCG, jacDef, epsilonR, epsilonA));
}

void CppADCGOperationTest::test0(const std::string& test,
                                 ADFun<double>* (*func1)(const std::vector<AD<double> >&),
                                 ADFun<CG<double> >* (*func2)(const std::vector<AD<CG<double> > >&),
                                 const std::vector<std::vector<double> >& indV,
                                 int& comparisons,
                                 double epsilonR, double epsilonA) {

    using std::vector;

    assert(!indV.empty());

    /**
     * Determine the values using the default CppAD
     */
    // independent variable vector, indices, values, and declaration
    std::vector< AD<double> > U1(indV[0].size());
    for (size_t i = 0; i < U1.size(); i++) {
        U1[i] = indV[0][i];
    }
    Independent(U1);

    // create f: U -> Z and vectors used for derivative calculations
    CppAD::ADFun<double>* f1 = (*func1)(U1);

    vector<vector<double> > depsDef = runDefault0(*f1, indV);

    delete f1;

    /**
     * Determine the values using the compiled version
     */
    vector<AD<CG<double> > > u2(indV[0].size());
    // values must be given during tapping in order to avoid NaN
    for (size_t i = 0; i < u2.size(); i++) {
        u2[i] = indV[0][i];
    }
    Independent(u2);

    std::unique_ptr<CppAD::ADFun<CG<double> > > f2((*func2)(u2));

    vector<vector<double> > depsCG;

    std::string library = "./tmp/test_" + test + ".so";
    std::string function = "test_" + test;
    depsCG = run0(*f2.get(), library, function, indV, comparisons);

    /**
     * compare results
     */
    // Forward 0
    ASSERT_TRUE(compareValues(depsCG, depsDef, epsilonR, epsilonA));
}

} // END cg namespace
} // END CppAD namespace

#endif