#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

#include <cppad/cppad_code_gen.hpp>

#include "gcc_load_dynamic.hpp"
#include "cppad/local/code_gen/ad_fun_code_gen.hpp"
#include "cppad/local/ad_fun.hpp"
#include "cppad/vector.hpp"

using namespace std;
using namespace CppAD;

void* loadLibrary(const string& library) throw (CppAD::TestException) {
    void * libHandle = dlopen(library.c_str(), RTLD_NOW);
    if (!libHandle) {
        throw CppAD::TestException("Failed to dynamically load library");
    }
    return libHandle;
}

void closeLibrary(void* libHandle) {
    dlclose(libHandle);
}

void* getFunction(void * libHandle, const string& functionName) throw (CppAD::TestException) {
    void* functor = dlsym(libHandle, functionName.c_str());
    char *error;
    if ((error = dlerror()) != NULL) {
        throw CppAD::TestException(error);
    }
    return functor;
}

void compile(const string& source, const string& library) throw (CppAD::TestException) {
    int fd[2];
    //Create pipe for piping source to gcc
    if (pipe(fd) < 0) {
        throw CppAD::TestException("Failed to create pipe");
    }

    //Fork a gcc, pipe source to it, wait for gcc to exit
    pid_t pid = fork();
    if (pid < 0) {
        throw CppAD::TestException("Failed to fork program");
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
         *   -O2                   Optimization level
         *   -x c                  C source
         *   -pipe                 Use pipes between gcc stages
         *   -fPIC -shared         Make shared object
         *   -Wl,-soname, library  Pass suitable options to linker
         * 
         */
        string linker = "-Wl,-soname," + library;
        execl("/usr/bin/gcc", "gcc", "-x", "c", "-O2", "-pipe", "-", "-fPIC", "-shared",
                linker.c_str(), "-o", library.c_str(), (char *) NULL);

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
        throw CppAD::TestException("Failed while waiting for gcc");
    }
}

bool runTest0(ADFunCodeGen<double>& f, const string& library, const string& function,
        const std::vector<double>& ind, const CPPAD_TEST_VECTOR< AD<double> >& dep,
        double epsilonR, double epsilonA) {

    assert(dep.size() == f.Range());
    std::vector<double> depCGen(dep.size());

    int comparisons;

    if (!run0(f, library, function, ind, comparisons, depCGen)) {
        return false;
    }

    return compareValues(depCGen, dep, "depCGEN", "depTape", epsilonR, epsilonA);
}

bool run0(ADFunCodeGen<double>& f, const string& library, const string& function,
        const std::vector<double>& ind,
        int& comparisons, std::vector<double>& depCGen) {

    stringstream code;
    f.ForwardCodeGen(0, code);

    CodeGenNameProvider<double>* n = f.getCodeGenNameProvider();

    string source = "#include <math.h>\n\n"
            "int " + function + "(const double* ind, double* dep) {\n";

    // declare variables
    source += n->baseTypeName() + " " + n->tempBaseVarName() + n->endl();
    source += "int " + n->compareChangeCounter() + " = 0 " + n->endl();

    const std::vector<VarID>& vars = n->getUsedVariables();
    source += n->baseTypeName() + " " + n->generateVarName(vars[0].order, vars[0].taddr);
    for (size_t i = 1; i < vars.size(); i++) {
        source += ", " + n->generateVarName(vars[i].order, vars[i].taddr);
    }
    source += n->endl();
    // set independent variable values
    const CppAD::vector<size_t>& indepAdd = f.IndependentTapeAddr();
    assert(ind.size() == indepAdd.size());

    for (size_t i = 0; i < indepAdd.size(); i++) {
        source += n->generateVarName(0, indepAdd[i]) + " = ind[" + n->toString(i) + "]" + n->endl();
    }

    source += code.str();

    // get dependent variable values
    const CppAD::vector<size_t>& depAdd = f.DependentTapeAddr();

    for (size_t i = 0; i < depAdd.size(); i++) {
        source += "dep[" + n->toString(i) + "] = " + n->generateVarName(0, depAdd[i]) + n->endl();
    }

    source += "return " + n->compareChangeCounter() + n->endl();
    source += "}";

    void* libHandle;
    try {
        compile(source, library);

        libHandle = loadLibrary(library);
    } catch (const exception& ex) {
        cerr << ex.what() << std::endl;
        return false;
    }

    int (*fn)(const double*, double*) = NULL;

    try {
        *(void **) (&fn) = getFunction(libHandle, function);
    } catch (const exception& ex) {
        cerr << ex.what() << std::endl;
        closeLibrary(libHandle);
        return false;
    }

    comparisons = (*fn)(&ind[0], &depCGen[0]);

    closeLibrary(libHandle);

    return true;
}

bool runTestSparseJac(ADFunCodeGen<double>& f, const string& library, const string& function,
        const std::vector<double>& ind, const std::vector<double>& jac,
        double epsilonR, double epsilonA) {

    stringstream code;
    f.SparseJacobianCodeGen(code);

    CodeGenNameProvider<double>* n = f.getCodeGenNameProvider();

    string source = "#include <math.h>\n\n"
            "int " + function + "(const double* ind, double* jac) {\n";

    // declare variables
    source += n->baseTypeName() + " " + n->tempBaseVarName() + n->endl();
    source += "int " + n->compareChangeCounter() + " = 0 " + n->endl();

    const std::vector<VarID>& vars = n->getUsedVariables();
    source += n->baseTypeName() + " " + n->generateVarName(vars[0].order, vars[0].taddr);
    for (size_t i = 1; i < vars.size(); i++) {
        source += ", " + n->generateVarName(vars[i].order, vars[i].taddr);
    }
    source += n->endl();
    // set independent variable values
    const CppAD::vector<size_t>& indepAdd = f.IndependentTapeAddr();
    assert(ind.size() == indepAdd.size());

    for (size_t i = 0; i < indepAdd.size(); i++) {
        source += n->generateVarName(0, indepAdd[i]) + " = ind[" + n->toString(i) + "]" + n->endl();
    }

    source += code.str();

    // get dependent variable values
    //
    // currently the jacobian variable name is hardcoded
    source += "return " + n->compareChangeCounter() + n->endl();
    source += "}";


    void* libHandle;
    try {
        compile(source, library);

        libHandle = loadLibrary(library);
    } catch (const exception& ex) {
        cerr << ex.what();
        return false;
    }

    int (*fn)(const double*, double*) = NULL;

    try {
        *(void **) (&fn) = getFunction(libHandle, function);
    } catch (const exception& ex) {
        cerr << ex.what();
        closeLibrary(libHandle);
        return false;
    }

    std::vector<double> jacOut(jac.size());
    (*fn)(&ind[0], &jacOut[0]);

    closeLibrary(libHandle);

    return compareValues(jacOut, jac, "jacCGEN", "jacTape", epsilonR, epsilonA);
}

bool test0nJac(const string& test, ADFunCodeGen<double>& f,
        const std::vector<double>& ind, const CPPAD_TEST_VECTOR< AD<double> >& dep,
        double epsilonR, double epsilonA) {

    bool ok = true;
    /**
     * forward zero order mode
     */
    // forward computation of partials w.r.t. s
    // generate the code
    string library = "./tmp/test_" + test + ".so";
    string function = "test_" + test;

    ok &= runTest0(f, library, function, ind, dep, epsilonR, epsilonA);

    /**
     * Jacobian
     */
    CodeGenNameProvider<double>* n = f.getCodeGenNameProvider();
    n->clearUsedVariables();

    std::vector<double> jac = f.SparseJacobian(ind);

    library = "./tmp/test_" + test + "_jac.so";
    function = "test_" + test + "_jac";
    ok &= runTestSparseJac(f, library, function, ind, jac, epsilonR, epsilonA);

    return ok;
}