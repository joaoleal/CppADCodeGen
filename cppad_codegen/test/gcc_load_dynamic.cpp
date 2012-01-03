#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>

#include <cppad_codegen/cppad_codegen.hpp>
#include <cppad_codegen/local/ad_fun_code_gen.hpp>

#include "gcc_load_dynamic.hpp"
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
              const std::vector<std::vector<double> >& ind,
              double epsilonR, double epsilonA) {

    int comparisons;

    return run0(f, library, function, ind, comparisons, epsilonR, epsilonA);
}

bool run0(ADFunCodeGen<double>& f, const string& library, const string& function,
          const std::vector<std::vector<double> >& indV,
          int& comparisons, double epsilonR, double epsilonA) {

    stringstream code;
    f.ForwardCodeGen(0, code);

    CodeGenNameProvider<double>* n = f.getCodeGenNameProvider();

    string source = "#include <math.h>\n\n"
            "int " + function + "(const double* ind, double* dep) {\n";

    // declare variables
    source += n->baseTypeName() + " " + n->tempBaseVarName() + n->endl();
    source += "int " + n->tempIntegerVarName() + n->endl();
    source += "int " + n->compareChangeCounter() + " = 0 " + n->endl();

    const std::vector<VarID>& vars = n->getUsedVariables();
    source += n->baseTypeName() + " " + n->generateVarName(vars[0].order, vars[0].taddr);
    for (size_t i = 1; i < vars.size(); i++) {
        source += ", " + n->generateVarName(vars[i].order, vars[i].taddr);
    }
    source += n->endl();
    // set independent variable values
    const CppAD::vector<size_t>& indepAdd = f.IndependentTapeAddr();

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

    bool ok = true;

    std::vector<double> depCGen(depAdd.size());
    for (size_t i = 0; i < indV.size(); i++) {
        const std::vector<double>& ind = indV[i];
        assert(ind.size() == indepAdd.size());

        std::vector<double> dep = f.Forward(0, ind);

        comparisons = (*fn)(&ind[0], &depCGen[0]);

        if (!compareValues(depCGen, dep, "depCGEN", "depTape", epsilonR, epsilonA)) {
            ok = false;
        }
    }

    closeLibrary(libHandle);

    return ok;
}

bool runTestSparseJac(ADFunCodeGen<double>& f, const string& library,
                      const string& functionFor, const string& functionRev,
                      const std::vector<std::vector<double> >& indV,
                      double epsilonR, double epsilonA) {

    CodeGenNameProvider<double>* n = f.getCodeGenNameProvider();

    stringstream code;

    string source = "#include <math.h>\n\n";
    DiffMode mode;
    /**
     * forward mode
     */
    mode = FORWARD;
    f.SparseJacobianCodeGen<std::vector<bool> >(code, mode);

    source += "int " + functionFor + "(const double* ind, double* jac) {\n";

    // declare variables
    source += n->baseTypeName() + " " + n->tempBaseVarName() + n->endl();
    source += "int " + n->tempIntegerVarName() + n->endl();
    source += "int " + n->compareChangeCounter() + " = 0 " + n->endl();

    const std::vector<VarID>& vars = n->getUsedVariables();
    source += n->baseTypeName() + " " + n->generateVarName(vars[0].order, vars[0].taddr);
    for (size_t i = 1; i < vars.size(); i++) {
        source += ", " + n->generateVarName(vars[i].order, vars[i].taddr);
    }
    source += n->endl();
    // set independent variable values
    const CppAD::vector<size_t>& indepAdd = f.IndependentTapeAddr();
    for (size_t i = 0; i < indepAdd.size(); i++) {
        source += n->generateVarName(0, indepAdd[i]) + " = ind[" + n->toString(i) + "]" + n->endl();
    }

    source += code.str();

    // get dependent variable values
    //
    // currently the jacobian variable name is hardcoded
    source += "return " + n->compareChangeCounter() + n->endl();
    source += "}\n\n";

    /**
     * reverse mode
     */
    mode = REVERSE;
    
    code.str("");
    n->clearUsedVariables();

    f.SparseJacobianCodeGen<std::vector<bool> >(code, mode);

    source += "int " + functionRev + "(const double* ind, double* jac) {\n";

    // declare variables
    source += n->baseTypeName() + " " + n->tempBaseVarName() + n->endl();
    source += "int " + n->tempIntegerVarName() + n->endl();
    source += "int " + n->compareChangeCounter() + " = 0 " + n->endl();

    // taylor
    const std::vector<VarID>& tvars = n->getUsedVariables();
    if (tvars.size() > 0) {
        source += n->baseTypeName() + " " + n->generateVarName(tvars[0].order, tvars[0].taddr);
        for (size_t i = 1; i < tvars.size(); i++) {
            source += ", " + n->generateVarName(tvars[i].order, tvars[i].taddr);
        }
        source += n->endl();
    }
    // partials
    const std::vector<VarID>& pvars = n->getUsedPartials();
    source += n->baseTypeName() + " " + n->generatePartialName(pvars[0].order, pvars[0].taddr);
    for (size_t i = 1; i < pvars.size(); i++) {
        source += ", " + n->generatePartialName(pvars[i].order, pvars[i].taddr);
    }
    source += n->endl();

    // set independent variable values
    for (size_t i = 0; i < indepAdd.size(); i++) {
        source += n->generateVarName(0, indepAdd[i]) + " = ind[" + n->toString(i) + "]" + n->endl();
    }

    source += code.str();

    // get dependent variable values
    //
    // currently the jacobian variable name is hardcoded
    source += "return " + n->compareChangeCounter() + n->endl();
    source += "}\n";

    /**
     * Compile
     */
    void* libHandle;
    try {
        compile(source, library);

        libHandle = loadLibrary(library);
    } catch (const exception& ex) {
        cerr << ex.what();
        return false;
    }

    int (*fnFor)(const double*, double*) = NULL;
    int (*fnRev)(const double*, double*) = NULL;

    try {
        *(void **) (&fnFor) = getFunction(libHandle, functionFor);
        *(void **) (&fnRev) = getFunction(libHandle, functionRev);
    } catch (const exception& ex) {
        cerr << ex.what();
        closeLibrary(libHandle);
        return false;
    }

    bool ok = true;
    for (size_t i = 0; i < indV.size(); i++) {
        const std::vector<double>& ind = indV[i];
        assert(ind.size() == indepAdd.size());

        /**
         * Jacobian
         */
        CodeGenNameProvider<double>* n = f.getCodeGenNameProvider();
        n->clearUsedVariables();

        const std::vector<double> jac = f.SparseJacobian(ind);

        /**
         * test forward jacobian
         */
        std::vector<double> jacOut(jac.size());
        (*fnFor)(&ind[0], &jacOut[0]);

        if (!compareValues(jacOut, jac, "jacCGEN", "jacTape", epsilonR, epsilonA)) {
            cerr << "Forward mode failed (dataset: " << i << ")" << endl;
            ok = false;
        }

        /**
         * test reverse jacobian
         */
        (*fnRev)(&ind[0], &jacOut[0]);

        if (!compareValues(jacOut, jac, "jacCGEN", "jacTape", epsilonR, epsilonA)) {
            cerr << "Reverse mode failed (dataset: " << i << ")" << endl;
            ok = false;
        }
    }

    closeLibrary(libHandle);

    return ok;
}

bool test0nJac(const std::string& test,
               CppAD::ADFunCodeGen<double>& f,
               const std::vector<double>& ind,
               CPPAD_TEST_VECTOR<CppAD::AD<double> > w,
               double epsilonR, double epsilonA) {

    std::vector<std::vector<double> > indV;
    indV.push_back(ind);

    return test0nJac(test, f, indV, epsilonR, epsilonA);
}

bool test0nJac(const string& test, ADFunCodeGen<double>& f,
               const std::vector<std::vector<double> >& indV,
               double epsilonR, double epsilonA) {

    bool ok = true;
    /**
     * forward zero order mode
     */
    // forward computation of partials w.r.t. s
    // generate the code
    string library = "./tmp/test_" + test + ".so";
    string function = "test_" + test;

    ok &= runTest0(f, library, function, indV, epsilonR, epsilonA);

    library = "./tmp/test_" + test + "_jac.so";
    string functionFor = "test_" + test + "_jac_for";
    string functionRev = "test_" + test + "_jac_rev";
    ok &= runTestSparseJac(f, library, functionFor, functionRev, indV, epsilonR, epsilonA);

    return ok;
}