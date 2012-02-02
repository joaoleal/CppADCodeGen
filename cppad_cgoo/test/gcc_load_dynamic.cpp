#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <vector>

#include "gcc_load_dynamic.hpp"
#include "cppad/local/ad_fun.hpp"
#include "cppad/vector.hpp"
#include "cppad_cgoo/cg.hpp"

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
         *   -O0                   Optimization level
         *   -x c                  C source
         *   -pipe                 Use pipes between gcc stages
         *   -fPIC -shared         Make shared object
         *   -Wl,-soname, library  Pass suitable options to linker
         * 
         */
        string linker = "-Wl,-soname," + library;
        execl("/usr/bin/gcc", "gcc", "-x", "c", "-O0", "-pipe", "-", "-fPIC", "-shared",
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

std::vector<std::vector<double> > runDefault0(ADFun<double>& f,
                                              const std::vector<std::vector<double> >& ind) {
    int comparisons;

    return runDefault0(f, ind, comparisons);
}

std::vector<std::vector<double> > runDefault0(ADFun<double>& f,
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

std::vector<std::vector<double> > run0(ADFun<CG<double> >& f,
                                       const string& library, const string& function,
                                       const std::vector<std::vector<double> >& ind) {

    int comparisons;

    return run0(f, library, function, ind, comparisons);
}

std::vector<std::vector<double> > run0(ADFun<CG<double> >& f,
                                       const string& library, const string& function,
                                       const std::vector<std::vector<double> >& indV,
                                       int& comparisons) {
    using std::vector;

    ostringstream code;
    CodeHandler<double> handler(code);

    vector<CG<double> > indVars(indV.begin()->size());
    handler.makeVariables(indVars);

    vector<CG<double> > dep = f.Forward(0, indVars);

    string spaces = "   ";
    string source = "#include <math.h>\n\n"
            "int " + function + "(const double* ind, double* dep) {\n";

    // declare variables
    const size_t maxID = handler.getMaximumVariableID();
    source += spaces + "double " + handler.createVariableName(1);
    for (size_t i = 2; i <= maxID; i++) {
        source += ", " + handler.createVariableName(i);
    }
    source += ";\n";
    // set independent variable values
    for (size_t i = 0; i < indVars.size(); i++) {
        source += spaces + handler.createVariableName(indVars[i]) + " = ind[" + handler.toString(i) + "];\n";
    }

    source += code.str();

    // get dependent variable values
    for (size_t i = 0; i < dep.size(); i++) {
        source += spaces + "dep[" + handler.toString(i) + "] = " + handler.createVariableName(dep[i].getVariableID()) + ";\n";
    }

    //source += "return " + n->compareChangeCounter() + n->endl();
    source += spaces + "return 0;\n";
    source += "}";

    cout << endl << source << endl;

    compile(source, library);

    void* libHandle = loadLibrary(library);

    int (*fn)(const double*, double*) = NULL;

    try {
        *(void **) (&fn) = getFunction(libHandle, function);
    } catch (const exception& ex) {
        closeLibrary(libHandle);
        throw;
    }

    vector<vector<double> >depCGen(indV.size());
    for (size_t i = 0; i < indV.size(); i++) {

        vector<double>& depi = depCGen[i];
        depi.resize(dep.size());

        const vector<double>& ind = indV[i];

        // the compiled version
        comparisons = (*fn)(&ind[0], &depi[0]);

        for (size_t j = 0; j < ind.size(); j++) {
            cout << " ind[" << j << "] = " << ind[j] << "\n";
        }

        for (size_t j = 0; j < depi.size(); j++) {
            cout << " dep[" << j << "] = " << depi[j] << "\n";
        }
    }

    closeLibrary(libHandle);

    return depCGen;
}

std::vector<std::vector<double> > runSparseJacDefault(CppAD::ADFun<double>& f,
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

std::vector<std::vector<double> > runSparseJac(ADFun<CG<double> >& f,
                                               const string& library,
                                               const string& functionJac,
                                               const std::vector<std::vector<double> >& indV) {
    assert(!indV.empty());

    using std::vector;

    ostringstream code;
    CodeHandler<double> handler(code);

    vector<CG<double> > indVars(indV[0].size());
    handler.makeVariables(indVars);

    string source = "#include <math.h>\n\n";
    string spaces = "   ";

    vector<CG<double> > jacCG = f.SparseJacobian(indVars);

    source += "int " + functionJac + "(const double* ind, double* jac) {\n";

    // declare variables
    const size_t maxID = handler.getMaximumVariableID();
    source += spaces + "double " + handler.createVariableName(1);
    for (size_t i = 2; i <= maxID; i++) {
        source += ", " + handler.createVariableName(i);
    }
    source += ";\n";

    // set independent variable values
    for (size_t i = 0; i < indVars.size(); i++) {
        source += spaces + handler.createVariableName(indVars[i]) + " = ind[" + handler.toString(i) + "];\n";
    }

    source += code.str();

    for (size_t i = 0; i < jacCG.size(); i++) {
        source += spaces + "jac[" + handler.toString(i) + "] = " + handler.operations(jacCG[i]) + ";\n";
    }


    // get dependent variable values
    //
    // currently the jacobian variable name is hardcoded
    //source += "return " + n->compareChangeCounter() + n->endl();
    source += spaces + "return 0;\n";
    source += "}\n\n";

    cout << endl << source << endl;

    /**
     * Compile
     */
    compile(source, library);

    void* libHandle = loadLibrary(library);

    int (*fn)(const double*, double*) = NULL;

    try {
        *(void **) (&fn) = getFunction(libHandle, functionJac);
    } catch (const exception& ex) {
        closeLibrary(libHandle);
        throw;
    }

    vector<vector<double> >jac(indV.size());

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

bool test0nJac(const std::string& test,
               ADFun<double>* (*func1)(const std::vector<AD<double> >&),
               ADFun<CG<double> >* (*func2)(const std::vector<AD<CG<double> > >&),
               const std::vector<double>& ind,
               double epsilonR, double epsilonA) {

    std::vector<std::vector<double> > indV;
    indV.push_back(ind);

    return test0nJac(test, func1, func2, indV, epsilonR, epsilonA);
}

bool test0nJac(const string& test,
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
    Independent(u2);

    CppAD::ADFun<CG<double> >* f2 = (*func2)(u2);

    vector<vector<double> > depsCG;
    vector<vector<double> > jacCG;
    try {
        string library = "./tmp/test_" + test + ".so";
        string function = "test_" + test;
        depsCG = run0(*f2, library, function, indV);

        library = "./tmp/test_" + test + "_jac.so";
        string functionJac = "test_" + test + "_jac_for";
        jacCG = runSparseJac(*f2, library, functionJac, indV);
    } catch (const exception& ex) {
        cerr << ex.what() << std::endl;
        delete f2;
        return false;
    }

    delete f2;

    /**
     * compare results
     */
    bool ok = true;

    ok &= compareValues("Forward 0", depsCG, depsDef, epsilonR, epsilonA);

    ok &= compareValues("Jacobian", jacCG, jacDef, epsilonR, epsilonA);


    return ok;
}

bool test0(const string& test,
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
    Independent(u2);

    CppAD::ADFun<CG<double> >* f2 = (*func2)(u2);

    vector<vector<double> > depsCG;
    try {
        string library = "./tmp/test_" + test + ".so";
        string function = "test_" + test;
        depsCG = run0(*f2, library, function, indV, comparisons);
    } catch (const exception& ex) {
        cerr << ex.what() << std::endl;
        delete f2;
        return false;
    }

    delete f2;

    /**
     * compare results
     */
    bool ok = true;

    ok &= compareValues("Forward 0", depsCG, depsDef, epsilonR, epsilonA);

    return ok;
}

bool compareValues(const std::string& testType,
                   const std::vector<std::vector<double> >& depCGen,
                   const std::vector<std::vector<double> >& dep,
                   double epsilonR, double epsilonA) {

    assert(depCGen.size() == dep.size());

    bool ok = true;

    for (size_t i = 0; i < depCGen.size(); i++) {
        if (!compareValues(depCGen[i], dep[i], "depCGEN", "depTape", epsilonR, epsilonA)) {
            cerr << testType << " failed (dataset: " << i << ")" << endl;
            ok = false;
        }
    }

    return ok;
}