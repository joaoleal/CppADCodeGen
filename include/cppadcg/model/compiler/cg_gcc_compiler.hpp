#ifndef CPPAD_CG_GCC_COMPILER_INCLUDED
#define CPPAD_CG_GCC_COMPILER_INCLUDED
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

namespace CppAD {

    /**
     * C compiler class used to create a dynamic library
     * 
     * @author Joao Leal
     */
    template<class Base>
    class GccCompiler : public AbstractCCompiler<Base> {
    public:

        GccCompiler() :
            AbstractCCompiler<Base>("/usr/bin/gcc") {

            this->_compileFlags.push_back("-O2"); // Optimization level
            this->_compileLibFlags.push_back("-O2"); // Optimization level
            this->_compileLibFlags.push_back("-shared"); // Make shared object
            this->_compileLibFlags.push_back("-rdynamic"); // add all symbols to the dynamic symbol table

        }

        GccCompiler(const std::string& gccPath) :
            AbstractCCompiler<Base>(gccPath) {

            this->_compileFlags.push_back("-O2"); // Optimization level
            this->_compileLibFlags.push_back("-O2"); // Optimization level
            this->_compileLibFlags.push_back("-shared"); // Make shared object
            this->_compileLibFlags.push_back("-rdynamic"); // add all symbols to the dynamic symbol table
        }

        /**
         * Creates a dynamic library from a set of object files
         * 
         * @param library the path to the dynamic library to be created
         */
        virtual void buildDynamic(const std::string& library,
                                  JobTimer* timer = nullptr) override {

            std::string linkerFlags = "-Wl,-soname," + system::filenameFromPath(library);
            for (size_t i = 0; i < this->_linkFlags.size(); i++)
                linkerFlags += "," + this->_linkFlags[i];

            std::vector<std::string> args;
            args.push_back("gcc");
            args.insert(args.end(), this->_compileLibFlags.begin(), this->_compileLibFlags.end());
            args.push_back(linkerFlags); // Pass suitable options to linker
            args.push_back("-o"); // Output file name
            args.push_back(library); // Output file name
            std::set<std::string>::const_iterator it;
            for (it = this->_ofiles.begin(); it != this->_ofiles.end(); ++it) {
                args.push_back(*it);
            }

            if (timer != nullptr) {
                timer->startingJob("'" + library + "'", JobTimer::COMPILING_DYNAMIC_LIBRARY);
            } else if (this->_verbose) {
                std::cout << "building library '" << library << "'" << std::endl;
            }

            system::callExecutable(this->_path, args);

            if (timer != nullptr) {
                timer->finishedJob();
            }
        }

        virtual ~GccCompiler() {
        }

    protected:

        /**
         * Compiles a single source file into an object file
         * 
         * @param source the content of the source file
         * @param output the compiled output file name (the object file path)
         */
        virtual void compile(const std::string& source, const std::string& output, bool posIndepCode) override {
            std::vector<std::string> args;
            args.push_back("gcc");
            args.push_back("-x");
            args.push_back("c"); // C source files
            args.insert(args.end(), this->_compileFlags.begin(), this->_compileFlags.end());
            args.push_back("-c");
            args.push_back("-");
            if (posIndepCode) {
                args.push_back("-fPIC"); // position-independent code for dynamic linking
            }
            args.push_back("-o");
            args.push_back(output);

            system::callExecutable(this->_path, args, true, source);
        }

    private:

        GccCompiler(const GccCompiler& orig); // not implemented
        GccCompiler& operator=(const GccCompiler& rhs); // not implemented
    };

}

#endif
