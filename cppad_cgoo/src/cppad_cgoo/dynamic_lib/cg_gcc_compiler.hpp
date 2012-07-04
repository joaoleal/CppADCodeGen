#ifndef CPPAD_CG_GCC_COMPILER_INCLUDED
#define	CPPAD_CG_GCC_COMPILER_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iosfwd>
#include <errno.h>
#include <string.h>

namespace CppAD {

    /**
     * C compiler class used to create a dynamic library
     * 
     * \author Joao Leal
     */
    template<class Base>
    class GccCompiler : public CLangCompiler<Base> {
    protected:
        std::string _gccPath;
        std::vector<std::string> _ofiles;
    public:

        GccCompiler() :
            _gccPath("/usr/bin/gcc") {
        }

        GccCompiler(const std::string& gccPath) :
            _gccPath(gccPath) {
        }

        std::string getGccPath() const {
            return _gccPath;
        }

        void setGccPath(const std::string& gccPath) {
            _gccPath = gccPath;
        }

        /**
         * creates a dynamic library with the provided C source code
         * 
         * \param library the path of the dynamic library to be created
         * \param sources maps the names to the content of the source files
         * \param savefiles whether or not to save the content of the source 
         *                  files in the sources folder
         */
        virtual void compileDynamic(const std::string& library,
                                    const std::map<std::string, std::string>& sources,
                                    bool savefiles) {

            if (savefiles) {
                CLangCompiler<Base>::createFolder(this->_sourcesFolder);

                std::map<std::string, std::string>::const_iterator it;
                for (it = sources.begin(); it != sources.end(); ++it) {
                    // for debugging purposes only
                    std::ofstream sourceFile;
                    std::string file = createPath(this->_sourcesFolder, it->first);
                    sourceFile.open(file.c_str());
                    sourceFile << it->second;
                    sourceFile.close();
                }
            }

            try {
                CLangCompiler<Base>::createFolder(this->_tmpFolder);

                // compile each source code file into a different object file
                _ofiles.clear();

                std::map<std::string, std::string>::const_iterator it;
                size_t count = 0;
                std::cout << std::endl;
                for (it = sources.begin(); it != sources.end(); ++it) {
                    count++;
                    _ofiles.push_back(createPath(this->_tmpFolder, it->first + ".o"));
                    std::cout << "[" << count << "/" << sources.size() << "] compiling '" << _ofiles.back() << "' ... ";
                    std::cout.flush();
                    compile(it->second, _ofiles.back());
                    std::cout << "done" << std::endl;
                }

                // make dynamic library
                buildDynamic(library);

            } catch (...) {
                cleanup();
                throw;
            }
            cleanup();
        }

        virtual ~GccCompiler() {
        }

    protected:

        /**
         * Compiles a single source file into an object file
         * 
         * \param source the content of the source file
         * \param output the compiled output file name (the object file path)
         */
        virtual void compile(const std::string& source, const std::string& output) {
            std::vector<std::string> args;
            args.push_back("gcc");
            args.push_back("-x");
            args.push_back("c");
            args.push_back("-O2");
            args.push_back("-c");
            args.push_back("-");
            args.push_back("-fPIC");
            args.push_back("-o");
            args.push_back(output);

            CLangCompiler<Base>::callExecutable(_gccPath, args, true, source);
        }

        /**
         * Creates a dynamic library from a set of object files
         * 
         * \param library the path to the dynamic library to be created
         * \param files the object files to include in the dynamic library
         */
        virtual void buildDynamic(const std::string& library) {

            std::vector<std::string> args;
            args.push_back("gcc");
            args.push_back("-O2"); // Optimization level
            args.push_back("-shared"); // Make shared object
            args.push_back("-Wl,-soname," + library); // Pass suitable options to linker
            args.push_back("-o"); // Output file name
            args.push_back(library); // Output file name
            for (size_t i = 0; i < _ofiles.size(); i++) {
                args.push_back(_ofiles[i]);
            }

            std::cout << "building library" << std::endl;
            CLangCompiler<Base>::callExecutable(_gccPath, args);
        }

        virtual void cleanup() {
            // clean up
            for (size_t i = 0; i < _ofiles.size(); i++) {
                if (remove(_ofiles[i].c_str()) != 0)
                    std::cerr << "Failed to delete temporary file '" << _ofiles[i] << "'" << std::endl;
            }
            _ofiles.clear();

            remove(this->_tmpFolder.c_str());
        }

    private:

        GccCompiler(const GccCompiler& orig); // not implemented
        GccCompiler& operator=(const GccCompiler& rhs); // not implemented
    };

}

#endif
