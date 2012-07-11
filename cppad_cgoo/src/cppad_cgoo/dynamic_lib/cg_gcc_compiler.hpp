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
        std::set<std::string> _ofiles;
        std::set<std::string> _sfiles;
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

        virtual const std::set<std::string>& getObjectFiles() const {
            return _ofiles;
        }

        virtual const std::set<std::string>& getSourceFiles() const {
            return _sfiles;
        }

        /**
         * creates a dynamic library with the provided C source code
         * 
         * \param library the path of the dynamic library to be created
         * \param sources maps the names to the content of the source files
         * \param savefiles whether or not to save the content of the source 
         *                  files in the sources folder
         */
        virtual void compileSources(const std::map<std::string, std::string>& sources,
                                    bool savefiles) {

            if (savefiles) {
                system::createFolder(this->_sourcesFolder);

                std::map<std::string, std::string>::const_iterator it;
                for (it = sources.begin(); it != sources.end(); ++it) {
                    // for debugging purposes only
                    std::ofstream sourceFile;
                    std::string file = system::createPath(this->_sourcesFolder, it->first);
                    sourceFile.open(file.c_str());
                    sourceFile << it->second;
                    sourceFile.close();
                }
            }

            system::createFolder(this->_tmpFolder);

            // compile each source code file into a different object file
            size_t maxsize = 0;
            std::map<std::string, std::string>::const_iterator it;
            for (it = sources.begin(); it != sources.end(); ++it) {
                _sfiles.insert(it->first);
                std::string file = system::createPath(this->_tmpFolder, it->first + ".o");
                _ofiles.insert(file);
                maxsize = std::max(maxsize, file.size());
            }

            size_t count = 0;
            std::cout << std::endl;
            for (it = sources.begin(); it != sources.end(); ++it) {
                count++;
                std::string file = system::createPath(this->_tmpFolder, it->first + ".o");

                double beginTime = system::currentTime();
                std::cout << "[" << count << "/" << sources.size() << "] compiling '" << std::setw(maxsize + 8) << (file + "' ...  ");
                std::cout.flush();

                compile(it->second, file);

                double endTime = system::currentTime();
                std::cout << "done [" << (endTime - beginTime) << "]" << std::endl;
            }
        }

        /**
         * Creates a dynamic library from a set of object files
         * 
         * \param library the path to the dynamic library to be created
         */
        virtual void buildDynamic(const std::string& library) {

            std::vector<std::string> args;
            args.push_back("gcc");
            args.push_back("-O2"); // Optimization level
            args.push_back("-shared"); // Make shared object
            args.push_back("-Wl,-soname," + library); // Pass suitable options to linker
            args.push_back("-o"); // Output file name
            args.push_back(library); // Output file name
            std::set<std::string>::const_iterator it;
            for (it = _ofiles.begin(); it != _ofiles.end(); ++it) {
                args.push_back(*it);
            }

            std::cout << "building library" << std::endl;
            system::callExecutable(_gccPath, args);
        }

        virtual void cleanup() {
            // clean up
            std::set<std::string>::const_iterator it;
            for (it = _ofiles.begin(); it != _ofiles.end(); ++it) {
                if (remove(it->c_str()) != 0)
                    std::cerr << "Failed to delete temporary file '" << *it << "'" << std::endl;
            }
            _ofiles.clear();
            _sfiles.clear();

            remove(this->_tmpFolder.c_str());
        }

        virtual ~GccCompiler() {
            cleanup();
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

            system::callExecutable(_gccPath, args, true, source);
        }

    private:

        GccCompiler(const GccCompiler& orig); // not implemented
        GccCompiler& operator=(const GccCompiler& rhs); // not implemented
    };

}

#endif
