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
    class GccCompiler : public CLangCompiler<Base> {
    protected:
        std::string _gccPath; // the path to the gcc executable
        std::set<std::string> _ofiles; // compiled object files
        std::set<std::string> _sfiles; // compiled source files
        std::vector<std::string> _compileFlags;
        std::vector<std::string> _compileLibFlags;
        std::vector<std::string> _linkFlags;
        bool _verbose;
    public:

        GccCompiler() :
            _gccPath("/usr/bin/gcc"),
            _verbose(false) {

            _compileFlags.push_back("-O2"); // Optimization level
            _compileLibFlags.push_back("-O2"); // Optimization level
            _compileLibFlags.push_back("-shared"); // Make shared object
            _compileLibFlags.push_back("-rdynamic"); // add all symbols to the dynamic symbol table

        }

        GccCompiler(const std::string& gccPath) :
            _gccPath(gccPath),
            _verbose(false) {

            _compileFlags.push_back("-O2"); // Optimization level
            _compileLibFlags.push_back("-O2"); // Optimization level
            _compileLibFlags.push_back("-shared"); // Make shared object
            _compileLibFlags.push_back("-rdynamic"); // add all symbols to the dynamic symbol table
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

        const std::vector<std::string>& getCompileFlags() const {
            return _compileFlags;
        }

        void setCompileFlags(const std::vector<std::string>& compileFlags) {
            _compileFlags = compileFlags;
        }

        const std::vector<std::string>& getLinkFlags() const {
            return _linkFlags;
        }

        void setLinkFlags(const std::vector<std::string>& linkFlags) {
            _linkFlags = linkFlags;
        }

        const std::vector<std::string>& getCompileLibFlags() const {
            return _compileLibFlags;
        }

        void setCompileLibFlags(const std::vector<std::string>& compileLibFlags) {
            _compileLibFlags = compileLibFlags;
        }

        virtual bool isVerbose() const {
            return _verbose;
        }

        virtual void setVerbose(bool verbose) {
            _verbose = verbose;
        }

        /**
         * creates a dynamic library with the provided C source code
         * 
         * @param library the path of the dynamic library to be created
         * @param sources maps the names to the content of the source files
         * @param posIndepCode whether or not to create position-independent
         *                     code for dynamic linking
         * @param savefiles whether or not to save the content of the source 
         *                  files in the sources folder
         */
        virtual void compileSources(const std::map<std::string, std::string>& sources,
                                    bool posIndepCode,
                                    bool savefiles,
                                    JobTimer* timer = NULL) {

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

            size_t countWidth = std::ceil(std::log10(sources.size()));

            size_t count = 0;
            if (timer != NULL) {
                size_t ms = 3 + 2 * countWidth + 1 + JobTypeHolder<>::COMPILING.getActionName().size() + 2 + maxsize + 5;
                ms += timer->getJobCount() * 2;
                if (timer->getMaxLineWidth() < ms)
                    timer->setMaxLineWidth(ms);
            } else if (_verbose) {
                std::cout << std::endl;
            }

            std::ostringstream os;


            for (it = sources.begin(); it != sources.end(); ++it) {
                count++;
                std::string file = system::createPath(this->_tmpFolder, it->first + ".o");

                double beginTime = 0.0;

                if (timer != NULL || _verbose) {
                    os << "[" << std::setw(countWidth) << std::setfill(' ') << std::right << count
                            << "/" << sources.size() << "]";
                }

                if (timer != NULL) {
                    timer->startingJob("'" + file + "'", JobTypeHolder<>::COMPILING, os.str());
                    os.str("");
                } else if (_verbose) {
                    beginTime = system::currentTime();
                    char f = std::cout.fill();
                    std::cout << os.str() << " compiling "
                            << std::setw(maxsize + 9) << std::setfill('.') << std::left
                            << ("'" + file + "' ") << " ";
                    os.str("");
                    std::cout.flush();
                    std::cout.fill(f); // restore fill character
                }

                compile(it->second, file, posIndepCode);


                if (timer != NULL) {
                    timer->finishedJob();
                } else if (_verbose) {
                    double endTime = system::currentTime();
                    std::cout << "done [" << std::fixed << std::setprecision(3)
                            << (endTime - beginTime) << "]" << std::endl;
                }

            }

        }

        /**
         * Creates a dynamic library from a set of object files
         * 
         * @param library the path to the dynamic library to be created
         */
        virtual void buildDynamic(const std::string& library,
                                  JobTimer* timer = NULL) {

            std::string linkerFlags = "-Wl,-soname," + system::filenameFromPath(library);
            for (size_t i = 0; i < _linkFlags.size(); i++)
                linkerFlags += "," + _linkFlags[i];

            std::vector<std::string> args;
            args.push_back("gcc");
            args.insert(args.end(), _compileLibFlags.begin(), _compileLibFlags.end());
            args.push_back(linkerFlags); // Pass suitable options to linker
            args.push_back("-o"); // Output file name
            args.push_back(library); // Output file name
            std::set<std::string>::const_iterator it;
            for (it = _ofiles.begin(); it != _ofiles.end(); ++it) {
                args.push_back(*it);
            }

            if (timer != NULL) {
                timer->startingJob("'" + library + "'", JobTimer::COMPILING_DYNAMIC_LIBRARY);
            } else if (_verbose) {
                std::cout << "building library '" << library << "'" << std::endl;
            }

            system::callExecutable(_gccPath, args);

            if (timer != NULL) {
                timer->finishedJob();
            }
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
         * @param source the content of the source file
         * @param output the compiled output file name (the object file path)
         */
        virtual void compile(const std::string& source, const std::string& output, bool posIndepCode) {
            std::vector<std::string> args;
            args.push_back("gcc");
            args.push_back("-x");
            args.push_back("c"); // C source files
            args.insert(args.end(), _compileFlags.begin(), _compileFlags.end());
            args.push_back("-c");
            args.push_back("-");
            if (posIndepCode) {
                args.push_back("-fPIC"); // position-independent code for dynamic linking
            }
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
