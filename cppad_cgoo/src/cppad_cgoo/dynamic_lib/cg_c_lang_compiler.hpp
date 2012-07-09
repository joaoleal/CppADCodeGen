#ifndef CPPAD_CG_C_LANG_COMPILER_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILER_INCLUDED
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
    class CLangCompiler {
    protected:
        std::string _sourcesFolder;
        std::string _tmpFolder;
    public:

        CLangCompiler() :
            _sourcesFolder("cppadcg_sources"),
            _tmpFolder("cppadcg_tmp") {
        }

        /**
         * Provides the path to the folder where source files are saved
         * (if requested).
         * 
         * \return path to the folder with the generated source files 
         */
        inline const std::string& getSourcesFolder() const {
            return _sourcesFolder;
        }

        /**
         * Defines the path to the folder where source files are saved
         * (if requested).
         * 
         * \param sourcesFolder path to the folder where the generated source 
         *                      files are saved
         */
        inline void setSourcesFolder(const std::string& sourcesFolder) {
            _sourcesFolder = sourcesFolder;
        }

        /**
         * Provides the path to a temporary folder that should not exist
         * (it will be deleted after the dynamic library is created)
         * 
         * \return path to a temporary folder.
         */
        inline const std::string& getTemporaryFolder() const {
            return _tmpFolder;
        }

        /**
         * Defines the path to a temporary folder that should not exist
         * (it will be deleted after the dynamic library is created)
         * 
         * \param tmpFolder path to a temporary folder.
         */
        inline void setTemporaryFolder_(const std::string& tmpFolder) {
            _tmpFolder = tmpFolder;
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
                                    bool savefiles) = 0;

    protected:

        /**
         * creates a new folder (system dependent)
         * 
         * \param folder the path to the folder
         */
        static void createFolder(const std::string& folder);

        /**
         * Creates a new path (system dependent)
         * 
         * \param baseFolder the path to the base folder
         * \param file the file or folder name inside the base folder
         * \return the new path
         */
        static std::string createPath(const std::string& baseFolder, const std::string& file);

        /**
         * Calls an external executable (system dependent)
         * 
         * @param executable the executable path
         * @param args the command line arguments to the executable
         * @param pipe whether or not to create a pipe to the executable
         * @param message the information to pass in the pipe
         */
        static void callExecutable(const std::string& executable,
                                   const std::vector<std::string>& args,
                                   bool pipe = false,
                                   const std::string& message = "");
        
        /**
         * Escapes a file or folder path (system dependent)
         * 
         * \param path the file/folder path
         * \return the escaped file/folder path
         */
        static std::string escapePath(const std::string& path);
    };

}

#endif
