#ifndef CPPAD_CG_C_LANG_COMPILER_INCLUDED
#define CPPAD_CG_C_LANG_COMPILER_INCLUDED
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
         * @return path to the folder with the generated source files 
         */
        inline const std::string& getSourcesFolder() const {
            return _sourcesFolder;
        }

        /**
         * Defines the path to the folder where source files are saved
         * (if requested).
         * 
         * @param sourcesFolder path to the folder where the generated source 
         *                      files are saved
         */
        inline void setSourcesFolder(const std::string& sourcesFolder) {
            _sourcesFolder = sourcesFolder;
        }

        /**
         * Provides the path to a temporary folder that should not exist
         * (it will be deleted after the dynamic library is created)
         * 
         * @return path to a temporary folder.
         */
        inline const std::string& getTemporaryFolder() const {
            return _tmpFolder;
        }

        /**
         * Defines the path to a temporary folder that should not exist
         * (it will be deleted after the dynamic library is created)
         * 
         * @param tmpFolder path to a temporary folder.
         */
        inline void setTemporaryFolder(const std::string& tmpFolder) {
            _tmpFolder = tmpFolder;
        }

        virtual const std::set<std::string>& getObjectFiles() const = 0;

        virtual const std::set<std::string>& getSourceFiles() const = 0;

        virtual bool isVerbose() const = 0;

        virtual void setVerbose(bool verbose) = 0;

        /**
         * creates a dynamic library with the provided C source code
         * 
         * @param sources maps the names to the content of the source files
         * @param posIndepCode whether or not to create position-independent
         *                     code for dynamic linking
         * @param savefiles whether or not to save the content of the source 
         *                  files in the sources folder
         */
        virtual void compileSources(const std::map<std::string, std::string>& sources,
                                    bool posIndepCode,
                                    bool savefiles,
                                    JobTimer* timer = NULL) = 0;

        /**
         * Creates a dynamic library from the previously compiled object files
         * 
         * @param library the path to the dynamic library to be created
         */
        virtual void buildDynamic(const std::string& library,
                                  JobTimer* timer = NULL) = 0;

        /**
         * Deletes the previously compiled object files and clears of files
         * to include in a dynamic library
         */
        virtual void cleanup() = 0;

        inline virtual ~CLangCompiler() {
        }

    };

}

#endif
