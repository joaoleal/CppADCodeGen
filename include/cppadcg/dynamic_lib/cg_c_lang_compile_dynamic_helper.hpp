#ifndef CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    /**
     * Useful class for generating source code for the creation of a dynamic
     * library.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangCompileDynamicHelper {
    public:
        static const std::string FUNCTION_VERSION;
        static const std::string FUNCTION_MODELS;
        static const unsigned long int API_VERSION;
    protected:
        static const std::string CONST;
    protected:
        std::map<std::string, CLangCompileModelHelper<Base>*> _models; // holds all models
        std::map<std::string, std::string> _customSource; // custom functions to be compiled in the dynamic library
        std::string _libraryName; // the path of the dynamic library to be created
        std::ostringstream _cache;
        bool _saveSourceFiles;
        //
        bool _verbose;
    public:

        CLangCompileDynamicHelper(CLangCompileModelHelper<Base>* model, bool saveSourceFiles = true) :
            _libraryName("cppad_cg_model"),
            _saveSourceFiles(saveSourceFiles),
            _verbose(false) {

            CPPADCG_ASSERT_KNOWN(model != NULL, "The model cannot be null");
            _models[model->getName()] = model;
        }

        inline const std::string& getLibraryName() const {
            return _libraryName;
        }

        inline void setLibraryName(const std::string& libraryName) {
            CPPADCG_ASSERT_KNOWN(!libraryName.empty(), "Library name cannot be empty");

            _libraryName = libraryName;
        }

        inline bool isVerbose() const {
            return _verbose;
        }

        inline void setVerbose(bool verbose) {
            _verbose = verbose;
        }

        inline void addModel(CLangCompileModelHelper<Base>* model) {
            CPPADCG_ASSERT_KNOWN(model != NULL, "The model cannot be null");
            CPPADCG_ASSERT_KNOWN(_models.find(model->getName()) == _models.end(),
                                 "Another model with the same name was already registered");

            _models[model->getName()] = model;
        }

        const std::map<std::string, CLangCompileModelHelper<Base>*>& getModels() const {
            return _models;
        }

        void addCustomFunctionSource(const std::string& filename, const std::string& source) {
            CPPADCG_ASSERT_KNOWN(!filename.empty(), "The filename name cannot be empty");

            _customSource[filename] = source;
        }

        DynamicLib<Base>* createDynamicLibrary(CLangCompiler<Base>& compiler);
        
        void createStaticLibrary(CLangCompiler<Base>& compiler, Archiver& ar, bool posIndepCode = false);

        inline virtual ~CLangCompileDynamicHelper() {
        };

    protected:

        virtual void generateVerionSource(std::map<std::string, std::string>& sources);

        virtual void generateModelsSource(std::map<std::string, std::string>& sources);

        virtual DynamicLib<Base>* loadDynamicLibrary();

    private:
        CLangCompileDynamicHelper(const CLangCompileDynamicHelper&); // not implemented

        CLangCompileDynamicHelper& operator=(const CLangCompileDynamicHelper&); // not implemented
    };

}

#endif
