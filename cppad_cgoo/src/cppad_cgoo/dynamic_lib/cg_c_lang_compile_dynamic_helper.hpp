#ifndef CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_INCLUDED

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
        bool _savedSourceFiles;
        //
        bool _verbose;
    public:

        CLangCompileDynamicHelper(CLangCompileModelHelper<Base>* model, bool savedSourceFiles = true) :
            _libraryName("cppad_cg_model.so"),
            _savedSourceFiles(savedSourceFiles),
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
