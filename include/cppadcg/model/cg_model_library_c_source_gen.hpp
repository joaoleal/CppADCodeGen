#ifndef CPPAD_CG_MODEL_LIBRARY_C_SOURCE_GEN_INCLUDED
#define CPPAD_CG_MODEL_LIBRARY_C_SOURCE_GEN_INCLUDED
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
     * Generates C source code for a bundle of models.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class ModelLibraryCSourceGen : public JobTimer {
    public:
        static const std::string FUNCTION_VERSION;
        static const std::string FUNCTION_MODELS;
        static const unsigned long API_VERSION;
    protected:
        static const std::string CONST;
    protected:
        std::map<std::string, ModelCSourceGen<Base>*> _models; // holds all models
        std::map<std::string, std::string> _customSource; // custom functions to be compiled in the dynamic library
        std::ostringstream _cache;
        /**
         * Library level generated source files
         */
        std::map<std::string, std::string> _libSources;
        /**
         * 
         */
        std::set<JobListener*> _listeners;
    public:

        /**
         * Creates a new helper class for the generation of dynamic libraries
         * using the C language.
         * 
         * @param model A model compilation helper (must only be deleted after
         *              this object)
         */
        ModelLibraryCSourceGen(ModelCSourceGen<Base>& model) {
            _models[model.getName()] = &model;
        }

        /**
         * Adds additional models to be compiled into the created library.
         * 
         * @param model a model compilation helper (must only be deleted after
         *              this object)
         */
        inline void addModel(ModelCSourceGen<Base>& model) {
            CPPADCG_ASSERT_KNOWN(_models.find(model.getName()) == _models.end(),
                                 "Another model with the same name was already registered");

            _models[model.getName()] = &model;
        }

        const std::map<std::string, ModelCSourceGen<Base>*>& getModels() const {
            return _models;
        }

        void addCustomFunctionSource(const std::string& filename, const std::string& source) {
            CPPADCG_ASSERT_KNOWN(!filename.empty(), "The filename name cannot be empty");

            _customSource[filename] = source;
        }

        /**
         * Provides the user defined custom sources. 
         * 
         * @return maps filenames to the file content for the user defined
         *         sources.
         */
        const std::map<std::string, std::string>& getCustomSources() const {
            return _customSource;
        }

        /**
         * Saves the generated C source code into several files.
         * 
         * @param sourcesFolder A directory path where the files should be
         *                      created (any existing files with the same names
         *                      will be overridden).
         */
        void saveSources(const std::string& sourcesFolder);

        inline virtual ~ModelLibraryCSourceGen() {
        };

    protected:
        virtual const std::map<std::string, std::string>& getLibrarySources();

        virtual void generateVerionSource(std::map<std::string, std::string>& sources);

        virtual void generateModelsSource(std::map<std::string, std::string>& sources);

        static void saveSources(const std::string& sourcesFolder,
                                const std::map<std::string, std::string>& sources);

    private:
        ModelLibraryCSourceGen(const ModelLibraryCSourceGen&); // not implemented

        ModelLibraryCSourceGen& operator=(const ModelLibraryCSourceGen&); // not implemented

        friend class ModelLibraryProcessor<Base>;
    };

}

#endif
