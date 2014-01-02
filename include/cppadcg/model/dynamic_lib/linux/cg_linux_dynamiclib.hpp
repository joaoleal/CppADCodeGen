#ifndef CPPAD_CG_LINUX_DYNAMICLIB_INCLUDED
#define CPPAD_CG_LINUX_DYNAMICLIB_INCLUDED
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
#ifdef __linux__

#include <typeinfo>
#include <dlfcn.h>

namespace CppAD {

    /**
     * Useful class to call the compiled source code in a dynamic library.
     * For the Linux Operating System only.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LinuxDynamicLib : public DynamicLib<Base> {
    protected:
        const std::string _dynLibName;
        /// the dynamic library handler
        void* _dynLibHandle;
        unsigned long _version; // API version
        std::set<std::string> _modelNames;
        std::set<LinuxDynamicLibModel<Base>*> _models;
    public:

        LinuxDynamicLib(const std::string& dynLibName) :
            _dynLibName(dynLibName),
            _dynLibHandle(NULL) {

            std::string path;
            if (dynLibName[0] == '/') {
                path = dynLibName; // absolute path
            } else if (!(dynLibName[0] == '.' && dynLibName[1] == '/') &&
                    !(dynLibName[0] == '.' && dynLibName[1] == '.')) {
                path = "./" + dynLibName; // relative path
            } else {
                path = dynLibName;
            }

            // load the dynamic library
            _dynLibHandle = dlopen(path.c_str(), RTLD_NOW);
            CPPADCG_ASSERT_KNOWN(_dynLibHandle != NULL, ("Failed to dynamically load library '" + dynLibName + "': " + dlerror()).c_str());

            // validate the dynamic library
            validate();
        }

        virtual std::set<std::string> getModelNames() {
            return _modelNames;
        }

        virtual GenericModel<Base>* model(const std::string& modelName) {
            typename std::set<std::string>::const_iterator it = _modelNames.find(modelName);
            if (it == _modelNames.end()) {
                return NULL;
            }
            LinuxDynamicLibModel<Base>* m = new LinuxDynamicLibModel<Base> (this, modelName);
            _models.insert(m);
            return m;
        }

        virtual unsigned long getAPIVersion() {
            return _version;
        }

        virtual void* loadFunction(const std::string& functionName, bool required = true) throw (CGException) {
            void* functor = dlsym(_dynLibHandle, functionName.c_str());

            if (required) {
                char *err = dlerror();
                if (err != NULL)
                    throw CGException("Failed to load function '" + functionName + "': " + err);
            }

            return functor;
        }

        virtual ~LinuxDynamicLib() {
            typename std::set<LinuxDynamicLibModel<Base>*>::const_iterator it;
            for (it = _models.begin(); it != _models.end(); ++it) {
                LinuxDynamicLibModel<Base>* model = *it;
                model->modelLibraryClosed();
            }

            if (_dynLibHandle != NULL) {
                dlclose(_dynLibHandle);
                _dynLibHandle = NULL;
            }
        }

    protected:

        inline void validate() throw (CGException) {
            /**
             * Check the version
             */
            unsigned long (*versionFunc)();
            *(void **) (&versionFunc) = loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_VERSION);

            _version = (*versionFunc)();
            if (ModelLibraryCSourceGen<Base>::API_VERSION != _version)
                throw CGException("The API version of the dynamic library is incompatible with the current version");

            /**
             * Load the list of models
             */
            void (*modelsFunc)(char const *const**, int*);
            *(void **) (&modelsFunc) = loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_MODELS);

            char const*const* model_names = NULL;
            int model_count;
            (*modelsFunc)(&model_names, &model_count);

            for (int i = 0; i < model_count; i++) {
                _modelNames.insert(model_names[i]);
            }
        }

        virtual void destroyed(LinuxDynamicLibModel<Base>* model) {
            _models.erase(model);
        }

    private:
        LinuxDynamicLib(const LinuxDynamicLib&); // not implemented

        LinuxDynamicLib& operator=(const LinuxDynamicLib&); // not implemented

        friend class LinuxDynamicLibModel<Base>;

    };

}

#endif

#endif