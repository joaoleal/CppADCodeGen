#ifndef CPPAD_CG_LINUX_DYNAMICLIB_INCLUDED
#define	CPPAD_CG_LINUX_DYNAMICLIB_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#ifdef __linux__

#include <typeinfo>
#include <dlfcn.h>

namespace CppAD {

    /**
     * Useful class to call the compiled source code in a dynamic library.
     * For the Linux Operating System only.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class LinuxDynamicLib : public DynamicLib<Base> {
    protected:
        const std::string _dynLibName;
        /// the dynamic library handler
        void* _dynLibHandle;
        unsigned long int _version; // API version
        std::map<std::string, DynamicLibModel<Base>*> _models;
    public:

        LinuxDynamicLib(const std::string& dynLibName) :
            _dynLibName(dynLibName),
            _dynLibHandle(NULL) {

            std::string path;
            if (dynLibName[0] == '/') {
                path = dynLibName; // absolute path
            } else {
                path = "./" + dynLibName; // relative path
            }

            // load the dynamic library
            _dynLibHandle = dlopen(path.c_str(), RTLD_NOW);
            CPPADCG_ASSERT_KNOWN(_dynLibHandle != NULL, ("Failed to dynamically load library '" + dynLibName + "'").c_str());

            // validate the dynamic library
            validate();
        }

        virtual std::map<std::string, DynamicLibModel<Base>*> getModels() {
            return _models;
        }

        virtual DynamicLibModel<Base>* model(const std::string& modelName) {
            typename std::map<std::string, DynamicLibModel<Base>*>::iterator it = _models.find(modelName);
            if (it != _models.end())
                return it->second;
            else
                return NULL;
        }

        virtual size_t getAPIVersion() {
            return _version;
        }
        
        virtual void* loadFunction(const std::string& functionName, std::string& error) {
            void* functor = dlsym(_dynLibHandle, functionName.c_str());
            char *err;
            error.clear();
            if ((err = dlerror()) != NULL) {
                error += "Failed to load function '";
                error += functionName + ": " + err;
                return NULL;
            }
            return functor;
        }

        virtual ~LinuxDynamicLib() {
            typename std::map<std::string, DynamicLibModel<Base>*>::const_iterator it;
            for (it != _models.begin(); it != _models.end(); ++it) {
                DynamicLibModel<Base>* model = it->second;
                delete model;
            }

            if (_dynLibHandle != NULL) {
                dlclose(_dynLibHandle);
                _dynLibHandle = NULL;
            }
        }

    protected:

        inline void validate() {
            std::string error;

            /**
             * Check the version
             */
            unsigned long int (*versionFunc)();
            *(void **) (&versionFunc) = loadFunction(CLangCompileDynamicHelper<Base>::FUNCTION_VERSION, error);
            CPPADCG_ASSERT_KNOWN(error.empty(), error.c_str());

            _version = (*versionFunc)();
            CPPADCG_ASSERT_KNOWN(CLangCompileDynamicHelper<Base>::API_VERSION == _version,
                                 "The API version of the dynamic library is incompatible with the current version");

            /**
             * Load the models
             */
            void (*modelsFunc)(char const *const**, int*);
            *(void **) (&modelsFunc) = loadFunction(CLangCompileDynamicHelper<Base>::FUNCTION_MODELS, error);
            CPPADCG_ASSERT_KNOWN(error.empty(), error.c_str());

            char const*const* model_names = NULL;
            int model_count;
            (*modelsFunc)(&model_names, &model_count);

            for (int i = 0; i < model_count; i++) {
                _models[model_names[i]] = new LinuxDynamicLibModel<Base>(this, model_names[i]);
            }
        }

    private:
        LinuxDynamicLib(const LinuxDynamicLib&); // not implemented

        LinuxDynamicLib& operator=(const LinuxDynamicLib&); // not implemented

        friend class LinuxDynamicLibModel<Base>;

    };

}

#endif

#endif