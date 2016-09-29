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
#if CPPAD_CG_SYSTEM_LINUX

#include <typeinfo>
#include <dlfcn.h>

namespace CppAD {
namespace cg {

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
    void (*_onClose)();
    void (*_setThreadPoolDisabled)(int);
    void (*_setThreads)(unsigned int);
    unsigned int (*_getThreads)();
    void (*_setSchedulerStrategy)(int);
    int (*_getSchedulerStrategy)();
    void (*_setThreadPoolVerbose)(int v);
    int (*_isThreadPoolVerbose)();
    void (*_setThreadPoolMultiJobMaxWork)(float v);
    float (*_getThreadPoolMultiJobMaxWork)();
    std::set<std::string> _modelNames;
    std::set<LinuxDynamicLibModel<Base>*> _models;
public:

    LinuxDynamicLib(const std::string& dynLibName,
                    int dlOpenMode = RTLD_NOW) :
        _dynLibName(dynLibName),
        _dynLibHandle(nullptr),
        _onClose(nullptr),
        _setThreadPoolDisabled(nullptr),
        _setThreads(nullptr),
        _getThreads(nullptr),
        _setSchedulerStrategy(nullptr),
        _getSchedulerStrategy(nullptr),
        _setThreadPoolVerbose(nullptr),
        _isThreadPoolVerbose(nullptr),
        _setThreadPoolMultiJobMaxWork(nullptr),
        _getThreadPoolMultiJobMaxWork(nullptr) {

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
        _dynLibHandle = dlopen(path.c_str(), dlOpenMode);
        CPPADCG_ASSERT_KNOWN(_dynLibHandle != nullptr, ("Failed to dynamically load library '" + dynLibName + "': " + dlerror()).c_str());

        // validate the dynamic library
        validate();
    }

    LinuxDynamicLib(const LinuxDynamicLib&) = delete;
    LinuxDynamicLib& operator=(const LinuxDynamicLib&) = delete;

    virtual std::set<std::string> getModelNames() override {
        return _modelNames;
    }

    virtual LinuxDynamicLibModel<Base>* model(const std::string& modelName) override {
        std::set<std::string>::const_iterator it = _modelNames.find(modelName);
        if (it == _modelNames.end()) {
            return nullptr;
        }
        LinuxDynamicLibModel<Base>* m = new LinuxDynamicLibModel<Base> (this, modelName);
        _models.insert(m);
        return m;
    }

    virtual unsigned long getAPIVersion() override {
        return _version;
    }

    virtual void setThreadPoolDisabled(bool disabled) override {
        if(_setThreadPoolDisabled != nullptr) {
            (*_setThreadPoolDisabled)(disabled);
        }
    }

    virtual unsigned int getThreadNumber() const override {
        if (_getThreads != nullptr) {
            return (*_getThreads)();
        }
        return 1;
    }

    virtual void setThreadNumber(unsigned int n) override {
        if (_setThreads != nullptr) {
            (*_setThreads)(n);
        }
    }

    virtual ThreadPoolScheduleStrategy getThreadPoolSchedulerStrategy() const override {
        if (_getSchedulerStrategy != nullptr) {
            return ThreadPoolScheduleStrategy((*_getSchedulerStrategy)());
        }
        return ThreadPoolScheduleStrategy::SINGLE_JOB;
    }

    virtual void setThreadPoolSchedulerStrategy(ThreadPoolScheduleStrategy s) override {
        if (_setSchedulerStrategy != nullptr) {
            (*_setSchedulerStrategy)(int(s));
        }
    }

    virtual void setThreadPoolVerbose(bool v) override {
        if (_setThreadPoolVerbose != nullptr) {
            (*_setThreadPoolVerbose)(int(v));
        }
    }

    virtual bool isThreadPoolVerbose() override {
        if (_isThreadPoolVerbose != nullptr) {
            return bool((*_isThreadPoolVerbose)());
        }
        return false;
    }

    virtual void setThreadPoolMultiJobMaxWork(float v) override {
        if (_setThreadPoolMultiJobMaxWork != nullptr) {
            (*_setThreadPoolMultiJobMaxWork)(v);
        }
    }

    virtual float getThreadPoolMultiJobMaxWork() override {
        if (_getThreadPoolMultiJobMaxWork != nullptr) {
            return (*_getThreadPoolMultiJobMaxWork)();
        }
        return 1.0;
    }

    virtual void* loadFunction(const std::string& functionName, bool required = true) override {
        void* functor = dlsym(_dynLibHandle, functionName.c_str());

        if (required) {
            char *err = dlerror();
            if (err != nullptr)
                throw CGException("Failed to load function '", functionName, "': ", err);
        }

        return functor;
    }

    virtual ~LinuxDynamicLib() {
        for (LinuxDynamicLibModel<Base>* model : _models) {
            model->modelLibraryClosed();
        }

        if (_dynLibHandle != nullptr) {
            if(_onClose != nullptr) {
                (*_onClose)();
            }

            dlclose(_dynLibHandle);
            _dynLibHandle = nullptr;
        }
    }

protected:

    inline void validate() {
        /**
         * Check the version
         */
        unsigned long (*versionFunc)();
        versionFunc = reinterpret_cast<decltype(versionFunc)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_VERSION));

        _version = (*versionFunc)();
        if (ModelLibraryCSourceGen<Base>::API_VERSION != _version)
            throw CGException("The API version of the dynamic library (", _version,
                              ") is incompatible with the current version (",
                              ModelLibraryCSourceGen<Base>::API_VERSION, ")");

        /**
         * Load the list of models
         */
        void (*modelsFunc)(char const *const**, int*);
        modelsFunc = reinterpret_cast<decltype(modelsFunc)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_MODELS));

        char const*const* model_names = nullptr;
        int model_count;
        (*modelsFunc)(&model_names, &model_count);

        for (int i = 0; i < model_count; i++) {
            _modelNames.insert(model_names[i]);
        }

        /**
         * Load the the on close function
         */
        _onClose = reinterpret_cast<decltype(_onClose)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_ONCLOSE, false));

        /**
         * Thread pool related functions
         */
        _setThreadPoolDisabled = reinterpret_cast<decltype(_setThreadPoolDisabled)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADPOOLDISABLED, false));
        _setThreads = reinterpret_cast<decltype(_setThreads)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADS, false));
        _getThreads = reinterpret_cast<decltype(_getThreads)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_GETTHREADS, false));
        _setSchedulerStrategy = reinterpret_cast<decltype(_setSchedulerStrategy)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADSCHEDULERSTRAT, false));
        _getSchedulerStrategy = reinterpret_cast<decltype(_getSchedulerStrategy)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_GETTHREADSCHEDULERSTRAT, false));
        _setThreadPoolVerbose = reinterpret_cast<decltype(_setThreadPoolVerbose)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADPOOLVERBOSE, false));
        _isThreadPoolVerbose = reinterpret_cast<decltype(_isThreadPoolVerbose)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_ISTHREADPOOLVERBOSE, false));
        _setThreadPoolMultiJobMaxWork = reinterpret_cast<decltype(_setThreadPoolMultiJobMaxWork)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADPOOLMULTIJOBMAXGROUPWORK, false));
        _getThreadPoolMultiJobMaxWork = reinterpret_cast<decltype(_getThreadPoolMultiJobMaxWork)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_GETTHREADPOOLMULTIJOBMAXGROUPWORK, false));

        if(_setThreads != nullptr) {
            (*_setThreads)(std::thread::hardware_concurrency());
        }
    }

    virtual void destroyed(LinuxDynamicLibModel<Base>* model) {
        _models.erase(model);
    }

    friend class LinuxDynamicLibModel<Base>;

};

} // END cg namespace
} // END CppAD namespace

#endif
#endif