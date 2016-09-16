#ifndef CPPAD_CG_LLVM_MODEL_LIBRARY_INCLUDED
#define CPPAD_CG_LLVM_MODEL_LIBRARY_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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
namespace cg {

template<class Base>
class LlvmModel;

/**
 * Abstract class used to load JIT'ed models by LLVM
 * 
 * @author Joao Leal
 */
template<class Base>
class LlvmModelLibrary : public FunctorModelLibrary<Base> {
protected:
    unsigned long _version; // API version
    std::set<std::string> _modelNames;
    std::set<LlvmModel<Base>*> _models;
    void (*_onClose)();
    void (*_setThreadPoolDisabled)(int);
    void (*_setThreads)(unsigned int);
    unsigned int (*_getThreads)();
public:
    inline LlvmModelLibrary() :
            _onClose(nullptr),
            _setThreadPoolDisabled(nullptr),
            _setThreads(nullptr),
            _getThreads(nullptr) {
    }

    virtual std::set<std::string> getModelNames() override {
        return _modelNames;
    }

    virtual LlvmModel<Base>* model(const std::string& modelName) override {
        typename std::set<std::string>::const_iterator it = _modelNames.find(modelName);
        if (it == _modelNames.end()) {
            return nullptr;
        }
        LlvmModel<Base>* m = new LlvmModel<Base> (this, modelName);
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

    inline virtual ~LlvmModelLibrary() {
        // do not call clean-up here
        // cleanUp() must be called by the subclass (before destruction of the execution engine...)
    }

protected:

    inline void cleanUp() {
        for (LlvmModel<Base>* model : _models) {
            model->modelLibraryClosed();
        }

        if(_onClose != nullptr) {
            (*_onClose)();
            _onClose = nullptr;
        }
    }

    virtual void destroyed(LlvmModel<Base>* model) {
        _models.erase(model);
    }

    inline void validate() {
        /**
         * Check the version
         */
        unsigned long (*versionFunc)();
        versionFunc = reinterpret_cast<decltype(versionFunc)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_VERSION));

        this->_version = (*versionFunc)();
        if (ModelLibraryCSourceGen<Base>::API_VERSION != this->_version)
            throw CGException("The API version of the dynamic library (", this->_version,
                              ") is incompatible with the current version (",
                              ModelLibraryCSourceGen<Base>::API_VERSION, ")");

        /**
         * Load the list of models
         */
        void (*modelsFunc)(char const *const**, int*);
        modelsFunc = reinterpret_cast<decltype(modelsFunc)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_MODELS));

        char const*const* model_names = nullptr;
        int model_count;
        (*modelsFunc)(&model_names, &model_count);

        for (int i = 0; i < model_count; i++) {
            this->_modelNames.insert(model_names[i]);
        }

        /**
         * Load the the on close function
         */
        _onClose = reinterpret_cast<decltype(_onClose)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_ONCLOSE));

        /**
         * Thread pool related functions
         */
        _setThreadPoolDisabled = reinterpret_cast<decltype(_setThreadPoolDisabled)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADPOOLDISABLED, false));
        _setThreads = reinterpret_cast<decltype(_setThreads)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_SETTHREADS, false));
        _getThreads = reinterpret_cast<decltype(_getThreads)> (this->loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_GETTHREADS, false));

        if(_setThreads != nullptr) {
            (*_setThreads)(std::thread::hardware_concurrency());
        }
    }

    friend class LlvmModel<Base>;
};

} // END cg namespace
} // END CppAD namespace

#endif