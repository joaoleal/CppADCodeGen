#ifndef CPPAD_CG_LLVM_MODEL_LIBRARY_3_6_INCLUDED
#define CPPAD_CG_LLVM_MODEL_LIBRARY_3_6_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2015 Ciengis
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

template<class Base> class LlvmModel;

/**
 * Class used to load JIT'ed models by LLVM 3.2
 * 
 * @author Joao Leal
 */
template<class Base>
class LlvmModelLibrary3_6 : public LlvmModelLibrary<Base> {
protected:
    llvm::Module* _module; // owned by _executionEngine
    std::unique_ptr<llvm::LLVMContext> _context;
    std::unique_ptr<llvm::ExecutionEngine> _executionEngine;
    std::unique_ptr<llvm::FunctionPassManager> _fpm;
public:

    LlvmModelLibrary3_6(llvm::Module* module,
                        llvm::LLVMContext* context) :
        _module(module),
        _context(context) {
        using namespace llvm;

        // Create the JIT.  This takes ownership of the module.
        std::unique_ptr<llvm::Module> m(_module);
        std::string errStr;
        _executionEngine.reset(EngineBuilder(std::move(m))
                               .setErrorStr(&errStr)
                               .setEngineKind(EngineKind::JIT)
#ifndef NDEBUG
                .setVerifyModules(true)
#endif
                // .setMCJITMemoryManager(llvm::make_unique<llvm::SectionMemoryManager>())
                               .create());
        if (!_executionEngine.get()) {
            throw CGException("Could not create ExecutionEngine: " + errStr);
        }

        _fpm.reset(new llvm::FunctionPassManager(_module));

        preparePassManager();

        _fpm->doInitialization();

        /**
         * 
         */
        validate();
    }

    LlvmModelLibrary3_6(const LlvmModelLibrary3_6&) = delete;
    LlvmModelLibrary3_6& operator=(const LlvmModelLibrary3_6&) = delete;

    /**
     * Set up the optimizer pipeline
     */
    virtual void preparePassManager() {
        llvm::PassManagerBuilder builder;
        builder.OptLevel = 2;
        builder.populateFunctionPassManager(*_fpm);
        //_fpm.add(new DataLayoutPass());
    }

    virtual void* loadFunction(const std::string& functionName, bool required = true) throw (CGException) override {
        llvm::Function* func = _module->getFunction(functionName);
        if (func == nullptr) {
            if (required)
                throw CGException("Unable to find function '" + functionName + "' in LLVM module");
            return nullptr;
        }

#ifndef NDEBUG
        // Validate the generated code, checking for consistency.
        llvm::raw_os_ostream os(std::cerr);
        bool failed = llvm::verifyFunction(*func, &os);
        if (failed)
            throw CGException("Function '" + functionName + "' verification failed");
#endif

        // Optimize the function.
        _fpm->run(*func);

        // JIT the function, returning a function pointer.
        uint64_t fPtr = _executionEngine->getFunctionAddress(functionName);
        if (fPtr == 0 && required) {
            throw CGException("Unable to find function '" + functionName + "' in LLVM module");
        }
        return (void*) fPtr;
    }

    inline virtual ~LlvmModelLibrary3_6() {
    }

protected:

    inline void validate() throw (CGException) {
        /**
         * Check the version
         */
        unsigned long (*versionFunc)();
        versionFunc = reinterpret_cast<decltype(versionFunc)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_VERSION));

        this->_version = (*versionFunc)();
        if (ModelLibraryCSourceGen<Base>::API_VERSION != this->_version)
            throw CGException("The API version of the dynamic library is incompatible with the current version");

        /**
         * Load the list of models
         */
        void (*modelsFunc)(char const *const**, int*);
        modelsFunc = reinterpret_cast<decltype(modelsFunc)> (loadFunction(ModelLibraryCSourceGen<Base>::FUNCTION_MODELS));

        char const*const* model_names = nullptr;
        int model_count;
        (*modelsFunc)(&model_names, &model_count);

        for (int i = 0; i < model_count; i++) {
            this->_modelNames.insert(model_names[i]);
        }
    }

    friend class LlvmModel<Base>;

};

} // END cg namespace
} // END CppAD namespace

#endif