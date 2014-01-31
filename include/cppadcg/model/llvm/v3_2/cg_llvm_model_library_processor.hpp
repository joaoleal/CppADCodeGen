#ifndef CPPAD_CG_LLVM_MODEL_LIBRARY_PROCESSOR_INCLUDED
#define CPPAD_CG_LLVM_MODEL_LIBRARY_PROCESSOR_INCLUDED
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

    /**
     * Useful class for generating a JIT evaluated model library.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LlvmModelLibraryProcessor : public ModelLibraryProcessor<Base> {
    protected:
        std::vector<std::string> _includePaths;
        std::auto_ptr<llvm::Linker> _linker;
        std::auto_ptr<llvm::LLVMContext> _context;
    public:

        /**
         * 
         * @param modelLibraryHelper
         */
        LlvmModelLibraryProcessor(ModelLibraryCSourceGen<Base>& modelLibraryHelper) :
            ModelLibraryProcessor<Base>(modelLibraryHelper) {
        }

        inline void setIncludePaths(const std::vector<std::string>& includePaths) {
            _includePaths = includePaths;
        }

        inline const std::vector<std::string>& getIncludePaths() const {
            return _includePaths;
        }

        LlvmModelLibrary<Base>* create() throw (CGException) {
            _linker.release();

            this->modelLibraryHelper_->startingJob("", JobTimer::JIT_MODEL_LIBRARY);

            llvm::InitializeAllTargets();
            llvm::InitializeAllAsmPrinters();

            _context.reset(new llvm::LLVMContext());

            const std::map<std::string, ModelCSourceGen<Base>*>& models = this->modelLibraryHelper_->getModels();
            typename std::map<std::string, ModelCSourceGen<Base>*>::const_iterator it;
            for (it = models.begin(); it != models.end(); ++it) {
                const std::map<std::string, std::string>& modelSources = this->getSources(*it->second);
                createLlvmModules(modelSources);
            }

            const std::map<std::string, std::string>& sources = this->getLibrarySources();
            createLlvmModules(sources);

            const std::map<std::string, std::string>& customSource = this->modelLibraryHelper_->getCustomSources();
            createLlvmModules(customSource);

            llvm::InitializeNativeTarget();

            LlvmModelLibrary3_2<Base>* lib = new LlvmModelLibrary3_2<Base>(_linker->releaseModule(), _context.release());

            this->modelLibraryHelper_->finishedJob();

            return lib;
        }

        static inline LlvmModelLibrary<Base>* create(ModelLibraryCSourceGen<Base>& modelLibraryHelper) throw (CGException) {
            LlvmModelLibraryProcessor<Base> p(modelLibraryHelper);
            return p.create();
        }

        virtual ~LlvmModelLibraryProcessor() {
        }

    protected:

        virtual void createLlvmModules(const std::map<std::string, std::string>& sources) throw (CGException) {
            std::map<std::string, std::string>::const_iterator its;
            for (its = sources.begin(); its != sources.end(); ++its) {
                createLlvmModule(its->first, its->second);
            }
        }

        virtual void createLlvmModule(const std::string& filename,
                                      const std::string& source) throw (CGException) {
            using namespace llvm;
            using namespace clang;

            static const char* argv [] = {"program", "-x", "c", "string-input"};
            static const int argc = sizeof (argv) / sizeof (argv[0]);

            IntrusiveRefCntPtr<DiagnosticOptions> diagOpts = new DiagnosticOptions();
            TextDiagnosticPrinter *diagClient = new TextDiagnosticPrinter(llvm::errs(), &*diagOpts); // will be owned by diags
            IntrusiveRefCntPtr<DiagnosticIDs> diagID(new DiagnosticIDs());
            IntrusiveRefCntPtr<DiagnosticsEngine> diags(new DiagnosticsEngine(diagID, &*diagOpts, diagClient));

            ArrayRef<const char *> args(argv + 1, // skip program name
                                        argc - 1);
            std::auto_ptr<CompilerInvocation> invocation(createInvocationFromCommandLine(args, diags));
            if (invocation.get() == NULL)
                throw CGException("Failed to create compiler invocation");
            CompilerInvocation::setLangDefaults(*invocation->getLangOpts(), IK_C,
                                                LangStandard::lang_unspecified);
            invocation->getFrontendOpts().DisableFree = false; // make sure we free memory (by default it does not)

            // Create a compiler instance to handle the actual work.
            CompilerInstance compiler;
            compiler.setInvocation(invocation.release());

            // Create the compilers actual diagnostics engine.
            compiler.createDiagnostics(argc, const_cast<char**> (argv));
            if (!compiler.hasDiagnostics())
                throw CGException("No diagnostics");

            // Create memory buffer with source text
            llvm::MemoryBuffer * buffer = llvm::MemoryBuffer::getMemBufferCopy(source, "SIMPLE_BUFFER");
            if (buffer == NULL)
                throw CGException("Failed to create memory buffer");

            // Remap auxiliary name "string-input" to memory buffer
            PreprocessorOptions& po = compiler.getInvocation().getPreprocessorOpts();
            po.addRemappedFile("string-input", buffer);

            HeaderSearchOptions& hso = compiler.getInvocation().getHeaderSearchOpts();
            for (size_t s = 0; s < _includePaths.size(); s++)
                hso.AddPath(llvm::StringRef(_includePaths[s]), clang::frontend::Angled, true, false, false);

            // Create and execute the frontend to generate an LLVM bitcode module.
            OwningPtr<CodeGenAction> action(new clang::EmitLLVMOnlyAction(_context.get()));
            if (!compiler.ExecuteAction(*action))
                throw CGException("Failed to emit LLVM bitcode");

            llvm::Module* module = action->takeModule();
            if (module == NULL)
                throw CGException("No module");

            if (_linker.get() == NULL) {
                _linker.reset(new llvm::Linker(std::string("MyLinker"), module));
            } else {
            std::string errorMsg;
            if (_linker->LinkInModule(module, &errorMsg)) {
                throw CGException(errorMsg);
            }
            }

            // NO delete module;
            // NO delete invocation;
            //llvm::llvm_shutdown();
        }

        inline llvm::Module* mergeModules(const std::vector<llvm::Module*>& modules) throw (CGException) {
            if (modules.empty())
                return NULL;

            std::string progName("MyLinker");
            std::auto_ptr<llvm::Linker> ld(new llvm::Linker(progName, modules[0]));

            for (size_t m = 1; m < modules.size(); m++) {
                std::string errorMsg;
                if (ld->LinkInModule(modules[m], &errorMsg)) {
                    throw CGException(errorMsg);
                }
            }

            return ld->releaseModule();
        }

    };
}

#endif