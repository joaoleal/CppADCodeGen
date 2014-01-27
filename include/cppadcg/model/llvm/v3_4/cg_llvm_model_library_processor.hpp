#ifndef CPPAD_CG_LLVM_MODEL_LIBRARY_PROCESSOR_INCLUDED
#define CPPAD_CG_LLVM_MODEL_LIBRARY_PROCESSOR_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2014 Ciengis
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

#include <llvm/ADT/OwningPtr.h>
#include <llvm/Bitcode/ReaderWriter.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>
#include <llvm/Support/ManagedStatic.h>
#include <llvm/Support/MemoryBuffer.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/Support/system_error.h>
#include <llvm/Linker.h>

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
            ClangCompiler<Base> clang;
            return create(clang);
        }

        LlvmModelLibrary<Base>* create(ClangCompiler<Base>& clang) throw (CGException) {
            using namespace llvm;

            _linker.release();

            LlvmModelLibrary3_4<Base>* lib = NULL;

            this->modelLibraryHelper_->startingJob("", JobTimer::JIT_MODEL_LIBRARY);

            const std::map<std::string, ModelCSourceGen<Base>*>& models = this->modelLibraryHelper_->getModels();
            try {
                /**
                 * generate bit code
                 */
                typename std::map<std::string, ModelCSourceGen<Base>*>::const_iterator it;
                for (it = models.begin(); it != models.end(); ++it) {
                    const std::map<std::string, std::string>& modelSources = this->getSources(*it->second);

                    this->modelLibraryHelper_->startingJob("", JobTimer::COMPILING_FOR_MODEL);
                    clang.generateLLVMBitCode(modelSources, this->modelLibraryHelper_);
                    this->modelLibraryHelper_->finishedJob();
                }

                const std::map<std::string, std::string>& sources = this->getLibrarySources();
                clang.generateLLVMBitCode(sources, this->modelLibraryHelper_);

                const std::map<std::string, std::string>& customSource = this->modelLibraryHelper_->getCustomSources();
                clang.generateLLVMBitCode(customSource, this->modelLibraryHelper_);

                /**
                 * Load bit code and create a single module
                 */
                llvm::InitializeAllTargets();
                llvm::InitializeAllAsmPrinters();

                _context.reset(new llvm::LLVMContext());

                const std::set<std::string>& bcFiles = clang.getBitCodeFiles();
                typename std::set<std::string>::const_iterator itbc;
                for (itbc = bcFiles.begin(); itbc != bcFiles.end(); ++itbc) {
                    // load bitcode file
                    OwningPtr<MemoryBuffer> buffer;

                    std::string errMsg;
                    error_code ec = MemoryBuffer::getFile(*itbc, buffer);
                    if (buffer.get() == NULL)
                        throw CGException(ec.message());

                    // create the module
                    Module* module = llvm::ParseBitcodeFile(buffer.get(), *_context.get(), &errMsg);

                    // link modules together
                    if (_linker.get() == NULL) {
                        _linker.reset(new llvm::Linker(module)); // module not destroyed
                    } else {
                        if (_linker->linkInModule(module, &errMsg)) { // module destroyed
                            throw CGException(errMsg);
                        }
                    }
                }

                llvm::InitializeNativeTarget();

                // voila
                lib = new LlvmModelLibrary3_4<Base>(_linker->getModule(), _context.release());

            } catch (...) {
                clang.cleanup();
                throw;
            }
            clang.cleanup();

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
#if 0

        /**
         * Creates LLVM modules in a separate process.
         */
        static void createnPrintModule() {
            using namespace llvm;
            using namespace clang;

            /**
             * load source file
             */
            std::string source;
            std::cin >> source;

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
            compiler.createDiagnostics();
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
                hso.AddPath(llvm::StringRef(_includePaths[s]), clang::frontend::Angled, true, false);

            // Create and execute the frontend to generate an LLVM bitcode module.
            OwningPtr<CodeGenAction> action(new clang::EmitLLVMOnlyAction(_context.get()));

            if (!compiler.ExecuteAction(*action))
                throw CGException("Failed to emit LLVM bitcode");

            llvm::Module* module = action->takeModule();
            if (module == NULL)
                throw CGException("No module");

            /**
             * Print out the IR
             */
            //std::cout << *module;
            raw_fd_ostream os(STDOUT_FILENO, true);
            module->print(os);
            delete module;
        }
#endif
    };
}

#endif