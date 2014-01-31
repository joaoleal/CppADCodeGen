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
    public:

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

        inline virtual ~LlvmModelLibrary() {
            typename std::set<LlvmModel<Base>*>::const_iterator it;
            for (it = _models.begin(); it != _models.end(); ++it) {
                LlvmModel<Base>* model = *it;
                model->modelLibraryClosed();
            }
        }

    protected:

        virtual void destroyed(LlvmModel<Base>* model) {
            _models.erase(model);
        }

        friend class LlvmModel<Base>;
    };

}

#endif


