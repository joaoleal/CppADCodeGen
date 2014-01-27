#ifndef CPPAD_CG_LLVM_MODEL_INCLUDED
#define CPPAD_CG_LLVM_MODEL_INCLUDED
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
     * Useful class to call the JIT'ed models with LLVM.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LlvmModel : public FunctorGenericModel<Base> {
    protected:
        /// the dynamic library
        LlvmModelLibrary<Base>* _dynLib;

    public:

        virtual ~LlvmModel() {
            if (_dynLib != NULL) {
                _dynLib->destroyed(this);
            }
        }

    protected:

        /**
         * Creates a new model 
         * 
         * @param name The model name
         */
        LlvmModel(LlvmModelLibrary<Base>* dynLib,
                  const std::string& name) :
            FunctorGenericModel<Base>(name),
            _dynLib(dynLib) {

            CPPADCG_ASSERT_UNKNOWN(_dynLib != NULL);

            this->init();
        }

        virtual void* loadFunction(const std::string& functionName, bool required = true) throw (CGException) {
            return _dynLib->loadFunction(functionName, required);
        }

        virtual void modelLibraryClosed() {
            _dynLib = NULL;
            FunctorGenericModel<Base>::modelLibraryClosed();
        }

    private:

        LlvmModel(const LlvmModel&); // not implemented

        LlvmModel& operator=(const LlvmModel&); // not implemented

        friend class LlvmModelLibrary<Base>;
    };

}

#endif