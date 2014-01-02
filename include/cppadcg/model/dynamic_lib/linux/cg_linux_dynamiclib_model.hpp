#ifndef CPPAD_CG_LINUX_DYNAMICLIB_MODEL_INCLUDED
#define CPPAD_CG_LINUX_DYNAMICLIB_MODEL_INCLUDED
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
     * Useful class to call the compiled model in a dynamic library.
     * For the Linux Operating System only.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LinuxDynamicLibModel : public FunctorGenericModel<Base> {
    protected:
        /// the dynamic library
        LinuxDynamicLib<Base>* _dynLib;

    public:

        virtual ~LinuxDynamicLibModel() {
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
        LinuxDynamicLibModel(LinuxDynamicLib<Base>* dynLib, const std::string& name) :
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

        LinuxDynamicLibModel(const LinuxDynamicLibModel&); // not implemented

        LinuxDynamicLibModel& operator=(const LinuxDynamicLibModel&); // not implemented

        friend class LinuxDynamicLib<Base>;
    };

}

#endif

#endif