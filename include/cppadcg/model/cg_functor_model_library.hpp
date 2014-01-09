#ifndef CPPAD_CG_FUNCTOR_MODEL_LIBRARY_INCLUDED
#define CPPAD_CG_FUNCTOR_MODEL_LIBRARY_INCLUDED
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
     * Abstract class used to load models
     * 
     * @author Joao Leal
     */
    template<class Base>
    class FunctorModelLibrary : public ModelLibrary<Base> {
    public:

        /**
         * Creates a new FunctorGenericModel object that can be used to evaluate
         * the model.
         * This object must be released by the user!
         * 
         * @param modelName The model name.
         * @return The model object (must be released by the user) or NULL if 
         *         no model exists with the provided name 
         */
        virtual FunctorGenericModel<Base>* model(const std::string& modelName) = 0;

        /**
         * Provides the API version used to create the model library.
         * 
         * @return the API version
         */
        virtual unsigned long getAPIVersion() = 0;

        /**
         * Provides a pointer to a function in the model library.
         * 
         * @param functionName The name of the function in the dynamic library
         * @param required Whether or not the function symbol must exist in the
         *                 library. If the function is required and does not
         *                 exist then the CppAD error handler is called, if it 
         *                 is not required and it does not exist then NULL is
         *                 return.
         * @return A pointer to the function symbol in the dynamic library if it
         *         exists, NULL otherwise.
         * @throws CGException If there is a problem loading the function symbol
         */
        virtual void* loadFunction(const std::string& functionName,
                                   bool required = true) throw (CGException) = 0;

        inline virtual ~FunctorModelLibrary() {
        }

    };

}

#endif


