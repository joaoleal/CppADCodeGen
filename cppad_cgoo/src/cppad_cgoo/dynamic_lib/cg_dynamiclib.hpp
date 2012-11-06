#ifndef CPPAD_CG_DYNAMICLIB_INCLUDED
#define	CPPAD_CG_DYNAMICLIB_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    /**
     * Abstract class used to load compiled models in a dynamic library
     * 
     * \author Joao Leal
     */
    template<class Base>
    class DynamicLib {
    public:
        /**
         * Provides the model names in the dynamic library.
         * 
         * \return the model names
         */
        virtual std::set<std::string> getModelNames() = 0;

        /**
         * Creates a new DynamicLibModel object that can be used to evaluate the
         * model. This object must be released by the user.
         * 
         * \param modelName The model name.
         * \return The model object or NULL if no model exists with the provided
         *         name
         */
        virtual DynamicLibModel<Base>* model(const std::string& modelName) = 0;

        /**
         * Provides the API version used to create the dynamic library.
         * 
         * \return the API version
         */
        virtual unsigned long int getAPIVersion() = 0;

        /**
         * Provides a pointer to a function in the dynamic library.
         * 
         * \param functionName The name of the function in the dynamic library
         * \param required Whether or not the function symbol must exist in the
         *                 library. If the function is required and does not
         *                 exist then the CppAD error handler is called, if it 
         *                 is not required and it does not exist then NULL is
         *                 return.
         * \return A pointer to the function symbol in the dynamic library if it
         *         exists, NULL otherwise.
         */
        virtual void* loadFunction(const std::string& functionName, bool required = true) = 0;

        /**
         * Provides a pointer to a function in the dynamic library.
         * 
         * \param functionName The name of the function in the dynamic library
         * \param error If there is a problem loading the function symbol this
         *              string will contain a non-empty error message.
         * \return A pointer to the function symbol in the dynamic library if it
         *         exists, NULL otherwise.
         */
        virtual void* loadFunction(const std::string& functionName, std::string& error) = 0;

        inline virtual ~DynamicLib() {
        }

    };

}

#endif


