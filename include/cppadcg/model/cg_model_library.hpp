#ifndef CPPAD_CG_MODEL_LIBRARY_INCLUDED
#define CPPAD_CG_MODEL_LIBRARY_INCLUDED
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
    class ModelLibrary {
    public:
        /**
         * Provides the model names in the dynamic library.
         * 
         * @return the model names
         */
        virtual std::set<std::string> getModelNames() = 0;

        /**
         * Creates a new GenericModel object that can be used to evaluate the
         * model.
         * This object must be released by the user!
         * 
         * @param modelName The model name.
         * @return The model object (must be released by the user) or nullptr if 
         *         no model exists with the provided name 
         */
        virtual GenericModel<Base>* model(const std::string& modelName) = 0;

        inline virtual ~ModelLibrary() {
        }

    };

}

#endif


