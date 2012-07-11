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
        virtual std::map<std::string, DynamicLibModel<Base>*> getModels() = 0;

        virtual DynamicLibModel<Base>* model(const std::string& modelName) = 0;

        virtual size_t getAPIVersion() = 0;

        virtual void* loadFunction(const std::string& functionName, std::string& error) = 0;

    };

}

#endif


