#ifndef CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_IMPL_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_IMPL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <typeinfo>

#include "cg_c_lang_compile_model_helper.hpp"

namespace CppAD {

    template<class Base>
    const unsigned long int CLangCompileDynamicHelper<Base>::API_VERSION = 2;

    template<class Base>
    const std::string CLangCompileDynamicHelper<Base>::FUNCTION_VERSION = "cppad_cg_version";

    template<class Base>
    const std::string CLangCompileDynamicHelper<Base>::FUNCTION_MODELS = "cppad_cg_models";

    template<class Base>
    const std::string CLangCompileDynamicHelper<Base>::CONST = "const";

    template<class Base>
    DynamicLib<Base>* CLangCompileDynamicHelper<Base>::createDynamicLibrary(CLangCompiler<Base>& compiler) {
        try {
            typename std::map<std::string, CLangCompileModelHelper<Base>*>::const_iterator it;
            for (it = _models.begin(); it != _models.end(); ++it) {
                it->second->setVerbose(_verbose);
                it->second->compileSources(compiler);
            }
            
            std::map<std::string, std::string> sources;
            generateVerionSource(sources);
            generateModelsSource(sources);

            sources.insert(_customSource.begin(), _customSource.end());

            compiler.compileSources(sources, _savedSourceFiles);

            compiler.buildDynamic(_libraryName);
        } catch (...) {
            compiler.cleanup();
            throw;
        }
        compiler.cleanup();

        return loadDynamicLibrary();
    }

    template<class Base>
    void CLangCompileDynamicHelper<Base>::generateVerionSource(std::map<std::string, std::string>& sources) {
        _cache.str("");
        _cache << "unsigned long int " << FUNCTION_VERSION << "() {\n"
                << "return " << API_VERSION << "u;\n"
                << "}\n\n";

        sources[FUNCTION_VERSION + ".c"] = _cache.str();
    }

    template<class Base>
    void CLangCompileDynamicHelper<Base>::generateModelsSource(std::map<std::string, std::string>& sources) {
        _cache.str("");
        _cache << "void " << FUNCTION_MODELS << "(char const *const** names, int* count) {\n"
                "const static char* models[] = {\n";

        typename std::map<std::string, CLangCompileModelHelper<Base>*>::const_iterator it;
        for (it = _models.begin(); it != _models.end(); ++it) {
            if (it != _models.begin()) {
                _cache << ",\n";
            }
            _cache << "\"" << it->first << "\"";
        }
        _cache << "};\n"
                "*names = models;\n"
                "*count = " << _models.size() << ";\n"
                "}\n\n";

        sources[FUNCTION_MODELS + ".c"] = _cache.str();
    }

}

#endif
