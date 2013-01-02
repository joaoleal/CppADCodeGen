#ifndef CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_IMPL_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_DYNAMIC_HELPER_IMPL_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

#include <typeinfo>

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
        
        compiler.setVerbose(_verbose);
        
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

            compiler.compileSources(sources, _saveSourceFiles);

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
                << "   return " << API_VERSION << "u;\n"
                << "}\n\n";

        sources[FUNCTION_VERSION + ".c"] = _cache.str();
    }

    template<class Base>
    void CLangCompileDynamicHelper<Base>::generateModelsSource(std::map<std::string, std::string>& sources) {
        _cache.str("");
        _cache << "void " << FUNCTION_MODELS << "(char const *const** names, int* count) {\n"
                "   static const char* const models[] = {\n";

        typename std::map<std::string, CLangCompileModelHelper<Base>*>::const_iterator it;
        for (it = _models.begin(); it != _models.end(); ++it) {
            if (it != _models.begin()) {
                _cache << ",\n";
            }
            _cache << "      \"" << it->first << "\"";
        }
        _cache << "};\n"
                "   *names = models;\n"
                "   *count = " << _models.size() << ";\n"
                "}\n\n";

        sources[FUNCTION_MODELS + ".c"] = _cache.str();
    }

}

#endif
