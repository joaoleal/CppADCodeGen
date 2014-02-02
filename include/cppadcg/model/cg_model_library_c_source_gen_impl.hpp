#ifndef CPPAD_CG_MODEL_LIBRARY_C_SOURCE_GEN_IMPL_INCLUDED
#define CPPAD_CG_MODEL_LIBRARY_C_SOURCE_GEN_IMPL_INCLUDED
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

#include <typeinfo>

namespace CppAD {

    template<class Base>
    const unsigned long ModelLibraryCSourceGen<Base>::API_VERSION = 5;

    template<class Base>
    const std::string ModelLibraryCSourceGen<Base>::FUNCTION_VERSION = "cppad_cg_version";

    template<class Base>
    const std::string ModelLibraryCSourceGen<Base>::FUNCTION_MODELS = "cppad_cg_models";

    template<class Base>
    const std::string ModelLibraryCSourceGen<Base>::CONST = "const";

    template<class Base>
    void ModelLibraryCSourceGen<Base>::saveSources(const std::string& sourcesFolder) {

        // create the folder if it does not exist
        system::createFolder(sourcesFolder);

        // save/generate model sources
        for (const auto& it : _models) {
            saveSources(sourcesFolder, it.second->getSources());
        }

        // save/generate library sources
        saveSources(sourcesFolder, getLibrarySources());

        // save custom user sources
        saveSources(sourcesFolder, getCustomSources());
    }

    template<class Base>
    void ModelLibraryCSourceGen<Base>::saveSources(const std::string& sourcesFolder,
                                                   const std::map<std::string, std::string>& sources) {
        for (const auto& it : sources) {
            // for debugging purposes only
            std::ofstream sourceFile;
            std::string file = system::createPath(sourcesFolder, it.first);
            sourceFile.open(file.c_str());
            sourceFile << it.second;
            sourceFile.close();
        }
    }

    template<class Base>
    const std::map<std::string, std::string>& ModelLibraryCSourceGen<Base>::getLibrarySources() {
        if (_libSources.empty()) {
            generateVerionSource(_libSources);
            generateModelsSource(_libSources);
        }

        return _libSources;
    }

    template<class Base>
    void ModelLibraryCSourceGen<Base>::generateVerionSource(std::map<std::string, std::string>& sources) {
        _cache.str("");
        _cache << "unsigned long " << FUNCTION_VERSION << "() {\n"
                << "   return " << API_VERSION << "u;\n"
                << "}\n\n";

        sources[FUNCTION_VERSION + ".c"] = _cache.str();
    }

    template<class Base>
    void ModelLibraryCSourceGen<Base>::generateModelsSource(std::map<std::string, std::string>& sources) {
        _cache.str("");
        _cache << "void " << FUNCTION_MODELS << "(char const *const** names, int* count) {\n"
                "   static const char* const models[] = {\n";

        for (auto it = _models.begin(); it != _models.end(); ++it) {
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
