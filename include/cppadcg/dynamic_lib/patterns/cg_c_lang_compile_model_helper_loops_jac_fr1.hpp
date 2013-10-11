#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_JAC_FR1_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_JAC_FR1_INCLUDED
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
    void CLangCompileModelHelper<Base>::generateSparseJacobianWithLoopsSourceFromForRev(std::map<std::string, std::string>& sources,
                                                                                        const std::map<size_t, std::vector<std::set<size_t> > >& userJacElLocation,
                                                                                        const std::map<size_t, bool>& keyOrdered,
                                                                                        size_t maxCompressedSize,
                                                                                        const std::string& localFunctionTypeName,
                                                                                        const std::string& suffix,
                                                                                        const std::string& keyName,
                                                                                        const std::map<size_t, std::set<size_t> >& nonLoopElements,
                                                                                        const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                                                                        void (*generateLocalFunctionName)(std::ostringstream& cache, const std::string& modelName, const LoopModel<Base>& loop, size_t g)) throw (CGException) {
        using namespace std;
        using namespace CppAD::loops;

        /**
         * Generate the source code
         */
        CLanguage<Base> langC(_baseTypeName);
        string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        string model_function = _name + "_" + FUNCTION_SPARSE_JACOBIAN;
        string localFunction = _name + "_" + localFunctionTypeName;
        string nlSuffix = "noloop_" + suffix;
        IndexDclrOperationNode<Base> indexIt("it");

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";

        generateFunctionDeclarationSource(_cache, localFunction, nlSuffix, nonLoopElements, argsDcl);
        generateFunctionDeclarationSourceLoopForRev(_cache, langC, _name, keyName, loopGroups, generateLocalFunctionName);

        _cache << "\n";
        printForRevUsageFunction<Base>(_cache, _baseTypeName, _name,
                model_function, 2,
                localFunction, suffix,
                keyName, indexIt, "jac",
                loopGroups,
                nonLoopElements,
                userJacElLocation, keyOrdered,
                generateLocalFunctionName,
                _jacSparsity.rows.size(), maxCompressedSize);

        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

}

#endif