#ifndef CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_R2_INCLUDED
#define CPPAD_CG_C_LANG_COMPILE_MODEL_HELPER_LOOPS_HESS_R2_INCLUDED
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
    void CLangCompileModelHelper<Base>::generateSparseHessianWithLoopsSourceFromRev2(std::map<std::string, std::string>& sources,
                                                                                     const std::map<size_t, std::vector<std::set<size_t> > >& userHessElLocation,
                                                                                     const std::map<size_t, bool>& jrowOrdered,
                                                                                     size_t maxCompressedSize) throw (CGException) {
        using namespace std;
        using namespace CppAD::loops;

        /**
         * Generate the source code
         */
        CLanguage<Base> langC(_baseTypeName);
        string argsDcl = langC.generateDefaultFunctionArgumentsDcl();

        string model_function = _name + "_" + FUNCTION_SPARSE_HESSIAN;
        string functionRev2 = _name + "_" + FUNCTION_SPARSE_REVERSE_TWO;
        string suffix = "indep";
        string nlRev2Suffix = "noloop_" + suffix;
        IndexDclrOperationNode<Base> indexIt("it");

        _cache.str("");
        _cache << "#include <stdlib.h>\n"
                << CLanguage<Base>::ATOMICFUN_STRUCT_DEFINITION << "\n\n";
        generateFunctionDeclarationSource(_cache, functionRev2, nlRev2Suffix, _nonLoopRev2Elements, argsDcl);
        generateFunctionDeclarationSourceLoopForRev(_cache, langC, _name, "jrow", _loopRev2Groups, generateFunctionNameLoopRev2);

        _cache << "\n";

        printForRevUsageFunction(_cache, _baseTypeName, _name,
                                 model_function, 3,
                                 functionRev2, suffix,
                                 "jrow", indexIt, "hess",
                                 _loopRev2Groups,
                                 _nonLoopRev2Elements,
                                 userHessElLocation, jrowOrdered,
                                 generateFunctionNameLoopRev2,
                                 _hessSparsity.rows.size(), maxCompressedSize);

        sources[model_function + ".c"] = _cache.str();
        _cache.str("");
    }

}

#endif