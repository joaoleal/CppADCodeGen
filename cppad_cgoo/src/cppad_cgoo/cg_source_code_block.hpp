#ifndef CPPAD_CG_SOURCE_CODE_BLOCK_INCLUDED
#define	CPPAD_CG_SOURCE_CODE_BLOCK_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <string>
#include <set>

namespace CppAD {

    typedef struct SourceCodeBlock {
        // variable ID that was altered/assigned in this source code
        size_t id;
        // source code (it can be an empty string for independent variables)
        std::string code;
        // flag indicating if the this block should be printed
        bool print;
        // the code blocks this block depends upon (empty for independent 
        // variables and possibly for the 1st assignment of a dependent variable)
        std::set<SourceCodeBlock*> depends;

        SourceCodeBlock(size_t varID) :
            id(varID),
            print(false) {
        }

        template<class Base>
        inline void addDependency(const CG<Base> &var) {
            if (var.isVariable()) {
                assert(var.getSourceCodeBlock() != NULL);
                depends.insert(var.getSourceCodeBlock());
            } else if (var.isTemporaryVariable()) {
                std::set<SourceCodeBlock*>& deps1 = var.getSourceCodeFragment()->depends;
                if (depends.empty()) {
                    depends = deps1;
                } else {
                    depends.insert(deps1.begin(), deps1.end());
                }

                for (std::set<SourceCodeBlock*>::const_iterator it = deps1.begin(); it != deps1.end(); ++it) {
                    assert(*it != NULL);
                }
            }
        }
    } CodeBlock;

}

#endif

