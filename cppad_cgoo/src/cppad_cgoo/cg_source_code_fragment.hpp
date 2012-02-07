#ifndef CPPAD_CG_SOURCE_CODE_FRAGMENT_INCLUDED
#define	CPPAD_CG_SOURCE_CODE_FRAGMENT_INCLUDED
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

    typedef struct SourceCodeFragment {
        // the operations used to create this variable (temporary variables only)
        std::string operations;
        // status of the operations
        OpContainement opTypes;
        // the code blocks this block depends upon
        std::set<SourceCodeBlock*> depends;

        SourceCodeFragment(const std::string& ops, OpContainement opContaint) :
            operations(ops),
            opTypes(opContaint) {
        }

        SourceCodeFragment(const std::string& ops, OpContainement opContaint, const std::set<SourceCodeBlock*>& deps) :
            operations(ops),
            opTypes(opContaint),
            depends(deps) {
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

    } CodeFragment;

}

#endif

