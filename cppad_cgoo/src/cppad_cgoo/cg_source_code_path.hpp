#ifndef CPPAD_CG_SOURCE_CODE_PATH_INCLUDED
#define	CPPAD_CG_SOURCE_CODE_PATH_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    template<class Base>
    struct SourceCodePathNode {
        size_t arg_index;
        SourceCodeFragment<Base>* node;

        inline SourceCodePathNode(SourceCodeFragment<Base>* node_, size_t arg_index_) :
            arg_index(arg_index_),
            node(node_) {
        }

    };

    template<class Base>
    inline std::vector<std::vector<SourceCodePathNode<Base> > > findPaths(SourceCodeFragment<Base>* root,
                                                                          SourceCodeFragment<Base>* code,
                                                                          size_t max) {

        std::vector<std::vector<SourceCodePathNode<Base> > > found;

        if (max > 0) {
            std::vector<SourceCodePathNode<Base> > path2node;
            path2node.reserve(30);
            path2node.push_back(SourceCodePathNode<Base > (root, 0));

            if (root == code) {
                found.push_back(path2node);
            } else {
                findPaths(path2node, code, found, max);
            }
        }

        return found;
    }

    template<class Base>
    inline void findPaths(std::vector<SourceCodePathNode<Base> >& currPath,
                          SourceCodeFragment<Base>* code,
                          std::vector<std::vector<SourceCodePathNode<Base> > >& found,
                          size_t max) {

        if (found.size() == max) {
            return;
        }

        if (code == currPath.back().node) {
            found.push_back(currPath);
            return;
        }

        const std::vector<Argument<Base> >& args = currPath.back().node->arguments();
        for (size_t i = 0; i < args.size(); ++i) {
            SourceCodeFragment<Base>* a = args[i].operation();
            if (a != NULL) {
                currPath.push_back(SourceCodePathNode<Base > (a, i));
                findPaths(currPath, code, found, max);
                currPath.pop_back();
            }
        }
    }
}

#endif
