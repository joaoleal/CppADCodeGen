#ifndef CPPAD_CG_SOURCE_CODE_PATH_INCLUDED
#define CPPAD_CG_SOURCE_CODE_PATH_INCLUDED
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

    /**
     * Finds occurences of a source code fragment in an operation graph.
     * 
     * \param root the operation graph where to search
     * \param code the source code fragment to find in root
     * \param max the maximum number of occurences of code to find in root
     * \return the paths from root to code
     */
    template<class Base>
    inline std::vector<std::vector<SourceCodePathNode<Base> > > findPaths(SourceCodeFragment<Base>& root,
                                                                          SourceCodeFragment<Base>& code,
                                                                          size_t max) {

        std::vector<std::vector<SourceCodePathNode<Base> > > found;

        if (max > 0) {
            std::vector<SourceCodePathNode<Base> > path2node;
            path2node.reserve(30);
            path2node.push_back(SourceCodePathNode<Base > (&root, 0));

            if (&root == &code) {
                found.push_back(path2node);
            } else {
                findPaths(path2node, code, found, max);
            }
        }

        return found;
    }

    template<class Base>
    inline void findPaths(std::vector<SourceCodePathNode<Base> >& currPath,
                          SourceCodeFragment<Base>& code,
                          std::vector<std::vector<SourceCodePathNode<Base> > >& found,
                          size_t max) {

        if (found.size() == max) {
            return;
        }

        if (&code == currPath.back().node) {
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
