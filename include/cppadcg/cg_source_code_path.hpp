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

        inline SourceCodePathNode() :
            arg_index(0),
            node(NULL) {
        }

        inline SourceCodePathNode(SourceCodeFragment<Base>* node_, size_t arg_index_) :
            arg_index(arg_index_),
            node(node_) {
        }

    };

    template<class Base>
    inline std::vector<std::vector<SourceCodePathNode<Base> > > CodeHandler<Base>::findPaths(SourceCodeFragment<Base>& root,
                                                                                             SourceCodeFragment<Base>& code,
                                                                                             size_t max) {
        resetCounters();

        std::vector<std::vector<SourceCodePathNode<Base> > > found;

        if (max > 0) {
            std::vector<SourceCodePathNode<Base> > path2node;
            path2node.reserve(30);
            path2node.push_back(SourceCodePathNode<Base> (&root, 0));

            if (&root == &code) {
                found.push_back(path2node);
            } else {
                findPaths(path2node, code, found, max);
            }
        }

        return found;
    }

    template<class Base>
    inline void CodeHandler<Base>::findPaths(std::vector<SourceCodePathNode<Base> >& currPath,
                                             SourceCodeFragment<Base>& code,
                                             std::vector<std::vector<SourceCodePathNode<Base> > >& found,
                                             size_t max) {

        SourceCodeFragment<Base>* currNode = currPath.back().node;
        if (&code == currNode) {
            found.push_back(currPath);
            return;
        }

        const std::vector<Argument<Base> >& args = currNode->arguments();
        if (args.empty())
            return; // nothing to look in

        if (currNode->usageCount() > 0) {
            // already searched inside this node
            // any match would have been saved in found
            std::vector<SourceCodePath> pathsFromNode = findPathsFromNode(found, *currNode);
            typename std::vector<std::vector<SourceCodePathNode<Base> > >::const_iterator it;
            for (it = pathsFromNode.begin(); it != pathsFromNode.end(); ++it) {
                const SourceCodePath& pathFromNode = *it;
                SourceCodePath newPath(currPath.size() + pathFromNode.size());
                std::copy(currPath.begin(), currPath.end(), newPath.begin());
                std::copy(pathFromNode.begin(), pathFromNode.end(), newPath.begin() + currPath.size());
                found.push_back(newPath);
            }

        } else {
            // not visited yet
            currNode->increaseUsageCount(); // mark node as visited

            size_t size = args.size();
            for (size_t i = 0; i < size; ++i) {
                SourceCodeFragment<Base>* a = args[i].operation();
                if (a != NULL) {
                    currPath.push_back(SourceCodePathNode<Base> (a, i));
                    findPaths(currPath, code, found, max);
                    currPath.pop_back();
                    if (found.size() == max) {
                        return;
                    }
                }
            }
        }
    }

    template<class Base>
    inline std::vector<std::vector<SourceCodePathNode<Base> > > CodeHandler<Base>::findPathsFromNode(const std::vector<std::vector<SourceCodePathNode<Base> > > nodePaths,
                                                                                                     SourceCodeFragment<Base>& node) {

        std::vector<SourceCodePath> foundPaths;
        std::set<size_t> argsFound;

        typename std::vector<SourceCodePath>::const_iterator it;
        for (it = nodePaths.begin(); it != nodePaths.end(); ++it) {
            const SourceCodePath& path = *it;
            size_t size = path.size();
            for (size_t i = 0; i < size - 1; i++) {
                const SourceCodePathNode<Base>& pnode = path[i];
                if (pnode.node == &node) {
                    if (argsFound.find(path[i + 1].arg_index) == argsFound.end()) {
                        foundPaths.push_back(SourceCodePath(path.begin() + i + 1, path.end()));
                        argsFound.insert(path[i + 1].arg_index);
                    }
                }
            }
        }

        return foundPaths;

    }
}

#endif
