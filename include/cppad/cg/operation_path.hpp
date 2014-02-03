#ifndef CPPAD_CG_OPERATION_PATH_INCLUDED
#define CPPAD_CG_OPERATION_PATH_INCLUDED
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

namespace CppAD {
namespace cg {

template<class Base>
struct OperationPathNode {
    size_t arg_index;
    OperationNode<Base>* node;

    inline OperationPathNode() :
        arg_index(0),
        node(nullptr) {
    }

    inline OperationPathNode(OperationNode<Base>* node_, size_t arg_index_) :
        arg_index(arg_index_),
        node(node_) {
    }

};

template<class Base>
inline std::vector<std::vector<OperationPathNode<Base> > > CodeHandler<Base>::findPaths(OperationNode<Base>& root,
                                                                                        OperationNode<Base>& code,
                                                                                        size_t max) {
    resetManagedNodes();

    std::vector<std::vector<OperationPathNode<Base> > > found;

    if (max > 0) {
        std::vector<OperationPathNode<Base> > path2node;
        path2node.reserve(30);
        path2node.push_back(OperationPathNode<Base> (&root, 0));

        if (&root == &code) {
            found.push_back(path2node);
        } else {
            findPaths(path2node, code, found, max);
        }
    }

    return found;
}

template<class Base>
inline void CodeHandler<Base>::findPaths(std::vector<OperationPathNode<Base> >& currPath,
                                         OperationNode<Base>& code,
                                         std::vector<std::vector<OperationPathNode<Base> > >& found,
                                         size_t max) {

    OperationNode<Base>* currNode = currPath.back().node;
    if (&code == currNode) {
        found.push_back(currPath);
        return;
    }

    const std::vector<Argument<Base> >& args = currNode->getArguments();
    if (args.empty())
        return; // nothing to look in

    if (currNode->getUsageCount() > 0) {
        // already searched inside this node
        // any match would have been saved in found
        std::vector<SourceCodePath> pathsFromNode = findPathsFromNode(found, *currNode);
        for (const SourceCodePath& pathFromNode : pathsFromNode) {
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
            OperationNode<Base>* a = args[i].getOperation();
            if (a != nullptr) {
                currPath.push_back(OperationPathNode<Base> (a, i));
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
inline std::vector<std::vector<OperationPathNode<Base> > > CodeHandler<Base>::findPathsFromNode(const std::vector<std::vector<OperationPathNode<Base> > > nodePaths,
                                                                                                OperationNode<Base>& node) {

    std::vector<SourceCodePath> foundPaths;
    std::set<size_t> argsFound;

    for (const SourceCodePath& path : nodePaths) {
        size_t size = path.size();
        for (size_t i = 0; i < size - 1; i++) {
            const OperationPathNode<Base>& pnode = path[i];
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

} // END cg namespace
} // END CppAD namespace

#endif