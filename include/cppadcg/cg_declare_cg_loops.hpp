#ifndef CPPAD_CG_DECLARE_CG_LOOPS_INCLUDED
#define CPPAD_CG_DECLARE_CG_LOOPS_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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

// forward declarations
namespace CppAD {

    template<class Base>
    class vector;

    namespace loops {

        typedef std::pair<size_t, size_t> SizeN1stIt;

        typedef std::pair<size_t, size_t> pairss;

        class JacobianWithLoopsRowInfo;

        class HessianElement;

        template<class Base>
        class IfBranchInfo;

        template <class Base>
        class IfElseInfo;

        template<class Base>
        class JacobianTermContrib;

        template<class Base>
        class JacobianColGroup;

        template<class Base>
        class HessianWithLoopsInfo;

        template<class Base>
        class HessianRowGroup;

        class ArrayGroup;

        template<class Base>
        inline vector<CG<Base> > createIndexedIndependents(CodeHandler<Base>& handler,
                                                           LoopModel<Base>& loop,
                                                           IndexOperationNode<Base>& iterationIndexOp);

        template<class Base>
        inline vector<CG<Base> > createLoopIndependentVector(CodeHandler<Base>& handler,
                                                             LoopModel<Base>& loop,
                                                             const vector<CG<Base> >& indexedIndeps,
                                                             const vector<CG<Base> >& nonIndexedIndeps,
                                                             const vector<CG<Base> >& nonIndexedTmps);

        template<class Base>
        inline vector<CG<Base> > createLoopDependentVector(CodeHandler<Base>& handler,
                                                           LoopModel<Base>& loop,
                                                           IndexOperationNode<Base>& iterationIndexOp);

        template<class Base>
        inline LoopEndOperationNode<Base>* createLoopEnd(CodeHandler<Base>& handler,
                                                         LoopStartOperationNode<Base>& loopStart,
                                                         const vector<std::pair<CG<Base>, IndexPattern*> >& indexedLoopResults,
                                                         const std::set<IndexOperationNode<Base>*>& indexesOps,
                                                         size_t assignOrAdd);

        template<class Base>
        inline void moveNonIndexedOutsideLoop(LoopStartOperationNode<Base>& loopStart,
                                              LoopEndOperationNode<Base>& loopEnd);

        template<class Base>
        inline bool findNonIndexedNodes(OperationNode<Base>& node,
                                        std::set<OperationNode<Base>*>& nonIndexed,
                                        const IndexDclrOperationNode<Base>& loopIndex);

        template<class Base>
        inline IfElseInfo<Base>* findExistingIfElse(vector<IfElseInfo<Base> >& ifElses,
                                                    const std::map<SizeN1stIt, std::pair<size_t, std::set<size_t> > >& first2Iterations);

        template<class Base>
        OperationNode<Base>* createIndexConditionExpression(const std::set<size_t>& iterations,
                                                            const std::set<size_t>& usedIter,
                                                            size_t maxIter,
                                                            IndexOperationNode<Base>& iterationIndexOp);

        template<class Base>
        inline void determineForRevUsagePatterns(const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& loopGroups,
                                                 const std::map<size_t, std::vector<std::set<size_t> > >& userElLocation,
                                                 const std::map<size_t, bool>& ordered,
                                                 std::map<size_t, std::map<LoopModel<Base>*, std::map<size_t, ArrayGroup*> > >& loopCalls,
                                                 SmartVectorPointer<ArrayGroup>& garbage);

        template<class Base>
        void generateFunctionDeclarationSourceLoopForRev(std::ostringstream& cache,
                                                         CLanguage<Base>& langC,
                                                         const std::string& modelName,
                                                         const std::string& keyName,
                                                         const std::map<LoopModel<Base>*, std::map<size_t, std::map<size_t, std::set<size_t> > > >& _loopRev2Groups,
                                                         void (*generateFunctionNameLoopRev2)(std::ostringstream& cache, const std::string& modelName, const LoopModel<Base>& loop, size_t g));
    }

}

#endif
