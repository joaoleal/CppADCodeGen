#ifndef CPPAD_CG_CODE_HANDLER_INCLUDED
#define CPPAD_CG_CODE_HANDLER_INCLUDED
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

    template<class Base>
    class CG;

    /**
     * Helper class to analyze the operation graph and generate source code
     * for several languages
     * 
     * @author Joao Leal
     */
    template<class Base>
    class CodeHandler {
    public:
        typedef std::vector<OperationPathNode<Base> > SourceCodePath;
        typedef std::vector<ScopePathElement<Base> > ScopePath;
    protected:
        // counter used to generate variable IDs
        size_t _idCount;
        // counter used to generate array variable IDs
        size_t _idArrayCount;
        // counter used to generate IDs for atomic functions
        size_t _idAtomicCount;
        // the independent variables
        std::vector<OperationNode<Base> *> _independentVariables;
        // the current dependent variables
        vector<CG<Base> >* _dependents;
        // all the source code blocks created with the CG<Base> objects
        std::vector<OperationNode<Base> *> _codeBlocks;
        // the order for the variable creation in the source code
        std::vector<OperationNode<Base> *> _variableOrder;
        // maps the ids of the atomic functions
        std::map<size_t, CGAbstractAtomicFun<Base>*> _atomicFunctions;
        // maps the loop ids of the loop atomic functions
        std::map<size_t, LoopModel<Base>*> _loops;
        // the used indexes
        std::set<const IndexDclrOperationNode<Base>*> _indexes;
        //
        std::vector<const IndexPattern*> _loopDependentIndexPatterns;
        std::vector<const IndexPattern*> _loopDependentIndexPatternManaged; // garbage collection
        std::vector<const IndexPattern*> _loopIndependentIndexPatterns;
        /**
         * already used atomic function names (may contain names which were 
         * used by previous calls to this/other CondeHandlers)
         */
        std::set<std::string> _atomicFunctionsSet;
        /**
         * the order of the atomic functions(may contain names which were 
         * used by previous calls to this/other CondeHandlers)
         */
        std::vector<std::string>* _atomicFunctionsOrder;
        // a flag indicating if this handler was previously used to generate code
        bool _used;
        // a flag indicating whether or not to reuse the IDs of destroyed variables
        bool _reuseIDs;
        // the used variables inside a loop from the outside (for different loop depths)
        std::vector<std::set<OperationNode<Base>*> > _loopOuterVars;
        // the current loop depth (-1 means no loop)
        int _loopDepth;
        // the evaluation order of the loop start for each loop depth
        std::vector<size_t> _loopStartEvalOrder;
        // scope color/index counter
        size_t _scopeColorCount;
        // the current scope color/index counter
        size_t _currentScopeColor;
        // all the scopes
        std::vector<ScopePath> _scopes;
        //
        std::map<size_t, std::set<OperationNode<Base>*> > _scopeExtraDependecy;
        // the language used for source code generation
        Language<Base>* _lang;
        // the lowest ID used for temporary variables
        size_t _minTemporaryVarID;
        /**
         * whether or not the dependent variables should be zeroed before 
         * executing the operation graph
         */
        bool _zeroDependents;
        //
        bool _verbose;
    public:

        CodeHandler(size_t varCount = 50) :
            _idCount(1),
            _idArrayCount(1),
            _idAtomicCount(1),
            _dependents(NULL),
            _atomicFunctionsOrder(NULL),
            _used(false),
            _reuseIDs(true),
            _loopDepth(-1),
            _scopeColorCount(0),
            _currentScopeColor(0),
            _lang(NULL),
            _minTemporaryVarID(0),
            _zeroDependents(false),
            _verbose(false) {
            _codeBlocks.reserve(varCount);
            _variableOrder.reserve(1 + varCount / 3);
        }

        inline void setReuseVariableIDs(bool reuse) {
            _reuseIDs = reuse;
        }

        inline bool isReuseVariableIDs() const {
            return _reuseIDs;
        }

        template<class VectorCG>
        inline void makeVariables(VectorCG& variables) {
            for (size_t i = 0; i < variables.size(); i++) {
                makeVariable(variables[i]);
            }
        }

        inline void makeVariables(std::vector<AD<CG<Base> > >& variables) {
            for (typename std::vector<AD<CG<Base> > >::iterator it = variables.begin(); it != variables.end(); ++it) {
                makeVariable(*it);
            }
        }

        inline void makeVariable(AD<CG<Base> >& variable) {
            CG<Base> v;
            makeVariable(v); // make it a codegen variable
            variable = v; // variable id now the same as v
        }

        inline void makeVariable(CG<Base>& variable) {
            _independentVariables.push_back(new OperationNode<Base> (CGInvOp));
            variable.makeVariable(*this, _independentVariables.back());
        }

        size_t getIndependentVariableSize() const {
            return _independentVariables.size();
        }

        size_t getIndependentVariableIndex(const OperationNode<Base>& var) const throw (CGException) {
            assert(var.getOperationType() == CGInvOp);

            typename std::vector<OperationNode<Base> *>::const_iterator it =
                    std::find(_independentVariables.begin(), _independentVariables.end(), &var);
            if (it == _independentVariables.end()) {
                throw CGException("Variable not found in the independent variable vector");
            }

            return it - _independentVariables.begin();
        }

        inline size_t getMaximumVariableID() const {
            return _idCount;
        }

        inline bool isVerbose() const {
            return _verbose;
        }

        inline void setVerbose(bool verbose) {
            _verbose = verbose;
        }

        /**
         * Determines whether or not the dependent variables will be set to zero
         * before executing the operation graph
         * 
         * @return true if the dependents will be zeroed
         */
        inline bool isZeroDependents() const {
            return _zeroDependents;
        }

        /**
         * Defines whether or not the dependent variables should be set to zero
         * executing the operation graph
         * 
         * @param true if the dependents should be zeroed
         */
        inline void setZeroDependents(bool zeroDependents) {
            _zeroDependents = zeroDependents;
        }

        /**
         * Provides the name used by an atomic function with a given ID.
         * 
         * @param id the atomic function ID.
         * @return a pointer to the atomic function name if it was registered
         *         or NULL otherwise
         */
        inline const std::string* getAtomicFunctionName(size_t id) const {
            typename std::map<size_t, CGAbstractAtomicFun<Base>*>::const_iterator it;
            it = _atomicFunctions.find(id);
            if (it != _atomicFunctions.end())
                return &(it->second->afun_name());
            else
                return NULL;
        }

        /**
         * Provides a map with all the currently registered atomic functions.
         * 
         * @return a map with the atomic function ID as key and the atomic 
         *         function as value
         */
        inline const std::map<size_t, CGAbstractAtomicFun<Base>* >& getAtomicFunctions() const {
            return _atomicFunctions;
        }

        /**
         * Provides the name used by a loop atomic function with a given ID.
         * 
         * @param id the atomic function ID.
         * @return a pointer to the atomic loop function name if it was
         *         registered or NULL otherwise
         */
        inline const std::string* getLoopName(size_t id) const {
            typename std::map<size_t, LoopModel<Base>*>::const_iterator it;
            it = _loops.find(id);
            if (it != _loops.end())
                return &(it->second->afun_name());
            else
                return NULL;
        }

        inline const std::vector<ScopePath>& getScopes() const {
            return _scopes;
        }

        /***********************************************************************
         *                   Graph management functions
         **********************************************************************/
        /**
         * Finds occurences of a source code fragment in an operation graph.
         * 
         * @param root the operation graph where to search
         * @param code the source code fragment to find in root
         * @param max the maximum number of occurences of code to find in root
         * @return the paths from root to code
         */
        inline std::vector<SourceCodePath> findPaths(OperationNode<Base>& root,
                                                     OperationNode<Base>& code,
                                                     size_t max);

        inline bool isSolvable(const SourceCodePath& path) throw (CGException);

        /***********************************************************************
         *                   Source code generation
         **********************************************************************/

        /**
         * Creates the source code from the operations registered so far.
         * 
         * @param out The output stream where the source code is to be printed.
         * @param lang The targeted language.
         * @param dependent The dependent variables for which the source code
         *                  should be generated. By defining this vector the 
         *                  number of operations in the source code can be 
         *                  reduced and thus providing a more optimized code.
         * @param nameGen Provides the rules for variable name creation.
         */
        virtual void generateCode(std::ostream& out,
                                  CppAD::Language<Base>& lang,
                                  vector<CG<Base> >& dependent,
                                  VariableNameGenerator<Base>& nameGen,
                                  const std::string& jobName = "source") {
            std::vector<std::string> atomicFunctions;
            generateCode(out, lang, dependent, nameGen, atomicFunctions, jobName);
        }

        /**
         * Creates the source code from the operations registered so far.
         * 
         * @param out The output stream where the source code is to be printed.
         * @param lang The targeted language.
         * @param dependent The dependent variables for which the source code
         *                  should be generated. By defining this vector the 
         *                  number of operations in the source code can be 
         *                  reduced and thus providing a more optimized code.
         * @param nameGen Provides the rules for variable name creation.
         * @param atomicFunctions The order of the atomic functions.
         */
        virtual void generateCode(std::ostream& out,
                                  CppAD::Language<Base>& lang,
                                  vector<CG<Base> >& dependent,
                                  VariableNameGenerator<Base>& nameGen,
                                  std::vector<std::string>& atomicFunctions,
                                  const std::string& jobName = "source") {
            double beginTime;
            if (_verbose) {
                std::cout << "generating source for '" << jobName << "' ... ";
                std::cout.flush();
                beginTime = system::currentTime();
            }

            _lang = &lang;
            _idCount = 1;
            _idArrayCount = 1;
            _idAtomicCount = 1;
            _dependents = &dependent;
            _atomicFunctionsOrder = &atomicFunctions;
            _atomicFunctionsSet.clear();
            _indexes.clear();
            _loopOuterVars.clear();
            _loopDepth = -1;
            _scopeColorCount = 0;
            _currentScopeColor = 0;
            _scopes.reserve(4);
            _scopeExtraDependecy.clear();
            _loopStartEvalOrder.clear();
            for (size_t i = 0; i < atomicFunctions.size(); i++) {
                _atomicFunctionsSet.insert(atomicFunctions[i]);
            }

            if (_used) {
                resetManagedNodes();
            }
            _used = true;

            /**
             * the first variable IDs are for the independent variables
             */
            size_t n = _independentVariables.size();
            for (size_t j = 0; j < n; j++) {
                _independentVariables[j]->setVariableID(_idCount++);
            }

            size_t m = dependent.size();
            for (size_t i = 0; i < m; i++) {
                OperationNode<Base>* node = dependent[i].getOperationNode();
                if (node != NULL && node->getVariableID() == 0) {
                    node->setVariableID(_idCount++);
                }
            }

            _minTemporaryVarID = _idCount;

            /**
             * determine the number of times each variable is used
             */
            for (size_t i = 0; i < m; i++) {
                OperationNode<Base>* node = dependent[i].getOperationNode();
                if (node != NULL) {
                    markCodeBlockUsed(*node);
                }
            }

            /**
             * add some additional arguments to some nodes related with other 
             * nodes being used in different scopes
             */
            typename std::map<size_t, std::set<OperationNode<Base>*> >::const_iterator itse;
            for (itse = _scopeExtraDependecy.begin(); itse != _scopeExtraDependecy.end(); ++itse) {
                assert(itse->first < _scopes.size());
                assert(!_scopes[itse->first].empty());
                assert(_scopes[itse->first].back().beginning != NULL);

                OperationNode<Base>& div = *_scopes[itse->first].back().beginning;
                const std::set<OperationNode<Base>*>& extraDeps = itse->second;

                typename std::set<OperationNode<Base>*>::const_iterator itn;
                for (itn = extraDeps.begin(); itn != extraDeps.end(); ++itn) {
                    OperationNode<Base>& code = **itn;
                    if (!containsArgument(div, code)) {
                        div.getArguments().push_back(Argument<Base>(code));
                    }
                }
            }
            _scopeExtraDependecy.clear();

            /**
             * determine the variable creation order
             */
            for (size_t i = 0; i < m; i++) {
                CG<Base>& var = dependent[i];
                if (var.getOperationNode() != NULL) {
                    OperationNode<Base>& code = *var.getOperationNode();
                    if (code.getUsageCount() == 0) {
                        // dependencies not visited yet
                        checkVariableCreation(code);

                        // make sure new temporary variables are NOT created for
                        // the independent variables and that a dependency did
                        // not use it first
                        if ((code.getVariableID() == 0 || !isIndependent(code)) && code.getUsageCount() == 0) {
                            addToEvaluationQueue(code);
                        }
                    }
                    code.increaseUsageCount();
                }
            }

            if (_reuseIDs) {
                reduceTemporaryVariables(dependent);
            }

            nameGen.setTemporaryVariableID(_minTemporaryVarID, _idCount - 1, _idArrayCount - 1);

            std::map<std::string, size_t> atomicFunctionName2Id;
            typename std::map<size_t, CGAbstractAtomicFun<Base>*>::iterator itA;
            for (itA = _atomicFunctions.begin(); itA != _atomicFunctions.end(); ++itA) {
                atomicFunctionName2Id[itA->second->afun_name()] = itA->first;
            }

            std::map<size_t, size_t> atomicFunctionId2Index;
            std::map<size_t, std::string> atomicFunctionId2Name;
            for (size_t i = 0; i < _atomicFunctionsOrder->size(); i++) {
                const std::string& atomicName = (*_atomicFunctionsOrder)[i];
                std::map<std::string, size_t>::const_iterator it = atomicFunctionName2Id.find(atomicName);
                if (it != atomicFunctionName2Id.end()) {
                    atomicFunctionId2Index[it->second] = i;
                    atomicFunctionId2Name[it->second] = atomicName;
                }
            }

            /**
             * Creates the source code for a specific language
             */
            LanguageGenerationData<Base> info(_independentVariables, dependent,
                                              _minTemporaryVarID, _variableOrder,
                                              nameGen,
                                              atomicFunctionId2Index, atomicFunctionId2Name,
                                              _reuseIDs,
                                              _indexes,
                                              _loopDependentIndexPatterns, _loopIndependentIndexPatterns,
                                              _zeroDependents);
            lang.generateSourceCode(out, info);

            _atomicFunctionsSet.clear();

            if (_verbose) {
                double endTime = system::currentTime();
                std::cout << "done [" << std::fixed << std::setprecision(3)
                        << (endTime - beginTime) << "]" << std::endl;
            }
        }

        size_t getTemporaryVariableCount() const {
            if (_idCount == 1)
                return 0; // no code generated
            else
                return _idCount - _minTemporaryVarID;
        }

        size_t getTemporaryArraySize() const {
            return _idArrayCount - 1;
        }

        /***********************************************************************
         *                   Reusing handler and nodes
         **********************************************************************/

        /**
         * Resets this handler for a usage with completely different nodes.
         * @warning all managed memory will be deleted
         */
        virtual void reset() {
            typename std::vector<OperationNode<Base> *>::iterator itc;
            for (itc = _codeBlocks.begin(); itc != _codeBlocks.end(); ++itc) {
                delete *itc;
            }
            _codeBlocks.clear();
            _independentVariables.clear();
            _idCount = 1;
            _idArrayCount = 1;
            _idAtomicCount = 1;

            _loops.clear();
            _indexes.clear();
            _loopDependentIndexPatterns.clear();
            _loopIndependentIndexPatterns.clear();

            std::vector<const IndexPattern*>::const_iterator itip;
            for (itip = _loopDependentIndexPatternManaged.begin(); itip != _loopDependentIndexPatternManaged.end(); ++itip) {
                delete *itip;
            }
            _loopDependentIndexPatternManaged.clear();

            _used = false;
        }

        /**
         * Resets the previously used dependents and their children so that they
         * can be reused again by this handler.
         */
        inline void resetNodes() {
            if (_dependents != NULL)
                resetNodes(*_dependents);
        }

        /**
         * Resets the nodes and their children so that they can be reused again
         * by a handler
         * @param dependents the nodes to be reset
         */
        static inline void resetNodes(vector<CG<Base> >& dependents) {
            for (size_t i = 0; i < dependents.size(); i++) {
                resetNodes(dependents[i].getOperationNode());
            }
        }

        /**
         * Resets a node and its children so that they can be reused again by a
         * handler
         * @param node the node to be reset
         */
        static inline void resetNodes(OperationNode<Base>* node) {
            if (node == NULL || node->getTotalUsageCount() == 0)
                return;

            node->resetHandlerCounters();
            node->setColor(0);

            const std::vector<Argument<Base> >& args = node->getArguments();
            for (size_t a = 0; a < args.size(); a++) {
                resetNodes(args[a].getOperation());
            }
        }

        /***********************************************************************
         *                        Value generation
         **********************************************************************/
        CG<Base> createCG(const Argument<Base>& arg) {
            return CG<Base>(*this, arg);
        }

        CG<Base> createCG(OperationNode<Base>* node) {
            return CG<Base>(*this, node);
        }

        /***********************************************************************
         *                        Loop management
         **********************************************************************/

        const std::map<size_t, LoopModel<Base>*>& getLoops() const;

        LoopModel<Base>* getLoop(size_t loopId) const;

        size_t addLoopDependentIndexPattern(const IndexPattern& jacPattern);

        void manageLoopDependentIndexPattern(const IndexPattern* pattern);

        size_t addLoopIndependentIndexPattern(const IndexPattern& pattern, size_t hint);

        /***********************************************************************
         *                   Operation graph manipulation
         **********************************************************************/

        /**
         * Solves an expression (e.g. f(x, y) == 0) for a given variable (e.g. x)
         * The variable can appear only once in the expression.
         * 
         * @param expression  The original expression (f(x, y))
         * @param code  The variable to solve for
         * @return  The expression for variable
         */
        inline CG<Base> solveFor(OperationNode<Base>& expression,
                                 OperationNode<Base>& code) throw (CGException);

        inline CG<Base> solveFor(const std::vector<OperationPathNode<Base> >& path) throw (CGException);

        /**
         * Eliminates an independent variable by substitution using the provided
         * dependent variable which is assumed to be a residual of an equation.
         * If successful the model will contain one less independent variable.
         * 
         * @param indep The independent variable to eliminate.
         * @param dep The dependent variable representing a residual
         * @param removeFromIndeps Whether or not to immediatelly remove the
         *                         independent variable from the list of
         *                         independents in the model. The subtitution
         *                         operation can only be reversed if the 
         *                         variable is not removed.
         */
        inline void substituteIndependent(const CG<Base>& indep,
                                          const CG<Base>& dep,
                                          bool removeFromIndeps = true) throw (CGException);

        inline void substituteIndependent(OperationNode<Base>& indep,
                                          OperationNode<Base>& dep,
                                          bool removeFromIndeps = true) throw (CGException);

        /**
         * Reverts a subtitution of an independent variable that has not been 
         * removed from the list of independents yet.
         * Warning: it does not recover any custom name assigned to the variable.
         * 
         * @param indep The independent variable
         */
        inline void undoSubstituteIndependent(OperationNode<Base>& indep) throw (CGException);

        /**
         * Finallizes the subtitution of an independent variable by eliminating
         * it from the list of independents. After this operation the variable
         * subtitution cannot be undone.
         * 
         * @param indep The independent variable
         */
        inline void removeIndependent(OperationNode<Base>& indep) throw (CGException);

        /**
         * Adds an operation node to the list of nodes to be deleted when this
         * handler is destroyed.
         * 
         * @param code The operation node to be managed.
         * @return true if the node was successfuly added to the list or
         *         false if it had already been previously added.
         */
        inline bool manageOperationNodeMemory(OperationNode<Base>* code) {
            if (std::find(_codeBlocks.begin(), _codeBlocks.end(), code) == _codeBlocks.end()) {
                manageOperationNode(code);
                return false;
            }
            return true;
        }

        /**
         * Destructor
         */
        inline virtual ~CodeHandler() {
            reset();
        }

    protected:

        virtual void manageOperationNode(OperationNode<Base>* code) {
            //assert(std::find(_codeBlocks.begin(), _codeBlocks.end(), code) == _codeBlocks.end()); // <<< too great of an impact in performance
            if (_codeBlocks.capacity() == _codeBlocks.size()) {
                _codeBlocks.reserve((_codeBlocks.size()*3) / 2 + 1);
            }

            _codeBlocks.push_back(code);
        }

        virtual void markCodeBlockUsed(OperationNode<Base>& code) {
            code.total_use_count_++;

            CGOpCode op = code.getOperationType();

            if (op == CGAliasOp) {
                /**
                 * Alias operations are always followed so that there is a 
                 * correct usage count at the operation that it points to
                 */
                assert(code.getArguments().size() == 1);
                OperationNode<Base>* arg = code.getArguments()[0].getOperation();
                if (arg != NULL) {
                    markCodeBlockUsed(*arg);
                }

            } else if (code.getTotalUsageCount() == 1) {
                // first time this operation is visited

                size_t previousScope = _currentScopeColor;

                code.setColor(_currentScopeColor);

                // check if there is a scope change
                if (op == CGLoopStartOp || op == CGStartIfOp || op == CGElseIfOp || op == CGElseOp) {
                    // leaving a scope
                    ScopePath& sPath = _scopes[_currentScopeColor];
                    assert(sPath.back().beginning == NULL);
                    if (op == CGLoopStartOp || op == CGStartIfOp) {
                        sPath.back().beginning = &code; // save the initial node
                    } else {
                        assert(!code.getArguments().empty() &&
                               code.getArguments()[0].getOperation() != NULL &&
                               code.getArguments()[0].getOperation()->getOperationType() == CGStartIfOp);
                        sPath.back().beginning = code.getArguments()[0].getOperation(); // save the initial node
                    }
                    _currentScopeColor = sPath.size() > 1 ? sPath[sPath.size() - 2].color : 0;
                }

                if (op == CGLoopEndOp || op == CGEndIfOp || op == CGElseIfOp || op == CGElseOp) {
                    // entering a new scope
                    _currentScopeColor = ++_scopeColorCount;

                    _scopes.resize(_currentScopeColor + 1);
                    _scopes[_currentScopeColor] = _scopes[previousScope];

                    // change current scope
                    if (op == CGLoopEndOp || op == CGEndIfOp) {
                        // one more scope level
                        _scopes[_currentScopeColor].push_back(ScopePathElement<Base>(_currentScopeColor, &code));
                    } else {
                        // same level but different scope
                        _scopes[_currentScopeColor].back() = ScopePathElement<Base>(_currentScopeColor, &code);
                    }
                }

                /**
                 * loop arguments
                 */
                const std::vector<Argument<Base> >& args = code.arguments_;

                typename std::vector<Argument<Base> >::const_iterator it;
                for (it = args.begin(); it != args.end(); ++it) {
                    if (it->getOperation() != NULL) {
                        OperationNode<Base>& arg = *it->getOperation();
                        markCodeBlockUsed(arg);
                    }
                }

                if (op == CGIndexOp) {
                    const IndexOperationNode<Base>& inode = static_cast<const IndexOperationNode<Base>&> (code);
                    // indexes that don't depend on a loop start or an index assignment are declared elsewhere
                    if (inode.getArguments().size() > 0) {
                        _indexes.insert(&inode.getIndex());
                    }
                } else if (op == CGDependentRefRhsOp) {
                    assert(code.getInfo().size() == 1);
                    size_t depIndex = code.getInfo()[0];

                    assert(_dependents->size() > depIndex);
                    OperationNode<Base>* depNode = (*_dependents)[depIndex].getOperationNode();
                    assert(depNode != NULL && depNode->getOperationType() != CGInvOp);

                    code.setVariableID(depNode->getVariableID());
                }

                /**
                 * reset scope
                 */
                if (previousScope != _currentScopeColor) {
                    _currentScopeColor = previousScope;
                }

            } else {
                // been to this node before

                if (code.getColor() != _currentScopeColor) {
                    /**
                     * node previously used in a different scope
                     * must make sure it is defined before being used in both
                     * scopes
                     */
                    size_t depth;
                    std::set<size_t> divergence = findFirstDifferentScopeNodes(code.getColor(), _currentScopeColor, depth);
                    std::set<size_t>::const_iterator itScope;
                    for (itScope = divergence.begin(); itScope != divergence.end(); ++itScope) {
                        _scopeExtraDependecy[*itScope].insert(&code);
                    }

                    // update the scope where it should be defined
                    if (depth == 0)
                        code.setColor(0);
                    else
                        code.setColor(_scopes[_currentScopeColor][depth - 1].color);
                }
            }
        }

        inline std::set<size_t> findFirstDifferentScopeNodes(size_t color1, size_t color2, size_t& depth) {
            assert(color1 < _scopes.size());
            assert(color2 < _scopes.size());

            ScopePath& scopePath1 = _scopes[color1];
            ScopePath& scopePath2 = _scopes[color2];

            size_t s1 = scopePath1.size();
            size_t s2 = scopePath2.size();

            std::set<size_t> divergence;
            for (depth = 0; depth < s1 && depth < s2; depth++) {
                if (scopePath1[depth].color != scopePath2[depth].color) {
                    divergence.insert(scopePath1[depth].color);
                    divergence.insert(scopePath2[depth].color);
                    return divergence;
                }
            }

            if (s1 < s2) {
                divergence.insert(scopePath2[s1].color);
            } else {
                divergence.insert(scopePath1[s2].color);
            }

            return divergence;
        }

        inline static bool containsArgument(const OperationNode<Base>& node, const OperationNode<Base>& arg) {
            const std::vector<Argument<Base> >& args = node.getArguments();
            for (size_t a = 0; a < args.size(); a++) {
                if (args[a].getOperation() == &arg) {
                    return true;
                }
            }
            return false;
        }

        virtual void registerAtomicFunction(CGAbstractAtomicFun<Base>& atomic) {
            _atomicFunctions[atomic.getId()] = &atomic;
        }

        /***********************************************************************
         *                        Loop management
         **********************************************************************/
        virtual void registerLoop(LoopModel<Base>& loop);

        /***********************************************************************
         * 
         **********************************************************************/
        virtual void checkVariableCreation(OperationNode<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            if (code.getOperationType() == CGAliasOp) {
                /**
                 * avoid creating temporary variables for alias operations,
                 * the temporary should be the variable where it points to
                 */
                assert(args.size() == 1);
                if (args[0].getOperation() != NULL) {
                    checkVariableCreation(*args[0].getOperation());
                }
                return;
            }

            typename std::vector<Argument<Base> >::const_iterator it;

            for (it = args.begin(); it != args.end(); ++it) {
                if (it->getOperation() != NULL) {
                    OperationNode<Base>& arg = *it->getOperation();

                    if (arg.getUsageCount() == 0) {
                        // dependencies not visited yet
                        checkVariableCreation(arg);

                        CGOpCode type = arg.getOperationType();
                        if (type == CGLoopEndOp || type == CGElseIfOp || type == CGElseOp || type == CGEndIfOp) {
                            /**
                             * Some types of operations must be added immediatelly 
                             * after its arguments
                             * in order to avoid having other arguments inside
                             * that stack frame (or scope)
                             */
                            if (arg.getVariableID() == 0) {
                                addToEvaluationQueue(arg);
                                // ID value is not really used but must be non-zero
                                arg.setVariableID(std::numeric_limits<size_t>::max());
                            }
                        } else if (type == CGAtomicForwardOp || type == CGAtomicReverseOp) {
                            /**
                             * Save atomic function related information
                             */
                            assert(arg.getArguments().size() > 1);
                            assert(arg.getInfo().size() > 1);
                            size_t id = arg.getInfo()[0];
                            const std::string& atomicName = _atomicFunctions.at(id)->afun_name();
                            if (_atomicFunctionsSet.find(atomicName) == _atomicFunctionsSet.end()) {
                                _atomicFunctionsSet.insert(atomicName);
                                _atomicFunctionsOrder->push_back(atomicName);
                            }
                        }
                    }
                }
            }

            for (it = args.begin(); it != args.end(); ++it) {
                if (it->getOperation() != NULL) {
                    OperationNode<Base>& arg = *it->getOperation();
                    CGOpCode aType = arg.getOperationType();
                    /**
                     * make sure new temporary variables are NOT created for
                     * the independent variables and that a dependency did
                     * not use it first
                     */
                    if ((arg.getVariableID() == 0 || !isIndependent(arg)) && arg.getUsageCount() == 0) {
                        if (aType == CGLoopIndexedIndepOp) {
                            // ID value not really used but must be non-zero
                            arg.setVariableID(std::numeric_limits<size_t>::max());
                        } else if (aType == CGLoopEndOp || aType == CGElseIfOp ||
                                aType == CGElseOp || aType == CGEndIfOp) {
                            continue; // already added
                        } else if (aType == CGAliasOp) {
                            continue; // should never be added to the evaluation queue
                        } else if (aType == CGLoopStartOp) {
                            if (arg.getVariableID() == 0) {
                                addToEvaluationQueue(arg);
                                // ID value is not really used but must be non-zero
                                arg.setVariableID(std::numeric_limits<size_t>::max());
                            }
                        } else {
                            size_t argIndex = it - args.begin();
                            if (_lang->createsNewVariable(arg) ||
                                    _lang->requiresVariableArgument(code.getOperationType(), argIndex)) {

                                addToEvaluationQueue(arg);

                                if (arg.getVariableID() == 0) {
                                    if (aType == CGAtomicForwardOp || aType == CGAtomicReverseOp) {
                                        arg.setVariableID(_idAtomicCount);
                                        _idAtomicCount++;
                                    } else if (aType == CGLoopIndexedDepOp) {
                                        // ID value not really used but must be non-zero
                                        arg.setVariableID(std::numeric_limits<size_t>::max());
                                    } else if (aType == CGArrayCreationOp) {
                                        // a temporary array
                                        size_t arraySize = arg.getArguments().size();
                                        arg.setVariableID(_idArrayCount);
                                        _idArrayCount += arraySize;
                                    } else {
                                        // a single temporary variable
                                        arg.setVariableID(_idCount);
                                        _idCount++;
                                    }
                                }
                            }
                        }
                    }

                    arg.increaseUsageCount();
                }
            }

        }

        inline void addToEvaluationQueue(OperationNode<Base>& arg) {
            if (_variableOrder.size() == _variableOrder.capacity()) {
                _variableOrder.reserve((_variableOrder.size()*3) / 2 + 1);
            }

            _variableOrder.push_back(&arg);
            arg.setEvaluationOrder(_variableOrder.size());

            dependentAdded2EvaluationQueue(arg);
        }

        inline void reduceTemporaryVariables(vector<CG<Base> >& dependent) {

            /**
             * determine the last line where each temporary variable is used
             */
            resetUsageCount();

            for (size_t i = 0; i < dependent.size(); i++) {
                OperationNode<Base>* node = dependent[i].getOperationNode();
                if (node != NULL) {
                    if (node->use_count_ == 0) {
                        // dependencies not visited yet
                        determineLastTempVarUsage(*node);
                    }
                    node->increaseUsageCount();
                }
            }

            // where temporary variables can be released
            vector<std::vector<OperationNode<Base>* > > tempVarRelease(_variableOrder.size());
            for (size_t i = 0; i < _variableOrder.size(); i++) {
                OperationNode<Base>* var = _variableOrder[i];
                if (isTemporary(*var) || isTemporaryArray(*var)) {
                    size_t releaseLocation = var->getLastUsageEvaluationOrder() - 1;
                    tempVarRelease[releaseLocation].push_back(var);
                }
            }


            /**
             * Redefine temporary variable IDs
             */
            std::vector<size_t> freedVariables; // variable IDs no longer in use
            std::vector<const Argument<Base>*> tmpArrayValues(_idArrayCount, NULL); //likely values in temporary array
            std::map<size_t, size_t> freeArrayStartSpace; // [start] = end
            std::map<size_t, size_t> freeArrayEndSpace; // [end] = start
            _idCount = _minTemporaryVarID;
            _idArrayCount = 1;

            for (size_t i = 0; i < _variableOrder.size(); i++) {
                OperationNode<Base>& var = *_variableOrder[i];

                const std::vector<OperationNode<Base>* >& released = tempVarRelease[i];
                for (size_t r = 0; r < released.size(); r++) {
                    if (isTemporary(*released[r])) {
                        freedVariables.push_back(released[r]->getVariableID());
                    } else if (isTemporaryArray(*released[r])) {
                        addFreeArraySpace(*released[r], freeArrayStartSpace, freeArrayEndSpace);
                        assert(freeArrayStartSpace.size() == freeArrayEndSpace.size());
                    }
                }

                if (isTemporary(var)) {
                    // a single temporary variable
                    if (freedVariables.empty()) {
                        var.setVariableID(_idCount);
                        _idCount++;
                    } else {
                        size_t id = freedVariables.back();
                        freedVariables.pop_back();
                        var.setVariableID(id);
                    }
                } else if (isTemporaryArray(var)) {
                    // a temporary array
                    size_t arrayStart = reserveArraySpace(var, freeArrayStartSpace, freeArrayEndSpace, tmpArrayValues);
                    assert(freeArrayStartSpace.size() == freeArrayEndSpace.size());
                    var.setVariableID(arrayStart + 1);
                }

            }
        }

        inline static void addFreeArraySpace(const OperationNode<Base>& released,
                                             std::map<size_t, size_t>& freeArrayStartSpace,
                                             std::map<size_t, size_t>& freeArrayEndSpace) {
            size_t arrayStart = released.getVariableID() - 1;
            const size_t arraySize = released.getArguments().size();
            size_t arrayEnd = arrayStart + arraySize - 1;

            std::map<size_t, size_t>::iterator it;
            if (arrayStart > 0) {
                it = freeArrayEndSpace.find(arrayStart - 1); // previous
                if (it != freeArrayEndSpace.end()) {
                    arrayStart = it->second; // merge space
                    freeArrayEndSpace.erase(it);
                    freeArrayStartSpace.erase(arrayStart);
                }
            }
            it = freeArrayStartSpace.find(arrayEnd + 1); // next
            if (it != freeArrayStartSpace.end()) {
                arrayEnd = it->second; // merge space 
                freeArrayStartSpace.erase(it);
                freeArrayEndSpace.erase(arrayEnd);
            }

            freeArrayStartSpace[arrayStart] = arrayEnd;
            freeArrayEndSpace[arrayEnd] = arrayStart;
        }

        inline size_t reserveArraySpace(const OperationNode<Base>& newArray,
                                        std::map<size_t, size_t>& freeArrayStartSpace,
                                        std::map<size_t, size_t>& freeArrayEndSpace,
                                        std::vector<const Argument<Base>*>& tmpArrayValues) {
            size_t arraySize = newArray.getArguments().size();

            std::set<size_t> blackList;
            const std::vector<Argument<Base> >& args = newArray.getArguments();
            for (size_t i = 0; i < args.size(); i++) {
                const OperationNode<Base>* argOp = args[i].getOperation();
                if (argOp != NULL && argOp->getOperationType() == CGArrayElementOp) {
                    const OperationNode<Base>& otherArray = *argOp->getArguments()[0].getOperation();
                    assert(otherArray.getVariableID() > 0); // make sure it had already been assigned space
                    size_t otherArrayStart = otherArray.getVariableID() - 1;
                    size_t index = argOp->getInfo()[0];
                    blackList.insert(otherArrayStart + index);
                }
            }

            /**
             * Find the best location for the new array
             */
            std::map<size_t, size_t>::reverse_iterator it;
            std::map<size_t, size_t>::reverse_iterator itBestFit = freeArrayStartSpace.rend();
            size_t bestCommonValues = 0; // the number of values likely to be the same
            for (it = freeArrayStartSpace.rbegin(); it != freeArrayStartSpace.rend(); ++it) {
                size_t start = it->first;
                size_t end = it->second;
                size_t space = end - start + 1;
                if (space < arraySize) {
                    continue;
                }

                std::set<size_t>::const_iterator itBlack = blackList.lower_bound(start);
                if (itBlack != blackList.end() && *itBlack <= end) {
                    continue; // cannot use this space
                }

                //possible candidate
                if (itBestFit == freeArrayStartSpace.rend()) {
                    itBestFit = it;
                } else {
                    size_t bestSpace = itBestFit->second - itBestFit->first + 1;

                    size_t commonVals = 0;
                    for (size_t i = 0; i < arraySize; i++) {
                        if (isSameArrayElement(tmpArrayValues[start + i], args[i])) {
                            commonVals++;
                        }
                    }

                    if (space < bestSpace || commonVals > bestCommonValues) {
                        // better fit
                        itBestFit = it;
                        bestCommonValues = commonVals;
                        if (bestCommonValues == arraySize) {
                            break; // jackpot
                        }
                    }
                }
            }

            size_t bestStart = std::numeric_limits<size_t>::max();
            if (itBestFit != freeArrayStartSpace.rend()) {
                /**
                 * Use available space
                 */
                bestStart = itBestFit->first;
                size_t bestEnd = itBestFit->second;
                size_t bestSpace = bestEnd - bestStart + 1;
                freeArrayStartSpace.erase(bestStart);
                if (bestSpace == arraySize) {
                    // entire space 
                    freeArrayEndSpace.erase(bestEnd);
                } else {
                    // some space left
                    size_t newFreeStart = bestStart + arraySize;
                    freeArrayStartSpace[newFreeStart] = bestEnd;
                    freeArrayEndSpace.at(bestEnd) = newFreeStart;
                }

            } else {
                /**
                 * no space available, need more
                 */
                // check if there is some free space at the end
                std::map<size_t, size_t>::iterator itEnd;
                itEnd = freeArrayEndSpace.find(_idArrayCount - 1);
                if (itEnd != freeArrayEndSpace.end()) {
                    // check if it can be used
                    size_t lastSpotStart = itEnd->second;
                    size_t lastSpotEnd = itEnd->first;
                    size_t lastSpotSize = lastSpotEnd - lastSpotStart + 1;
                    std::set<size_t>::const_iterator itBlack = blackList.lower_bound(lastSpotStart);
                    if (itBlack == blackList.end()) {
                        // can use this space
                        size_t newEnd = lastSpotStart + arraySize - 1;

                        freeArrayEndSpace.erase(itEnd);
                        freeArrayEndSpace[newEnd] = lastSpotStart;
                        freeArrayStartSpace[lastSpotStart] = newEnd;

                        _idArrayCount += arraySize - lastSpotSize;
                        bestStart = lastSpotStart;
                    }
                }

                if (bestStart == std::numeric_limits<size_t>::max()) {
                    // brand new space
                    size_t id = _idArrayCount;
                    _idArrayCount += arraySize;
                    bestStart = id - 1;
                }
            }

            for (size_t i = 0; i < arraySize; i++) {
                tmpArrayValues[bestStart + i] = &args[i];
            }

            return bestStart;
        }

        inline static bool isSameArrayElement(const Argument<Base>* oldArg, const Argument<Base>& arg) {
            if (oldArg != NULL) {
                if (oldArg->getParameter() != NULL) {
                    if (arg.getParameter() != NULL) {
                        return (*arg.getParameter() == *oldArg->getParameter());
                    }
                } else {
                    return (arg.getOperation() == oldArg->getOperation());
                }
            }
            return false;
        }

        /**
         * Determines when each temporary variable is last used in the
         * evaluation order
         * 
         * @param code The current node to determine the number of usages
         */
        inline void determineLastTempVarUsage(OperationNode<Base>& code) {
            if (code.getOperationType() == CGLoopEndOp) {
                LoopEndOperationNode<Base>& loopEnd = static_cast<LoopEndOperationNode<Base>&> (code);
                _loopDepth++;
                _loopOuterVars.resize(_loopDepth + 1);
                _loopStartEvalOrder.push_back(loopEnd.getLoopStart().getEvaluationOrder());

            } else if (code.getOperationType() == CGLoopStartOp) {
                _loopDepth--; // leaving the currrent loop
            }

            /**
             * count variable usage
             */
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;
            for (it = args.begin(); it != args.end(); ++it) {
                if (it->getOperation() != NULL) {
                    OperationNode<Base>& arg = *it->getOperation();

                    if (arg.use_count_ == 0) {
                        // dependencies not visited yet
                        determineLastTempVarUsage(arg);
                    }

                    arg.increaseUsageCount();

                    size_t order = code.getEvaluationOrder();
                    OperationNode<Base>* aa = getOperationFromAlias(arg); // follow alias!
                    if (aa != NULL) {
                        if (aa->getLastUsageEvaluationOrder() < order) {
                            aa->setLastUsageEvaluationOrder(order);
                        }

                        if (_loopDepth >= 0 &&
                                aa->getEvaluationOrder() < _loopStartEvalOrder[_loopDepth] &&
                                isTemporary(*aa)) {
                            // outer variable used inside the loop
                            _loopOuterVars[_loopDepth].insert(aa);
                        }
                    }
                }
            }

            if (code.getOperationType() == CGLoopEndOp) {
                /**
                 * temporary variables from outside the loop which are used
                 * whithin the loop cannot be overwritten inside that loop
                 */
                const std::set<OperationNode<Base>*>& outerLoopUsages = _loopOuterVars.back();
                typename std::set<OperationNode<Base>*>::const_iterator it;
                for (it = outerLoopUsages.begin(); it != outerLoopUsages.end(); ++it) {
                    OperationNode<Base>* outerVar = *it;
                    size_t order = code.getEvaluationOrder();

                    OperationNode<Base>* aa = getOperationFromAlias(*outerVar); // follow alias!
                    if (aa != NULL && aa->getLastUsageEvaluationOrder() < order)
                        aa->setLastUsageEvaluationOrder(order);
                }

                _loopDepth--;
                _loopOuterVars.pop_back();
                _loopStartEvalOrder.pop_back();

            } else if (code.getOperationType() == CGLoopStartOp) {
                _loopDepth++; // comming back to the loop
            }

        }

        inline void resetUsageCount() {
            typename std::vector<OperationNode<Base> *>::const_iterator it;
            for (it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                OperationNode<Base>* block = *it;
                block->use_count_ = 0;
            }
        }

        /**
         * Defines the evaluation order for the code fragments that do not
         * create variables
         * @param code The operation just added to the evaluation order
         */
        inline void dependentAdded2EvaluationQueue(OperationNode<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;

            for (it = args.begin(); it != args.end(); ++it) {
                if (it->getOperation() != NULL) {
                    OperationNode<Base>& arg = *it->getOperation();
                    if (arg.getEvaluationOrder() == 0) {
                        arg.setEvaluationOrder(code.getEvaluationOrder());
                        dependentAdded2EvaluationQueue(arg);
                    }
                }
            }
        }

        inline bool isIndependent(const OperationNode<Base>& arg) const {
            if (arg.getOperationType() == CGArrayCreationOp ||
                    arg.getOperationType() == CGAtomicForwardOp ||
                    arg.getOperationType() == CGAtomicReverseOp)
                return false;

            size_t id = arg.getVariableID();
            return id > 0 && id <= _independentVariables.size();
        }

        inline bool isTemporary(const OperationNode<Base>& arg) const {
            CGOpCode op = arg.getOperationType();
            return op != CGArrayCreationOp &&
                    op != CGAtomicForwardOp &&
                    op != CGAtomicReverseOp &&
                    op != CGLoopStartOp &&
                    op != CGLoopEndOp &&
                    op != CGLoopIndexedDepOp &&
                    op != CGLoopIndexedIndepOp &&
                    op != CGIndexOp &&
                    op != CGIndexAssignOp &&
                    arg.getVariableID() >= _minTemporaryVarID;
        }

        inline bool isTemporaryArray(const OperationNode<Base>& arg) const {
            return arg.getOperationType() == CGArrayCreationOp;
        }

        inline static OperationNode<Base>* getOperationFromAlias(OperationNode<Base>& alias) {
            if (alias.getOperationType() != CGAliasOp) {
                return &alias;
            } else {
                OperationNode<Base>* aa = &alias;
                do {
                    assert(aa->getArguments().size() == 1);
                    aa = aa->getArguments()[0].getOperation();
                } while (aa != NULL && aa->getOperationType() == CGAliasOp);
                return aa;
            }
        }

        inline void resetManagedNodes() {
            _variableOrder.clear();

            for (typename std::vector<OperationNode<Base> *>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                OperationNode<Base>* block = *it;
                block->resetHandlerCounters();
                block->setColor(0);
            }
        }

        /***********************************************************************
         *                   Graph management functions
         **********************************************************************/

        inline void findPaths(SourceCodePath& path2node,
                              OperationNode<Base>& code,
                              std::vector<SourceCodePath>& found,
                              size_t max);

        static inline std::vector<SourceCodePath> findPathsFromNode(const std::vector<SourceCodePath> nodePaths,
                                                                    OperationNode<Base>& node);

    private:

        CodeHandler(const CodeHandler&); // not implemented

        CodeHandler& operator=(const CodeHandler&); // not implemented

        friend class CG<Base>;
        friend class CGAbstractAtomicFun<Base>;
        friend class BaseAbstractAtomicFun<Base>;
        friend class LoopModel<Base>;

    };

}

#endif
