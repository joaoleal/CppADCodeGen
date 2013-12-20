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
        /**
         * the order for the variable creation in the source code 
         */
        std::vector<OperationNode<Base> *> _variableOrder;
        /**
         * the order for the variable creation in the source code 
         * (each level represents a different variable scope)
         */
        std::vector<std::vector<OperationNode<Base> *> > _scopedVariableOrder;
        // maps the ids of the atomic functions
        std::map<size_t, CGAbstractAtomicFun<Base>*> _atomicFunctions;
        // maps the loop ids of the loop atomic functions
        std::map<size_t, LoopModel<Base>*> _loops;
        // the used indexes
        std::set<const IndexDclrOperationNode<Base>*> _indexes;
        // the used random index patterns
        std::set<RandomIndexPattern*> _indexRandomPatterns;
        //
        std::vector<IndexPattern*> _loopDependentIndexPatterns;
        std::vector<const IndexPattern*> _loopDependentIndexPatternManaged; // garbage collection
        std::vector<IndexPattern*> _loopIndependentIndexPatterns;
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
        // possible altered nodes due to scope conditionals (altered node <-> clone of original)
        std::list<std::pair<OperationNode<Base>*, OperationNode<Base>* > > _alteredNodes;
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
        //
        JobTimer* _jobTimer;
    public:

        CodeHandler(size_t varCount = 50) :
            _idCount(1),
            _idArrayCount(1),
            _idAtomicCount(1),
            _dependents(NULL),
            _scopedVariableOrder(1),
            _atomicFunctionsOrder(NULL),
            _used(false),
            _reuseIDs(true),
            _loopDepth(-1),
            _scopeColorCount(0),
            _currentScopeColor(0),
            _lang(NULL),
            _minTemporaryVarID(0),
            _zeroDependents(false),
            _verbose(false),
            _jobTimer(NULL) {
            _codeBlocks.reserve(varCount);
            //_variableOrder.reserve(1 + varCount / 3);
            _scopedVariableOrder[0].reserve(1 + varCount / 3);

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
            CPPADCG_ASSERT_UNKNOWN(var.getOperationType() == CGInvOp);

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

        inline JobTimer* getJobTimer() const {
            return _jobTimer;
        }

        inline void setJobTimer(JobTimer* jobTimer) {
            _jobTimer = jobTimer;
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
            double beginTime = 0.0;

            if (_jobTimer != NULL) {
                _jobTimer->startingJob("source for '" + jobName + "'");
            } else if (_verbose) {
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
            _indexRandomPatterns.clear();
            _loopOuterVars.clear();
            _loopDepth = -1;
            _scopeColorCount = 0;
            _currentScopeColor = 0;
            _scopes.reserve(4);
            _scopes.resize(1);
            _alteredNodes.clear();
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
             * determine the variable creation order
             */
            _scopedVariableOrder.reserve(std::max(size_t(1), _scopes.size()) + 10); // some additional scopes might still be added

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

            /**
             * Generate flat variable order (without scopes)
             */
            if (_scopedVariableOrder.size() == 1) {
                _variableOrder.swap(_scopedVariableOrder[0]); // most common situation
            } else {
                optimizeIfs(); // reduce the number of adjoining ifs

                size_t vosize = 0;
                for (size_t s = 0; s < _scopedVariableOrder.size(); s++) {
                    vosize += _scopedVariableOrder[s].size();
                }
                _variableOrder.resize(vosize);

                size_t e = 0;
                addScopeToVarOrder(0, e);

                // if e > vosize then some nodes (marking the beginning of scopes)
                // must have been added more than once
                CPPADCG_ASSERT_UNKNOWN(_variableOrder.size() == e);
            }

            for (size_t p = 0; p < _variableOrder.size(); p++) {
                OperationNode<Base>& arg = *_variableOrder[p];
                arg.setEvaluationOrder(p + 1);
                dependentAdded2EvaluationQueue(arg);
            }

            /**
             * Reuse temporary variables
             */
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
                                              _indexes, _indexRandomPatterns,
                                              _loopDependentIndexPatterns, _loopIndependentIndexPatterns,
                                              _zeroDependents);
            lang.generateSourceCode(out, info);

            /**
             * clean-up
             */
            _atomicFunctionsSet.clear();

            // restore altered nodes
            typename std::list<std::pair<OperationNode<Base>*, OperationNode<Base>* > >::const_iterator itAlt;
            for (itAlt = _alteredNodes.begin(); itAlt != _alteredNodes.end(); ++itAlt) {
                OperationNode<Base>* tmp = itAlt->first;
                OperationNode<Base>* opClone = itAlt->second;
                if (tmp->getOperationType() == CGTmpOp && !tmp->getInfo().empty()) { // some might have already been restored
                    tmp->setOperation(opClone->getOperationType(), opClone->getArguments());
                    tmp->getInfo() = opClone->getInfo();
                }
            }
            _alteredNodes.clear();

            if (_jobTimer != NULL) {
                _jobTimer->finishedJob();
            } else if (_verbose) {
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
            _indexRandomPatterns.clear();
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

        size_t addLoopDependentIndexPattern(IndexPattern& jacPattern);

        void manageLoopDependentIndexPattern(const IndexPattern* pattern);

        size_t addLoopIndependentIndexPattern(IndexPattern& pattern, size_t hint);

        /***********************************************************************
         *                           Index patterns
         **********************************************************************/
        static inline void findRandomIndexPatterns(IndexPattern* ip,
                                                   std::set<RandomIndexPattern*>& found) {
            if (ip == NULL)
                return;

            if (ip->getType() == RANDOM1D || ip->getType() == RANDOM2D) {
                found.insert(static_cast<RandomIndexPattern*> (ip));
            } else {
                std::set<IndexPattern*> indexes;
                ip->getSubIndexes(indexes);
                for (std::set<IndexPattern*>::const_iterator itIp = indexes.begin(); itIp != indexes.end(); ++itIp) {
                    IndexPattern* sip = *itIp;
                    if (sip->getType() == RANDOM1D || sip->getType() == RANDOM2D)
                        found.insert(static_cast<RandomIndexPattern*> (sip));
                }
            }
        }

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
         * @param removeFromIndeps Whether or not to immediately remove the
         *                         independent variable from the list of
         *                         independents in the model. The substitution
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
         * Reverts a substitution of an independent variable that has not been 
         * removed from the list of independents yet.
         * Warning: it does not recover any custom name assigned to the variable.
         * 
         * @param indep The independent variable
         */
        inline void undoSubstituteIndependent(OperationNode<Base>& indep) throw (CGException);

        /**
         * Finalizes the substitution of an independent variable by eliminating
         * it from the list of independents. After this operation the variable
         * substitution cannot be undone.
         * 
         * @param indep The independent variable
         */
        inline void removeIndependent(OperationNode<Base>& indep) throw (CGException);

        /**
         * Adds an operation node to the list of nodes to be deleted when this
         * handler is destroyed.
         * 
         * @param code The operation node to be managed.
         * @return true if the node was successfully added to the list or
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
            //CPPADCG_ASSERT_UNKNOWN(std::find(_codeBlocks.begin(), _codeBlocks.end(), code) == _codeBlocks.end()); // <<< too great of an impact in performance
            if (_codeBlocks.capacity() == _codeBlocks.size()) {
                _codeBlocks.reserve((_codeBlocks.size()*3) / 2 + 1);
            }

            _codeBlocks.push_back(code);
        }

        virtual void markCodeBlockUsed(OperationNode<Base>& code) {
            code.increaseTotalUsageCount();

            CGOpCode op = code.getOperationType();
            if (isIndependent(code)) {
                return; // nothing to do
            } else if (op == CGAliasOp) {
                /**
                 * Alias operations are always followed so that there is a 
                 * correct usage count at the operation that it points to
                 */
                CPPADCG_ASSERT_UNKNOWN(code.getArguments().size() == 1);
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
                    CPPADCG_ASSERT_UNKNOWN(sPath.back().beginning == NULL);
                    if (op == CGLoopStartOp || op == CGStartIfOp) {
                        sPath.back().beginning = &code; // save the initial node
                    } else {
                        CPPADCG_ASSERT_UNKNOWN(!code.getArguments().empty() &&
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
                    if (inode.isDefinedLocally()) {
                        _indexes.insert(&inode.getIndex());
                    }
                } else if (op == CGLoopIndexedIndepOp || op == CGLoopIndexedDepOp || op == CGIndexAssignOp) {
                    IndexPattern* ip;
                    if (op == CGLoopIndexedDepOp) {
                        size_t pos = code.getInfo()[0];
                        ip = _loopDependentIndexPatterns[pos];
                    } else if (op == CGLoopIndexedIndepOp) {
                        size_t pos = code.getInfo()[1];
                        ip = _loopIndependentIndexPatterns[pos];
                    } else {
                        ip = &static_cast<IndexAssignOperationNode<Base>&> (code).getIndexPattern();
                    }

                    findRandomIndexPatterns(ip, _indexRandomPatterns);

                } else if (op == CGDependentRefRhsOp) {
                    CPPADCG_ASSERT_UNKNOWN(code.getInfo().size() == 1);
                    size_t depIndex = code.getInfo()[0];

                    CPPADCG_ASSERT_UNKNOWN(_dependents->size() > depIndex);
                    OperationNode<Base>* depNode = (*_dependents)[depIndex].getOperationNode();
                    CPPADCG_ASSERT_UNKNOWN(depNode != NULL && depNode->getOperationType() != CGInvOp);

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

                if (op == CGTmpOp && !code.getInfo().empty()) {
                    /**
                     * this node was previously altered to ensure that the 
                     * evaluation of the expression is only performed for the 
                     * required iterations
                     */
                    if (code.getColor() == _currentScopeColor) {
                        // outside an if (defined for all iterations)
                        restoreTemporaryVar(code);
                    } else {
                        updateTemporaryVarInDiffScopes(code);
                    }

                } else if (code.getColor() != _currentScopeColor && op != CGLoopIndexedIndepOp) {
                    size_t oldScope = code.getColor();
                    /**
                     * node previously used in a different scope
                     * must make sure it is defined before being used in both
                     * scopes
                     */
                    size_t depth = findFirstDifferentScope(oldScope, _currentScopeColor);

                    // update the scope where it should be defined
                    size_t newScope;
                    if (depth == 0)
                        newScope = 0;
                    else
                        newScope = _scopes[_currentScopeColor][depth - 1].color;

                    if (oldScope != newScope) {
                        /**
                         * does this variable require a condition based on indexes?
                         */
                        bool addedIf = handleTemporaryVarInDiffScopes(code, oldScope, newScope);

                        if (!addedIf) {
                            code.setColor(newScope);

                            /**
                             * Must also update the scope of the arguments used by this operation
                             */
                            const std::vector<Argument<Base> >& args = code.getArguments();
                            size_t aSize = args.size();
                            for (size_t a = 0; a < aSize; a++) {
                                updateVarScopeUsage(args[a].getOperation(), newScope, oldScope);
                            }
                        }
                    }
                }
            }
        }

        inline bool handleTemporaryVarInDiffScopes(OperationNode<Base>& code,
                                                   size_t oldScope, size_t newScope) {
            if (_currentScopeColor == 0)
                return false;

            /**
             * @TODO allow Array elements to use a CGTmp instead of a CGArrayCreationOp
             */
            CPPADCG_ASSERT_KNOWN(code.getOperationType() != CGArrayCreationOp, "Not supported yet");

            /**
             * does this variable require a condition based on indexes?
             */
            std::vector<size_t> iterationRegions;
            OperationNode<Base>* bScopeNewEnd = _scopes[_currentScopeColor].back().end;
            OperationNode<Base>* bScopeOldEnd = _scopes[oldScope].back().end;

            CGOpCode bNewOp = bScopeNewEnd->getOperationType();
            CGOpCode bOldOp = bScopeOldEnd->getOperationType();

            if ((bNewOp == CGEndIfOp || bNewOp == CGElseOp || bNewOp == CGElseIfOp) &&
                    (bOldOp == CGEndIfOp || bOldOp == CGElseOp || bOldOp == CGElseIfOp)) {
                // used in 2 different if/else branches

                /**
                 * determine the iterations which use this temporary variable
                 */
                OperationNode<Base>* bScopeNew = bScopeNewEnd->getArguments()[0].getOperation();
                OperationNode<Base>* bScopeOld = bScopeOldEnd->getArguments()[0].getOperation();

                IndexOperationNode<Base>* newIterIndexOp = NULL;
                iterationRegions = ifBranchIterationRanges(bScopeNew, newIterIndexOp);
                CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);

                IndexOperationNode<Base>* oldIterIndexOp = NULL;
                std::vector<size_t> oldIterRegions = ifBranchIterationRanges(bScopeOld, oldIterIndexOp);
                combineOverlapingIterationRanges(iterationRegions, oldIterRegions);
                CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);
                CPPADCG_ASSERT_UNKNOWN(newIterIndexOp != NULL && newIterIndexOp == oldIterIndexOp);

                if (iterationRegions.size() > 2 ||
                        iterationRegions[0] != 0 ||
                        iterationRegions[1] != std::numeric_limits<size_t>::max()) {
                    // this temporary variable is not used by all iterations

                    replaceWithConditionalTempVar(code, *newIterIndexOp, iterationRegions, oldScope, newScope);
                    return true;
                }
            }

            return false;
        }

        inline void replaceWithConditionalTempVar(OperationNode<Base>& tmp,
                                                  IndexOperationNode<Base>& iterationIndexOp,
                                                  const std::vector<size_t>& iterationRegions,
                                                  size_t oldScope,
                                                  size_t commonScopeColor) {
            /**
             * clone
             */
            OperationNode<Base>* opClone = new OperationNode<Base>(tmp);
            manageOperationNode(opClone);

            /**
             * Create condition
             */
            OperationNode<Base>* tmpDclVar = new OperationNode<Base>(CGTmpDclOp);
            manageOperationNode(tmpDclVar);
            Argument<Base> tmpArg(*tmpDclVar);

            std::vector<Argument<Base> > args(1);
            args[0] = Argument<Base>(iterationIndexOp);
            OperationNode<Base>* cond = new OperationNode<Base>(CGIndexCondExprOp, iterationRegions, args);
            manageOperationNode(cond);

            // if
            OperationNode<Base>* ifStart = new OperationNode<Base>(CGStartIfOp, Argument<Base>(*cond));
            manageOperationNode(ifStart);

            OperationNode<Base>* tmpAssign = new OperationNode<Base>(CGLoopIndexedTmpOp, tmpArg, Argument<Base>(*opClone));
            manageOperationNode(tmpAssign);
            OperationNode<Base>* ifAssign = new OperationNode<Base>(CGCondResultOp, Argument<Base>(*ifStart), Argument<Base>(*tmpAssign));
            manageOperationNode(ifAssign);

            // end if
            OperationNode<Base>* endIf = new OperationNode<Base>(CGEndIfOp, Argument<Base>(*ifStart), Argument<Base>(*ifAssign));
            manageOperationNode(endIf);

            /**
             * Change original variable
             */
            std::vector<Argument<Base> > arguments(2);
            arguments[0] = tmpArg;
            arguments[1] = Argument<Base>(*endIf);
            tmp.setOperation(CGTmpOp, arguments);
            tmp.getInfo().resize(1); // used to mark that this node was altered here

            /**
             * add the new scope
             */
            size_t newScopeColor = ++_scopeColorCount;
            _scopes.resize(newScopeColor + 1);
            _scopes[newScopeColor] = _scopes[commonScopeColor];

            // one more scope level
            _scopes[newScopeColor].push_back(ScopePathElement<Base>(newScopeColor, endIf, ifStart));

            // apply scope colors
            tmpDclVar->setColor(commonScopeColor);
            ifStart->setColor(newScopeColor);
            cond->setColor(newScopeColor);
            opClone->setColor(newScopeColor);
            ifAssign->setColor(newScopeColor);
            tmpAssign->setColor(newScopeColor);
            endIf->setColor(commonScopeColor);
            tmp.setColor(commonScopeColor);

            // total usage count
            tmpDclVar->setTotalUsageCount(1);
            ifStart->setTotalUsageCount(1);
            cond->setTotalUsageCount(1);
            opClone->setTotalUsageCount(1);
            ifAssign->setTotalUsageCount(1);
            tmpAssign->setTotalUsageCount(1);
            endIf-> setTotalUsageCount(1);

            /**
             * Must also update the scope of the arguments used by this operation
             */
            const std::vector<Argument<Base> >& cargs = opClone->getArguments();
            size_t aSize = cargs.size();
            for (size_t a = 0; a < aSize; a++) {
                updateVarScopeUsage(cargs[a].getOperation(), newScopeColor, oldScope);
            }

            _alteredNodes.push_back(std::make_pair(&tmp, opClone));
        }

        inline void updateTemporaryVarInDiffScopes(OperationNode<Base>& code) {
            if (code.getColor() != _currentScopeColor) {
                return; //nothing to change
            }

            /**
             * does this variable require a condition based on indexes?
             */
            if (_currentScopeColor == 0)
                restoreTemporaryVar(code);

            /**
             * Determine if it should be moved into a different scope
             * so that it is defined before being used in both
             * scopes
             */
            size_t oldScope = code.getColor();

            size_t depth = findFirstDifferentScope(oldScope, _currentScopeColor);

            // update the scope where it should be defined
            size_t newScope = depth == 0 ? 0 : _scopes[_currentScopeColor][depth - 1].color;

            /**
             * does this variable require a condition based on indexes?
             */
            std::vector<size_t> iterationRegions;
            OperationNode<Base>* bScopeNewEnd = _scopes[_currentScopeColor].back().end;
            OperationNode<Base>* endif = code.getArguments()[0].getOperation();
            CPPADCG_ASSERT_UNKNOWN(endif->getOperationType() == CGEndIfOp);
            OperationNode<Base>* bScopeOldEnd = _scopes[endif->getColor()].back().end;

            CGOpCode bNewOp = bScopeNewEnd->getOperationType();

            if (bNewOp == CGEndIfOp || bNewOp == CGElseOp || bNewOp == CGElseIfOp) {
                // used in 2 different if/else branches

                /**
                 * determine the iterations which use this temporary variable
                 */
                OperationNode<Base>* bScopeNew = bScopeNewEnd->getArguments()[0].getOperation();
                OperationNode<Base>* bScopeOld = bScopeOldEnd->getArguments()[0].getOperation();

                IndexOperationNode<Base>* newIterIndexOp = NULL;
                iterationRegions = ifBranchIterationRanges(bScopeNew, newIterIndexOp);
                CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);

                IndexOperationNode<Base>* oldIterIndexOp = NULL;
                const std::vector<size_t> oldIterRegions = ifBranchIterationRanges(bScopeOld, oldIterIndexOp);
                combineOverlapingIterationRanges(iterationRegions, oldIterRegions);
                CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);
                CPPADCG_ASSERT_UNKNOWN(newIterIndexOp != NULL && newIterIndexOp == oldIterIndexOp);

                if (iterationRegions.size() == 2 &&
                        (iterationRegions[0] == 0 ||
                        iterationRegions[1] == std::numeric_limits<size_t>::max())) {
                    // this temporary variable is used by all iterations
                    // there is no need for an 'if'
                    restoreTemporaryVar(code);

                } else if (oldIterRegions != iterationRegions) {
                    OperationNode<Base>* cond = bScopeOld->getArguments()[0].getOperation();
                    CPPADCG_ASSERT_UNKNOWN(cond->getOperationType() == CGIndexCondExprOp);
                    cond->getInfo() = iterationRegions;
                }

            }

            if (oldScope != newScope) {
                code.setColor(newScope);
                /**
                 * Must also update the scope of the arguments used by this operation
                 */
                const std::vector<Argument<Base> >& cargs = code.getArguments();
                size_t aSize = cargs.size();
                for (size_t a = 0; a < aSize; a++) {
                    updateVarScopeUsage(cargs[a].getOperation(), newScope, oldScope);
                }
            }

        }

        inline void restoreTemporaryVar(OperationNode<Base>& tmp) {
            CPPADCG_ASSERT_UNKNOWN(tmp.getOperationType() == CGTmpOp && !tmp.getInfo().empty());

            OperationNode<Base>* endIf = tmp.getArguments()[1].getOperation();
            OperationNode<Base>* ifAssign = endIf->getArguments()[1].getOperation();
            OperationNode<Base>* tmpAssign = ifAssign->getArguments()[1].getOperation();
            OperationNode<Base>* opClone = tmpAssign->getArguments()[1].getOperation();
            tmp.setOperation(opClone->getOperationType(), opClone->getArguments());
            tmp.getInfo() = opClone->getInfo();

            tmp.setColor(_currentScopeColor);

            /**
             * Must also update the scope of the arguments used by this operation
             */
            const std::vector<Argument<Base> >& args = tmp.getArguments();
            size_t aSize = args.size();
            for (size_t a = 0; a < aSize; a++) {
                updateVarScopeUsage(args[a].getOperation(), _currentScopeColor, opClone->getColor());
            }
        }

        inline void restoreTemporaryVar(OperationNode<Base>* tmp,
                                        OperationNode<Base>* opClone) {
            CPPADCG_ASSERT_UNKNOWN(tmp.getOperationType() == CGTmpOp && !tmp.getInfo().empty());

            tmp.setOperation(opClone->getOperationType(), opClone->getArguments());
            tmp.getInfo() = opClone->getInfo();

            tmp.setColor(_currentScopeColor);

            /**
             * Must also update the scope of the arguments used by this operation
             */
            const std::vector<Argument<Base> >& args = tmp.getArguments();
            size_t aSize = args.size();
            for (size_t a = 0; a < aSize; a++) {
                updateVarScopeUsage(args[a].getOperation(), _currentScopeColor, opClone->getColor());
            }
        }

        inline void updateVarScopeUsage(OperationNode<Base>* node,
                                        size_t usageScope,
                                        size_t oldUsageScope) {
            if (node == NULL || node->getColor() == usageScope)
                return;


            size_t oldScope = node->getColor();
            size_t newScope;

            if (oldScope == oldUsageScope) {
                newScope = usageScope;
            } else {
                size_t depth = findFirstDifferentScope(oldScope, usageScope);

                newScope = (depth == 0) ? 0 : _scopes[usageScope][depth - 1].color;
            }

            if (newScope == oldScope)
                return;

            node->setColor(newScope);

            const std::vector<Argument<Base> >& args = node->getArguments();
            size_t aSize = args.size();
            for (size_t a = 0; a < aSize; a++) {
                updateVarScopeUsage(args[a].getOperation(), newScope, oldScope);
            }
        }

        inline void addScopeToVarOrder(size_t scope, size_t& e) {
            std::vector<OperationNode<Base> *>& vorder = _scopedVariableOrder[scope];

            const size_t vsize = vorder.size();
            for (size_t p = 0; p < vsize; p++) {
                OperationNode<Base>* node = vorder[p];
                CGOpCode op = node->getOperationType();

                if (op == CGLoopEndOp || op == CGEndIfOp || op == CGElseIfOp || op == CGElseOp) {
                    CPPADCG_ASSERT_UNKNOWN(!node->getArguments().empty());

                    OperationNode<Base>* beginScopeNode = node->getArguments()[0].getOperation();
                    CPPADCG_ASSERT_UNKNOWN(beginScopeNode != NULL);

                    addScopeToVarOrder(beginScopeNode->getColor(), e);
                }

                //std::cout << "e:" << e << "  " << vorder[p] << "  scope:" << scope << "  p:" << p << "  " << *vorder[p] << std::endl;
                _variableOrder[e++] = vorder[p];
            }
        }

        /**
         * Determines the depth of the first different scope from scope paths of
         * two scopes
         * 
         * @param color1 scope color 1
         * @param color2 scope color 2
         * @return the depth of the first different scope
         */
        inline size_t findFirstDifferentScope(size_t color1, size_t color2) {
            CPPADCG_ASSERT_UNKNOWN(color1 < _scopes.size());
            CPPADCG_ASSERT_UNKNOWN(color2 < _scopes.size());

            ScopePath& scopePath1 = _scopes[color1];
            ScopePath& scopePath2 = _scopes[color2];

            size_t s1 = scopePath1.size();
            size_t s2 = scopePath2.size();
            size_t depth;
            for (depth = 0; depth < s1 && depth < s2; depth++) {
                if (scopePath1[depth].color != scopePath2[depth].color) {
                    break;
                }
            }

            return depth;
        }

        /**
         * Attempt to reduce the number of ifs when there consecutive ifs with
         * the same condition
         */
        inline void optimizeIfs() {
            if (_scopedVariableOrder.size() < 3)
                return; // there has to be at least 2 ifs

            for (size_t scope = 0; scope < _scopedVariableOrder.size(); scope++) {
                std::vector<OperationNode<Base> *>& vorder = _scopedVariableOrder[scope];

                for (long p = vorder.size() - 1; p > 0; p--) {
                    OperationNode<Base>* endIf = vorder[p];
                    if (endIf->getOperationType() != CGEndIfOp)
                        continue;

                    long p1 = p - 1;
                    while (p1 >= 0) {
                        if (vorder[p1]->getOperationType() == CGTmpDclOp) {
                            p1--;
                        } else {
                            break;
                        }
                    }
                    OperationNode<Base>* endIf1 = vorder[p1];
                    if (endIf1->getOperationType() != CGEndIfOp)
                        continue;

                    // 2 consecutive ifs
                    OperationNode<Base>* startIf = endIf->getArguments()[0].getOperation();
                    OperationNode<Base>* startIf1 = endIf1->getArguments()[0].getOperation();
                    if (startIf->getOperationType() != CGStartIfOp || startIf1->getOperationType() != CGStartIfOp)
                        continue;

                    OperationNode<Base>* cond = startIf->getArguments()[0].getOperation();
                    OperationNode<Base>* cond1 = startIf1->getArguments()[0].getOperation();

                    CPPADCG_ASSERT_UNKNOWN(cond->getOperationType() == CGIndexCondExprOp || cond1->getOperationType() == CGIndexCondExprOp);
                    if (cond->getInfo() == cond1->getInfo()) {
                        /**
                         * same condition -> combine the contents into a single if
                         */
                        const std::vector<Argument<Base> >& eArgs = endIf->getArguments();
                        std::vector<Argument<Base> >& eArgs1 = endIf1->getArguments();

                        size_t ifScope = startIf->getColor();
                        size_t ifScope1 = startIf1->getColor();
                        std::vector<OperationNode<Base> *>& vorderIf = _scopedVariableOrder[ifScope];
                        std::vector<OperationNode<Base> *>& vorderIf1 = _scopedVariableOrder[ifScope1];

                        // break cycles caused by dependencies on the previous if
                        for (size_t a = 1; a < eArgs.size(); a++) { // exclude the initial startIf
                            CPPADCG_ASSERT_UNKNOWN(eArgs[a].getOperation() != NULL && eArgs[a].getOperation()->getOperationType() == CGCondResultOp);
                            breakCyclicDependency(eArgs[a].getOperation(), ifScope, endIf1);
                            replaceScope(eArgs[a].getOperation(), ifScope, ifScope1); // update scope
                        }

                        vorderIf1.insert(vorderIf1.end(), vorderIf.begin() + 1, vorderIf.end()); // exclude the initial startIf

                        vorderIf.clear();

                        // update startIf
                        for (size_t a = 1; a < eArgs.size(); a++) { // exclude the initial startIf
                            CPPADCG_ASSERT_UNKNOWN(eArgs[a].getOperation() != NULL && eArgs[a].getOperation()->getOperationType() == CGCondResultOp);
                            eArgs[a].getOperation()->getArguments()[0] = Argument<Base>(*startIf1);
                        }

                        // update endIf
                        eArgs1.insert(eArgs1.end(), eArgs.begin() + 1, eArgs.end());

                        // replace endIf
                        std::vector<Argument<Base> > endIfArgs(1);
                        endIfArgs[0] = Argument<Base>(*endIf1);
                        endIf->setOperation(CGAliasOp, endIfArgs);

                        // remove one of the ifs
                        vorder.erase(vorder.begin() + p);

                        // move nodes in scope containing the ifs
                        for (long pp = p1; pp < p - 1; pp++) {
                            vorder[pp] = vorder[pp + 1];
                        }
                        vorder[p - 1] = endIf1;
                    }
                }
            }
        }

        inline void replaceScope(OperationNode<Base>* node, size_t oldScope, size_t newScope) {
            if (node == NULL || node->getColor() != oldScope)
                return;

            node->setColor(newScope);

            const std::vector<Argument<Base> >& args = node->getArguments();
            for (size_t a = 0; a < args.size(); a++) {
                replaceScope(args[a].getOperation(), oldScope, newScope);
            }
        }

        /**
         * Removes cyclic dependencies when 'ifs' are merged together.
         * Relative variable order must have already been defined.
         * 
         * @todo: avoid visiting the same node!
         * 
         * @param node the node being visited
         * @param scope the scope where the cyclic dependency could appear (or scopes inside it)
         * @param endIf the dependency to remove
         */
        inline void breakCyclicDependency(OperationNode<Base>* node,
                                          size_t scope,
                                          OperationNode<Base>* endIf) {
            if (node == NULL)
                return;

            CGOpCode op = node->getOperationType();
            std::vector<Argument<Base> >& args = node->getArguments();

            if (op == CGTmpOp && args.size() > 1) {
                OperationNode<Base>* arg = args[1].getOperation();
                if (arg == endIf) {
                    // a dependency on CGLoopIndexedTmpOp could be added but
                    // it is not required since variable order was already decided
                    args.erase(args.begin() + 1);
                }
            }

            if (!containedInScope(*node, scope)) {
                return;
            }

            for (size_t a = 0; a < args.size(); a++) {
                OperationNode<Base>* arg = args[a].getOperation();
                if (arg == endIf) {
                    if (op == CGStartIfOp || op == CGLoopStartOp) {
                        args.erase(args.begin() + a);
                        a--;
                    }
                } else {
                    breakCyclicDependency(arg, scope, endIf);
                }
            }
        }

        inline bool containedInScope(const OperationNode<Base>& node, size_t scope) {
            size_t nScope = node.getColor();
            if (nScope == scope)
                return true;

            return _scopes[nScope].size() >= _scopes[scope].size() &&
                    _scopes[nScope][_scopes[scope].size() - 1].color == scope;
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

            size_t aSize = args.size();
            for (size_t argIndex = 0; argIndex < aSize; argIndex++) {
                if (args[argIndex].getOperation() == NULL) {
                    continue;
                }

                OperationNode<Base>& arg = *args[argIndex].getOperation();
                CGOpCode aType = arg.getOperationType();

                if (arg.getUsageCount() == 0) {
                    // dependencies not visited yet
                    checkVariableCreation(arg);

                    if (aType == CGLoopEndOp || aType == CGElseIfOp || aType == CGElseOp || aType == CGEndIfOp) {
                        if (arg.getVariableID() == 0) {
                            // ID value is not really used but must be non-zero
                            arg.setVariableID(std::numeric_limits<size_t>::max());
                        }
                    } else if (aType == CGAtomicForwardOp || aType == CGAtomicReverseOp) {
                        /**
                         * Save atomic function related information
                         */
                        CPPADCG_ASSERT_UNKNOWN(arg.getArguments().size() > 1);
                        CPPADCG_ASSERT_UNKNOWN(arg.getInfo().size() > 1);
                        size_t id = arg.getInfo()[0];
                        const std::string& atomicName = _atomicFunctions.at(id)->afun_name();
                        if (_atomicFunctionsSet.find(atomicName) == _atomicFunctionsSet.end()) {
                            _atomicFunctionsSet.insert(atomicName);
                            _atomicFunctionsOrder->push_back(atomicName);
                        }
                    }

                    /**
                     * make sure new temporary variables are NOT created for
                     * the independent variables and that a dependency did
                     * not use it first
                     */
                    if (arg.getVariableID() == 0 || !isIndependent(arg)) {
                        if (aType == CGLoopIndexedIndepOp) {
                            // ID value not really used but must be non-zero
                            arg.setVariableID(std::numeric_limits<size_t>::max());
                        } else if (aType == CGAliasOp) {
                            continue; // should never be added to the evaluation queue
                        } else if (aType == CGTmpOp) {
                            arg.setVariableID(std::numeric_limits<size_t>::max());
                        } else if (aType == CGLoopStartOp ||
                                aType == CGLoopEndOp ||
                                aType == CGStartIfOp ||
                                aType == CGElseIfOp ||
                                aType == CGElseOp ||
                                aType == CGEndIfOp) {
                            /**
                             * Operation that mark a change in variable scope
                             * are always added
                             */
                            addToEvaluationQueue(arg);
                            if (arg.getVariableID() == 0) {
                                // ID value is not really used but must be non-zero
                                arg.setVariableID(std::numeric_limits<size_t>::max());
                            }
                        } else if (aType == CGPriOp) {
                            addToEvaluationQueue(arg);
                            if (arg.getVariableID() == 0) {
                                // ID value is not really used but must be non-zero
                                arg.setVariableID(std::numeric_limits<size_t>::max());
                            }
                        } else if (aType == CGTmpDclOp) {
                            addToEvaluationQueue(arg);

                            arg.setVariableID(_idCount);
                            _idCount++;

                        } else if (_lang->createsNewVariable(arg) ||
                                _lang->requiresVariableArgument(code.getOperationType(), argIndex)) {

                            addToEvaluationQueue(arg);

                            if (arg.getVariableID() == 0) {
                                if (aType == CGAtomicForwardOp || aType == CGAtomicReverseOp) {
                                    arg.setVariableID(_idAtomicCount);
                                    _idAtomicCount++;
                                } else if (aType == CGLoopIndexedDepOp || aType == CGLoopIndexedTmpOp) {
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

        inline void addToEvaluationQueue(OperationNode<Base>& arg) {
            size_t scope = arg.getColor();
            if (scope >= _scopedVariableOrder.size()) {
                _scopedVariableOrder.resize(scope + 1);
            }

            if (_scopedVariableOrder[scope].empty() &&
                    scope != 0 && // the upper most scope does not need any special node at the beginning
                    _scopes[scope].back().end->getArguments()[0].getOperation() != &arg) {
                // the first node must be a beginning of a scope
                checkVariableCreation(*_scopes[scope].back().end); // go inside a scope from the end 
            }

            // must be after checkVariableCreation() because _scopedVariableOrder might be resized
            std::vector<OperationNode<Base> *>& varOrder = _scopedVariableOrder[scope];

            if (varOrder.size() == varOrder.capacity()) {
                varOrder.reserve((varOrder.size() * 3) / 2 + 1);
            }

            varOrder.push_back(&arg);
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
                        CPPADCG_ASSERT_UNKNOWN(freeArrayStartSpace.size() == freeArrayEndSpace.size());
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
                    CPPADCG_ASSERT_UNKNOWN(freeArrayStartSpace.size() == freeArrayEndSpace.size());
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
                    CPPADCG_ASSERT_UNKNOWN(otherArray.getVariableID() > 0); // make sure it had already been assigned space
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
            CGOpCode op = code.getOperationType();

            if (op == CGLoopEndOp) {
                LoopEndOperationNode<Base>& loopEnd = static_cast<LoopEndOperationNode<Base>&> (code);
                _loopDepth++;
                _loopOuterVars.resize(_loopDepth + 1);
                _loopStartEvalOrder.push_back(loopEnd.getLoopStart().getEvaluationOrder());

            } else if (op == CGLoopStartOp) {
                _loopDepth--; // leaving the current loop
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

            if (op == CGLoopEndOp) {
                /**
                 * temporary variables from outside the loop which are used
                 * within the loop cannot be overwritten inside that loop
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

            } else if (op == CGLoopStartOp) {
                _loopDepth++; // coming back to the loop
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
                    op != CGStartIfOp &&
                    op != CGElseIfOp &&
                    op != CGElseOp &&
                    op != CGEndIfOp &&
                    op != CGLoopIndexedDepOp &&
                    op != CGLoopIndexedIndepOp &&
                    op != CGLoopIndexedTmpOp && // not considered as a temporary (the temporary is CGTmpDclOp)
                    op != CGIndexOp &&
                    op != CGIndexAssignOp &&
                    op != CGTmpOp && // not considered as a temporary (the temporary is CGTmpDclOp)
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
                    CPPADCG_ASSERT_UNKNOWN(aa->getArguments().size() == 1);
                    aa = aa->getArguments()[0].getOperation();
                } while (aa != NULL && aa->getOperationType() == CGAliasOp);
                return aa;
            }
        }

        inline void resetManagedNodes() {
            _variableOrder.clear();
            _scopedVariableOrder.resize(1);
            _scopedVariableOrder[0].clear();

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
