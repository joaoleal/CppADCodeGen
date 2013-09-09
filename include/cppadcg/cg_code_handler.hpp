#ifndef CPPAD_CG_CODE_HANDLER_INCLUDED
#define CPPAD_CG_CODE_HANDLER_INCLUDED
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
        std::vector<CG<Base> >* _dependents;
        // all the source code blocks created with the CG<Base> objects
        std::vector<OperationNode<Base> *> _codeBlocks;
        // the order for the variable creation in the source code
        std::vector<OperationNode<Base> *> _variableOrder;
        // maps the ids of the atomic functions
        std::map<size_t, CGAbstractAtomicFun<Base>*> _atomicFunctions;
        // maps the loop ids of the loop atomic functions
        std::map<size_t, LoopAtomicFun<Base>*> _loops;
        //
        std::set<const Index*> _indexes;
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
        // the language used for source code generation
        Language<Base>* _lang;
        // the lowest ID used for temporary variables
        size_t _minTemporaryVarID;
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
            _lang(NULL),
            _minTemporaryVarID(0),
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
            typename std::map<size_t, LoopAtomicFun<Base>*>::const_iterator it;
            it = _loops.find(id);
            if (it != _loops.end())
                return &(it->second->afun_name());
            else
                return NULL;
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
                                  std::vector<CG<Base> >& dependent,
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
                                  std::vector<CG<Base> >& dependent,
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
            for (size_t i = 0; i < atomicFunctions.size(); i++) {
                _atomicFunctionsSet.insert(atomicFunctions[i]);
            }

            if (_used) {
                resetCounters();
            }
            _used = true;

            /**
             * the first variable IDs are for the independent variables
             */
            for (typename std::vector<OperationNode<Base> *>::iterator it = _independentVariables.begin(); it != _independentVariables.end(); ++it) {
                (*it)->setVariableID(_idCount++);
            }

            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                if (it->getOperationNode() != NULL && it->getOperationNode()->getVariableID() == 0) {
                    it->getOperationNode()->setVariableID(_idCount++);
                }
            }

            _minTemporaryVarID = _idCount;

            // determine the number of times each variable is used
            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.getOperationNode() != NULL) {
                    OperationNode<Base>& code = *var.getOperationNode();
                    markCodeBlockUsed(code);
                }
            }

            // determine the variable creation order

            for (size_t i = 0; i < dependent.size(); i++) {
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

            //assert(_idCount - 1 + _idArrayCount == _variableOrder.size() + _independentVariables.size());

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
                                              _loopDependentIndexPatterns, _loopIndependentIndexPatterns);
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

        const std::map<size_t, LoopAtomicFun<Base>*>& getLoops() const;

        LoopAtomicFun<Base>* getLoop(size_t loopId) const;

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

            if (code.getTotalUsageCount() == 1) {
                // first time this operation is visited

                const std::vector<Argument<Base> >& args = code.arguments_;

                typename std::vector<Argument<Base> >::const_iterator it;
                for (it = args.begin(); it != args.end(); ++it) {
                    if (it->getOperation() != NULL) {
                        OperationNode<Base>& arg = *it->getOperation();
                        markCodeBlockUsed(arg);
                    }
                }

                if (code.getOperationType() == CGIndexOp) {
                    const IndexOperationNode<Base>& inode = static_cast<const IndexOperationNode<Base>&> (code);
                    // indexes that don't depend on a loop start or an index assignment are declared elsewhere
                    if (inode.getArguments().size() > 0) {
                        _indexes.insert(&inode.getIndex());
                    }
                } else if (code.getOperationType() == CGDependentRefOp) {
                    assert(code.getInfo().size() == 1);
                    size_t depIndex = code.getInfo()[0];

                    assert(_dependents->size() > depIndex);
                    OperationNode<Base>* depNode = (*_dependents)[depIndex].getOperationNode();
                    assert(depNode != NULL && depNode->getOperationType() != CGInvOp);

                    code.setVariableID(depNode->getVariableID());
                }
            }
        }

        virtual void registerAtomicFunction(CGAbstractAtomicFun<Base>& atomic) {
            _atomicFunctions[atomic.getId()] = &atomic;
        }

        /***********************************************************************
         *                        Loop management
         **********************************************************************/
        virtual void registerLoop(LoopAtomicFun<Base>& loop);

        /***********************************************************************
         * 
         **********************************************************************/
        virtual void checkVariableCreation(OperationNode<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;

            for (it = args.begin(); it != args.end(); ++it) {
                if (it->getOperation() != NULL) {
                    OperationNode<Base>& arg = *it->getOperation();

                    if (arg.getUsageCount() == 0) {
                        // dependencies not visited yet
                        checkVariableCreation(arg);

                        /**
                         * Save atomic function related information
                         */
                        if (arg.getOperationType() == CGAtomicForwardOp || arg.getOperationType() == CGAtomicReverseOp) {
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
                    // make sure new temporary variables are NOT created for
                    // the independent variables and that a dependency did
                    // not use it first
                    if ((arg.getVariableID() == 0 || !isIndependent(arg)) && arg.getUsageCount() == 0) {
                        if (arg.getOperationType() == CGLoopIndexedIndepOp) {
                            // ID value not really used but must be non-zero
                            arg.setVariableID(std::numeric_limits<size_t>::max());
                        } else {
                            size_t argIndex = it - args.begin();
                            if (arg.getOperationType() == CGLoopStartOp || arg.getOperationType() == CGLoopEndOp) {
                                if (arg.getVariableID() == 0) {
                                    addToEvaluationQueue(arg);
                                    // ID value not really used but must be non-zero
                                    arg.setVariableID(std::numeric_limits<size_t>::max());
                                }
                            } else if (_lang->createsNewVariable(arg) ||
                                    _lang->requiresVariableArgument(code.getOperationType(), argIndex)) {
                                addToEvaluationQueue(arg);
                                if (arg.getVariableID() == 0) {
                                    if (arg.getOperationType() == CGAtomicForwardOp || arg.getOperationType() == CGAtomicReverseOp) {
                                        arg.setVariableID(_idAtomicCount);
                                        _idAtomicCount++;
                                    } else if (arg.getOperationType() == CGLoopIndexedDepOp) {
                                        // ID value not really used but must be non-zero
                                        arg.setVariableID(std::numeric_limits<size_t>::max());
                                    } else if (arg.getOperationType() == CGArrayCreationOp) {
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

        inline void reduceTemporaryVariables(std::vector<CG<Base> >& dependent) {

            /**
             * determine the last line where each temporary variable is used
             */
            resetUsageCount();

            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.getOperationNode() != NULL) {
                    OperationNode<Base>& code = *var.getOperationNode();
                    if (code.use_count_ == 0) {
                        // dependencies not visited yet
                        determineLastTempVarUsage(code);
                    }
                    code.increaseUsageCount();
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

        inline void determineLastTempVarUsage(OperationNode<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;

            /**
             * count variable usage
             */
            for (it = args.begin(); it != args.end(); ++it) {
                if (it->getOperation() != NULL) {
                    OperationNode<Base>& arg = *it->getOperation();

                    if (arg.use_count_ == 0) {
                        // dependencies not visited yet
                        determineLastTempVarUsage(arg);
                    }

                    arg.increaseUsageCount();

                    if (arg.getLastUsageEvaluationOrder() < code.getEvaluationOrder()) {
                        arg.setLastUsageEvaluationOrder(code.getEvaluationOrder());
                    }
                }
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

        virtual void resetCounters() {
            _variableOrder.clear();

            for (typename std::vector<OperationNode<Base> *>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                OperationNode<Base>* block = *it;
                block->resetHandlerCounters();
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
        friend class LoopAtomicFun<Base>;

    };

}
#endif

