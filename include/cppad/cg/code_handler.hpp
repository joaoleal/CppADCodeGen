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
namespace cg {

/**
 * Helper class to analyze the operation graph and generate source code
 * for several languages
 * 
 * @author Joao Leal
 */
template<class Base>
class CodeHandler {
    friend class CodeHandlerVectorSync<Base>;
public:
    typedef std::vector<OperationPathNode<Base> > SourceCodePath;
    typedef std::vector<ScopePathElement<Base> > ScopePath;
    typedef unsigned short ScopeIDType;
protected:
    struct LoopData; // forward declaration

protected:
    // counter used to determine visitation IDs for the operation tree
    size_t _idVisit;
    // counter used to generate variable IDs
    size_t _idCount;
    // counter used to generate array variable IDs
    size_t _idArrayCount;
    // counter used to generate sparse array variable IDs
    size_t _idSparseArrayCount;
    // counter used to generate IDs for atomic functions
    size_t _idAtomicCount;
    // the independent variables
    std::vector<OperationNode<Base> *> _independentVariables;
    // the current dependent variables
    CppAD::vector<CG<Base> >* _dependents;
    /**
     * nodes managed by this code handler which include all
     * all OperationNodes created by CG<Base> objects
     */
    std::vector<OperationNode<Base> *> _codeBlocks;
    /**
     * All CodeHandlerVector associated with this code handler
     */
    std::set<CodeHandlerVectorSync<Base>*> _managedVectors;
    /**
     * the ID of the last visit to each managed node
     */
    CodeHandlerVector<Base, size_t> _lastVisit;
    /**
     * scope of each managed operation node
     */
    CodeHandlerVector<Base, ScopeIDType> _scope;
    /**
     * evaluation order of each managed node
     * (zero means that an evaluation position was never assigned)
     */
    CodeHandlerVector<Base, size_t> _evaluationOrder;
    /**
     * the last index in the evaluation order for which an operation node 
     * is taken as an argument of another operation node.
     * (zero means that the node was never used)
     */
    CodeHandlerVector<Base, size_t> _lastUsageOrder;
    /**
     * the total number of times the result of an operation node  is used
     */
    CodeHandlerVector<Base, size_t> _totalUseCount;
    /**
     * Provides the variable ID that was altered/assigned to operation nodes.
     * Zero means that no variable is assigned.
     */
    CodeHandlerVector<Base, size_t> _varId;
    /**
     * the order for the variable creation in the source code 
     */
    std::vector<OperationNode<Base> *> _variableOrder;
    /**
     * the order for the variable creation in the source code 
     * (each level represents a different variable scope)
     */
    std::vector<std::vector<OperationNode<Base> *> > _scopedVariableOrder;
    /**
     *
     */
    LoopData _loops;
    /**
     * maps the IDs of the atomic functions
     */
    std::map<size_t, CGAbstractAtomicFun<Base>*> _atomicFunctions;
    /**
     * already used atomic function names (may contain names which were 
     * used by previous calls to this/other CondeHandlers)
     */
    std::map<std::string, size_t> _atomicFunctionName2Index;
    /**
     * the order of the atomic functions (may contain names which were 
     * used by previous calls to this/other CondeHandlers)
     */
    std::vector<std::string>* _atomicFunctionsOrder;
    /**
     * 
     */
    std::map<size_t, size_t> _atomicFunctionId2Index;
    /**
     * the maximum forward mode order each atomic function is called
     * (-1 means forward mode not used)
     */
    std::vector<int> _atomicFunctionsMaxForward;
    /**
     * the maximum reverse mode order each atomic function is called
     * (-1 means reverse mode not used)
     */
    std::vector<int> _atomicFunctionsMaxReverse;
    // a flag indicating if this handler was previously used to generate code
    bool _used;
    // a flag indicating whether or not to reuse the IDs of destroyed variables
    bool _reuseIDs;
    // scope color/index counter
    ScopeIDType _scopeColorCount;
    // the current scope color/index counter
    ScopeIDType _currentScopeColor;
    // all scopes
    std::vector<ScopePath> _scopes;
    // possible altered nodes due to scope conditionals (altered node <-> clone of original)
    std::list<std::pair<OperationNode<Base>*, OperationNode<Base>* > > _alteredNodes;
    // the language used for source code generation
    Language<Base>* _lang;
    /**
     * information sent to the language
     */
    std::unique_ptr<LanguageGenerationData<Base> > _info;
    // the lowest ID used for temporary variables
    size_t _minTemporaryVarID;
    /**
     * whether or not the dependent variables should be zeroed before 
     * executing the operation graph
     */
    bool _zeroDependents;
    //
    bool _verbose;
    /**
     * used to track evaluation times and print out messages
     */
    JobTimer* _jobTimer;
    /**
     * Auxiliary index declaration (might not be used)
     */
    OperationNode<Base>* _auxIndexI;
    /**
     * Auxiliary index (might not be used)
     */
    IndexOperationNode<Base>* _auxIterationIndexOp;
public:

    CodeHandler(size_t varCount = 50) :
        _idVisit(1),
        _idCount(1),
        _idArrayCount(1),
        _idSparseArrayCount(1),
        _idAtomicCount(1),
        _dependents(nullptr),
        _lastVisit(*this),
        _scope(*this),
        _evaluationOrder(*this),
        _lastUsageOrder(*this),
        _totalUseCount(*this),
        _varId(*this),
        _scopedVariableOrder(1),
        _atomicFunctionsOrder(nullptr),
        _used(false),
        _reuseIDs(true),
        _scopeColorCount(0),
        _currentScopeColor(0),
        _lang(nullptr),
        _minTemporaryVarID(0),
        _zeroDependents(false),
        _verbose(false),
        _jobTimer(nullptr) {
        _codeBlocks.reserve(varCount);
        //_variableOrder.reserve(1 + varCount / 3);
        _scopedVariableOrder[0].reserve(1 + varCount / 3);

        _auxIndexI = makeIndexDclrNode("i");
        _auxIterationIndexOp = makeIndexNode(*_auxIndexI);
    }

    CodeHandler(const CodeHandler&) = delete;
    CodeHandler& operator=(const CodeHandler&) = delete;

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
        _independentVariables.push_back(makeNode(CGOpCode::Inv));
        variable.makeVariable(*_independentVariables.back());
    }

    size_t getIndependentVariableSize() const {
        return _independentVariables.size();
    }

    /**
     * @throws CGException if a variable is not found in the independent vector
     */
    size_t getIndependentVariableIndex(const OperationNode<Base>& var) const {
        CPPADCG_ASSERT_UNKNOWN(var.getOperationType() == CGOpCode::Inv);

        typename std::vector<OperationNode<Base> *>::const_iterator it =
                std::find(_independentVariables.begin(), _independentVariables.end(), &var);
        if (it == _independentVariables.end()) {
            throw CGException("Variable not found in the independent variable vector");
        }

        return it - _independentVariables.begin();
    }

    /**
     * Provides variable IDs that were assigned to operation nodes.
     * Zero means that no variable is assigned.
     * The first IDs are reserved for the independent variables.
     * It can be an empty vector if IDs have not yet been assigned.
     */
    inline const CodeHandlerVector<Base, size_t>& getVariablesIDs() const {
      return _varId;
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
     * before executing the operation graph
     * 
     * @param true if the dependents should be zeroed
     */
    inline void setZeroDependents(bool zeroDependents) {
        _zeroDependents = zeroDependents;
    }

    inline size_t getOperationTreeVisitId() const {
        return _idVisit;
    }

    inline void startNewOperationTreeVisit() {
        assert(_idVisit < std::numeric_limits<size_t>::max());

        _lastVisit.adjustSize();
        _idVisit++;
    }

    inline bool isVisited(const OperationNode<Base>& node) const {
        size_t p = node.getHandlerPosition();
        return p < _lastVisit.size() && _lastVisit[node] == _idVisit;
    }

    inline void markVisited(OperationNode<Base>& node) {
        _lastVisit.adjustSize(node);
        _lastVisit[node] = _idVisit;
    }

    /**
     * Provides the name used by an atomic function with a given ID.
     * 
     * @param id the atomic function ID.
     * @return a pointer to the atomic function name if it was registered
     *         or nullptr otherwise
     */
    inline const std::string* getAtomicFunctionName(size_t id) const {
        typename std::map<size_t, CGAbstractAtomicFun<Base>*>::const_iterator it;
        it = _atomicFunctions.find(id);
        if (it != _atomicFunctions.end())
            return &(it->second->afun_name());
        else
            return nullptr;
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
     * Provides the maximum forward mode order used by all atomic functions
     * in the last call to ::generateCode 
     * (-1 means forward mode not used).
     */
    const std::vector<int>& getExternalFuncMaxForwardOrder() const {
        return _atomicFunctionsMaxForward;
    }

    /**
     * Provides the maximum reverse mode order used by all atomic functions
     * in the last call to ::generateCode
     * (-1 means forward mode not used).
     */
    const std::vector<int>& getExternalFuncMaxReverseOrder() const {
        return _atomicFunctionsMaxReverse;
    }

    /**
     * Provides the name used by a loop atomic function with a given ID.
     * 
     * @param id the atomic function ID.
     * @return a pointer to the atomic loop function name if it was
     *         registered or nullptr otherwise
     */
    inline const std::string* getLoopName(size_t id) const {
        return _loops.getLoopName(id);
    }

    inline const std::vector<ScopePath>& getScopes() const {
        return _scopes;
    }

    /**************************************************************************
     *                       Graph management functions
     *************************************************************************/
    /**
     * Finds occurrences of a source code fragment in an operation graph.
     * 
     * @param root the operation graph where to search
     * @param target the source code fragment to find in root
     * @param max the maximum number of occurrences of code to find in root
     * @return the paths from root to code
     */
    inline std::vector<SourceCodePath> findPaths(OperationNode<Base>& root,
                                                 OperationNode<Base>& target,
                                                 size_t max);

    inline BidirGraph<Base> findPathGraph(OperationNode<Base>& root,
                                          OperationNode<Base>& target) ;

    inline BidirGraph<Base> findPathGraph(OperationNode<Base>& root,
                                          OperationNode<Base>& target,
                                          size_t& bifurcations,
                                          size_t maxBifurcations = std::numeric_limits<size_t>::max());

    /**************************************************************************
     *                       Source code generation
     *************************************************************************/

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
                              Language<Base>& lang,
                              CppAD::vector<CG<Base> >& dependent,
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
                              Language<Base>& lang,
                              CppAD::vector<CG<Base> >& dependent,
                              VariableNameGenerator<Base>& nameGen,
                              std::vector<std::string>& atomicFunctions,
                              const std::string& jobName = "source") {
        using namespace std::chrono;
        steady_clock::time_point beginTime;

        if (_jobTimer != nullptr) {
            _jobTimer->startingJob("source for '" + jobName + "'");
        } else if (_verbose) {
            std::cout << "generating source for '" << jobName << "' ... ";
            std::cout.flush();
            beginTime = steady_clock::now();
        }

        _lang = &lang;
        _idCount = 1;
        _idArrayCount = 1;
        _idSparseArrayCount = 1;
        _idAtomicCount = 1;
        _dependents = &dependent;
        _atomicFunctionsOrder = &atomicFunctions;
        _atomicFunctionsMaxForward.resize(atomicFunctions.size());
        _atomicFunctionsMaxReverse.resize(atomicFunctions.size());
        _atomicFunctionName2Index.clear();
        _loops.prepare4NewSourceGen();
        _scopeColorCount = 0;
        _currentScopeColor = 0;
        _scopes.reserve(4);
        _scopes.resize(1);
        _alteredNodes.clear();
        _evaluationOrder.adjustSize();
        _lastUsageOrder.adjustSize();
        _totalUseCount.adjustSize();
        _varId.adjustSize();
        _scope.adjustSize();

        for (size_t i = 0; i < atomicFunctions.size(); i++) {
            _atomicFunctionName2Index[atomicFunctions[i]] = i;
        }
        std::fill(_atomicFunctionsMaxForward.begin(), _atomicFunctionsMaxForward.end(), -1);
        std::fill(_atomicFunctionsMaxReverse.begin(), _atomicFunctionsMaxReverse.end(), -1);

        if (_used) {
            resetManagedNodes();
        }
        _used = true;

        /**
         * the first variable IDs are for the independent variables
         */
        size_t n = _independentVariables.size();
        for (size_t j = 0; j < n; j++) {
            _varId[*_independentVariables[j]] = _idCount++;
        }

        size_t m = dependent.size();
        for (size_t i = 0; i < m; i++) {
            OperationNode<Base>* node = dependent[i].getOperationNode();
            if (node != nullptr && _varId[*node] == 0) {
                _varId[*node] = _idCount++;
            }
        }

        _minTemporaryVarID = _idCount;

        /**
         * determine the number of times each variable is used
         */
        for (size_t i = 0; i < m; i++) {
            OperationNode<Base>* node = dependent[i].getOperationNode();
            if (node != nullptr) {
                markCodeBlockUsed(*node);
            }
        }

        /**
         * determine the variable creation order
         */
        _scopedVariableOrder.reserve(std::max(size_t(1), _scopes.size()) + 10); // some additional scopes might still be added

        startNewOperationTreeVisit();

        for (size_t i = 0; i < m; i++) {
            CG<Base>& var = dependent[i];
            if (var.getOperationNode() != nullptr) {
                OperationNode<Base>& code = *var.getOperationNode();
                if (!isVisited(code)) {
                    // dependencies not visited yet
                    checkVariableCreation(code);

                    // make sure new temporary variables are NOT created for
                    // the independent variables and that a dependency did
                    // not use it first
                    if ((_varId[code] == 0 || !isIndependent(code)) && !isVisited(code)) {
                        addToEvaluationQueue(code);
                    }
                }
                markVisited(code);
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
            setEvaluationOrder(arg, p + 1);
            dependentAdded2EvaluationQueue(arg);
        }

        /**
         * Reuse temporary variables
         */
        if (_reuseIDs) {
            reduceTemporaryVariables(dependent);
        }

        nameGen.setTemporaryVariableID(_minTemporaryVarID, _idCount - 1, _idArrayCount - 1, _idSparseArrayCount - 1);

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
        _info.reset(new LanguageGenerationData<Base>(_independentVariables, dependent,
                    _minTemporaryVarID, _varId, _variableOrder,
                    nameGen,
                    atomicFunctionId2Index, atomicFunctionId2Name,
                    _atomicFunctionsMaxForward, _atomicFunctionsMaxReverse,
                    _reuseIDs,
                    _loops.indexes, _loops.indexRandomPatterns,
                    _loops.dependentIndexPatterns, _loops.independentIndexPatterns,
                    _totalUseCount, _scope, *_auxIterationIndexOp,
                    _zeroDependents));
        lang.generateSourceCode(out, _info);

        /**
         * clean-up
         */
        _atomicFunctionName2Index.clear();

        // restore altered nodes
        for (const auto& itAlt : _alteredNodes) {
            OperationNode<Base>* tmp = itAlt.first;
            OperationNode<Base>* opClone = itAlt.second;
            if (tmp->getOperationType() == CGOpCode::Tmp && !tmp->getInfo().empty()) { // some might have already been restored
                tmp->setOperation(opClone->getOperationType(), opClone->getArguments());
                tmp->getInfo() = opClone->getInfo();
            }
        }
        _alteredNodes.clear();

        if (_jobTimer != nullptr) {
            _jobTimer->finishedJob();
        } else if (_verbose) {
            duration<float> dt = steady_clock::now() - beginTime;
            std::cout << "done [" << std::fixed << std::setprecision(3)
                    << dt.count() << "]" << std::endl;
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

    size_t getTemporarySparseArraySize() const {
        return _idSparseArrayCount - 1;
    }

    /**************************************************************************
     *                       Reusing handler and nodes
     *************************************************************************/

    /**
     * Resets this handler for a usage with completely different nodes.
     * @warning all managed memory will be deleted
     */
    virtual void reset() {
        for (OperationNode<Base>* n : _codeBlocks) {
            delete n;
        }
        _codeBlocks.clear();
        _independentVariables.clear();
        _idCount = 1;
        _idArrayCount = 1;
        _idSparseArrayCount = 1;
        _idAtomicCount = 1;

        _loops.reset();

        _used = false;
    }

    /**
     * Resets the previously used dependents and their children so that they
     * can be reused again by this handler.
     */
    inline void resetNodes() {
        _scope.fill(0);
        _evaluationOrder.fill(0);
        _lastUsageOrder.fill(0);
        _totalUseCount.fill(0);
        _varId.fill(0);
    }

    /**************************************************************************
     *                         access to managed memory
     *************************************************************************/

    /**
     * Creates a shallow clone of an operation node
     */
    inline OperationNode<Base>* cloneNode(const OperationNode<Base>& n) {
        return manageOperationNode(new OperationNode<Base>(n));
    }

    inline OperationNode<Base>* makeNode(CGOpCode op) {
        return manageOperationNode(new OperationNode<Base>(this, op));
    }

    inline OperationNode<Base>* makeNode(CGOpCode op,
                                         const Argument<Base>& arg) {
        return manageOperationNode(new OperationNode<Base>(this, op, arg));
    }

    inline OperationNode<Base>* makeNode(CGOpCode op,
                                         std::vector<Argument<Base> >&& args) {
        return manageOperationNode(new OperationNode<Base>(this, op, std::move(args)));
    }

    inline OperationNode<Base>* makeNode(CGOpCode op,
                                         std::vector<size_t>&& info,
                                         std::vector<Argument<Base> >&& args) {
        return manageOperationNode(new OperationNode<Base>(this, op, std::move(info), std::move(args)));
    }

    inline OperationNode<Base>* makeNode(CGOpCode op,
                                         const std::vector<size_t>& info,
                                         const std::vector<Argument<Base> >& args) {
        return manageOperationNode(new OperationNode<Base>(this, op, info, args));
    }

    inline LoopStartOperationNode<Base>* makeLoopStartNode(OperationNode<Base>& indexDcl,
                                                           size_t iterationCount) {
        auto* n = new LoopStartOperationNode<Base>(this, indexDcl, iterationCount);
        manageOperationNode(n);
        return n;
    }

    inline LoopStartOperationNode<Base>* makeLoopStartNode(OperationNode<Base>& indexDcl,
                                                           IndexOperationNode<Base>& iterCount) {
        auto* n = new LoopStartOperationNode<Base>(this, indexDcl, iterCount);
        manageOperationNode(n);
        return n;
    }

    inline LoopEndOperationNode<Base>* makeLoopEndNode(LoopStartOperationNode<Base>& loopStart,
                                                       const std::vector<Argument<Base> >& endArgs) {
        auto* n = new LoopEndOperationNode<Base>(this, loopStart, endArgs);
        manageOperationNode(n);
        return n;
    }

    inline PrintOperationNode<Base>* makePrintNode(const std::string& before,
                                                   const Argument<Base>& arg,
                                                   const std::string& after) {
        auto* n = new PrintOperationNode<Base>(this, before, arg, after);
        manageOperationNode(n);
        return n;
    }

    inline IndexOperationNode<Base>* makeIndexNode(OperationNode<Base>& indexDcl) {
        auto* n = new IndexOperationNode<Base>(this, indexDcl);
        manageOperationNode(n);
        return n;
    }

    inline IndexOperationNode<Base>* makeIndexNode(LoopStartOperationNode<Base>& loopStart) {
        auto* n = new IndexOperationNode<Base>(this, loopStart);
        manageOperationNode(n);
        return n;
    }

    inline IndexOperationNode<Base>* makeIndexNode(IndexAssignOperationNode<Base>& indexAssign) {
        auto* n = new IndexOperationNode<Base>(this, indexAssign);
        manageOperationNode(n);
        return n;
    }

    inline IndexAssignOperationNode<Base>* makeIndexAssignNode(OperationNode<Base>& index,
                                                               IndexPattern& indexPattern,
                                                               IndexOperationNode<Base>& index1) {
        auto* n = new IndexAssignOperationNode<Base>(this, index, indexPattern, index1);
        manageOperationNode(n);
        return n;
    }

    inline IndexAssignOperationNode<Base>* makeIndexAssignNode(OperationNode<Base>& index,
                                                               IndexPattern& indexPattern,
                                                               IndexOperationNode<Base>* index1,
                                                               IndexOperationNode<Base>* index2) {
        auto* n = new IndexAssignOperationNode<Base>(this, index, indexPattern, index1, index2);
        manageOperationNode(n);
        return n;
    }

    inline OperationNode<Base>* makeIndexDclrNode(const std::string& name) {
        CPPADCG_ASSERT_KNOWN(!name.empty(), "index name cannot be empty");
        auto* n = manageOperationNode(new OperationNode<Base>(this, CGOpCode::IndexDeclaration));
        n->setName(name);
        return n;
    }

    /**
     * Provides the current number of OperationNodes created by the model.
     * This number is not the total number of operations in the final
     * model since it also contains Operations nodes marking
     * independent variables and there could be unused operations by
     * the model (dead-code).
     * @return The number of OperationNodes created by the model.
     */
    inline size_t getManagedNodesCount() const {
        return _codeBlocks.size();
    }

    /**
     * Provides the OperationNodes created by the model.
     */
    inline const std::vector<OperationNode<Base> *>& getManagedNodes() const {
        return _codeBlocks;
    }

    /**
     * Allows to delete OperationNodes that are managed internally.
     * @warning: This is a dangerous method, make sure these nodes are not used
     *           anywhere else!
     * @param start The index of the first OperationNode to be deleted
     * @param end The index after the last OperationNode to be deleted
     */
    inline void deleteManagedNodes(size_t start, size_t end) {
        if (start >= end)
            return;

        start = std::min<size_t>(start, _codeBlocks.size());
        end = std::min<size_t>(end, _codeBlocks.size());

        for (size_t i = start; i < end; ++i) {
            delete _codeBlocks[i];
        }
        _codeBlocks.erase(_codeBlocks.begin() + start, _codeBlocks.begin() + end);

        // update positions
        for (size_t i = start; i < _codeBlocks.size(); ++i) {
            _codeBlocks[i]->setHandlerPosition(i);
        }

        for (auto* v : _managedVectors) {
            v->nodesErased(start, end);
        }
    }

    /**************************************************************************
     *                           Value generation
     *************************************************************************/
    CG<Base> createCG(const Argument<Base>& arg) {
        return CG<Base>(arg);
    }

    /**************************************************************************
     *                           Loop management
     *************************************************************************/

    const std::map<size_t, LoopModel<Base>*>& getLoops() const;

    inline LoopModel<Base>* getLoop(size_t loopId) const {
        return _loops.getLoop(loopId);
    }

    inline size_t addLoopDependentIndexPattern(IndexPattern& jacPattern) {
        return _loops.addDependentIndexPattern(jacPattern);
    }

    inline void manageLoopDependentIndexPattern(const IndexPattern* pattern) {
        _loops.manageDependentIndexPattern(pattern);
    }

    inline size_t addLoopIndependentIndexPattern(IndexPattern& pattern, size_t hint) {
        return _loops.addIndependentIndexPattern(pattern, hint);
    }

    /***********************************************************************
     *                           Index patterns
     **********************************************************************/
    static inline void findRandomIndexPatterns(IndexPattern* ip,
                                               std::set<RandomIndexPattern*>& found) {
        if (ip == nullptr)
            return;

        if (ip->getType() == IndexPatternType::Random1D || ip->getType() == IndexPatternType::Random2D) {
            found.insert(static_cast<RandomIndexPattern*> (ip));
        } else {
            std::set<IndexPattern*> indexes;
            ip->getSubIndexes(indexes);
            for (IndexPattern* sip : indexes) {
                if (sip->getType() == IndexPatternType::Random1D || sip->getType() == IndexPatternType::Random2D)
                    found.insert(static_cast<RandomIndexPattern*> (sip));
            }
        }
    }

    /**************************************************************************
     *                      Operation graph manipulation
     *************************************************************************/

    /**
     * Solves an expression (e.g. f(x, y) == 0) for a given variable (e.g. x).
     *
     * @param expression  The original expression (f(x, y))
     * @param var  The variable to solve for
     * @return  The expression for the variable
     * @throws CGException if it is not possible to solve the expression
     */
    inline CG<Base> solveFor(OperationNode<Base>& expression,
                             OperationNode<Base>& var);

    inline bool isSolvable(OperationNode<Base>& expression,
                           OperationNode<Base>& var);

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
     * @throws CGException if the dependent variable does not belong to this handler
     */
    inline void substituteIndependent(const CG<Base>& indep,
                                      const CG<Base>& dep,
                                      bool removeFromIndeps = true);

    inline void substituteIndependent(OperationNode<Base>& indep,
                                      OperationNode<Base>& dep,
                                      bool removeFromIndeps = true);

    /**
     * Reverts a substitution of an independent variable that has not been 
     * removed from the list of independents yet.
     * Warning: it does not recover any custom name assigned to the variable.
     * 
     * @param indep The independent variable
     * @throws CGException if the dependent variable does not belong to this handler
     */
    inline void undoSubstituteIndependent(OperationNode<Base>& indep);

    /**
     * Finalizes the substitution of an independent variable by eliminating
     * it from the list of independents. After this operation the variable
     * substitution cannot be undone.
     * 
     * @param indep The independent variable
     * @throws CGException if the dependent variable is not an not an alias or it does not belong to this handler
     */
    inline void removeIndependent(OperationNode<Base>& indep);

    /**
     * Adds an operation node to the list of nodes to be deleted when this
     * handler is destroyed.
     * 
     * @param code The operation node to be managed.
     * @return true if the node was successfully added to the list or
     *         false if it had already been previously added.
     */
    inline bool manageOperationNodeMemory(OperationNode<Base>* code) {
        size_t pos = code->getHandlerPosition();
        if (pos >= _codeBlocks.size() || code != _codeBlocks[pos]) {
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
        
        for (auto* v : _managedVectors) {
            v->handler_ = nullptr;
        }
    }

protected:

    virtual OperationNode<Base>* manageOperationNode(OperationNode<Base>* code) {
        //CPPADCG_ASSERT_UNKNOWN(std::find(_codeBlocks.begin(), _codeBlocks.end(), code) == _codeBlocks.end()); // <<< too great of an impact in performance
        if (_codeBlocks.capacity() == _codeBlocks.size()) {
            _codeBlocks.reserve((_codeBlocks.size()*3) / 2 + 1);
        }

        code->setHandlerPosition(_codeBlocks.size());
        _codeBlocks.push_back(code);
        return code;
    }

    inline void addVector(CodeHandlerVectorSync<Base>* v) {
        _managedVectors.insert(v);
    }

    inline void removeVector(CodeHandlerVectorSync<Base>* v) {
        _managedVectors.erase(v);
    }

    virtual void markCodeBlockUsed(OperationNode<Base>& code) {
        increaseTotalUsageCount(code);

        CGOpCode op = code.getOperationType();
        if (isIndependent(code)) {
            return; // nothing to do
        } else if (op == CGOpCode::Alias) {
            /**
             * Alias operations are always followed so that there is a 
             * correct usage count at the operation that it points to
             */
            CPPADCG_ASSERT_UNKNOWN(code.getArguments().size() == 1);
            OperationNode<Base>* arg = code.getArguments()[0].getOperation();
            if (arg != nullptr) {
                markCodeBlockUsed(*arg);
            }

        } else if (getTotalUsageCount(code) == 1) {
            // first time this operation is visited

            size_t previousScope = _currentScopeColor;

            _scope[code] = _currentScopeColor;

            // check if there is a scope change
            if (op == CGOpCode::LoopStart || op == CGOpCode::StartIf || op == CGOpCode::ElseIf || op == CGOpCode::Else) {
                // leaving a scope
                ScopePath& sPath = _scopes[_currentScopeColor];
                CPPADCG_ASSERT_UNKNOWN(sPath.back().beginning == nullptr);
                if (op == CGOpCode::LoopStart || op == CGOpCode::StartIf) {
                    sPath.back().beginning = &code; // save the initial node
                } else {
                    CPPADCG_ASSERT_UNKNOWN(!code.getArguments().empty() &&
                                           code.getArguments()[0].getOperation() != nullptr &&
                                           code.getArguments()[0].getOperation()->getOperationType() == CGOpCode::StartIf);
                    sPath.back().beginning = code.getArguments()[0].getOperation(); // save the initial node
                }
                _currentScopeColor = sPath.size() > 1 ? sPath[sPath.size() - 2].color : 0;
            }

            if (op == CGOpCode::LoopEnd || op == CGOpCode::EndIf || op == CGOpCode::ElseIf || op == CGOpCode::Else) {
                // entering a new scope
                _currentScopeColor = ++_scopeColorCount;

                _scopes.resize(_currentScopeColor + 1);
                _scopes[_currentScopeColor] = _scopes[previousScope];

                // change current scope
                if (op == CGOpCode::LoopEnd || op == CGOpCode::EndIf) {
                    // one more scope level
                    _scopes[_currentScopeColor].push_back(ScopePathElement<Base>(_currentScopeColor, &code));
                } else {
                    // same level but different scope
                    _scopes[_currentScopeColor].back() = ScopePathElement<Base>(_currentScopeColor, &code);
                }

                if (op == CGOpCode::LoopEnd) {
                    _loops.addLoopEndNode(code);
                }
            }

            /**
             * loop arguments
             */
            for (const Argument<Base>& it : code.getArguments()) {
                if (it.getOperation() != nullptr) {
                    markCodeBlockUsed(*it.getOperation());
                }
            }

            if (op == CGOpCode::Index) {
                const IndexOperationNode<Base>& inode = static_cast<const IndexOperationNode<Base>&> (code);
                // indexes that don't depend on a loop start or an index assignment are declared elsewhere
                if (inode.isDefinedLocally()) {
                    _loops.indexes.insert(&inode.getIndex());
                }
            } else if (op == CGOpCode::LoopIndexedIndep || op == CGOpCode::LoopIndexedDep || op == CGOpCode::IndexAssign) {
                IndexPattern* ip;
                if (op == CGOpCode::LoopIndexedDep) {
                    size_t pos = code.getInfo()[0];
                    ip = _loops.dependentIndexPatterns[pos];
                } else if (op == CGOpCode::LoopIndexedIndep) {
                    size_t pos = code.getInfo()[1];
                    ip = _loops.independentIndexPatterns[pos];
                } else {
                    ip = &static_cast<IndexAssignOperationNode<Base>&> (code).getIndexPattern();
                }

                findRandomIndexPatterns(ip, _loops.indexRandomPatterns);

            } else if (op == CGOpCode::DependentRefRhs) {
                CPPADCG_ASSERT_UNKNOWN(code.getInfo().size() == 1);
                size_t depIndex = code.getInfo()[0];

                CPPADCG_ASSERT_UNKNOWN(_dependents->size() > depIndex);
                OperationNode<Base>* depNode = (*_dependents)[depIndex].getOperationNode();
                CPPADCG_ASSERT_UNKNOWN(depNode != nullptr && depNode->getOperationType() != CGOpCode::Inv);

                _varId[code] = _varId[*depNode];
            }

            /**
             * reset scope
             */
            if (previousScope != _currentScopeColor) {
                _currentScopeColor = previousScope;
            }

        } else {
            // been to this node before

            if (op == CGOpCode::Tmp && !code.getInfo().empty()) {
                /**
                 * this node was previously altered to ensure that the 
                 * evaluation of the expression is only performed for the 
                 * required iterations
                 */
                if (_scope[code] == _currentScopeColor) {
                    // outside an if (defined for all iterations)
                    restoreTemporaryVar(code);
                } else {
                    updateTemporaryVarInDiffScopes(code);
                }

            } else if (_scope[code] != _currentScopeColor && op != CGOpCode::LoopIndexedIndep) {
                ScopeIDType oldScope = _scope[code];
                /**
                 * node previously used in a different scope
                 * must make sure it is defined before being used in both
                 * scopes
                 */
                size_t depth = findFirstDifferentScope(oldScope, _currentScopeColor);

                // update the scope where it should be defined
                ScopeIDType newScope;
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
                        _scope[code] = newScope;

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
        CPPADCG_ASSERT_KNOWN(code.getOperationType() != CGOpCode::ArrayCreation, "Not supported yet");
        CPPADCG_ASSERT_KNOWN(code.getOperationType() != CGOpCode::SparseArrayCreation, "Not supported yet");

        /**
         * does this variable require a condition based on indexes?
         */
        std::vector<size_t> iterationRegions;
        OperationNode<Base>* bScopeNewEnd = _scopes[_currentScopeColor].back().end;
        OperationNode<Base>* bScopeOldEnd = _scopes[oldScope].back().end;

        CGOpCode bNewOp = bScopeNewEnd->getOperationType();
        CGOpCode bOldOp = bScopeOldEnd->getOperationType();

        if ((bNewOp == CGOpCode::EndIf || bNewOp == CGOpCode::Else || bNewOp == CGOpCode::ElseIf) &&
                (bOldOp == CGOpCode::EndIf || bOldOp == CGOpCode::Else || bOldOp == CGOpCode::ElseIf)) {
            // used in 2 different if/else branches

            /**
             * determine the iterations which use this temporary variable
             */
            OperationNode<Base>* bScopeNew = bScopeNewEnd->getArguments()[0].getOperation();
            OperationNode<Base>* bScopeOld = bScopeOldEnd->getArguments()[0].getOperation();

            IndexOperationNode<Base>* newIterIndexOp = nullptr;
            iterationRegions = ifBranchIterationRanges(bScopeNew, newIterIndexOp);
            CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);

            IndexOperationNode<Base>* oldIterIndexOp = nullptr;
            std::vector<size_t> oldIterRegions = ifBranchIterationRanges(bScopeOld, oldIterIndexOp);
            combineOverlapingIterationRanges(iterationRegions, oldIterRegions);
            CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);
            CPPADCG_ASSERT_UNKNOWN(newIterIndexOp != nullptr && newIterIndexOp == oldIterIndexOp);

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
                                              ScopeIDType oldScope,
                                              ScopeIDType commonScopeColor) {
        /**
         * clone
         */
        OperationNode<Base>* opClone = cloneNode(tmp);

        /**
         * Create condition
         */
        OperationNode<Base>* tmpDclVar = makeNode(CGOpCode::TmpDcl);
        Argument<Base> tmpArg(*tmpDclVar);

        OperationNode<Base>* cond = makeNode(CGOpCode::IndexCondExpr, iterationRegions,{iterationIndexOp});

        // if
        OperationNode<Base>* ifStart = makeNode(CGOpCode::StartIf, *cond);

        OperationNode<Base>* tmpAssign = makeNode(CGOpCode::LoopIndexedTmp,{tmpArg, *opClone});
        OperationNode<Base>* ifAssign = makeNode(CGOpCode::CondResult,{*ifStart, *tmpAssign});

        // end if
        OperationNode<Base>* endIf = makeNode(CGOpCode::EndIf,{*ifStart, *ifAssign});

        /**
         * Change original variable
         */
        tmp.setOperation(CGOpCode::Tmp,{tmpArg, *endIf});
        tmp.getInfo().resize(1); // used to mark that this node was altered here

        // created new nodes, must adjust vector sizes
        _scope.adjustSize();
        _lastVisit.adjustSize();
        _scope.adjustSize();
        _evaluationOrder.adjustSize();
        _lastUsageOrder.adjustSize();
        _totalUseCount.adjustSize();
        _varId.adjustSize();

        /**
         * add the new scope
         */
        size_t newScopeColor = ++_scopeColorCount;
        _scopes.resize(newScopeColor + 1);
        _scopes[newScopeColor] = _scopes[commonScopeColor];

        // one more scope level
        _scopes[newScopeColor].push_back(ScopePathElement<Base>(newScopeColor, endIf, ifStart));

        // apply scope colors
        _scope[*tmpDclVar] = commonScopeColor;
        _scope[*ifStart] = newScopeColor;
        _scope[*cond] = newScopeColor;
        _scope[*opClone] = newScopeColor;
        _scope[*ifAssign] = newScopeColor;
        _scope[*tmpAssign] = newScopeColor;
        _scope[*endIf] = commonScopeColor;
        _scope[tmp] = commonScopeColor;

        // total usage count
        setTotalUsageCount(*tmpDclVar, 1);
        setTotalUsageCount(*ifStart, 1);
        setTotalUsageCount(*cond, 1);
        setTotalUsageCount(*opClone, 1);
        setTotalUsageCount(*ifAssign, 1);
        setTotalUsageCount(*tmpAssign, 1);
        setTotalUsageCount(*endIf, 1);

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
        if (_scope[code] != _currentScopeColor) {
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
        ScopeIDType oldScope = _scope[code];

        size_t depth = findFirstDifferentScope(oldScope, _currentScopeColor);

        // update the scope where it should be defined
        ScopeIDType newScope = depth == 0 ? 0 : _scopes[_currentScopeColor][depth - 1].color;

        /**
         * does this variable require a condition based on indexes?
         */
        std::vector<size_t> iterationRegions;
        OperationNode<Base>* bScopeNewEnd = _scopes[_currentScopeColor].back().end;
        OperationNode<Base>* endif = code.getArguments()[0].getOperation();
        CPPADCG_ASSERT_UNKNOWN(endif->getOperationType() == CGOpCode::EndIf);
        OperationNode<Base>* bScopeOldEnd = _scopes[_scope[*endif]].back().end;

        CGOpCode bNewOp = bScopeNewEnd->getOperationType();

        if (bNewOp == CGOpCode::EndIf || bNewOp == CGOpCode::Else || bNewOp == CGOpCode::ElseIf) {
            // used in 2 different if/else branches

            /**
             * determine the iterations which use this temporary variable
             */
            OperationNode<Base>* bScopeNew = bScopeNewEnd->getArguments()[0].getOperation();
            OperationNode<Base>* bScopeOld = bScopeOldEnd->getArguments()[0].getOperation();

            IndexOperationNode<Base>* newIterIndexOp = nullptr;
            iterationRegions = ifBranchIterationRanges(bScopeNew, newIterIndexOp);
            CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);

            IndexOperationNode<Base>* oldIterIndexOp = nullptr;
            const std::vector<size_t> oldIterRegions = ifBranchIterationRanges(bScopeOld, oldIterIndexOp);
            combineOverlapingIterationRanges(iterationRegions, oldIterRegions);
            CPPADCG_ASSERT_UNKNOWN(iterationRegions.size() >= 2);
            CPPADCG_ASSERT_UNKNOWN(newIterIndexOp != nullptr && newIterIndexOp == oldIterIndexOp);

            if (iterationRegions.size() == 2 &&
                    (iterationRegions[0] == 0 ||
                    iterationRegions[1] == std::numeric_limits<size_t>::max())) {
                // this temporary variable is used by all iterations
                // there is no need for an 'if'
                restoreTemporaryVar(code);

            } else if (oldIterRegions != iterationRegions) {
                OperationNode<Base>* cond = bScopeOld->getArguments()[0].getOperation();
                CPPADCG_ASSERT_UNKNOWN(cond->getOperationType() == CGOpCode::IndexCondExpr);
                cond->getInfo() = iterationRegions;
            }

        }

        if (oldScope != newScope) {
            _scope[code] = newScope;
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
        CPPADCG_ASSERT_UNKNOWN(tmp.getOperationType() == CGOpCode::Tmp && !tmp.getInfo().empty());

        OperationNode<Base>* endIf = tmp.getArguments()[1].getOperation();
        OperationNode<Base>* ifAssign = endIf->getArguments()[1].getOperation();
        OperationNode<Base>* tmpAssign = ifAssign->getArguments()[1].getOperation();
        OperationNode<Base>* opClone = tmpAssign->getArguments()[1].getOperation();
        tmp.setOperation(opClone->getOperationType(), opClone->getArguments());
        tmp.getInfo() = opClone->getInfo();

        _scope[tmp] = _currentScopeColor;

        /**
         * Must also update the scope of the arguments used by this operation
         */
        const std::vector<Argument<Base> >& args = tmp.getArguments();
        size_t aSize = args.size();
        for (size_t a = 0; a < aSize; a++) {
            updateVarScopeUsage(args[a].getOperation(), _currentScopeColor, _scope[*opClone]);
        }
    }

    inline void restoreTemporaryVar(OperationNode<Base>* tmp,
                                    OperationNode<Base>* opClone) {
        CPPADCG_ASSERT_UNKNOWN(tmp->getOperationType() == CGOpCode::Tmp && !tmp->getInfo().empty());

        tmp->setOperation(opClone->getOperationType(), opClone->getArguments());
        tmp->getInfo() = opClone->getInfo();

        _scope[*tmp] = _currentScopeColor;

        /**
         * Must also update the scope of the arguments used by this operation
         */
        const std::vector<Argument<Base> >& args = tmp->getArguments();
        size_t aSize = args.size();
        for (size_t a = 0; a < aSize; a++) {
            updateVarScopeUsage(args[a].getOperation(), _currentScopeColor, _scope[*opClone]);
        }
    }

    inline void updateVarScopeUsage(OperationNode<Base>* node,
                                    ScopeIDType usageScope,
                                    ScopeIDType oldUsageScope) {
        if (node == nullptr || _scope[*node] == usageScope)
            return;


        ScopeIDType oldScope = _scope[*node];
        ScopeIDType newScope;

        if (oldScope == oldUsageScope) {
            newScope = usageScope;
        } else {
            size_t depth = findFirstDifferentScope(oldScope, usageScope);

            newScope = (depth == 0) ? 0 : _scopes[usageScope][depth - 1].color;
        }

        if (newScope == oldScope)
            return;

        _scope[*node] = newScope;

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

            if (op == CGOpCode::LoopEnd || op == CGOpCode::EndIf || op == CGOpCode::ElseIf || op == CGOpCode::Else) {
                CPPADCG_ASSERT_UNKNOWN(!node->getArguments().empty());

                OperationNode<Base>* beginScopeNode = node->getArguments()[0].getOperation();
                CPPADCG_ASSERT_UNKNOWN(beginScopeNode != nullptr);

                addScopeToVarOrder(_scope[*beginScopeNode], e);
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
     * Attempt to reduce the number of ifs when there are consecutive ifs with
     * the same condition
     */
    inline void optimizeIfs() {
        if (_scopedVariableOrder.size() < 3)
            return; // there has to be at least 2 ifs

        for (size_t scope = 0; scope < _scopedVariableOrder.size(); scope++) {
            std::vector<OperationNode<Base> *>& vorder = _scopedVariableOrder[scope];

            for (long p = vorder.size() - 1; p > 0; p--) {
                OperationNode<Base>* endIf = vorder[p];
                if (endIf->getOperationType() != CGOpCode::EndIf)
                    continue;

                long p1 = p - 1;
                while (p1 >= 0) {
                    if (vorder[p1]->getOperationType() == CGOpCode::TmpDcl) {
                        p1--;
                    } else {
                        break;
                    }
                }
                OperationNode<Base>* endIf1 = vorder[p1];
                if (endIf1->getOperationType() != CGOpCode::EndIf)
                    continue;

                // 2 consecutive ifs
                OperationNode<Base>* startIf = endIf->getArguments()[0].getOperation();
                OperationNode<Base>* startIf1 = endIf1->getArguments()[0].getOperation();
                if (startIf->getOperationType() != CGOpCode::StartIf || startIf1->getOperationType() != CGOpCode::StartIf)
                    continue;

                OperationNode<Base>* cond = startIf->getArguments()[0].getOperation();
                OperationNode<Base>* cond1 = startIf1->getArguments()[0].getOperation();

                CPPADCG_ASSERT_UNKNOWN(cond->getOperationType() == CGOpCode::IndexCondExpr || cond1->getOperationType() == CGOpCode::IndexCondExpr);
                if (cond->getInfo() == cond1->getInfo()) {
                    /**
                     * same condition -> combine the contents into a single if
                     */
                    const std::vector<Argument<Base> >& eArgs = endIf->getArguments();
                    std::vector<Argument<Base> >& eArgs1 = endIf1->getArguments();

                    ScopeIDType ifScope = _scope[*startIf];
                    ScopeIDType ifScope1 = _scope[*startIf1];
                    std::vector<OperationNode<Base> *>& vorderIf = _scopedVariableOrder[ifScope];
                    std::vector<OperationNode<Base> *>& vorderIf1 = _scopedVariableOrder[ifScope1];

                    startNewOperationTreeVisit();

                    // break cycles caused by dependencies on the previous if
                    for (size_t a = 1; a < eArgs.size(); a++) { // exclude the initial startIf
                        CPPADCG_ASSERT_UNKNOWN(eArgs[a].getOperation() != nullptr && eArgs[a].getOperation()->getOperationType() == CGOpCode::CondResult);
                        breakCyclicDependency(eArgs[a].getOperation(), ifScope, endIf1);
                        replaceScope(eArgs[a].getOperation(), ifScope, ifScope1); // update scope
                    }

                    vorderIf1.insert(vorderIf1.end(), vorderIf.begin() + 1, vorderIf.end()); // exclude the initial startIf

                    vorderIf.clear();

                    // update startIf
                    for (size_t a = 1; a < eArgs.size(); a++) { // exclude the initial startIf
                        CPPADCG_ASSERT_UNKNOWN(eArgs[a].getOperation() != nullptr && eArgs[a].getOperation()->getOperationType() == CGOpCode::CondResult);
                        eArgs[a].getOperation()->getArguments()[0] = Argument<Base>(*startIf1);
                    }

                    // update endIf
                    eArgs1.insert(eArgs1.end(), eArgs.begin() + 1, eArgs.end());

                    // replace endIf
                    endIf->setOperation(CGOpCode::Alias,{*endIf1});

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

    inline void replaceScope(OperationNode<Base>* node,
                             ScopeIDType oldScope,
                             ScopeIDType newScope) {
        if (node == nullptr || _scope[*node] != oldScope)
            return;

        _scope[*node] = newScope;

        const std::vector<Argument<Base> >& args = node->getArguments();
        for (size_t a = 0; a < args.size(); a++) {
            replaceScope(args[a].getOperation(), oldScope, newScope);
        }
    }

    /**
     * Removes cyclic dependencies when 'ifs' are merged together.
     * Relative variable order must have already been defined.
     * 
     * @param node the node being visited
     * @param scope the scope where the cyclic dependency could appear (or scopes inside it)
     * @param endIf the dependency to remove
     */
    inline void breakCyclicDependency(OperationNode<Base>* node,
                                      size_t scope,
                                      OperationNode<Base>* endIf) {
        if (node == nullptr || isVisited(*node))
            return;

        markVisited(*node);

        CGOpCode op = node->getOperationType();
        std::vector<Argument<Base> >& args = node->getArguments();

        if (op == CGOpCode::Tmp && args.size() > 1) {
            OperationNode<Base>* arg = args[1].getOperation();
            if (arg == endIf) {
                // a dependency on LoopIndexedTmp could be added but
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
                if (op == CGOpCode::StartIf || op == CGOpCode::LoopStart) {
                    args.erase(args.begin() + a);
                    a--;
                }
            } else {
                breakCyclicDependency(arg, scope, endIf);
            }
        }
    }

    inline bool containedInScope(const OperationNode<Base>& node,
                                 ScopeIDType scope) {
        ScopeIDType nScope = _scope[node];
        if (nScope == scope)
            return true;

        return _scopes[nScope].size() >= _scopes[scope].size() &&
                _scopes[nScope][_scopes[scope].size() - 1].color == scope;
    }

    inline static bool containsArgument(const OperationNode<Base>& node,
                                        const OperationNode<Base>& arg) {
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
     * 
     **********************************************************************/
    virtual void checkVariableCreation(OperationNode<Base>& code) {
        const std::vector<Argument<Base> >& args = code.getArguments();

        size_t aSize = args.size();
        for (size_t argIndex = 0; argIndex < aSize; argIndex++) {
            if (args[argIndex].getOperation() == nullptr) {
                continue;
            }

            OperationNode<Base>& arg = *args[argIndex].getOperation();
            CGOpCode aType = arg.getOperationType();

            if (!isVisited(arg)) {
                // dependencies not visited yet
                checkVariableCreation(arg);

                if (aType == CGOpCode::LoopEnd || aType == CGOpCode::ElseIf || aType == CGOpCode::Else || aType == CGOpCode::EndIf) {
                    if (_varId[arg] == 0) {
                        // ID value is not really used but must be non-zero
                        _varId[arg] = std::numeric_limits<size_t>::max();
                    }
                } else if (aType == CGOpCode::AtomicForward || aType == CGOpCode::AtomicReverse) {
                    /**
                     * Save atomic function related information
                     */
                    CPPADCG_ASSERT_UNKNOWN(arg.getArguments().size() > 1);
                    CPPADCG_ASSERT_UNKNOWN(arg.getInfo().size() > 1);
                    size_t id = arg.getInfo()[0];

                    size_t pos;
                    const std::string& atomicName = _atomicFunctions.at(id)->afun_name();
                    std::map<std::string, size_t>::const_iterator itName2Idx;
                    itName2Idx = _atomicFunctionName2Index.find(atomicName);

                    if (itName2Idx == _atomicFunctionName2Index.end()) {
                        pos = _atomicFunctionsOrder->size();
                        _atomicFunctionsOrder->push_back(atomicName);
                        _atomicFunctionName2Index[atomicName] = pos;
                        _atomicFunctionsMaxForward.push_back(-1);
                        _atomicFunctionsMaxReverse.push_back(-1);
                    } else {
                        pos = itName2Idx->second;
                    }

                    if (aType == CGOpCode::AtomicForward) {
                        int p = arg.getInfo()[2];
                        _atomicFunctionsMaxForward[pos] = std::max(_atomicFunctionsMaxForward[pos], p);
                    } else {
                        int p = arg.getInfo()[1];
                        _atomicFunctionsMaxReverse[pos] = std::max(_atomicFunctionsMaxReverse[pos], p);
                    }
                }

                /**
                 * make sure new temporary variables are NOT created for
                 * the independent variables and that a dependency did
                 * not use it first
                 */
                if (_varId[arg] == 0 || !isIndependent(arg)) {
                    if (aType == CGOpCode::LoopIndexedIndep) {
                        // ID value not really used but must be non-zero
                        _varId[arg] = std::numeric_limits<size_t>::max();
                    } else if (aType == CGOpCode::Alias) {
                        continue; // should never be added to the evaluation queue
                    } else if (aType == CGOpCode::Tmp) {
                        _varId[arg] = std::numeric_limits<size_t>::max();
                    } else if (aType == CGOpCode::LoopStart ||
                            aType == CGOpCode::LoopEnd ||
                            aType == CGOpCode::StartIf ||
                            aType == CGOpCode::ElseIf ||
                            aType == CGOpCode::Else ||
                            aType == CGOpCode::EndIf) {
                        /**
                         * Operation that mark a change in variable scope
                         * are always added
                         */
                        addToEvaluationQueue(arg);
                        if (_varId[arg] == 0) {
                            // ID value is not really used but must be non-zero
                            _varId[arg] = std::numeric_limits<size_t>::max();
                        }
                    } else if (aType == CGOpCode::Pri) {
                        addToEvaluationQueue(arg);
                        if (_varId[arg] == 0) {
                            // ID value is not really used but must be non-zero
                            _varId[arg] = std::numeric_limits<size_t>::max();
                        }
                    } else if (aType == CGOpCode::TmpDcl) {
                        addToEvaluationQueue(arg);

                        _varId[arg] = _idCount;
                        _idCount++;

                    } else if (_lang->createsNewVariable(arg, getTotalUsageCount(arg)) ||
                            _lang->requiresVariableArgument(code.getOperationType(), argIndex)) {

                        addToEvaluationQueue(arg);

                        if (_varId[arg] == 0) {
                            if (aType == CGOpCode::AtomicForward || aType == CGOpCode::AtomicReverse) {
                                _varId[arg] = _idAtomicCount;
                                _idAtomicCount++;
                            } else if (aType == CGOpCode::LoopIndexedDep || aType == CGOpCode::LoopIndexedTmp) {
                                // ID value not really used but must be non-zero
                                _varId[arg] = std::numeric_limits<size_t>::max();
                            } else if (aType == CGOpCode::ArrayCreation) {
                                // a temporary array
                                size_t arraySize = arg.getArguments().size();
                                _varId[arg] = _idArrayCount;
                                _idArrayCount += arraySize;
                            } else if (aType == CGOpCode::SparseArrayCreation) {
                                // a temporary array
                                size_t nnz = arg.getArguments().size();
                                _varId[arg] = _idSparseArrayCount;
                                _idSparseArrayCount += nnz;
                            } else {
                                // a single temporary variable
                                _varId[arg] = _idCount;
                                _idCount++;
                            }
                        }
                    }

                }
            }

            markVisited(arg);

        }

    }

    inline void addToEvaluationQueue(OperationNode<Base>& arg) {
        ScopeIDType scope = _scope[arg];
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

    inline void reduceTemporaryVariables(CppAD::vector<CG<Base> >& dependent) {
        using CppAD::vector;

        reorderOperations(dependent);

        /**
         * determine the last line where each temporary variable is used
         */
        startNewOperationTreeVisit();

        for (size_t i = 0; i < dependent.size(); i++) {
            OperationNode<Base>* node = dependent[i].getOperationNode();
            if (node != nullptr) {
                if (!isVisited(*node)) {
                    // dependencies not visited yet
                    determineLastTempVarUsage(*node);
                }
                markVisited(*node);
            }
        }

        // where temporary variables can be released
        vector<std::vector<OperationNode<Base>* > > tempVarRelease(_variableOrder.size());
        for (size_t i = 0; i < _variableOrder.size(); i++) {
            OperationNode<Base>* var = _variableOrder[i];
            if (isTemporary(*var) || isTemporaryArray(*var) || isTemporarySparseArray(*var)) {
                size_t releaseLocation = getLastUsageEvaluationOrder(*var) - 1;
                tempVarRelease[releaseLocation].push_back(var);
            }
        }


        /**
         * Redefine temporary variable IDs
         */
        std::vector<size_t> freedVariables; // variable IDs no longer in use
        _idCount = _minTemporaryVarID;
        ArrayIdCompresser<Base> arrayComp(_varId, _idArrayCount);
        ArrayIdCompresser<Base> sparseArrayComp(_varId, _idSparseArrayCount);

        for (size_t i = 0; i < _variableOrder.size(); i++) {
            OperationNode<Base>& var = *_variableOrder[i];

            const std::vector<OperationNode<Base>* >& released = tempVarRelease[i];
            for (size_t r = 0; r < released.size(); r++) {
                if (isTemporary(*released[r])) {
                    freedVariables.push_back(_varId[*released[r]]);
                } else if (isTemporaryArray(*released[r])) {
                    arrayComp.addFreeArraySpace(*released[r]);
                } else if (isTemporarySparseArray(*released[r])) {
                    sparseArrayComp.addFreeArraySpace(*released[r]);
                }
            }

            if (isTemporary(var)) {
                // a single temporary variable
                if (freedVariables.empty()) {
                    _varId[var] = _idCount;
                    _idCount++;
                } else {
                    size_t id = freedVariables.back();
                    freedVariables.pop_back();
                    _varId[var] = id;
                }
            } else if (isTemporaryArray(var)) {
                // a temporary array
                size_t arrayStart = arrayComp.reserveArraySpace(var);
                _varId[var] = arrayStart + 1;
            } else if (isTemporarySparseArray(var)) {
                // a temporary array
                size_t arrayStart = sparseArrayComp.reserveArraySpace(var);
                _varId[var] = arrayStart + 1;
            }

        }

        _idArrayCount = arrayComp.getIdCount();
        _idSparseArrayCount = sparseArrayComp.getIdCount();
    }

    /**
     * Change operation order so that the total number of temporary variables is
     * reduced.
     * @param dependent The vector of dependent variable values
     */
    inline void reorderOperations(CppAD::vector<CG<Base> >& dependent) {
        // determine the location of the last temporary variable used for each dependent
        startNewOperationTreeVisit();

        // normal dependent nodes
        for (size_t i = 0; i < dependent.size(); ++i) {
            OperationNode<Base>* node = dependent[i].getOperationNode();
            if (node != nullptr) {
                reorderOperation(*node);
            }
        }

        // dependent nodes defined inside loops
        for (const LoopEndOperationNode<Base>* endNode : _loops.endNodes) {
            const std::vector<Argument<Base> >& args = endNode->getArguments();
            for (size_t i = 1; i < args.size(); ++i) {
                CPPADCG_ASSERT_UNKNOWN(args[i].getOperation() != nullptr);
                // TODO: also consider CGOpCode::LoopIndexedDep inside a CGOpCode::endIf
                if (args[i].getOperation()->getOperationType() == CGOpCode::LoopIndexedDep) {
                    reorderOperation(*args[i].getOperation());
                }
            }
        }
    }

    inline void reorderOperation(OperationNode<Base>& node) {
        /**
         * determine the location of the last temporary variable
         */
        size_t depPos = getEvaluationOrder(node);
        size_t lastTmpPos = depPos;
        if (!isVisited(node)) {
            // dependencies not visited yet
            lastTmpPos = findLastTemporaryLocation(node);
        }
        markVisited(node);

        /**
         * move dependent if beneficial
         */
        if (lastTmpPos == depPos || lastTmpPos + 1 == depPos) {
            return; // should not change location of the evaluation of this dependent
        }

        // should only move if there are temporaries which could use other temporaries released by the dependent
        bool foundTemporaries = false;
        size_t newPos;
        for (size_t l = lastTmpPos + 1; l < depPos; ++l) {
            const auto* n = _variableOrder[l - 1];
            if (isTemporary(*n) || isTemporaryArray(*n) || isTemporarySparseArray(*n)) {
                foundTemporaries = true;
                newPos = l;
                break;
            } else {
                /**
                 * Must be in same scope
                 */
                CGOpCode op = n->getOperationType();
                if (op == CGOpCode::StartIf) {
                    // must not change scope (find the end of this conditional statement)
                    ++l;
                    while (l < depPos) {
                        const auto* node2 = _variableOrder[l - 1];
                        if (node2->getOperationType() == CGOpCode::EndIf &&
                                node2->getArguments()[0].getOperation() == n) {
                            break; // found the end (returned to the same scope)
                        }
                        ++l;
                    }

                } else if (op == CGOpCode::LoopStart) {
                    // must not change scope (find the end of this loop statement)
                    ++l;
                    while (l < depPos) {
                        const auto* node2 = _variableOrder[l - 1];
                        if (node2->getOperationType() == CGOpCode::LoopEnd &&
                                &static_cast<const LoopEndOperationNode<Base>*> (node2)->getLoopStart() == n) {
                            break; // found the end (returned to the same scope)
                        }
                        ++l;
                    }
                }
            }
        }

        if (foundTemporaries) {
            // move variables
            repositionEvaluationQueue(depPos, newPos);
        }
    }

    /**
     * Determine the highest location in the evaluation queue of temporary
     * variables used by an operation node in the same scope.
     * @return the highest location of the temporary variables or 
     *         the location of node itself if it doesn't use any temporary 
     *         variable (in the same scope)
     */
    inline size_t findLastTemporaryLocation(OperationNode<Base>& node) {
        size_t depOrder = getEvaluationOrder(node);
        size_t maxTmpOrder = 0; // lowest possible value is 1
        for (const Argument<Base>& it : node) {
            if (it.getOperation() != nullptr) {
                OperationNode<Base>& arg = *it.getOperation();
                CGOpCode aOp = arg.getOperationType();
                if (aOp == CGOpCode::LoopEnd || aOp == CGOpCode::EndIf || aOp == CGOpCode::ElseIf || aOp == CGOpCode::Else) {
                    continue; //should not move variables to a different scope
                }

                if (aOp == CGOpCode::Index) {
                    size_t iorder = getEvaluationOrder(static_cast<IndexOperationNode<Base>&> (arg).getIndexCreationNode());
                    if (iorder > maxTmpOrder)
                        maxTmpOrder = iorder;
                } else if (getEvaluationOrder(arg) == depOrder) {
                    // dependencies not visited yet
                    size_t orderNew = findLastTemporaryLocation(arg);
                    if (orderNew > maxTmpOrder)
                        maxTmpOrder = orderNew;
                } else {
                    // no need to visit dependencies
                    if (getEvaluationOrder(arg) > maxTmpOrder)
                        maxTmpOrder = getEvaluationOrder(arg);
                }
            }
        }

        return maxTmpOrder == 0 ? depOrder : maxTmpOrder;
    }

    inline void repositionEvaluationQueue(size_t fromPos, size_t toPos) {
        // Warning: there is an offset of 1 between the evaluation order saved 
        // in the node and the actual location in the _variableOrder
        CPPADCG_ASSERT_UNKNOWN(fromPos > toPos);
        OperationNode<Base>* node = _variableOrder[fromPos - 1]; // node to be moved

        // move variables in between the order change
        for (size_t l = fromPos - 1; l > toPos - 1; --l) {
            _variableOrder[l] = _variableOrder[l - 1];
            updateEvaluationQueueOrder(*_variableOrder[l], l + 1);
        }

        _variableOrder[toPos - 1] = node;
        updateEvaluationQueueOrder(*node, toPos);
    }

    /**
     * Determines when each temporary variable is last used in the
     * evaluation order
     * 
     * @param node The current node for which the number of usages is to be to determined
     */
    inline void determineLastTempVarUsage(OperationNode<Base>& node) {
        CGOpCode op = node.getOperationType();

        if (op == CGOpCode::LoopEnd) {
            LoopEndOperationNode<Base>& loopEnd = static_cast<LoopEndOperationNode<Base>&> (node);
            _loops.depth++;
            _loops.outerVars.resize(_loops.depth + 1);
            _loops.startEvalOrder.push_back(getEvaluationOrder(loopEnd.getLoopStart()));

        } else if (op == CGOpCode::LoopStart) {
            _loops.depth--; // leaving the current loop
        }

        /**
         * count variable usage
         */
        for (const Argument<Base>& it : node.getArguments()) {
            if (it.getOperation() != nullptr) {
                OperationNode<Base>& arg = *it.getOperation();

                if (!isVisited(arg)) {
                    // dependencies not visited yet
                    determineLastTempVarUsage(arg);
                }

                markVisited(arg);

                size_t order = getEvaluationOrder(node);
                OperationNode<Base>* aa = getOperationFromAlias(arg); // follow alias!
                if (aa != nullptr) {
                    if (getLastUsageEvaluationOrder(*aa) < order) {
                        setLastUsageEvaluationOrder(*aa, order);
                    }

                    if (_loops.depth >= 0 &&
                            getEvaluationOrder(*aa) < _loops.startEvalOrder[_loops.depth] &&
                            isTemporary(*aa)) {
                        // outer variable used inside the loop
                        _loops.outerVars[_loops.depth].insert(aa);
                    }
                }
            }
        }

        if (op == CGOpCode::LoopEnd) {
            /**
             * temporary variables from outside the loop which are used
             * within the loop cannot be overwritten inside that loop
             */
            size_t order = getEvaluationOrder(node);

            const std::set<OperationNode<Base>*>& outerLoopUsages = _loops.outerVars.back();
            for (OperationNode<Base>* outerVar : outerLoopUsages) {
                OperationNode<Base>* aa = getOperationFromAlias(*outerVar); // follow alias!
                if (aa != nullptr && getLastUsageEvaluationOrder(*aa) < order)
                    setLastUsageEvaluationOrder(*aa, order);
            }

            _loops.depth--;
            _loops.outerVars.pop_back();
            _loops.startEvalOrder.pop_back();

        } else if (op == CGOpCode::LoopStart) {
            _loops.depth++; // coming back to the loop
        }

    }

    /**
     * Defines the evaluation order for the code fragments that do not
     * create variables (right hand side variables)
     * @param code The operation just added to the evaluation order
     */
    inline void dependentAdded2EvaluationQueue(OperationNode<Base>& node) {
        for (const Argument<Base>& a : node.getArguments()) {
            if (a.getOperation() != nullptr) {
                OperationNode<Base>& arg = *a.getOperation();
                if (getEvaluationOrder(arg) == 0) {
                    setEvaluationOrder(arg, getEvaluationOrder(node));
                    dependentAdded2EvaluationQueue(arg);
                }
            }
        }
    }

    inline void updateEvaluationQueueOrder(OperationNode<Base>& node,
                                           size_t newEvalOrder) {
        size_t oldEvalOrder = getEvaluationOrder(node);

        setEvaluationOrder(node, newEvalOrder);

        for (const Argument<Base>& a : node.getArguments()) {
            if (a.getOperation() != nullptr) {
                OperationNode<Base>& arg = *a.getOperation();
                if (getEvaluationOrder(arg) == oldEvalOrder)
                    updateEvaluationQueueOrder(arg, newEvalOrder);
            }
        }
    }

    inline bool isIndependent(const OperationNode<Base>& arg) const {
        return arg.getOperationType() == CGOpCode::Inv;
    }

    inline bool isTemporary(const OperationNode<Base>& arg) const {
        CGOpCode op = arg.getOperationType();
        return op != CGOpCode::ArrayCreation && // classified as TemporaryArray
                op != CGOpCode::SparseArrayCreation && // classified as TemporarySparseArray
                op != CGOpCode::AtomicForward &&
                op != CGOpCode::AtomicReverse &&
                op != CGOpCode::LoopStart &&
                op != CGOpCode::LoopEnd &&
                op != CGOpCode::StartIf &&
                op != CGOpCode::ElseIf &&
                op != CGOpCode::Else &&
                op != CGOpCode::EndIf &&
                op != CGOpCode::LoopIndexedDep &&
                op != CGOpCode::LoopIndexedIndep &&
                op != CGOpCode::LoopIndexedTmp && // not considered as a temporary (the temporary is CGTmpDclOp)
                op != CGOpCode::Index &&
                op != CGOpCode::IndexAssign &&
                op != CGOpCode::Tmp && // not considered as a temporary (the temporary is CGTmpDclOp)
                _varId[arg] >= _minTemporaryVarID;
    }

    inline static bool isTemporaryArray(const OperationNode<Base>& arg) {
        return arg.getOperationType() == CGOpCode::ArrayCreation;
    }

    inline static bool isTemporarySparseArray(const OperationNode<Base>& arg) {
        return arg.getOperationType() == CGOpCode::SparseArrayCreation;
    }

    inline static OperationNode<Base>* getOperationFromAlias(OperationNode<Base>& alias) {
        if (alias.getOperationType() != CGOpCode::Alias) {
            return &alias;
        } else {
            OperationNode<Base>* aa = &alias;
            do {
                CPPADCG_ASSERT_UNKNOWN(aa->getArguments().size() == 1);
                aa = aa->getArguments()[0].getOperation();
            } while (aa != nullptr && aa->getOperationType() == CGOpCode::Alias);
            return aa;
        }
    }

    inline size_t getEvaluationOrder(const OperationNode<Base>& node) const {
        return _evaluationOrder[node];
    }

    inline void setEvaluationOrder(OperationNode<Base>& node,
                                   size_t order) {
        CPPADCG_ASSERT_UNKNOWN(order <= _variableOrder.size());
        _evaluationOrder[node] = order;
    }

    inline size_t getLastUsageEvaluationOrder(const OperationNode<Base>& node) const {
        return _lastUsageOrder[node];
    }

    inline void setLastUsageEvaluationOrder(const OperationNode<Base>& node,
                                            size_t last) {
        CPPADCG_ASSERT_UNKNOWN(last <= _variableOrder.size()); // _lastUsageOrder[node] = 0  means that it was never used
        _lastUsageOrder[node] = last;

        CGOpCode op = node.getOperationType();
        if (op == CGOpCode::ArrayElement) {
            OperationNode<Base>* array = node.getArguments()[0].getOperation();
            CPPADCG_ASSERT_UNKNOWN(array->getOperationType() == CGOpCode::ArrayCreation);
            if (getLastUsageEvaluationOrder(*array) < last) {
                setLastUsageEvaluationOrder(*array, last);
            }
        } else if (op == CGOpCode::Tmp) {
            OperationNode<Base>* declr = node.getArguments()[0].getOperation();
            CPPADCG_ASSERT_UNKNOWN(declr->getOperationType() == CGOpCode::TmpDcl);
            if (getLastUsageEvaluationOrder(*declr) < last) {
                setLastUsageEvaluationOrder(*declr, last);
            }
        }
    }

    /**
     * Provides the total number of times the result of an operation node is
     * being used as an argument for another operation.
     * @return the total usage count
     */
    inline size_t getTotalUsageCount(const OperationNode<Base>& node) const {
        return _totalUseCount[node];
    }

    inline void setTotalUsageCount(const OperationNode<Base>& node,
                                   size_t cout) {
        _totalUseCount[node] = cout;
    }

    inline void increaseTotalUsageCount(const OperationNode<Base>& node) {
        _totalUseCount[node]++;
    }

    inline void resetManagedNodes() {
        _variableOrder.clear();
        _scopedVariableOrder.resize(1);
        _scopedVariableOrder[0].clear();
        _evaluationOrder.fill(0);
        _lastUsageOrder.fill(0);
        _totalUseCount.fill(0);
        _varId.fill(0);
        _scope.fill(0);
    }

    /**************************************************************************
     *                       Graph management functions
     *************************************************************************/

    inline void findPaths(SourceCodePath& path2node,
                          OperationNode<Base>& code,
                          std::vector<SourceCodePath>& found,
                          size_t max);

    static inline std::vector<SourceCodePath> findPathsFromNode(const std::vector<SourceCodePath> nodePaths,
                                                                OperationNode<Base>& node);

    /**************************************************************************
     *                        Operation graph manipulation
     *************************************************************************/
    /**
     * Solves an expression (e.g. f(x, y) == 0) for a given variable (e.g. x)
     * The variable can appear only once in the expression.
     * This is also known as isolation.
     *
     * @param path  The path from the equation residual to the variable
     * @return  The expression for the variable
     * @throws CGException if it is not possible to solve the expression
     */
    inline CG<Base> solveFor(const SourceCodePath& path);

    /**
     * Reduces the number of occurrences of a variable in an equation.
     * For instance:
     *  f(x,y) = x + y + x
     * could become
     *  f(x,y) = 2 * x + y
     *
     * @param expression  The original expression (f(x, y))
     * @param path1  A path from the equation residual to where the variable
     *               is used.
     * @param path2  A different path from the equation residual to where the
     *               variable is used.
     * @return  The new expression for the equation
     * @throws CGException if it is not possible to combine the multiple
     *                     occurrences of the variable
     */
    inline CG<Base> collectVariable(OperationNode<Base>& expression,
                                    const SourceCodePath& path1,
                                    const SourceCodePath& path2,
                                    size_t bifPos);

    inline CG<Base> collectVariableAddSub(const SourceCodePath& pathLeft,
                                          const SourceCodePath& pathRight);

    inline bool isCollectableVariableAddSub(const SourceCodePath& pathLeft,
                                            const SourceCodePath& pathRight,
                                            bool throwEx);

    inline bool isSolvable(const SourceCodePath& path) const;

    /**************************************************************************
     *                     Loop related structure/methods
     *************************************************************************/
    struct LoopData {
        // maps the loop ids of the loop atomic functions
        std::map<size_t, LoopModel<Base>*> loopModels;
        std::vector<LoopEndOperationNode<Base>*> endNodes;
        // the used indexes
        std::set<const OperationNode<Base>*> indexes;
        // the used random index patterns
        std::set<RandomIndexPattern*> indexRandomPatterns;
        //
        std::vector<IndexPattern*> dependentIndexPatterns;
        std::vector<const IndexPattern*> dependentIndexPatternManaged; // garbage collection
        std::vector<IndexPattern*> independentIndexPatterns;
        // variables used inside a loop which are assigned outside (for different loop depths)
        std::vector<std::set<OperationNode<Base>*> > outerVars;
        // the current loop depth (-1 means no loop)
        int depth;
        // the evaluation order of the loop start for each loop depth
        std::vector<size_t> startEvalOrder;

        inline LoopData() :
            depth(-1) {
        }

        inline void prepare4NewSourceGen();

        inline void reset();

        /**
         * Provides the name used by a loop atomic function with a given ID.
         * 
         * @param id the atomic function ID.
         * @return a pointer to the atomic loop function name if it was
         *         registered, nullptr otherwise
         */
        inline const std::string* getLoopName(size_t id) const;

        inline void registerModel(LoopModel<Base>& loop);

        inline LoopModel<Base>* getLoop(size_t loopId) const;

        size_t addDependentIndexPattern(IndexPattern& jacPattern);

        void manageDependentIndexPattern(const IndexPattern* pattern);

        size_t addIndependentIndexPattern(IndexPattern& pattern, size_t hint);

        void addLoopEndNode(OperationNode<Base>& node);
    };

    /**************************************************************************
     *                                friends
     *************************************************************************/
    friend class CG<Base>;
    friend class CGAbstractAtomicFun<Base>;
    friend class BaseAbstractAtomicFun<Base>;
    friend class LoopModel<Base>;

};

} // END cg namespace
} // END CppAD namespace

#endif
