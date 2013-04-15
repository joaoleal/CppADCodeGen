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
        typedef std::vector<SourceCodePathNode<Base> > SourceCodePath;
    protected:
        // counter used to generate variable IDs
        size_t _idCount;
        // the independent variables
        std::vector<SourceCodeFragment<Base> *> _independentVariables;
        // all the source code blocks created with the CG<Base> objects
        std::vector<SourceCodeFragment<Base> *> _codeBlocks;
        // the order for the variable creation in the source code
        std::vector<SourceCodeFragment<Base> *> _variableOrder;
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
            _idCount(0),
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

        inline void makeVariables(std::vector<CG<Base> >& variables) {
            for (typename std::vector<CG<Base> >::iterator it = variables.begin(); it != variables.end(); ++it) {
                makeVariable(*it);
            }
        }

        inline void makeVariables(std::vector<AD<CG<Base> > >& variables) {
            for (typename std::vector<AD<CG<Base> > >::iterator it = variables.begin(); it != variables.end(); ++it) {
                CG<Base> v;
                makeVariable(v); // make it a codegen variable
                *it = v; // variable[i] id now the same as v
            }
        }

        inline void makeVariable(CG<Base>& variable) {
            _independentVariables.push_back(new SourceCodeFragment<Base > (CGInvOp));
            variable.makeVariable(*this, _independentVariables.back());
        }

        size_t getIndependentVariableSize() const {
            return _independentVariables.size();
        }

        size_t getIndependentVariableIndex(const SourceCodeFragment<Base>& var) const throw (CGException) {
            assert(var.operation_ == CGInvOp);

            typename std::vector<SourceCodeFragment<Base> *>::const_iterator it =
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
        inline std::vector<SourceCodePath> findPaths(SourceCodeFragment<Base>& root,
                                                     SourceCodeFragment<Base>& code,
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
            double beginTime;
            if (_verbose) {
                std::cout << "generating source for '" << jobName << "' ... ";
                std::cout.flush();
                beginTime = system::currentTime();
            }

            _lang = &lang;
            _idCount = 0;

            if (_used) {
                resetCounters();
            }
            _used = true;

            /**
             * the first variable IDs are for the independent variables
             */
            for (typename std::vector<SourceCodeFragment<Base> *>::iterator it = _independentVariables.begin(); it != _independentVariables.end(); ++it) {
                (*it)->setVariableID(++_idCount);
            }

            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                if (it->getSourceCodeFragment() != NULL && it->getSourceCodeFragment()->variableID() == 0) {
                    it->getSourceCodeFragment()->setVariableID(++_idCount);
                }
            }

            _minTemporaryVarID = _idCount + 1;

            // determine the number of times each variable is used
            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.getSourceCodeFragment() != NULL) {
                    SourceCodeFragment<Base>& code = *var.getSourceCodeFragment();
                    markCodeBlockUsed(code);
                }
            }

            // determine the variable creation order
            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.getSourceCodeFragment() != NULL) {
                    SourceCodeFragment<Base>& code = *var.getSourceCodeFragment();
                    if (code.use_count_ == 0) {
                        // dependencies not visited yet
                        checkVariableCreation(code);

                        // make sure new temporary variables are NOT created for
                        // the independent variables and that a dependency did
                        // not use it first
                        if ((code.variableID() == 0 || !isIndependent(code)) && code.use_count_ == 0) {
                            addToEvaluationQueue(code);
                        }
                    }
                    code.increaseUsageCount();
                }
            }

            assert(_idCount == _variableOrder.size() + _independentVariables.size());

            if (_reuseIDs) {
                reduceTemporaryVariables(dependent);
            }

            nameGen.setTemporaryVariableID(_minTemporaryVarID, _idCount);

            /**
             * Creates the source code for a specific language
             */
            LanguageGenerationData<Base> info(_independentVariables, dependent,
                                              _minTemporaryVarID, _variableOrder,
                                              nameGen, _reuseIDs);
            lang.generateSourceCode(out, info);

            if (_verbose) {
                double endTime = system::currentTime();
                std::cout << "done [" << std::fixed << std::setprecision(3)
                        << (endTime - beginTime) << "]" << std::endl;
            }
        }

        virtual void reset() {
            typename std::vector<SourceCodeFragment<Base> *>::iterator itc;
            for (itc = _codeBlocks.begin(); itc != _codeBlocks.end(); ++itc) {
                delete *itc;
            }
            _codeBlocks.clear();
            _independentVariables.clear();
            _idCount = 0;
            _used = false;
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
        inline CG<Base> solveFor(SourceCodeFragment<Base>& expression,
                                 SourceCodeFragment<Base>& code) throw (CGException);

        inline CG<Base> solveFor(const std::vector<SourceCodePathNode<Base> >& path) throw (CGException);

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

        inline void substituteIndependent(SourceCodeFragment<Base>& indep,
                                          SourceCodeFragment<Base>& dep,
                                          bool removeFromIndeps = true) throw (CGException);

        /**
         * Reverts a subtitution of an independent variable that has not been 
         * removed from the list of independents yet.
         * Warning: it does not recover any custom name assigned to the variable.
         * 
         * @param indep The independent variable
         */
        inline void undoSubstituteIndependent(SourceCodeFragment<Base>& indep) throw (CGException);

        /**
         * Finallizes the subtitution of an independent variable by eliminating
         * it from the list of independents. After this operation the variable
         * subtitution cannot be undone.
         * 
         * @param indep The independent variable
         */
        inline void removeIndependent(SourceCodeFragment<Base>& indep) throw (CGException);

        inline virtual ~CodeHandler() {
            reset();
        }

    protected:

        virtual void manageSourceCodeBlock(SourceCodeFragment<Base>* code) {
            //assert(std::find(_codeBlocks.begin(), _codeBlocks.end(), code) == _codeBlocks.end()); // <<< too great of an impact in performance
            if (_codeBlocks.capacity() == _codeBlocks.size()) {
                _codeBlocks.reserve((_codeBlocks.size()*3) / 2 + 1);
            }

            _codeBlocks.push_back(code);
        }

        virtual void markCodeBlockUsed(SourceCodeFragment<Base>& code) {
            code.total_use_count_++;

            if (code.total_use_count_ == 1) {
                // first time this operation is visited

                const std::vector<Argument<Base> >& args = code.arguments_;

                typename std::vector<Argument<Base> >::const_iterator it;
                for (it = args.begin(); it != args.end(); ++it) {
                    if (it->operation() != NULL) {
                        SourceCodeFragment<Base>& arg = *it->operation();
                        markCodeBlockUsed(arg);
                    }
                }
            }
        }

        virtual void checkVariableCreation(SourceCodeFragment<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;

            for (it = args.begin(); it != args.end(); ++it) {
                if (it->operation() != NULL) {
                    SourceCodeFragment<Base>& arg = *it->operation();

                    if (arg.use_count_ == 0) {
                        // dependencies not visited yet
                        checkVariableCreation(arg);

                        // make sure new temporary variables are NOT created for
                        // the independent variables and that a dependency did
                        // not use it first
                        if ((arg.variableID() == 0 || !isIndependent(arg)) && arg.use_count_ == 0) {

                            size_t argIndex = it - args.begin();
                            if (_lang->createsNewVariable(arg) ||
                                    _lang->requiresVariableArgument(code.operation(), argIndex)) {
                                addToEvaluationQueue(arg);
                                if (arg.variableID() == 0) {
                                    arg.setVariableID(++_idCount);
                                }
                            }
                        }
                    }

                    arg.increaseUsageCount();
                }
            }

        }

        inline void addToEvaluationQueue(SourceCodeFragment<Base>& arg) {
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
                if (var.getSourceCodeFragment() != NULL) {
                    SourceCodeFragment<Base>& code = *var.getSourceCodeFragment();
                    if (code.use_count_ == 0) {
                        // dependencies not visited yet
                        determineLastTempVarUsage(code);
                    }
                    code.increaseUsageCount();
                }
            }

            // location of where temporary variables can be released
            std::vector<std::vector<SourceCodeFragment<Base>* > > tempVarRelease(_variableOrder.size());
            for (size_t i = 0; i < _variableOrder.size(); i++) {
                SourceCodeFragment<Base>* var = _variableOrder[i];
                if (isTemporary(*var)) {
                    tempVarRelease[var->getLastUsageEvaluationOrder() - 1].push_back(var);
                }
            }


            /**
             * Redefine temporary variable IDs
             */
            std::vector<size_t> freedVariables; // variable IDs no longer in use
            _idCount = _minTemporaryVarID - 1;

            for (size_t i = 0; i < _variableOrder.size(); i++) {
                SourceCodeFragment<Base>& var = *_variableOrder[i];

                const std::vector<SourceCodeFragment<Base>* >& released = tempVarRelease[i];
                for (size_t r = 0; r < released.size(); r++) {
                    freedVariables.push_back(released[r]->variableID());
                }

                if (isTemporary(var)) {
                    if (freedVariables.empty()) {
                        var.setVariableID(++_idCount);
                    } else {
                        size_t id = freedVariables.back();
                        freedVariables.pop_back();
                        var.setVariableID(id);
                    }
                }
            }
        }

        inline void determineLastTempVarUsage(SourceCodeFragment<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;

            /**
             * count variable usage
             */
            for (it = args.begin(); it != args.end(); ++it) {
                if (it->operation() != NULL) {
                    SourceCodeFragment<Base>& arg = *it->operation();

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
            typename std::vector<SourceCodeFragment<Base> *>::const_iterator it;
            for (it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                SourceCodeFragment<Base>* block = *it;
                block->use_count_ = 0;

            }
        }

        /**
         * Defines the evaluation order for the code fragments that do not
         * create variables
         * @param code The operation just added to the evaluation order
         */
        inline void dependentAdded2EvaluationQueue(SourceCodeFragment<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;

            for (it = args.begin(); it != args.end(); ++it) {
                if (it->operation() != NULL) {
                    SourceCodeFragment<Base>& arg = *it->operation();
                    if (arg.getEvaluationOrder() == 0) {
                        arg.setEvaluationOrder(code.getEvaluationOrder());
                        dependentAdded2EvaluationQueue(arg);
                    }
                }
            }
        }

        inline bool isIndependent(const SourceCodeFragment<Base>& arg) const {
            size_t id = arg.variableID();
            return id > 0 && id <= _independentVariables.size();
        }

        inline bool isTemporary(const SourceCodeFragment<Base>& arg) const {
            return arg.variableID() >= _minTemporaryVarID;
        }

        virtual void resetCounters() {
            _variableOrder.clear();

            for (typename std::vector<SourceCodeFragment<Base> *>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                SourceCodeFragment<Base>* block = *it;
                block->resetHandlerCounters();
            }
        }

        /***********************************************************************
         *                   Graph management functions
         **********************************************************************/

        inline void findPaths(SourceCodePath& path2node,
                              SourceCodeFragment<Base>& code,
                              std::vector<SourceCodePath>& found,
                              size_t max);

        static inline std::vector<SourceCodePath> findPathsFromNode(const std::vector<SourceCodePath> nodePaths,
                                                                    SourceCodeFragment<Base>& node);

    private:

        CodeHandler(const CodeHandler&); // not implemented

        CodeHandler& operator=(const CodeHandler&); // not implemented

        friend class CG<Base>;

    };

}
#endif

