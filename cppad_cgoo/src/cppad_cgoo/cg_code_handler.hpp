#ifndef CPPAD_CG_CODE_HANDLER_INCLUDED
#define	CPPAD_CG_CODE_HANDLER_INCLUDED
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
    class CG;

    /**
     * Helper class to analyze the operation graph and generate source code
     * for several languages
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CodeHandler {
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
    public:

        CodeHandler(size_t varCount = 50) :
            _idCount(0),
            _used(false),
            _reuseIDs(true),
            _lang(NULL),
            _minTemporaryVarID(0) {
            _codeBlocks.reserve(varCount);
            _variableOrder.reserve(1 + varCount / 3);
        }

        void setReuseVariableIDs(bool reuse) {
            _reuseIDs = reuse;
        }

        bool isReuseVariableIDs() const {
            return _reuseIDs;
        }

        void makeVariables(std::vector<CG<Base> >& variables) {
            for (typename std::vector<CG<Base> >::iterator it = variables.begin(); it != variables.end(); ++it) {
                makeVariable(*it);
            }
        }

        void makeVariable(CG<Base>& variable) {
            _independentVariables.push_back(new SourceCodeFragment<Base > (CGInvOp));
            variable.makeVariable(*this, _independentVariables.back());
        }

        virtual size_t getMaximumVariableID() const {
            return _idCount;
        }

        /**
         * Creates the source code from the operations registered so far.
         * 
         * \param out The output stream where the source code is to be printed.
         * \param lang The targeted language.
         * \param dependent The dependent variables for which the source code
         *                  should be generated. By defining this vector the 
         *                  number of operations in the source code can be 
         *                  reduced and thus providing a more optimized code.
         * \param nameGen Provides the rules for variable name creation.
         */
        virtual void generateCode(std::ostream& out, CppAD::Language<Base>& lang, std::vector<CG<Base> >& dependent, VariableNameGenerator<Base>& nameGen) {
            _lang = &lang;
            _idCount = 0;

            if (_used) {
                _variableOrder.clear();

                for (typename std::vector<SourceCodeFragment<Base> *>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                    SourceCodeFragment<Base>* block = *it;
                    block->total_use_count_ = 0;
                    block->var_id_ = 0;
                }
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
                    code.use_count_++;
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

        virtual ~CodeHandler() {
            reset();
        }

    protected:

        virtual void manageSourceCodeBlock(SourceCodeFragment<Base>* code) {
            assert(std::find(_codeBlocks.begin(), _codeBlocks.end(), code) == _codeBlocks.end());

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

                    arg.use_count_++;
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
                    code.use_count_++;
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

                    arg.use_count_++;

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
         * \param code The operation just added to the evaluation order
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

    private:

        CodeHandler(const CodeHandler&); // not implemented

        CodeHandler& operator=(const CodeHandler&); // not implemented

        friend class CG<Base>;

    };

}
#endif

