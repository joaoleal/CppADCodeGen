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

#include <assert.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>

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
        // variable IDs no longer in use
        std::set<size_t> _freedVariables;
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
        // the maximum number of times a code fragment can be visited during the
        // usage count
        size_t _maxCodeBlockVisit;
        // the language used for source code generation
        Language<Base>* _lang;
    public:

        CodeHandler(size_t varCount = 50) :
            _idCount(0),
            _used(false),
            _reuseIDs(true),
            _maxCodeBlockVisit(2),
            _lang(NULL) {
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
         */
        virtual void generateCode(std::ostream& out, CppAD::Language<Base>& lang, std::vector<CG<Base> >& dependent, VariableNameGenerator<Base>& nameGen) {
            _lang = &lang;
            _idCount = 0;

            if (_used) {
                _variableOrder.clear();

                for (typename std::vector<SourceCodeFragment<Base> *>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                    SourceCodeFragment<Base>* block = *it;
                    block->use_count_ = 0;
                    block->var_id_ = 0;
                }
            }
            _used = true;

            _maxCodeBlockVisit = lang.getMaximumCodeBlockVisit();

            /**
             * the first variable IDs are for the independent variables
             */
            for (typename std::vector<SourceCodeFragment<Base> *>::iterator it = _independentVariables.begin(); it != _independentVariables.end(); ++it) {
                (*it)->setVariableID(++_idCount);
            }

            // determine the number of times each variable is used
            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.getSourceCodeFragment() != NULL) {
                    SourceCodeFragment<Base>& code = *var.getSourceCodeFragment();
                    if (code.use_count_ < _maxCodeBlockVisit) {
                        markCodeBlockUsed(code);
                    }
                }
            }

            // determine the variable creation order
            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.getSourceCodeFragment() != NULL) {
                    SourceCodeFragment<Base>& code = *var.getSourceCodeFragment();
                    if (code.variableID() == 0) {
                        checkVariableCreation(code); // dependencies not visited yet

                        code.setVariableID(++_idCount);
                        _variableOrder.push_back(&code);
                    }
                }
            }

            assert(_idCount == _variableOrder.size() + _independentVariables.size());

            /**
             * Creates the source code for a specific language
             */
            lang.generateSourceCode(out, _independentVariables, dependent, _variableOrder, nameGen);
        }

        virtual void reset() {
            typename std::vector<SourceCodeFragment<Base> *>::iterator itc;
            for (itc = _codeBlocks.begin(); itc != _codeBlocks.end(); ++itc) {
                delete *itc;
            }
            _codeBlocks.clear();
            _independentVariables.clear();
            _freedVariables.clear();
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
            code.use_count_++;

            if (code.use_count_ == _maxCodeBlockVisit) {
                return;
            }

            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;
            for (it = args.begin(); it != args.end(); ++it) {
                if (it->operation() != NULL) {
                    SourceCodeFragment<Base>& arg2 = *it->operation();
                    if (arg2.use_count_ < _maxCodeBlockVisit) {
                        markCodeBlockUsed(arg2);
                    }
                }
            }
        }

        virtual void checkVariableCreation(SourceCodeFragment<Base>& code) {
            const std::vector<Argument<Base> >& args = code.arguments_;

            typename std::vector<Argument<Base> >::const_iterator it;
            for (it = args.begin(); it != args.end(); ++it) {
                if (it->operation() != NULL) {
                    SourceCodeFragment<Base>& arg2 = *it->operation();

                    if (arg2.variableID() == 0) {
                        checkVariableCreation(arg2); // dependencies not visited yet

                        size_t argIndex = it - args.begin();

                        if (_lang->createsNewVariable(arg2) ||
                                _lang->requiresVariableArgument(code.operation(), argIndex)) {
                            _variableOrder.reserve((_variableOrder.size()*3) / 2 + 1);
                            arg2.setVariableID(++_idCount);
                            _variableOrder.push_back(&arg2);
                        }
                    }
                }
            }
        }

    private:

        CodeHandler(); // not implemented

        CodeHandler(const CodeHandler&); // not implemented

        CodeHandler& operator=(const CodeHandler&); // not implemented

        friend class CG<Base>;

    };

}
#endif

