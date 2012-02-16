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

    template<class Base>
    class CodeHandler {
    protected:
        size_t _idCount;
        // the variable IDs that were used during source code generation
        std::set<size_t> _usedVariables;
        // variable IDs no longer in use
        std::set<size_t> _freedVariables;
        // all the source code blocks created with the CG<Base> objects
        std::vector<CodeBlock*> _codeBlocks;
        //
        std::string _spaces;
        // references to existing variables
        std::map<size_t, std::vector<CG<Base>* > > proxies_;
        // pure variables (not references to other variables)
        std::map<size_t, CG<Base>* > variables_;
        // a flag indicating if this handler was previously used to generate code
        bool _used;
        // a flag indicating whether or not to reuse the IDs of destroyed variables
        bool _reuseIDs;
    public:

        CodeHandler(size_t varCount = 50, size_t spaces = 3) :
            _idCount(0),
            _used(false),
            _reuseIDs(true),
            _spaces(spaces, ' ') {
            _codeBlocks.reserve(varCount);
        }

        void setReuseVariableIDs(bool reuse) {
            _reuseIDs = reuse;
        }

        bool isReuseVariableIDs() const {
            return _reuseIDs;
        }

        void makeVariables(std::vector<CG<Base> >& variables) {
            for (typename std::vector<CG<Base> >::iterator it = variables.begin(); it != variables.end(); ++it) {
                it->makeVariable(*this);
            }
        }

        virtual CodeBlock* createSourceCodeBlock() {
            size_t id;
            if (_reuseIDs && !_freedVariables.empty()) {
                // reuse an existing ID
                id = *_freedVariables.begin();
                _freedVariables.erase(_freedVariables.begin());
            } else {
                // create a new variable ID
                id = ++_idCount;
            }
            CPPADCG_DEBUG_VARIABLE_CHECKID(id);

            return createSourceCodeBlock(id);
        }

        virtual CodeBlock* createSourceCodeBlock(size_t varID) {
            if (_codeBlocks.capacity() == _codeBlocks.size()) {
                _codeBlocks.reserve((_codeBlocks.size()*3) / 2 + 1);
            }

            _codeBlocks.push_back(new CodeBlock(varID));

            return _codeBlocks.back();
        }

        virtual size_t getMaximumVariableID() const {
            return _idCount;
        }

        virtual std::string createVariableName(const CG<Base>& val) const {
            return createVariableName(val.getVariableID());
        }

        virtual std::string createVariableName(size_t id) const {
            return std::string("var") + toString(id);
        }

        virtual std::string operations(const CG<Base>& val) const {
            assert(val.getCodeHandler() == NULL || val.getCodeHandler() == this);

            return val.isParameter() ? baseToString(val.getParameterValue()) : val.operations();
        }

        virtual void printOperationAssign(CG<Base>& var, const CG<Base> &right, bool newCodeBlock = true) {
            std::string code = var.createVariableName() + " = " + operations(right) + ";\n";

            changeSourceCode(var, code, newCodeBlock, right);
        }

        virtual void printOperationPlusAssign(CG<Base>& var, const CG<Base> &right) {
            std::string code = var.createVariableName() + " += " + operations(right) + ";\n";

            changeSourceCode(var, code, true, var, right);
        }

        virtual void printOperationMinusAssign(CG<Base>& var, const CG<Base> &right) {
            std::string code = var.createVariableName() + " -= " + operations(right) + ";\n";

            changeSourceCode(var, code, true, var, right);
        }

        virtual void printOperationMultAssign(CG<Base>& var, const CG<Base> &right) {
            std::string code = var.createVariableName() + " *= " + operations(right) + ";\n";

            changeSourceCode(var, code, true, var, right);
        }

        virtual void printOperationDivAssign(CG<Base>& var, const CG<Base> &right) {
            std::string code = var.createVariableName() + " /= " + operations(right) + ";\n";

            changeSourceCode(var, code, true, var, right);
        }

        virtual void printConditionalAssignment(enum CompareOp cop,
                                                CG<Base> &result,
                                                const CG<Base> &left,
                                                const CG<Base> &right,
                                                const CG<Base> &trueCase,
                                                const CG<Base> &falseCase) {
            assert(result.isVariable());
            assert(!result.isReference());

            SourceCodeBlock& sc = *result.sourceCode_;
            std::string& code = result.sourceCode_->code;

            std::string strOp = getComparison(cop);

            code = "if(" + operations(left) + strOp + operations(right) + ") {\n";
            printOperationAssign(result, trueCase, false);
            code += "} else {\n";
            printOperationAssign(result, falseCase, false);
            code += "}\n";

            sc.addDependency(left);
            sc.addDependency(right);
        }

        /**
         * Creates the source code from the operations registered so far.
         * 
         * \param out The output stream where the source code is to be printed.
         * \param dependent The dependent variables for which the source code
         *                  should be generated. By defining this vector the 
         *                  number of operations in the source code can be 
         *                  reduced and thus providing a more optimized code.
         */
        virtual void generateCode(std::ostream& out, std::vector<CG<Base> >& dependent) {

            if (_used) {
                _usedVariables.clear();

                for (typename std::vector<CodeBlock*>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                    CodeBlock* block = *it;
                    block->print = false;
                }
            }
            _used = true;

            for (typename std::vector<CG<Base> >::iterator it = dependent.begin(); it != dependent.end(); ++it) {
                CG<Base>& var = *it;
                if (var.isVariable()) {
                    SourceCodeBlock* code = var.getSourceCodeBlock();
                    if (!code->print) {
                        markPrintCodeBlock(out, *code);
                    }
                }
            }

            for (typename std::vector<CodeBlock*>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                CodeBlock* block = *it;
                if (block->print && !block->code.empty()) {
                    out << _spaces << block->code;
                }
            }
        }

        virtual void generateCode(std::ostream& out) {
            if (_used) {
                _usedVariables.clear();
            }
            _used = true;

            for (typename std::vector<CodeBlock*>::const_iterator it = _codeBlocks.begin(); it != _codeBlocks.end(); ++it) {
                CodeBlock* block = *it;
                _usedVariables.insert(block->id);
                if (!block->code.empty()) {
                    out << _spaces << block->code;
                }
            }
        }

        virtual const std::set<size_t>& getUsedVariableIDs() const {
            return _usedVariables;
        }

        virtual std::string baseToString(const Base& value) const {
            std::stringstream str;
            // make sure all digits of floating point values are printed
            str << std::scientific << std::setprecision(std::numeric_limits< Base >::digits10 + 2) << value;
            std::string result;
            str >> result;

            // take out the extra zeros
            size_t end = result.rfind('e') - 1;
            size_t start = end;
            while (result[start] == '0') {
                start--;
            }

            return result.substr(0, start + 1) + result.substr(end + 1);
        }

        virtual std::string getComparison(enum CompareOp op) {
            switch (op) {
                case CompareLt:
                    return "<";

                case CompareLe:
                    return "<=";

                case CompareEq:
                    return "==";

                case CompareGe:
                    return ">=";

                case CompareGt:
                    return ">";

                case CompareNe:
                    return "!=";

                default:
                    CPPAD_ASSERT_UNKNOWN(0);
            }
        }

        inline std::string toString(size_t v) const {
            std::stringstream str;
            str << v;
            std::string result;
            str >> result;
            return result;
        }

        virtual void reset() {
            typename std::vector<CG<Base>* >::reverse_iterator rit;
            typename std::map<size_t, std::vector<CG<Base>* > >::iterator it;
            for (it = proxies_.begin(); it != proxies_.end(); ++it) {
                std::vector<CG<Base>* >& p = it->second;
                for (rit = p.rbegin(); rit != p.rend(); ++rit) {
                    CG<Base>* ref = *rit;
                    ref->makeParameterNoChecks(0);
                }
            }
            proxies_.clear();

            typename std::map<size_t, CG<Base>* >::reverse_iterator it2;
            for (it2 = variables_.rbegin(); it2 != variables_.rend();) {
                typename std::map<size_t, CG<Base>* >::reverse_iterator it3 = it2;
                ++it2;
                it3->second->makeParameterNoChecks(0);
            }
            variables_.clear();

            typename std::vector<CodeBlock*>::iterator itc;
            for (itc = _codeBlocks.begin(); itc != _codeBlocks.end(); ++itc) {
                delete *itc;
            }
            _codeBlocks.clear();

            _usedVariables.clear();
            _freedVariables.clear();

            _idCount = 0;
            _used = false;
        }

        virtual ~CodeHandler() {
            reset();
        }

    protected:

        virtual void markPrintCodeBlock(std::ostream& out, CodeBlock& codeBlock) {
            const std::set<SourceCodeBlock*>& deps = codeBlock.depends;
            if (!deps.empty()) {
                // take care of dependencies
                for (std::set<SourceCodeBlock*>::const_iterator it = deps.begin(); it != deps.end(); ++it) {
                    SourceCodeBlock* scb = *it;
                    if (!scb->print) {
                        markPrintCodeBlock(out, *scb);
                    }
                }
            }

            codeBlock.print = true;
            _usedVariables.insert(codeBlock.id);
        }

        virtual void putVariableReference(const CG<Base>& variable, CG<Base>& reference) {
            std::vector<CG<Base>* >& p = proxies_[variable.getVariableID()];
            assert(std::find(p.begin(), p.end(), &reference) == p.end());
            p.push_back(&reference);
        }

        virtual void removeVariableReference(const CG<Base>& reference) {
            typename std::map<size_t, std::vector<CG<Base>* > >::iterator it = proxies_.find(reference.getVariableID());
            assert(it != proxies_.end()); // must exist in the map
            std::vector<CG<Base>* >& p = it->second;

            if (p.size() == 1) {
                assert(p[0] == &reference);
                proxies_.erase(it);
            } else {
                typename std::vector<CG<Base>* >::iterator it2 = std::find(p.begin(), p.end(), &reference);
                assert(it2 != p.end()); // must exist in the list
                p.erase(it2);
            }
        }

        virtual void releaseReferences(const CG<Base>& variable) {
            typename std::map<size_t, std::vector<CG<Base>* > >::iterator it = proxies_.find(variable.getVariableID());
            if (it == proxies_.end()) {
                return;
            }
            const std::vector<CG<Base>* >& p = it->second;
            assert(!p.empty());

            const size_t refCount = p.size();
            CG<Base>* newVar = p[0];
            newVar->makeVariable(*this); // will reduce by 1 the number of references to variable (or delete p)

            printOperationAssign(*newVar, variable);

            if (refCount > 1) {
                size_t newID = newVar->getVariableID();
                std::vector<CG<Base>* >& p2 = proxies_[newID];
                p2 = p;

                // make the other references to this variable a reference to the 
                // 1st instead
                for (size_t i = 0; i < p2.size(); i++) {
                    p2[i]->referenceTo_ = newVar;
                    p2[i]->id_ = newID;
                }

                proxies_.erase(it); // variable will no longer have any references
            }
        }

        virtual void addPureVariable(CG<Base>& variable) {
            assert(variable.isVariable());
            variables_[variable.getVariableID()] = &variable;
        }

        virtual void removePureVariable(const CG<Base>& variable) {
            CPPADCG_DEBUG_VARIABLE_CHECKID(variable.getVariableID());

            // remove the variable from the variables list
            size_t erased = variables_.erase(variable.getVariableID());
            assert(erased > 0);

            /**
             * allow the references to this variable to keep using its ID
             */
            typename std::map<size_t, std::vector<CG<Base>* > >::iterator it = proxies_.find(variable.getVariableID());
            if (it == proxies_.end()) {
                if (_reuseIDs) {
                    _freedVariables.insert(variable.getVariableID());
                }
                return; // no references
            }
            std::vector<CG<Base>* >& p = it->second;
            assert(!p.empty());

            // make the 1st reference a variable
            CG<Base>* newVar = p[0];
            assert(newVar->getVariableID() == variable.getVariableID());
            newVar->referenceTo_ = NULL; //make it a variable that can keep the current id
            newVar->sourceCode_ = variable.getSourceCodeBlock();
            addPureVariable(*newVar);

            if (p.size() == 1) {
                proxies_.erase(it);
                return;
            } else {
                p.erase(p.begin());

                // make the other references a reference to the 1st instead
                for (size_t i = 0; i < p.size(); i++) {
                    p[i]->referenceTo_ = newVar;
                }
            }
        }

        inline void changeSourceCode(CG<Base>& var, const std::string& code, bool newCodeBlock, const CG<Base>& depend1, const CG<Base>& depend2) {
            changeSourceCode(var, code, newCodeBlock, depend1, &depend2);
        }

        inline void changeSourceCode(CG<Base>& var, const std::string& code, bool newCodeBlock, const CG<Base>& depend1, const CG<Base>* depend2 = NULL) {
            SourceCodeBlock* sc;
            if (newCodeBlock) {
                sc = createSourceCodeBlock(var.id_);
            } else {
                sc = var.sourceCode_;
            }

            sc->code += code;

            sc->addDependency(depend1);

            if (depend2 != NULL) {
                sc->addDependency(*depend2);
            }

            if (newCodeBlock) {
                var.sourceCode_ = sc;
            }
        }

    private:

        CodeHandler() {
            throw 1;
        }

        friend class CG<Base>;

    };

}
#endif

