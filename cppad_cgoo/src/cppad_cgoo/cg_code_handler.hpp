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

namespace CppAD {

    template<class Base>
    class CG;

    template<class Base>
    class CodeHandler {
    protected:
        std::ostream* _out;
        size_t _idCount;
        std::string _spaces;
        // references to existing variables
        std::map<size_t, std::vector<CG<Base>* > > proxies_;
        // pure variables (not references to other variables)
        std::vector<CG<Base>* > variables_;
    public:

        CodeHandler(std::ostream& out, size_t spaces = 3) {
            _out = &out;
            _idCount = 0;
            _spaces = std::string(spaces, ' ');
        }

        std::ostream* getOutputStream() const {
            return _out;
        }

        void makeVariables(std::vector<CG<Base> >& variables) {
            for (typename std::vector<CG<Base> >::iterator it = variables.begin(); it != variables.end(); ++it) {
                it->makeVariable(*this);
            }
        }

        virtual size_t createID() {
            return ++_idCount;
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

        virtual void printOperationAssign(const std::string& var, const std::string& operations) const {
            (*this->_out) << _spaces << var << " = " << operations << ";\n";
        }

        virtual void printOperationPlusAssign(const std::string& var, const std::string& operations) const {
            (*this->_out) << _spaces << var << " += " << operations << ";\n";
        }

        virtual void printOperationMinusAssign(const std::string& var, const std::string& operations) const {
            (*this->_out) << _spaces << var << " -= " << operations << ";\n";
        }

        virtual void printOperationMultAssign(const std::string& var, const std::string& operations) const {
            (*this->_out) << _spaces << var << " *= " << operations << ";\n";
        }

        virtual void printOperationDivAssign(const std::string& var, const std::string& operations) const {
            (*this->_out) << _spaces << var << " /= " << operations << ";\n";
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

        virtual void printComparison(std::string leftOps, enum CompareOp op, std::string rightOps) {
            (*_out) << leftOps << " ";
            switch (op) {
                case CompareLt:
                    (*_out) << "<";
                    break;

                case CompareLe:
                    (*_out) << "<=";
                    break;

                case CompareEq:
                    (*_out) << "==";
                    break;

                case CompareGe:
                    (*_out) << ">=";
                    break;

                case CompareGt:
                    (*_out) << ">";
                    break;

                case CompareNe:
                    (*_out) << "!=";
                    break;

                default:
                    CPPAD_ASSERT_UNKNOWN(0);
            }

            (*_out) << " " << rightOps;
        }

        inline std::string toString(size_t v) const {
            std::stringstream str;
            str << v;
            std::string result;
            str >> result;
            return result;
        }

        virtual ~CodeHandler() {
            typename std::vector<CG<Base>* >::reverse_iterator rit;
            typename std::map<size_t, std::vector<CG<Base>* > >::iterator it;
            for (it = proxies_.begin(); it != proxies_.end(); ++it) {
                std::vector<CG<Base>* >& p = it->second;
                for (rit = p.rbegin(); rit != p.rend(); ++rit) {
                    CG<Base>* ref = *rit;
                    ref->makeParameterNoChecks(0);
                }
            }
            for (rit = variables_.rbegin(); rit != variables_.rend(); ++rit) {
                CG<Base>* ref = *rit;
                ref->makeParameterNoChecks(0);
            }
        }

    protected:

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

            // variable :: print operation
            std::string name = newVar->createVariableName();
            newVar->printOperationAssig(name, operations(variable));

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
            variables_.push_back(&variable);
        }

        virtual void removePureVariable(const CG<Base>& variable) {
            // remove the variable from the variables list
            typename std::vector<CG<Base>* >::iterator it2 = std::find(variables_.begin(), variables_.end(), &variable);
            assert(it2 != variables_.end());
            variables_.erase(it2);

            /**
             * allow the references to this variable to keep using its ID
             */
            typename std::map<size_t, std::vector<CG<Base>* > >::iterator it = proxies_.find(variable.getVariableID());
            if (it == proxies_.end()) {
                return; // no references
            }
            std::vector<CG<Base>* >& p = it->second;
            assert(!p.empty());

            // make the 1st reference a variable
            CG<Base>* newVar = p[0];
            assert(newVar->getVariableID() == variable.getVariableID());
            newVar->referenceTo_ = NULL; //make it a variable that can keep the current id
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

    private:

        CodeHandler() {
            throw 1;
        }

        friend class CG<Base>;

    };

}
#endif

