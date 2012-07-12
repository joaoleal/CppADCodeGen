#ifndef CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
#define	CPPAD_CG_VARIABLE_NAME_GENERATOR_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

namespace CppAD {

    /**
     * Function arguments
     */
    typedef struct FuncArgument {
        std::string name;
        bool array;

        inline FuncArgument() {
        }

        inline FuncArgument(const std::string & nam) :
            name(nam),
            array(true) {
        }

        inline FuncArgument(const std::string& nam, bool a) :
            name(nam),
            array(a) {
        }
    } FuncArgument;

    /**
     * Creates variables names for the source code.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class VariableNameGenerator {
    protected:
        std::vector<FuncArgument> _dependent;
        std::vector<FuncArgument> _independent;
        std::vector<FuncArgument> _temporary;
    public:

        virtual const std::vector<FuncArgument>& getDependent() const {
            return _dependent;
        }

        virtual const std::vector<FuncArgument>& getIndependent() const {
            return _independent;
        }

        virtual const std::vector<FuncArgument>& getTemporaryVariables() const {
            return _temporary;
        }

        virtual size_t getMinTemporaryVariableID() const = 0;

        virtual size_t getMaxTemporaryVariableID() const = 0;

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) = 0;

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& variable) = 0;

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) = 0;

        virtual void setTemporaryVariableID(size_t minTempID, size_t maxTempID) = 0;
    };
}

#endif
