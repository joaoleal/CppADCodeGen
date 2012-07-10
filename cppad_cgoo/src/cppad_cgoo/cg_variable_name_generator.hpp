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
     * Creates variables names for the source code.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class VariableNameGenerator {
    protected:
        // the lowest variable ID used for the temporary variables
        size_t _minTemporaryID;
        // the highest variable ID used for the temporary variables
        size_t _maxTemporaryID;
    public:
        virtual std::string generateDependent(const CG<Base>& variable, size_t index) = 0;

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& variable) = 0;

        virtual std::string generateTemporary(const SourceCodeFragment<Base>& variable) = 0;

        virtual bool isTemporaryVariablesArray() const = 0;

        virtual const std::string& getTemporaryVariablesArrayName() const = 0;

        inline size_t getMinTemporaryVariableID() const {
            return _minTemporaryID;
        }

        inline size_t getMaxTemporaryVariableID() const {
            return _maxTemporaryID;
        }

    protected:

        inline void setTemporaryVariableID(size_t minTempID, size_t maxTempID) {
            _minTemporaryID = minTempID;
            _maxTemporaryID = maxTempID;
            
            // if
            //  _minTemporaryID == _maxTemporaryID + 1
            // then no temporary variables are being used
            assert(_minTemporaryID <= _maxTemporaryID + 1);
        }

        friend class CodeHandler<Base>;
    };
}

#endif
