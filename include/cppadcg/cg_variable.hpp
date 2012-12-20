#ifndef CPPAD_CG_VARIABLE_INCLUDED
#define CPPAD_CG_VARIABLE_INCLUDED
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
    inline bool CG<Base>::isVariable() const {
        return sourceCode_ != NULL;
    }

    template<class Base>
    inline bool CG<Base>::isParameter() const {
        return value_ != NULL;
    }

    template<class Base>
    inline const Base& CG<Base>::getParameterValue() const throw (CGException) {
        if (!isParameter()) {
            throw CGException("getParameterValue() can only be used for parameters");
        }

        return *value_;
    }

    template<class Base>
    inline bool CG<Base>::IdenticalZero() const throw (CGException) {
        return CppAD::IdenticalZero(getParameterValue());
    }

    template<class Base>
    inline bool CG<Base>::IdenticalOne() const throw (CGException) {
        return CppAD::IdenticalOne(getParameterValue());
    }

    template<class Base>
    inline void CG<Base>::makeParameter(const Base &b) {
        sourceCode_ = NULL;
        handler_ = NULL;
        if (value_ != NULL) {
            *value_ = b;
        } else {
            value_ = new Base(b);
        }
    }

    template<class Base>
    inline void CG<Base>::makeVariable(CodeHandler<Base>& handler, SourceCodeFragment<Base>* operation) {
        assert(operation != NULL);
        sourceCode_ = operation;
        handler_ = &handler;
        delete value_;
        value_ = NULL;
        
        handler.manageSourceCodeBlock(operation);
    }

    template<class Base>
    inline SourceCodeFragment<Base>* CG<Base>::getSourceCodeFragment() const {
        return sourceCode_;
    }

    template<class Base>
    inline Argument<Base> CG<Base>::argument() const {
        if (sourceCode_ != NULL)
            return Argument<Base > (*sourceCode_);
        else
            return Argument<Base > (*value_);
    }
}

#endif

