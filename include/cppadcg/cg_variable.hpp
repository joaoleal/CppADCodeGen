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
        return opNode_ != NULL;
    }

    template<class Base>
    inline bool CG<Base>::isParameter() const {
        return opNode_ == NULL;
    }

    template<class Base>
    inline bool CG<Base>::isValueDefined() const {
        return value_ != NULL;
    }

    template<class Base>
    inline const Base& CG<Base>::getValue() const throw (CGException) {
        if (!isValueDefined()) {
            throw CGException("No value defined for this variable");
        }

        return *value_;
    }

    template<class Base>
    inline void CG<Base>::setValue(const Base& b) {
        if (value_ != NULL) {
            *value_ = b;
        } else {
            value_ = new Base(b);
        }
    }

    template<class Base>
    inline bool CG<Base>::IdenticalZero() const throw (CGException) {
        return CppAD::IdenticalZero(getValue());
    }

    template<class Base>
    inline bool CG<Base>::IdenticalOne() const throw (CGException) {
        return CppAD::IdenticalOne(getValue());
    }

    template<class Base>
    inline void CG<Base>::makeParameter(const Base &b) {
        opNode_ = NULL;
        handler_ = NULL;
        setValue(b);
    }

    template<class Base>
    inline void CG<Base>::makeVariable(CodeHandler<Base>& handler, OperationNode<Base>* operation) {
        assert(operation != NULL);
        opNode_ = operation;
        handler_ = &handler;
        delete value_;
        value_ = NULL;

        handler.manageSourceCodeBlock(operation);
    }

    template<class Base>
    inline OperationNode<Base>* CG<Base>::getOperationNode() const {
        return opNode_;
    }

    template<class Base>
    inline Argument<Base> CG<Base>::argument() const {
        if (opNode_ != NULL)
            return Argument<Base > (*opNode_);
        else
            return Argument<Base > (*value_);
    }
}

#endif

