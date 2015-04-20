#ifndef CPPAD_CG_VARIABLE_INCLUDED
#define CPPAD_CG_VARIABLE_INCLUDED
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

template<class Base>
inline bool CG<Base>::isVariable() const {
    return opNode_ != nullptr;
}

template<class Base>
inline bool CG<Base>::isParameter() const {
    return opNode_ == nullptr;
}

template<class Base>
inline bool CG<Base>::isValueDefined() const {
    return value_ != nullptr;
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
    if (value_ != nullptr) {
        *value_ = b;
    } else {
        value_ = new Base(b);
    }
}

template<class Base>
inline bool CG<Base>::isIdenticalZero() const throw (CGException) {
    return CppAD::IdenticalZero(getValue());
}

template<class Base>
inline bool CG<Base>::isIdenticalOne() const throw (CGException) {
    return CppAD::IdenticalOne(getValue());
}

template<class Base>
inline void CG<Base>::makeParameter(const Base &b) {
    opNode_ = nullptr;
    handler_ = nullptr;
    setValue(b);
}

template<class Base>
inline void CG<Base>::makeVariable(CodeHandler<Base>& handler,
                                   OperationNode<Base>* operation) {
    CPPADCG_ASSERT_UNKNOWN(operation != nullptr);
    opNode_ = operation;
    handler_ = &handler;
    delete value_;
    value_ = nullptr;

    handler.manageOperationNode(operation);
}

template<class Base>
inline void CG<Base>::makeVariable(CodeHandler<Base>& handler,
                                   OperationNode<Base>* operation,
                                   std::unique_ptr<Base>& value) {
    CPPADCG_ASSERT_UNKNOWN(operation != nullptr);
    opNode_ = operation;
    handler_ = &handler;
    delete value_;
    value_ = value.release();

    handler.manageOperationNode(operation);
}

template<class Base>
inline OperationNode<Base>* CG<Base>::getOperationNode() const {
    return opNode_;
}

template<class Base>
inline Argument<Base> CG<Base>::argument() const {
    if (opNode_ != nullptr)
        return Argument<Base> (*opNode_);
    else
        return Argument<Base> (*value_);
}

} // END cg namespace
} // END CppAD namespace

#endif