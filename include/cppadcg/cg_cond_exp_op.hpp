#ifndef CPPAD_CG_COND_EXP_OP_INCLUDED
#define CPPAD_CG_COND_EXP_OP_INCLUDED
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

    template<class Base>
    inline CG<Base> CondExpOp(enum CompareOp cop,
                              const CG<Base> &left,
                              const CG<Base> &right,
                              const CG<Base> &trueCase,
                              const CG<Base> &falseCase) {

        if (left.isParameter() && right.isParameter()) {
            switch (cop) {
                case CompareLt:
                    if (left.getValue() < right.getValue())
                        return trueCase;
                    else
                        return falseCase;

                case CompareLe:
                    if (left.getValue() <= right.getValue())
                        return trueCase;
                    else
                        return falseCase;

                case CompareEq:
                    if (left.getValue() == right.getValue())
                        return trueCase;
                    else
                        return falseCase;

                case CompareGe:
                    if (left.getValue() >= right.getValue())
                        return trueCase;
                    else
                        return falseCase;

                case CompareGt:
                    if (left.getValue() > right.getValue())
                        return trueCase;
                    else
                        return falseCase;

                default:
                    CPPAD_ASSERT_UNKNOWN(0);
                    return trueCase;
            }

        } else if ((trueCase.isParameter() && falseCase.isParameter() &&
                trueCase.getValue() == falseCase.getValue()) ||
                (trueCase.isVariable() && falseCase.isVariable() &&
                trueCase.getOperationNode() == falseCase.getOperationNode())) {
            return trueCase;
        } else {

            CodeHandler<Base>* handler;

            if (!left.isParameter()) {
                handler = left.getCodeHandler();
            } else if (!right.isParameter()) {
                handler = right.getCodeHandler();
            } else if (!trueCase.isParameter()) {
                handler = trueCase.getCodeHandler();
            } else if (!falseCase.isParameter()) {
                handler = falseCase.getCodeHandler();
            } else {
                CPPAD_ASSERT_UNKNOWN(0);
                throw CGException("Unexpected error!");
            }

            if ((!right.isParameter() && right.getCodeHandler() != handler)
                    || (!trueCase.isParameter() && trueCase.getCodeHandler() != handler)
                    || (!falseCase.isParameter() && falseCase.getCodeHandler() != handler)) {
                throw CGException("Attempting to use different source code generation handlers in the same source code generation");
            }

            CGOpCode op;
            switch (cop) {
                case CompareLt:
                    op = CGComOpLt;
                    break;
                case CompareLe:
                    op = CGComOpLe;
                    break;
                case CompareEq:
                    op = CGComOpEq;
                    break;
                case CompareGe:
                    op = CGComOpGe;
                    break;
                case CompareGt:
                    op = CGComOpGt;
                    break;
                case CompareNe:
                    op = CGComOpNe;
                    break;

                default:
                    CPPAD_ASSERT_UNKNOWN(0);
                    throw CGException("Unexpected error!");
            }

            CG<Base> result(*handler, new OperationNode<Base> (op, {left.argument(), right.argument(), trueCase.argument(), falseCase.argument()}));

            if (left.isValueDefined() && right.isValueDefined()) {
                switch (cop) {
                    case CompareLt:
                        if (left.getValue() < right.getValue()) {
                            if (trueCase.isValueDefined()) {
                                result.setValue(trueCase.getValue());
                            }
                        } else {
                            if (falseCase.isValueDefined()) {
                                result.setValue(falseCase.getValue());
                            }
                        }
                        break;

                    case CompareLe:
                        if (left.getValue() <= right.getValue()) {
                            if (trueCase.isValueDefined()) {
                                result.setValue(trueCase.getValue());
                            }
                        } else {
                            if (falseCase.isValueDefined()) {
                                result.setValue(falseCase.getValue());
                            }
                        }
                        break;

                    case CompareEq:
                        if (left.getValue() == right.getValue()) {
                            if (trueCase.isValueDefined()) {
                                result.setValue(trueCase.getValue());
                            }
                        } else {
                            if (falseCase.isValueDefined()) {
                                result.setValue(falseCase.getValue());
                            }
                        }
                        break;

                    case CompareGe:
                        if (left.getValue() >= right.getValue()) {
                            if (trueCase.isValueDefined()) {
                                result.setValue(trueCase.getValue());
                            }
                        } else {
                            if (falseCase.isValueDefined()) {
                                result.setValue(falseCase.getValue());
                            }
                        }
                        break;

                    case CompareGt:
                        if (left.getValue() > right.getValue()) {
                            if (trueCase.isValueDefined()) {
                                result.setValue(trueCase.getValue());
                            }
                        } else {
                            if (falseCase.isValueDefined()) {
                                result.setValue(falseCase.getValue());
                            }
                        }
                        break;

                    default:
                        break;
                }
            }

            return result;
        }

    }

}
#endif

