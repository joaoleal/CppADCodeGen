#ifndef CPPAD_CG_COND_EXP_OP_INCLUDED
#define	CPPAD_CG_COND_EXP_OP_INCLUDED
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
    inline CG<Base> CondExpOp(enum CompareOp cop,
                              const CG<Base> &left,
                              const CG<Base> &right,
                              const CG<Base> &trueCase,
                              const CG<Base> &falseCase) {

        if (left.isParameter() && right.isParameter()) {
            switch (cop) {
                case CompareLt:
                    if (left.getParameterValue() < right.getParameterValue())
                        return trueCase;
                    else return falseCase;
                    break;

                case CompareLe:
                    if (left.getParameterValue() <= right.getParameterValue())
                        return trueCase;
                    else return falseCase;
                    break;

                case CompareEq:
                    if (left.getParameterValue() == right.getParameterValue())
                        return trueCase;
                    else return falseCase;
                    break;

                case CompareGe:
                    if (left.getParameterValue() >= right.getParameterValue())
                        return trueCase;
                    else return falseCase;
                    break;

                case CompareGt:
                    if (left.getParameterValue() > right.getParameterValue())
                        return trueCase;
                    else return falseCase;
                    break;

                default:
                    CPPAD_ASSERT_UNKNOWN(0);
                    return trueCase;
            }

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
            }

            if ((!right.isParameter() && right.getCodeHandler() != handler)
                    || (!trueCase.isParameter() && trueCase.getCodeHandler() != handler)
                    || (!falseCase.isParameter() && falseCase.getCodeHandler() != handler)) {
                throw CGException("Attempting to use different source code generation handlers in the same source code generation");
            }

            std::string leftStr = handler->operations(left);
            std::string rightStr = handler->operations(right);
            std::string trueCaseStr = handler->operations(trueCase);
            std::string falseCaseStr = handler->operations(falseCase);

            CG<Base> result;
            result.makeVariable(*handler);
            std::string resultStr = handler->createVariableName(result);

            std::ostream* out = handler->getOutputStream();
            (*out) << "if(";
            handler->printComparison(leftStr, cop, rightStr);
            (*out) << ") {\n";
            handler->printOperationAssign(resultStr, trueCaseStr);
            (*out) << "} else {\n";
            handler->printOperationAssign(resultStr, falseCaseStr);
            (*out) << "}\n";

            return result;
        }

    }

}
#endif

