#ifndef CPPAD_CG_EVALUATOR_INCLUDED
#define	CPPAD_CG_EVALUATOR_INCLUDED
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
     * Utility class used for some code transformation
     */
    template<class Base, class BaseOut>
    class Evaluator {
    protected:
        const CodeHandler<Base>* handler_;
        const std::vector<CG<Base> > dep_;
        const std::vector<AD<BaseOut > >* indep_;
        std::map<SourceCodeFragment<Base>*, AD<BaseOut> > evals_;
    public:

        /**
         * \param handler The source code handler
         * \param dep Dependent variable vector (all variables must belong to
         *             the same code handler)
         */
        Evaluator(const CodeHandler<Base>& handler, const std::vector<CG<Base> >& dep) :
            handler_(&handler),
            dep_(dep),
            indep_(NULL) {
        }

        /**
         * Performs all the operations required to calculate the dependent 
         * variables with a (potentially) new data type
         * 
         * \param indep The new independent variables.
         * \return The dependent variable values
         */
        inline std::vector<AD<BaseOut> > evaluate(const std::vector<AD<BaseOut> >& indep) throw (CGException) {
            if (handler_->getIndependentVariableSize() != indep.size()) {
                std::stringstream ss;
                ss << "Invalid independent variable size. Expected " << handler_->getIndependentVariableSize() << " but got " << indep.size() << ".";
                throw CGException(ss.str());
            }

            if (indep.empty()) {
                return std::vector<AD<BaseOut> >(0); // nothing in, nothing out
            }

            indep_ = &indep;

            evals_.clear(); // clean-up

            std::vector<AD<BaseOut> > newDep(dep_.size());

            for (size_t i = 0; i < dep_.size(); i++) {
                newDep[i] = evalCG(dep_[i]);
            }

            evals_.clear(); // clean-up

            return newDep;
        }

    private:

        inline AD<BaseOut> evalCG(const CG<Base>& dep) const throw (CGException) {
            if (dep.isParameter()) {
                // parameter
                return AD<BaseOut > (dep.getParameterValue());
            } else {
                return evalSourceCodeFragment(*dep.getSourceCodeFragment());
            }
        }

        inline AD<BaseOut> evalArg(const Argument<Base>& arg) const throw (CGException) {
            if (arg.operation() != NULL) {
                return evalSourceCodeFragment(*arg.operation());
            } else {
                // parameter
                return AD<BaseOut > (*arg.parameter());
            }
        }

        inline AD<BaseOut> evalSourceCodeFragment(SourceCodeFragment<Base>& node) const throw (CGException) {
            // check if this node was previously determined
            typename std::map<SourceCodeFragment<Base>*, AD<BaseOut> >::const_iterator it;
            it = evals_.find(&node);
            if (it != evals_.end()) {
                return it->second;
            }

            // first evaluation of this node
            const std::vector<Argument<Base> >& args = node.arguments();
            const CGOpCode code = node.operation();
            switch (code) {
                case CGAbsOp: //  abs(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for abs()");
                    return abs(evalArg(args[0]));
                case CGAcosOp: // acos(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for acos()");
                    return acos(evalArg(args[0]));
                case CGAddOp: //  a + b
                    CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for addition");
                    return evalArg(args[0]) + evalArg(args[1]);
                case CGAliasOp:
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for alias");
                    return evalArg(args[0]);
                case CGAsinOp: // asin(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for asin()");
                    return asin(evalArg(args[0]));
                case CGAtanOp: // atan(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for atan()");
                    return atan(evalArg(args[0]));
                case CGComOpLt: // result = left < right? trueCase: falseCase
                    CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareLt, )");
                    return CondExpOp(CompareLt, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                case CGComOpLe: // result = left <= right? trueCase: falseCase
                    CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareLe, )");
                    return CondExpOp(CompareLe, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                case CGComOpEq: // result = left == right? trueCase: falseCase
                    CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareEq, )");
                    return CondExpOp(CompareEq, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                case CGComOpGe: // result = left >= right? trueCase: falseCase
                    CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareGe, )");
                    return CondExpOp(CompareGe, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                case CGComOpGt: // result = left > right? trueCase: falseCase
                    CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareGt, )");
                    return CondExpOp(CompareGt, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                case CGComOpNe: // result = left != right? trueCase: falseCase
                    CPPADCG_ASSERT_KNOWN(args.size() == 4, "Invalid number of arguments for CondExpOp(CompareNe, )");
                    return CondExpOp(CompareNe, evalArg(args[0]), evalArg(args[1]), evalArg(args[2]), evalArg(args[3]));
                case CGCoshOp: // cosh(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for cosh()");
                    return cosh(evalArg(args[0]));
                case CGCosOp: //  cos(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for cos()");
                    return cos(evalArg(args[0]));
                case CGDivOp: // a / b
                    CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for division");
                    return evalArg(args[0]) / evalArg(args[1]);
                case CGExpOp: //  exp(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for exp()");
                    return exp(evalArg(args[0]));
                case CGInvOp: //                             independent variable
                {
                    size_t index = handler_->getIndependentVariableIndex(node);
                    return (*indep_)[index];
                }
                case CGLogOp: //  log(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for log()");
                    return log(evalArg(args[0]));
                case CGMulOp: // a * b
                    CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for multiplication");
                    return evalArg(args[0]) * evalArg(args[1]);
                case CGPowOp: //  pow(a,   b)
                    CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for pow()");
                    return pow(evalArg(args[0]), evalArg(args[1]));
                    //case PriOp: //  PrintFor(text, parameter or variable, parameter or variable)
                case CGSignOp: // result = (x > 0)? 1.0:((x == 0)? 0.0:-1)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sign()");
                    return sign(evalArg(args[0]));
                case CGSinhOp: // sinh(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sinh()");
                    return sinh(evalArg(args[0]));
                case CGSinOp: //  sin(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sin()");
                    return sin(evalArg(args[0]));
                case CGSqrtOp: // sqrt(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for sqrt()");
                    return sqrt(evalArg(args[0]));
                case CGSubOp: //  a - b
                    CPPADCG_ASSERT_KNOWN(args.size() == 2, "Invalid number of arguments for subtraction");
                    return evalArg(args[0]) - evalArg(args[1]);
                case CGTanhOp: //  tanh(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for tanh()");
                    return tanh(evalArg(args[0]));
                case CGTanOp: //  tan(variable)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for tan()");
                    return tan(evalArg(args[0]));
                case CGUnMinusOp: // -(a)
                    CPPADCG_ASSERT_KNOWN(args.size() == 1, "Invalid number of arguments for unary minus");
                    return -evalArg(args[0]);
                default:
                {
                    std::stringstream ss;
                    ss << "Unknown operation code '" << code << "'";
                    throw CGException(ss.str());
                }
            }
        }
    };
}

#endif
