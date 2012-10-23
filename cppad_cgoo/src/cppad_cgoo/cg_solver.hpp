#ifndef CG_SOLVER_INCLUDED
#define	CG_SOLVER_INCLUDED

#include "cg_source_code_path.hpp"
#include "cg_source_code_fragment.hpp"
#include "cg_argument.hpp"

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
     * Solves an expression (e.g. f(x, y) == 0) for a given variable (e.g. x)
     * The variable can appear only once in the expression.
     * 
     * \param expression  The original expression (f(x, y))
     * \param code  The variable to solve for
     * \return  The expression for variable
     */
    template<class Base>
    inline CG<Base> CodeHandler<Base>::solveFor(SourceCodeFragment<Base>* expression,
                                                SourceCodeFragment<Base>* code) throw (CGException) {

        assert(expression != NULL);
        assert(code != NULL);

        using std::vector;

        // find code in expression
        if (expression == code)
            return CG<Base > (*this, Argument<Base > (*code));

        typedef vector<SourceCodePathNode<Base> > SourceCodePath;

        vector<SourceCodePath> paths = findPaths(expression, code, 2);

        if (paths.empty()) {
            throw CGException("The provided variable is not present in the expression");
        } else if (paths.size() > 1) {
            // todo: support multiple variable locations
            throw CGException("Unable to determine expression for variable:"
                              " the provided variable was found in multiple locations (not yet supported)");
        }


        CG<Base> rightHs(0.0);

        SourceCodePath& path0 = paths[0];

        assert(path0.back().node == code);

        for (size_t n = 0; n < path0.size() - 1; ++n) {
            SourceCodePathNode<Base>& pnodeOp = path0[n];
            size_t argIndex = path0[n + 1].arg_index;
            const std::vector<Argument<Base> >& args = pnodeOp.node->arguments();

            CGOpCode op = pnodeOp.node->operation();
            switch (op) {
                case CGMulOp:
                {
                    const Argument<Base>& other = args[argIndex == 0 ? 1 : 0];
                    rightHs /= CG<Base > (*this, other);
                    break;
                }
                case CGDivOp:
                    if (argIndex == 0) {
                        const Argument<Base>& other = args[argIndex == 0 ? 1 : 0];
                        rightHs *= CG<Base > (*this, other);
                    } else {
                        const Argument<Base>& other = args[argIndex == 0 ? 1 : 0];
                        rightHs = CG<Base > (*this, other) / rightHs;
                    }
                    break;

                case CGUnMinusOp:
                    rightHs *= Base(-1.0);
                    break;
                case CGAddOp:
                {
                    const Argument<Base>& other = args[argIndex == 0 ? 1 : 0];
                    rightHs -= CG<Base > (*this, other);
                    break;
                }
                case CGSubOp:
                {
                    if (argIndex == 0) {
                        rightHs += CG<Base > (*this, args[1]);
                    } else {
                        rightHs = CG<Base > (*this, args[0]) - rightHs;
                    }
                    break;
                }
                case CGExpOp:
                    rightHs = log(rightHs);
                    break;
                case CGLogOp:
                    rightHs = exp(rightHs);
                    break;
                case CGPowOp:
                {
                    if (argIndex == 0) {
                        // base
                        const Argument<Base>& exponent = args[1];
                        if (exponent.parameter() != NULL && *exponent.parameter() == Base(2.0)) {
                            rightHs = sqrt(rightHs);
                        } else {
                            rightHs = pow(rightHs, Base(1.0) / CG<Base > (*this, exponent));
                        }
                    } else {
                        // 
                        const Argument<Base>& base = args[0];
                        rightHs = log(rightHs) / log(CG<Base > (*this, base));
                    }
                    break;
                }
                case CGSqrtOp:
                    rightHs *= rightHs;
                    break;
                    //case CGAcosOp: // asin(variable)
                    //case CGAsinOp: // asin(variable)
                    //case CGAtanOp: // atan(variable)
                case CGCoshOp: // cosh(variable)
                {
                    rightHs = log(rightHs + sqrt(rightHs * rightHs - Base(1.0))); // asinh
                    break;
                    //case CGCosOp: //  cos(variable)
                }
                case CGSinhOp: // sinh(variable)
                    rightHs = log(rightHs + sqrt(rightHs * rightHs + Base(1.0))); // asinh
                    break;
                    //case CGSinOp: //  sin(variable)
                case CGTanhOp: //  tanh(variable)
                    rightHs = Base(0.5) * (log(Base(1.0) + rightHs) - log(Base(1.0) - rightHs)); // atanh
                    break;
                    //case CGTanOp: //  tan(variable)
                default:
                    std::ostringstream ss;
                    ss << "Unable to invert operation '" << op << "'";
                    throw CGException(ss.str());
            };
        }

        return rightHs;
    }
}

#endif
