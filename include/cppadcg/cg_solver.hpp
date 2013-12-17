#ifndef CPPAD_CG_SOLVER_INCLUDED
#define CPPAD_CG_SOLVER_INCLUDED
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
    inline CG<Base> CodeHandler<Base>::solveFor(OperationNode<Base>& expression,
                                                OperationNode<Base>& code) throw (CGException) {
        using std::vector;

        // find code in expression
        if (&expression == &code)
            return CG<Base > (*this, Argument<Base> (code));

        typedef vector<OperationPathNode<Base> > SourceCodePath;

        vector<SourceCodePath> paths = findPaths(expression, code, 2);
        if (paths.empty()) {
            throw CGException("The provided variable is not present in the expression");
        } else if (paths.size() > 1) {
            // todo: support multiple variable locations
            throw CGException("Unable to determine expression for variable:"
                              " the provided variable was found in multiple locations (not yet supported)");
        }

        CPPADCG_ASSERT_UNKNOWN(paths[0].back().node == &code);

        return solveFor(paths[0]);
    }

    template<class Base>
    inline CG<Base> CodeHandler<Base>::solveFor(const std::vector<OperationPathNode<Base> >& path) throw (CGException) {

        CG<Base> rightHs(0.0);

        for (size_t n = 0; n < path.size() - 1; ++n) {
            const OperationPathNode<Base>& pnodeOp = path[n];
            size_t argIndex = path[n + 1].arg_index;
            const std::vector<Argument<Base> >& args = pnodeOp.node->getArguments();

            CGOpCode op = pnodeOp.node->getOperationType();
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
                case CGAliasOp:
                    // do nothing 
                    break;
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
                        if (exponent.getParameter() != NULL && *exponent.getParameter() == Base(0.0)) {
                            throw CGException("Invalid zero exponent");
                        } else if (exponent.getParameter() != NULL && *exponent.getParameter() == Base(1.0)) {
                            continue; // do nothing
                        } else {
                            std::ostringstream ss;
                            ss << "Unable to invert operation '" << op << "'";
                            throw CGException(ss.str());
                            /*
                            if (exponent.getParameter() != NULL && *exponent.getParameter() == Base(2.0)) {
                                rightHs = sqrt(rightHs); // TODO: should -sqrt(rightHs) somehow be considered???
                            } else {
                                rightHs = pow(rightHs, Base(1.0) / CG<Base > (*this, exponent));
                            }
                             */
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

    template<class Base>
    inline bool isSolvable(const std::vector<OperationPathNode<Base> >& path) throw (CGException) {
        for (size_t n = 0; n < path.size() - 1; ++n) {
            const OperationPathNode<Base>& pnodeOp = path[n];
            size_t argIndex = path[n + 1].arg_index;
            const std::vector<Argument<Base> >& args = pnodeOp.node->getArguments();

            CGOpCode op = pnodeOp.node->getOperationType();
            switch (op) {
                case CGMulOp:
                case CGDivOp:
                case CGUnMinusOp:
                case CGAddOp:
                case CGAliasOp:
                case CGSubOp:
                case CGExpOp:
                case CGLogOp:
                case CGSqrtOp:
                case CGCoshOp: // cosh(variable)
                case CGSinhOp: // sinh(variable)
                case CGTanhOp: //  tanh(variable)
                    break;
                case CGPowOp:
                {
                    if (argIndex == 0) {
                        // base
                        const Argument<Base>& exponent = args[1];
                        if (exponent.getParameter() != NULL && *exponent.getParameter() == Base(0.0)) {
                            return false;
                        } else if (exponent.getParameter() != NULL && *exponent.getParameter() == Base(1.0)) {
                            break;
                        } else {
                            return false;
                        }
                    } else {
                        break;
                    }
                    break;
                }

                default:
                    return false;
            };
        }
        return true;
    }
}

#endif
