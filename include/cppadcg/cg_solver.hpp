#ifndef CPPAD_CG_SOLVER_INCLUDED
#define CPPAD_CG_SOLVER_INCLUDED
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
    inline CG<Base> CodeHandler<Base>::solveFor(SourceCodeFragment<Base>* expression,
                                                SourceCodeFragment<Base>* code) throw (CGException) {

        assert(expression != NULL);
        assert(code != NULL);

        using std::vector;

        // find code in expression
        if (expression == code)
            return CG<Base > (*this, Argument<Base > (*code));

        typedef vector<SourceCodePathNode<Base> > SourceCodePath;

        // find the location of the variable (code) in the expression
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

            CGOpCode op = pnodeOp.node->operationCode();
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
                case CGDivMulOp:
                {
                    vector<Argument<Base> > rhsArgs(args.size());
                    size_t j = 0;
                    for (size_t i = 0; i < args.size(); i++) {
                        if(i != argIndex)
                            rhsArgs[j++] = args[i];
                    }
                    if(rightHs.isVariable())
                        rhsArgs.back() = Argument<Base > (*rightHs.getSourceCodeFragment());
                    else
                        rhsArgs.back() = Argument<Base > (rightHs.getParameterValue());
                    
                    const vector<CGOpCodeExtra>& opInfoOrig = code->operationInfo();
                    vector<CGOpCodeExtra> operationInfo(opInfoOrig.size());
                    if (argIndex == 0 || opInfoOrig[argIndex - 1] == CGExtraMulOp) {
                        j = 0;
                        for (size_t i = 0; i < opInfoOrig.size(); i++) {
                            if (i + 1 != argIndex)
                                operationInfo[j++] = opInfoOrig[i] == CGExtraMulOp ? CGExtraDivOp : CGExtraMulOp;
                        }
                        operationInfo.back() = CGExtraMulOp;
                    } else {
                        j = 0;
                        for (size_t i = 0; i < args.size() - 1; i++) {
                            if (i + 1 != argIndex)
                                operationInfo[j++] = opInfoOrig[i];
                        }
                        operationInfo.back() = CGExtraDivOp;
                    }
                    rightHs = CG<Base > (*this, new SourceCodeFragment<Base>(CGDivMulOp, rhsArgs, operationInfo));
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
