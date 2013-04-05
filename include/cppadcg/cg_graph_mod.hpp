#ifndef CPPAD_CG_GRAPH_MOD_INCLUDED
#define CPPAD_CG_GRAPH_MOD_INCLUDED
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
    inline void CodeHandler<Base>::substituteIndependent(const CG<Base>& indep,
                                                         const CG<Base>& dep,
                                                         bool removeFromIndependents) throw (CGException) {
        substituteIndependent(*indep.getSourceCodeFragment(), *dep.getSourceCodeFragment(), removeFromIndependents);
    }

    template<class Base>
    inline void CodeHandler<Base>::substituteIndependent(SourceCodeFragment<Base>& indep,
                                                         SourceCodeFragment<Base>& dep,
                                                         bool removeFromIndependents) throw (CGException) {
        using std::vector;
        typedef CG<Base> CGBase;
        typedef AD<CGBase> ADCG;

        //check if the independent variable belongs to this handler
        size_t indepIndex = getIndependentVariableIndex(indep);

        //check if the dependent variable belongs to this handler
        typename std::vector<SourceCodeFragment<Base> *>::const_iterator it =
                std::find(_codeBlocks.begin(), _codeBlocks.end(), &dep);
        if (it == _codeBlocks.end()) {
            throw CGException("The dependent variable does not belong to this handler");
        }

        // determine the expression for the independent variable
        CGBase dummyExp = solveFor(dep, indep);

        Argument<Base> arg;
        // change the independent variable
        if (dummyExp.isVariable()) {
            arg = Argument<Base > (*dummyExp.getSourceCodeFragment());
        } else {
            // create a bogus variable to avoid searching for all occurrences of the independent variable
            arg = Argument<Base > (dummyExp.getParameterValue());
        }

        indep.makeAlias(arg);

        if (removeFromIndependents) {
            // remove the substituted variable from the independent variable vector
            _independentVariables.erase(_independentVariables.begin() + indepIndex);
        }
    }

    template<class Base>
    inline void CodeHandler<Base>::undoSubstituteIndependent(SourceCodeFragment<Base>& indep) throw (CGException) {
        typename std::vector<SourceCodeFragment<Base> *>::const_iterator it =
                std::find(_independentVariables.begin(), _independentVariables.end(), &indep);
        if (it == _independentVariables.end()) {
            throw CGException("Variable not found in the independent variable vector");
        }

        indep.setOperation(CGInvOp);
    }

    template<class Base>
    inline void CodeHandler<Base>::removeIndependent(SourceCodeFragment<Base>& indep) throw (CGException) {
        if (indep.operation() != CGAliasOp) {
            throw CGException("Cannot remove independent variable: not an alias");
        }

        typename std::vector<SourceCodeFragment<Base> *>::iterator it =
                std::find(_independentVariables.begin(), _independentVariables.end(), &indep);
        if (it == _independentVariables.end()) {
            throw CGException("Variable not found in the independent variable vector");
        }
        _independentVariables.erase(it);
    }

}

#endif
