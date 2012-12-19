#ifndef CPPAD_CG_LANGUAGE_INCLUDED
#define	CPPAD_CG_LANGUAGE_INCLUDED
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
     * Information required for the generation of source code for a language
     * 
     * \author Joao Leal
     */
    template<class Base>
    class LanguageGenerationData {
    public:
        // The independent variables
        const std::vector<SourceCodeFragment<Base> *>& independent;
        // The dependent variables
        const std::vector<CG<Base> >& dependent;
        // the lowest ID used for temporary variables
        size_t minTemporaryVarID;
        // The order of the assignment of the variables in the source code
        const std::vector<SourceCodeFragment<Base>*>& variableOrder;
        // Provides the rules for variable name creation
        VariableNameGenerator<Base>& nameGen;
        // a flag indicating whether or not temporary variable IDs have been recycled
        const bool reuseIDs;

    public:

        LanguageGenerationData(const std::vector<SourceCodeFragment<Base> *>& ind,
                               const std::vector<CG<Base> >& dep,
                               size_t minTempVID,
                               const std::vector<SourceCodeFragment<Base>*>& vo,
                               VariableNameGenerator<Base>& ng,
                               const bool ri) :
            independent(ind),
            dependent(dep),
            minTemporaryVarID(minTempVID),
            variableOrder(vo),
            nameGen(ng),
            reuseIDs(ri) {
        }
    };

    /**
     * Creates the source code for a specific language
     * 
     * \author Joao Leal
     */
    template<class Base>
    class Language {
    protected:
        virtual void generateSourceCode(std::ostream& out, LanguageGenerationData<Base>& info) = 0;

        virtual bool createsNewVariable(const SourceCodeFragment<Base>& op) = 0;

        virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) = 0;

        friend class CodeHandler<Base>;
    };

}

#endif

