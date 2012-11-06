#ifndef CPPAD_CG_C_LANG_CUSTOM_VAR_NAME_GEN_INCLUDED
#define	CPPAD_CG_C_LANG_CUSTOM_VAR_NAME_GEN_INCLUDED
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
     * Creates variables names for the source code using a list of provided
     * custom names.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangCustomVariableNameGenerator : public CLangDefaultVariableNameGenerator<Base> {
    protected:
        //
        const std::vector<std::string> depNames_;
        const std::vector<std::string> indepNames_;
    public:

        CLangCustomVariableNameGenerator(const std::vector<std::string>& depNames,
                                         const std::vector<std::string>& indepNames) :
            CLangDefaultVariableNameGenerator<Base>(),
            depNames_(depNames),
            indepNames_(indepNames) {
        }

        CLangCustomVariableNameGenerator(const std::vector<std::string>& depNames,
                                         const std::vector<std::string>& indepNames,
                                         const std::string& depName,
                                         const std::string& indepName,
                                         const std::string& tmpName) :
            CLangDefaultVariableNameGenerator<Base>(depName, indepName, tmpName),
            depNames_(depNames),
            indepNames_(indepNames) {
        }

        virtual std::string generateDependent(const CG<Base>& variable, size_t index) {
            if (index < depNames_.size() && !depNames_[index].empty()) {
                return depNames_[index];
            } else {
                return CLangDefaultVariableNameGenerator<Base>::generateDependent(variable, index);
            }
        }

        virtual std::string generateIndependent(const SourceCodeFragment<Base>& independent) {
            size_t index = independent.variableID() - 1;
            if (index < indepNames_.size() && !indepNames_[index].empty()) {
                return indepNames_[index];
            } else {
                return CLangDefaultVariableNameGenerator<Base>::generateIndependent(independent);
            }
        }

        inline virtual ~CLangCustomVariableNameGenerator() {
        }

    };
}

#endif
