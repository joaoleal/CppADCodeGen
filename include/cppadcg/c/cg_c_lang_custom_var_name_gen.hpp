#ifndef CPPAD_CG_C_LANG_CUSTOM_VAR_NAME_GEN_INCLUDED
#define CPPAD_CG_C_LANG_CUSTOM_VAR_NAME_GEN_INCLUDED
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

    /**
     * Creates variables names for the source code using a list of provided
     * custom names.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class CLangCustomVariableNameGenerator : public CLangDefaultVariableNameGenerator<Base> {
    protected:
        //
        const std::vector<std::string> depNames_;
        const std::vector<std::string> indepNames_;
    public:

        CLangCustomVariableNameGenerator(const std::vector<std::string>& depNames,
                                         const std::vector<std::string>& indepNames,
                                         const std::string& depName = "y",
                                         const std::string& indepName = "x",
                                         const std::string& tmpName = "v",
                                         const std::string& tmpArrayName = "array") :
            CLangDefaultVariableNameGenerator<Base>(depName, indepName, tmpName, tmpArrayName),
            depNames_(depNames),
            indepNames_(indepNames) {
        }

        virtual std::string generateDependent(size_t index) {
            if (index < depNames_.size() && !depNames_[index].empty()) {
                return depNames_[index];
            } else {
                return CLangDefaultVariableNameGenerator<Base>::generateDependent(index);
            }
        }

        virtual std::string generateIndependent(const OperationNode<Base>& independent) {
            size_t index = independent.getVariableID() - 1;
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
