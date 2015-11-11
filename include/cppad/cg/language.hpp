#ifndef CPPAD_CG_LANGUAGE_INCLUDED
#define CPPAD_CG_LANGUAGE_INCLUDED
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
namespace cg {

/**
 * Information required for the generation of source code for a language
 * 
 * @author Joao Leal
 */
template<class Base>
class LanguageGenerationData {
public:
    /**
     * The independent variables
     */
    const std::vector<OperationNode<Base> *>& independent;
    /**
     * The dependent variables
     */
    const CppAD::vector<CG<Base> >& dependent;
    /**
     * The lowest ID used for temporary variables
     */
    size_t minTemporaryVarID;
    /**
     * Variable assignment order in the source code
     */
    const std::vector<OperationNode<Base>*>& variableOrder;
    /**
     * Provides the rules for variable name creation
     */
    VariableNameGenerator<Base>& nameGen;
    /**
     * maps atomic function IDs to their internal index
     */
    const std::map<size_t, size_t>& atomicFunctionId2Index;
    /**
     * maps atomic function IDs to their names
     */
    const std::map<size_t, std::string>& atomicFunctionId2Name;
    /**
     * the maximum forward mode order each atomic function is called 
     * (-1 means forward mode not used)
     */
    const std::vector<int>& atomicFunctionsMaxForward;
    /**
     * the maximum reverse mode order each atomic function is called
     * (-1 means reverse mode not used)
     */
    const std::vector<int>& atomicFunctionsMaxReverse;
    /**
     * a flag indicating whether or not temporary variable IDs have been recycled
     */
    const bool reuseIDs;
    //
    const std::set<const IndexDclrOperationNode<Base>*>& indexes;
    //
    const std::set<RandomIndexPattern*>& indexRandomPatterns;
    //
    const std::vector<IndexPattern*>& loopDependentIndexPatterns;
    //
    const std::vector<IndexPattern*>& loopIndependentIndexPatterns;
    /**
     * whether or not the dependent variables should be zeroed before 
     * executing the operation graph
     */
    const bool zeroDependents;
public:

    LanguageGenerationData(const std::vector<OperationNode<Base> *>& ind,
                           const CppAD::vector<CG<Base> >& dep,
                           size_t minTempVID,
                           const std::vector<OperationNode<Base>*>& vo,
                           VariableNameGenerator<Base>& ng,
                           const std::map<size_t, size_t>& atomicId2Index,
                           const std::map<size_t, std::string>& atomicId2Name,
                           const std::vector<int>& atomicMaxForward,
                           const std::vector<int>& atomicMaxReverse,
                           const bool ri,
                           const std::set<const IndexDclrOperationNode<Base>*>& indexs,
                           const std::set<RandomIndexPattern*>& idxRandomPatterns,
                           const std::vector<IndexPattern*>& dependentIndexPatterns,
                           const std::vector<IndexPattern*>& independentIndexPatterns,
                           bool zero) :
        independent(ind),
        dependent(dep),
        minTemporaryVarID(minTempVID),
        variableOrder(vo),
        nameGen(ng),
        atomicFunctionId2Index(atomicId2Index),
        atomicFunctionId2Name(atomicId2Name),
        atomicFunctionsMaxForward(atomicMaxForward),
        atomicFunctionsMaxReverse(atomicMaxReverse),
        reuseIDs(ri),
        indexes(indexs),
        indexRandomPatterns(idxRandomPatterns),
        loopDependentIndexPatterns(dependentIndexPatterns),
        loopIndependentIndexPatterns(independentIndexPatterns),
        zeroDependents(zero) {
    }
};

/**
 * Creates the source code for a specific language
 * 
 * @author Joao Leal
 */
template<class Base>
class Language {
protected:
    virtual void generateSourceCode(std::ostream& out, const std::unique_ptr<LanguageGenerationData<Base> >& info) = 0;

    /**
     * Whether or not a new variable is created as a result of this operation
     * 
     * @param op Operation
     * @return true if a new variable is created
     */
    virtual bool createsNewVariable(const OperationNode<Base>& op) const = 0;

    virtual bool requiresVariableArgument(enum CGOpCode op, size_t argIndex) const = 0;

    friend class CodeHandler<Base>;
};

} // END cg namespace
} // END CppAD namespace

#endif