#ifndef CPPAD_CG_DAE_INDEX_REDUCTION_INCLUDED
#define CPPAD_CG_DAE_INDEX_REDUCTION_INCLUDED
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

#include <cppad/cg/cppadcg.hpp>
#include <cppad/cg/dae_index_reduction/dae_var_info.hpp>
#include <cppad/cg/dae_index_reduction/dae_equation_info.hpp>
#include <cppad/cg/dae_index_reduction/simple_logger.hpp>

namespace CppAD {
namespace cg {

/**
 * Base class for algorithms that perform automatic (differentiation) index
 * reduction of implicit DAEs.
 */
template<class Base>
class DaeIndexReduction : public SimpleLogger {
protected:
    /**
     * The original model
     */
    ADFun<CG<Base> > * const fun_;
public:

    /**
     * Creates a new DAE model index reduction algorithm.
     * 
     * @param fun  The original (high index) model
     * @param varInfo  DAE  system variable information (in the same order 
     *                 as in the tape)
     */
    DaeIndexReduction(ADFun<CG<Base> >* fun) :
        fun_(fun) {
    }

    inline virtual ~DaeIndexReduction() {
    }

    virtual ADFun<CG<Base> >* reduceIndex(std::vector<DaeVarInfo>& newVarInfo,
                                          std::vector<DaeEquationInfo>& equationInfo) = 0;

};

} // END cg namespace
} // END CppAD namespace

#endif	

