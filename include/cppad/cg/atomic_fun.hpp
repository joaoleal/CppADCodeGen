#ifndef CPPAD_CG_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_ATOMIC_FUN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
 *    Copyright (C) 2013 Ciengis
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
 * An atomic function for source code generation.
 * Dynamic parameters are expected to be at the end of the independent
 * variable vector.
 *
 * @warning This class is not thread-safe!
 *
 * @author Joao Leal
 */
template <class Base>
class CGAtomicFun : public CGAbstractAtomicFun<Base> {
protected:
    using CGB = CG<Base>;
protected:
    atomic_three<Base>& atomicFun_;
    const CppAD::vector<Base> xSparsity_; // independent vector used to determine sparsity patterns
    const size_t nPar_; // number of dynamic parameters
public:

    /**
     * Creates a new atomic function wrapper that is responsible for
     * defining the dependencies to calls of a user atomic function.
     *
     * @param atomicFun The atomic function to the called by the compiled
     *                  source.
     * @param xSparsity Default independent vector used to determine sparsity patterns
     *                  when the provided independent vector using the CG data type does
     *                  not have all values defined.
     * @param nPar The number of dynamic parameters.
     * @param standAlone Whether or not forward and reverse function calls
     *                   do not require the Taylor coefficients for the
     *                   dependent variables (ty) and the previous
     *                   evaluation of other forward/reverse modes.
     */
    CGAtomicFun(atomic_three<Base>& atomicFun,
                const CppAD::vector<Base>& xSparsity,
                size_t nPar,
                bool standAlone = false) :
        CGAbstractAtomicFun<Base>(atomicFun.atomic_name(), standAlone),
        atomicFun_(atomicFun),
        xSparsity_(xSparsity),
        nPar_(nPar) {
    }

    CGAtomicFun(atomic_three<Base>& atomicFun,
                ArrayView<const Base> xSparsity,
                size_t nPar,
                bool standAlone = false) :
            CGAtomicFun(atomicFun, values(xSparsity), nPar, standAlone) {
    }

    CGAtomicFun(atomic_three<Base>& atomicFun,
                ArrayView<const CppAD::AD<Base>> xSparsity,
                size_t nPar,
                bool standAlone = false) :
            CGAtomicFun(atomicFun, values(xSparsity), nPar, standAlone) {
    }

    virtual ~CGAtomicFun() = default;

    template <class ADVector>
    void operator()(const ADVector& ax,
                    ADVector& ay) {
        CPPADCG_ASSERT_KNOWN(ax.size() > nPar_, "ax has the wrong size. "
                                                "It must be the concatenation of independent variables and dynamic parameters.")

        CPPADCG_ASSERT_KNOWN(ax.size() == xSparsity_.size(), "ax has the wrong size. "
                                                             "It must be the concatenation of independent variables and dynamic parameters.")

        this->CGAbstractAtomicFun<Base>::operator()(ax, ay);
    }

    size_t size_dyn_ind() const override {
        return nPar_;
    }

    bool for_type(const CppAD::vector<CGB>& parameter_x,
                  const CppAD::vector<ad_type_enum>& type_x,
                  CppAD::vector<ad_type_enum>& type_y) override {

        this->validateTypeX(type_x);

        return atomicFun_.for_type(sparsityIndeps(parameter_x), type_x, type_y);
    }

    bool rev_depend(const CppAD::vector<CGB>& parameter_x,
                    const CppAD::vector<ad_type_enum>& type_x,
                    CppAD::vector<bool>& depend_x,
                    const CppAD::vector<bool>& depend_y) override {

        this->validateTypeX(type_x);

        return atomicFun_.rev_depend(sparsityIndeps(parameter_x), type_x, depend_x, depend_y);
    }

    bool jac_sparsity(const CppAD::vector<CGB>& parameter_x,
                      const CppAD::vector<ad_type_enum>& type_x,
                      bool dependency,
                      const CppAD::vector<bool>& select_x,
                      const CppAD::vector<bool>& select_y,
                      CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {

        this->validateTypeX(type_x);

        return atomicFun_.jac_sparsity(sparsityIndeps(parameter_x), type_x, dependency, select_x, select_y, pattern_out);
    }

    bool hes_sparsity(const CppAD::vector<CGB>& parameter_x,
                      const CppAD::vector<ad_type_enum>& type_x,
                      const CppAD::vector<bool>& select_x,
                      const CppAD::vector<bool>& select_y,
                      CppAD::sparse_rc<CppAD::vector<size_t> >& pattern_out) override {

        this->validateTypeX(type_x);

        return atomicFun_.hes_sparsity(sparsityIndeps(parameter_x), type_x, select_x, select_y, pattern_out);
    }

protected:

    bool atomicForward(const CppAD::vector<Base>& parameter_x,
                       const CppAD::vector<CppAD::ad_type_enum>& type_x,
                       size_t need_y,
                       size_t order_low,
                       size_t order_up,
                       const CppAD::vector<Base>& taylor_x,
                       CppAD::vector<Base>& taylor_y) override {

        this->validateTypeX(type_x);

        return atomicFun_.forward(parameter_x, type_x, need_y, order_low, order_up, taylor_x, taylor_y);
    }

    bool atomicReverse(const CppAD::vector<Base>& parameter_x,
                       const CppAD::vector<CppAD::ad_type_enum>& type_x,
                       size_t order_up,
                       const CppAD::vector<Base>& taylor_x,
                       const CppAD::vector<Base>& taylor_y,
                       CppAD::vector<Base>& partial_x,
                       const CppAD::vector<Base>& partial_y) override {

        this->validateTypeX(type_x);

        return atomicFun_.reverse(parameter_x, type_x, order_up, taylor_x, taylor_y, partial_x, partial_y);
    }

private:
    inline CppAD::vector<Base> sparsityIndeps(const CppAD::vector<CGB>& x) const {
        CPPADCG_ASSERT_UNKNOWN(x.size() == xSparsity_.size())

        size_t n = x.size();
        CppAD::vector<Base> out(n);
        for (size_t i = 0; i < n; ++i) {
            if (x[i].isValueDefined()) {
                out[i] = x[i].getValue();
            } else {
                out = xSparsity_;
                break;
            }
        }

        return out;
    }

    inline static CppAD::vector<Base> values(ArrayView<const CppAD::AD<Base>> x) {
        CppAD::vector<Base> out(x.size());
        for (size_t i = 0; i < out.size(); ++i) {
            out[i] = CppAD::Value(CppAD::Var2Par(x[i]));
        }
        return out;
    }

    inline static CppAD::vector<Base> values(ArrayView<const Base> x) {
        CppAD::vector<Base> out(x.size());
        for (size_t i = 0; i < out.size(); ++i) {
            out[i] = x[i];
        }
        return out;
    }
};

} // END cg namespace
} // END CppAD namespace

#endif
