#ifndef CPPAD_CG_ATOMIC_FUN_BRIDGE_INCLUDED
#define CPPAD_CG_ATOMIC_FUN_BRIDGE_INCLUDED
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
 * An atomic function wrapper for atomic functions using the type
 * ::CppAD::cg::CG.
 * This class can be useful when a CppAD::ADFun<CppAD::cg::CG> is going to
 * be used to create a compiled model library but has not been compiled yet.
 *
 * @warning This class is not thread-safe!
 *
 * @author Joao Leal
 */
template<class Base>
class CGAtomicFunBridge : public CGAbstractAtomicFun<Base> {
public:
    using CGB = CppAD::cg::CG<Base>;
    using ADCGD = CppAD::AD<CGB>;
protected:
    ADFun<CGB>& fun_;
    bool cacheSparsities_;
    CustomPosition custom_jac_;
    CustomPosition custom_hess_;
public:

    /**
     * Creates a new atomic function wrapper.
     *
     * @param name The atomic function name
     * @param fun The atomic function to be wrapped
     * @param standAlone Whether or not forward and reverse function calls
     *                   do not require the Taylor coefficients for the
     *                   dependent variables (ty) and the previous
     *                   evaluation of other forward/reverse modes.
     * @param cacheSparsities Whether or not to cache information related
     *                        with sparsity evaluation.
     */
    CGAtomicFunBridge(const std::string& name,
                      CppAD::ADFun<CGB>& fun,
                      bool standAlone = false,
                      bool cacheSparsities = true) :
            CGAbstractAtomicFun<Base>(name, standAlone),
            fun_(fun),
            cacheSparsities_(cacheSparsities) {
    }

    CGAtomicFunBridge(const CGAtomicFunBridge& orig) = delete;

    CGAtomicFunBridge& operator=(const CGAtomicFunBridge& rhs) = delete;

    virtual ~CGAtomicFunBridge() = default;

    using CGAbstractAtomicFun<Base>::operator();

    template<class VectorSize>
    inline void setCustomSparseJacobianElements(const VectorSize& row,
                                                const VectorSize& col) {
        custom_jac_ = CustomPosition(fun_.Range(), fun_.Domain(), row, col);
    }

    template<class VectorSet>
    inline void setCustomSparseJacobianElements(const VectorSet& elements) {
        custom_jac_ = CustomPosition(fun_.Range(), fun_.Domain(), elements);
    }

    template<class VectorSize>
    inline void setCustomSparseHessianElements(const VectorSize& row,
                                               const VectorSize& col) {
        size_t n = fun_.Domain();
        custom_hess_ = CustomPosition(n, n, row, col);
    }

    template<class VectorSet>
    inline void setCustomSparseHessianElements(const VectorSet& elements) {
        size_t n = fun_.Domain();
        custom_hess_ = CustomPosition(n, n, elements);
    }

    size_t size_dyn_ind() const override {
        return fun_.size_dyn_ind();
    }

    bool for_type(const CppAD::vector<CGB>& parameter_x,
                  const CppAD::vector<ad_type_enum>& type_x,
                  CppAD::vector<ad_type_enum>& type_y) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == fun_.Domain() + fun_.size_dyn_ind(), "type_x has the wrong size")
        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")

        this->validateTypeX(type_x);

        //size_t m = fun_.Range();
        size_t n = fun_.Domain();
        //size_t ndyn = type_x.size();

        ArrayView<const ad_type_enum> type_x_local(type_x.data(), n);

        if (cacheSparsities_ || custom_jac_.isFilterDefined() || custom_jac_.isFullDefined()) {
            if (!custom_jac_.isFullDefined()) {
                custom_jac_.setFullElements(CppAD::cg::jacobianSparsity<CGB>(fun_));
                fun_.size_forward_set(0);
            }

            CppAD::cg::for_type(custom_jac_.getFullElements(), type_x_local, type_y);

        } else {
            CppAD::cg::for_type(fun_, type_x_local, type_y);
        }

        return true;
    }

    bool rev_depend(const CppAD::vector<CGB>& parameter_x,
                    const CppAD::vector<ad_type_enum>& type_x,
                    CppAD::vector<bool>& depend_x,
                    const CppAD::vector<bool>& depend_y) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == fun_.Domain() + fun_.size_dyn_ind(), "type_x has the wrong size")
        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        this->validateTypeX(type_x);

        //size_t m = fun_.Range();
        size_t n = fun_.Domain();
        size_t ndyn = type_x.size();

        ArrayView<const ad_type_enum> type_x_local(type_x.data(), n);
        ArrayView<bool> depend_x_local(depend_x.data(), n);

        if (cacheSparsities_ || custom_jac_.isFilterDefined() || custom_jac_.isFullDefined()) {
            if (!custom_jac_.isFullDefined()) {
                custom_jac_.setFullElements(CppAD::cg::jacobianSparsity<CGB>(fun_));
                fun_.size_forward_set(0);
            }
            CppAD::cg::rev_depend(custom_jac_.getFullElements(), type_x_local, depend_x_local, depend_y);

        } else {
            CppAD::cg::rev_depend(fun_, type_x_local, depend_x_local, depend_y);
        }

        // model_ does not contain the parameters as independent variables
        for (size_t j = n; j < ndyn; ++j)
            depend_x[j] = false;

        return true;
    }

    using CGAbstractAtomicFun<Base>::jacobianSparsity;

    virtual sparse_rc<CppAD::vector<size_t>> jacobianSparsity(const CppAD::vector<CGB>& parameter_x) {

        CPPADCG_ASSERT_KNOWN(parameter_x.size() == fun_.Domain() + fun_.size_dyn_ind(), "parameter_x has the wrong size")

        size_t m = fun_.Range();

        CppAD::vector<ad_type_enum> type_x = makeDefaultTypeX();

        return jacobianSparsity(m, parameter_x, type_x);
    }

    bool jac_sparsity(const CppAD::vector<CGB>& parameter_x,
                      const CppAD::vector<ad_type_enum>& type_x,
                      bool dependency,
                      const CppAD::vector<bool>& select_x,
                      const CppAD::vector<bool>& select_y,
                      CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == fun_.Domain() + fun_.size_dyn_ind(), "type_x has the wrong size")
        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        this->validateTypeX(type_x);

        size_t n = fun_.Domain();
        size_t m = fun_.Range();
        size_t ndyn = type_x.size();

        ArrayView<const bool> select_x_local(select_x.data(), n);

        if (cacheSparsities_ || custom_jac_.isFilterDefined()) {
            if (!custom_jac_.isFullDefined()) {
                custom_jac_.setFullElements(CppAD::cg::jacobianSparsity<CGB>(fun_));
            }

            const auto& fullJac = custom_jac_.getFullElements();

            CppAD::cg::filter(fullJac, select_y, select_x_local, pattern_out);

        } else {
            bool transpose = false;
            CppAD::vector<bool> select_x_local_cppad(select_x_local.size());
            std::copy(select_x_local.begin(), select_x_local.end(), select_x_local_cppad.data());

            fun_.subgraph_sparsity(select_x_local_cppad, select_y, transpose, pattern_out);
        }

        pattern_out.resize(m, ndyn, pattern_out.nnz());

        return true;
    }

    using CGAbstractAtomicFun<Base>::hessianSparsity;

    virtual sparse_rc<CppAD::vector<size_t>> hessianSparsity(const CppAD::vector<CGB>& parameter_x) {

        CPPADCG_ASSERT_KNOWN(parameter_x.size() == fun_.Domain() + fun_.size_dyn_ind(), "parameter_x has the wrong size")

        size_t m = fun_.Range();

        CppAD::vector<ad_type_enum> type_x = makeDefaultTypeX();

        return hessianSparsity(m, parameter_x, type_x);
    }

    bool hes_sparsity(const CppAD::vector<CGB>& parameter_x,
                      const CppAD::vector<ad_type_enum>& type_x,
                      const CppAD::vector<bool>& select_x,
                      const CppAD::vector<bool>& select_y,
                      CppAD::sparse_rc<CppAD::vector<size_t> >& pattern_out) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == fun_.Domain() + fun_.size_dyn_ind(), "type_x has the wrong size")
        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        this->validateTypeX(type_x);

        size_t n = fun_.Domain();
        //size_t m = fun_.Range();
        size_t ndyn = type_x.size();

        ArrayView<const bool> select_x_local(select_x.data(), n);

        if (cacheSparsities_ || custom_hess_.isFilterDefined()) {
            if (!custom_hess_.isFullDefined()) {
                custom_hess_.setFullElements(CppAD::cg::hessianSparsity<CGB>(fun_));
            }

            const auto& fullHess = custom_hess_.getFullElements();

            CppAD::cg::filter(fullHess, select_x_local, pattern_out);

        } else {
            bool internal_bool = true;
            CppAD::vector<bool> select_x_local_cppad(select_x_local.size());
            std::copy(select_x_local.begin(), select_x_local.end(), select_x_local_cppad.data());

            fun_.for_hes_sparsity(select_x_local_cppad, select_y, internal_bool, pattern_out);
        }

        pattern_out.resize(ndyn, ndyn, pattern_out.nnz());

        return true;
    }

protected:

    bool atomicForward(const CppAD::vector<Base>& parameter_x,
                       const CppAD::vector<CppAD::ad_type_enum>& type_x,
                       size_t need_y,
                       size_t order_low,
                       size_t order_up,
                       const CppAD::vector<Base>& taylor_x,
                       CppAD::vector<Base>& taylor_y) override {
        using CppAD::vector;

        size_t n = fun_.Domain();
        size_t npar = fun_.size_dyn_ind();

        vector<CGB> par(npar);
        toCG(ArrayView<const Base>(parameter_x).tail(npar), par);

        vector<CGB> txcg(n * (order_up + 1));
        toCG(ArrayView<const Base>(taylor_x).head(txcg.size()), txcg);

        fun_.new_dynamic(par);

        vector<CGB> tycg = fun_.Forward(order_up, txcg);
        fromCG(tycg, taylor_y);

        fun_.capacity_order(0);

        return true;
    }

    bool atomicReverse(const CppAD::vector<Base>& parameter_x,
                       const CppAD::vector<CppAD::ad_type_enum>& type_x,
                       size_t order_up,
                       const CppAD::vector<Base>& taylor_x,
                       const CppAD::vector<Base>& taylor_y,
                       CppAD::vector<Base>& partial_x,
                       const CppAD::vector<Base>& partial_y) override {
        using CppAD::vector;

        size_t n = fun_.Domain();
        size_t npar = fun_.size_dyn_ind();

        vector<CGB> par(npar);
        toCG(ArrayView<const Base>(parameter_x).tail(npar), par);

        vector<CGB> txcg(n * (order_up + 1));
        toCG(ArrayView<const Base>(taylor_x).head(txcg.size()), txcg);

        vector<CGB> pycg(partial_y.size());

        toCG(partial_y, pycg);

        fun_.new_dynamic(par);

        fun_.Forward(order_up, txcg);

        vector<CGB> pxcg = fun_.Reverse(order_up + 1, pycg);

        fromCG(pxcg, ArrayView<Base>(partial_x).head(pxcg.size()));
        std::fill(partial_x.data() + pxcg.size(), partial_x.data() + partial_x.size(), Base(0));

        fun_.capacity_order(0);

        return true;
    }

private:

    inline CppAD::vector<ad_type_enum> makeDefaultTypeX() const {
        size_t n = fun_.Domain();
        //size_t m = fun_.Range();
        size_t np = fun_.size_dyn_ind();

        CppAD::vector<ad_type_enum> type_x(n + np);
        std::fill(type_x.data(), type_x.data() + n, ad_type_enum::variable_enum);
        std::fill(type_x.data() + n, type_x.data() + type_x.size(), ad_type_enum::dynamic_enum);

        return type_x;
    }

    inline void convertParams(ArrayView<const CGB> from,
                              ArrayView<CGB> to,
                              size_t order_up) {

        size_t p1 = order_up + 1;

        for (size_t j = 0; j < from.size(); j++) {
            to[j * p1] = from[j];

            for (size_t k = 1; k <= order_up; k++) {
                to[j * p1 + k] = Base(0);
            }
        }
    }

    static void toCG(ArrayView<const Base> from,
                     CppAD::vector<CGB>& to) {
        CPPAD_ASSERT_UNKNOWN(from.size() == to.size())

        for (size_t i = 0; i < from.size(); i++) {
            to[i] = from[i];
        }
    }

    static void fromCG(ArrayView<const CGB> from,
                       ArrayView<Base> to) {
        CPPAD_ASSERT_UNKNOWN(from.size() == to.size())

        for (size_t i = 0; i < from.size(); i++) {
            CPPADCG_ASSERT_KNOWN(from[i].isValueDefined(), "No value defined")
            to[i] = from[i].getValue();
        }
    }

};

} // END cg namespace
} // END CppAD namespace

#endif
