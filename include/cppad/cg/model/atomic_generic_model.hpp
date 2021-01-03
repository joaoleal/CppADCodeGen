#ifndef CPPAD_CG_ATOMIC_GENERIC_MODEL_INCLUDED
#define CPPAD_CG_ATOMIC_GENERIC_MODEL_INCLUDED
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
 * An atomic function that uses a compiled model created by CppADCodeGen.
 *
 * When the compiled models contains dynamic parameters, they must be
 * provided at the end of the independent variable vector.
 *
 * @author Joao Leal
 */
template <class Base>
class CGAtomicGenericModel : public atomic_three<Base> {
protected:
    GenericModel<Base>& model_;
public:

    /**
     * Creates a new atomic function wrapper that is responsible for
     * calling the appropriate methods of the compiled model.
     *
     * @param model The compiled model.
     */
    explicit CGAtomicGenericModel(GenericModel<Base>& model) :
        atomic_three<Base>(model.getName()),
        model_(model) {

    }

    virtual ~CGAtomicGenericModel() = default;

    template <class ADVector>
    void operator()(const ADVector& ax, ADVector& ay) {
        this->atomic_three<Base>::operator()(ax, ay);
    }

    bool for_type(const CppAD::vector<Base>& parameter_x,
                  const CppAD::vector<ad_type_enum>& type_x,
                  CppAD::vector<ad_type_enum>& type_y) override {

        if (!model_.isJacobianSparsityAvailable())
            return false;

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() == model_.Domain() + model_.Parameters(), "type_x has the wrong size")
        validateTypeX(type_x);

        CppAD::sparse_rc<CppAD::vector<size_t>> jac_pattern;
        model_.JacobianSparsity(jac_pattern);

        //size_t m = model_.Range();
        size_t n = model_.Domain();
        //size_t ndyn = type_x.size();

        ArrayView<const ad_type_enum> type_x_local(type_x.data(), n);

        CppAD::cg::for_type(jac_pattern, type_x_local, type_y);

        return true;
    }

    bool rev_depend(const CppAD::vector<Base>& parameter_x,
                    const CppAD::vector<ad_type_enum>& type_x,
                    CppAD::vector<bool>& depend_x,
                    const CppAD::vector<bool>& depend_y) override {

        if (!model_.isJacobianSparsityAvailable())
            return false;

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() == model_.Domain() + model_.Parameters(), "type_x has the wrong size")
        validateTypeX(type_x);

        CppAD::sparse_rc<CppAD::vector<size_t>> jac_pattern;
        model_.JacobianSparsity(jac_pattern);

        //size_t m = model_.Range();
        size_t n = model_.Domain();
        size_t ndyn = type_x.size();

        ArrayView<const ad_type_enum> type_x_local(type_x.data(), n);
        ArrayView<bool> depend_x_local(depend_x.data(), n);

        CppAD::cg::rev_depend(jac_pattern, type_x_local, depend_x, depend_y);

        // model_ does not contain the parameters as independent variables
        for (size_t j = n; j < ndyn; ++j)
            depend_x[j] = false;

        return true;
    }

    bool forward(const CppAD::vector<Base>& parameter_x,
                 const CppAD::vector<ad_type_enum>& type_x,
                 size_t need_y,
                 size_t order_low,
                 size_t order_up,
                 const CppAD::vector<Base>& taylor_x,
                 CppAD::vector<Base>& taylor_y) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() == model_.Domain() + model_.Parameters(), "type_x has the wrong size")
        validateTypeX(type_x);

        size_t n = model_.Domain();
        size_t np = model_.Parameters();

        ArrayView<const Base> par(parameter_x.data() + n, np);

        if (order_up == 0) {
            ArrayView<const Base> x(taylor_x.data(), n);

            model_.ForwardZero(x, par, taylor_y);
            return true;

        } else if (order_up == 1) {
            ArrayView<const Base> tx(taylor_x.data(), 2 * n);

            model_.ForwardOne(tx, par, taylor_y);
            return true;
        }

        return false;
    }

    bool reverse(const CppAD::vector<Base>& parameter_x,
                 const CppAD::vector<ad_type_enum>& type_x,
                 size_t order_up,
                 const CppAD::vector<Base>& taylor_x,
                 const CppAD::vector<Base>& taylor_y,
                 CppAD::vector<Base>& partial_x,
                 const CppAD::vector<Base>& partial_y) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() == model_.Domain() + model_.Parameters(), "type_x has the wrong size")
        validateTypeX(type_x);

        size_t n = model_.Domain();
        size_t np = model_.Parameters();

        ArrayView<const Base> par(parameter_x.data() + n, np);
        ArrayView<const Base> ty(taylor_y);
        ArrayView<const Base> py(partial_y);

        if (order_up == 0) {
            ArrayView<const Base> x(taylor_x.data(), n);
            ArrayView<Base> px(partial_x.data(), n);

            model_.ReverseOne(x, par, ty, px, py);

            std::fill(px.data() + n, px.data() + px.size(), Base(0));

            return true;

        } else if (order_up == 1) {
            ArrayView<const Base> tx(taylor_x.data(), 2 * n);
            ArrayView<Base> px(partial_x.data(), 2 * n);

            model_.ReverseTwo(tx, par, ty, px, py);

            std::fill(px.data() + 2 * n, px.data() + px.size(), Base(0));

            return true;
        }

        return false;
    }

    bool jac_sparsity(const CppAD::vector<Base>& parameter_x,
                      const CppAD::vector<ad_type_enum>& type_x,
                      bool dependency,
                      const CppAD::vector<bool>& select_x,
                      const CppAD::vector<bool>& select_y,
                      sparse_rc<CppAD::vector<size_t>>& pattern_out) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() == model_.Domain() + model_.Parameters(), "type_x has the wrong size")
        validateTypeX(type_x);

        model_.JacobianSparsity(pattern_out);

        size_t n = model_.Domain();
        size_t m = model_.Range();
        size_t ndyn = type_x.size();

        ArrayView<const bool> select_x_local(select_x.data(), n);

        filter(pattern_out, select_y, select_x_local);

        // model_ does not contain the parameters as independent variables
        pattern_out.resize(m, ndyn, pattern_out.nnz());

        return true;
    }

    bool hes_sparsity(const CppAD::vector<Base>& parameter_x,
                      const CppAD::vector<ad_type_enum>& type_x,
                      const CppAD::vector<bool>& select_x,
                      const CppAD::vector<bool>& select_y,
                      sparse_rc<CppAD::vector<size_t>>& pattern_out) override {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() == model_.Domain() + model_.Parameters(), "type_x has the wrong size")
        validateTypeX(type_x);

        size_t n = model_.Domain();
        size_t m = model_.Range();
        size_t ndyn = type_x.size();

        ArrayView<const bool> select_x_local(select_x.data(), n);

        if (isAllTrue(select_y)) {
            model_.HessianSparsity(pattern_out);
            filter(pattern_out, select_x_local);

        } else {
            std::vector<std::set<size_t> > sparsitySF2R(n);
            for (size_t i = 0; i < m; i++) {
                if (select_y[i]) {
                    CppAD::cg::multMatrixSparsity(model_.HessianSparsitySet(i), select_x_local, sparsitySF2R); // f''_i(x)
                }
            }

            // filter(sparsitySF2R, select_x, pattern_out);
        }

        // model_ does not contain the parameters as independent variables
        pattern_out.resize(ndyn, ndyn, pattern_out.nnz());

        return true;
    }

private:

    inline void validateTypeX(const CppAD::vector<ad_type_enum>& type_x) const {
        size_t ndyn = type_x.size();
        size_t n = model_.Domain();
        for (size_t j = n; j < ndyn; ++j) {
            if (type_x[j] != ad_type_enum::constant_enum && type_x[j] != ad_type_enum::dynamic_enum) {
                throw CGException("Wrong variable type at ", j, ". Expected either constant_enum or dynamic_enum.");
            }
        }
    }

};

} // END cg namespace
} // END CppAD namespace

#endif
