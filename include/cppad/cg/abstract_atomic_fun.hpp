#ifndef CPPAD_CG_ABSTRACT_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_ABSTRACT_ATOMIC_FUN_INCLUDED
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
class CGAbstractAtomicFun : public BaseAbstractAtomicFun<Base> {
public:
    using Super = BaseAbstractAtomicFun<Base>;
    using CGB = CppAD::cg::CG<Base>;
    using Arg = Argument<Base>;
protected:
    /**
     * A unique identifier for this atomic function type
     */
    const size_t id_;
    /**
     * Whether or not forward and reverse function calls do not require the
     * Taylor coefficients for the dependent variables (ty) and any previous
     * evaluation of other forward/reverse modes.
     */
    bool standAlone_;

protected:

    /**
     * Creates a new atomic function that is responsible for defining the
     * dependencies to calls of a user atomic function.
     *
     * @param name The atomic function name.
     * @param standAlone Whether or not forward and reverse function calls
     *                   do not require the Taylor coefficients for the
     *                   dependent variables (ty) and any previous
     *                   evaluation of other forward/reverse modes.
     */
    explicit CGAbstractAtomicFun(const std::string& name,
                                 bool standAlone = false) :
        Super(name),
        id_(createNewAtomicFunctionID()),
        standAlone_(standAlone) {
        CPPADCG_ASSERT_KNOWN(!name.empty(), "The atomic function name cannot be empty")
    }

public:
    virtual ~CGAbstractAtomicFun() = default;

    using Super::operator();

    /**
     * Provides a unique identifier for this atomic function type.
     *
     * @return a unique identifier ID
     */
    inline size_t getId() const {
        return id_;
    }

    /**
     * Whether or not forward and reverse function calls do not require the
     * Taylor coefficients for the dependent variables (taylor_y) and any
     * previous evaluation of other forward/reverse modes.
     */
    inline bool isStandAlone() const {
        return standAlone_;
    }

    /**
     * The number of dynamic parameters which are expected at the end of the
     * independent variable vector.
     *
     * @return The number of parameters.
     */
    virtual size_t size_dyn_ind() const = 0;

    bool forward(const CppAD::vector<CGB>& parameter_x,
                 const CppAD::vector<CppAD::ad_type_enum>& type_x,
                 size_t need_y,
                 size_t order_low,
                 size_t order_up,
                 const CppAD::vector<CGB>& taylor_x,
                 CppAD::vector<CGB>& taylor_y) override {
        using ::CppAD::vector;

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() > size_dyn_ind(), "type_x must have the size of the sum of the number of independent variables and independent dynamic parameters")
        validateTypeX(type_x);

        bool valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(taylor_x);

        bool allParameters = BaseAbstractAtomicFun<Base>::isParameters(taylor_x);
        if (allParameters) {
            vector<Base> taylor_y_b;
            if (!evalForwardValues(parameter_x, type_x, need_y, order_low, order_up, taylor_x, taylor_y_b, taylor_y.size()))
                return false;

            CPPADCG_ASSERT_UNKNOWN(taylor_y_b.size() == taylor_y.size())
            for (size_t i = 0; i < taylor_y.size(); i++) {
                taylor_y[i] = taylor_y_b[i];
            }
            return true;
        }

        size_t m = taylor_y.size() / (order_up + 1);
        size_t ndyn = taylor_x.size() / (order_up + 1);
        size_t npar = this->size_dyn_ind();
        size_t n = ndyn - npar;

        vector<bool> vyLocal;
        if (order_up >= 1) {
            /**
             * Use the Jacobian sparsity to determine which elements
             * will always be zero
             */

            vector<bool> select_x(ndyn);
            for (size_t j = 0; j < n; j++) {
                select_x[j]= !taylor_x[j * (order_up + 1) + 1].isIdenticalZero();
            }
            for (size_t j = n; j < ndyn; j++) {
                select_x[j]= false;
            }

            vector<bool> select_y(m);
            for (size_t i = 0; i < m; ++i) select_y[i] = true;

            sparse_rc<vector<size_t> > jac_pattern;
            bool good = this->jac_sparsity(parameter_x, type_x, false, select_x, select_y, jac_pattern);
            if (!good)
                return false;

            if (order_up == 1 && jac_pattern.nnz() == 0) {
                for (size_t i = 0; i < taylor_y.size(); i++) {
                    taylor_y[i] = Base(0.0);
                }
                return true;
            }

            vyLocal.resize(taylor_y.size());
            for (bool& i : vyLocal) {
                i = true;
            }

            size_t nnz = jac_pattern.nnz();
            for (size_t i = 0; i < m; i++) {
                vyLocal[i * (order_up + 1) + 1] = false;
            }

            for (size_t e = 0; e < nnz; e++) {
                size_t i = jac_pattern.row()[e];
                vyLocal[i * (order_up + 1) + 1] = true;
            }

        }

        vector<Base> tyb;
        if (valuesDefined) {
            if (!evalForwardValues(parameter_x, type_x, need_y, order_low, order_up, taylor_x, tyb, taylor_y.size()))
                return false;
        }

        size_t p1 = order_up + 1;

        ArrayView<const CGB> taylor_x_local(taylor_x.data(), n * p1);
        ArrayView<const CGB> pars(parameter_x.data() + n, npar);

        CodeHandler<Base>* handler = findHandler(taylor_x);
        CPPADCG_ASSERT_UNKNOWN(handler != nullptr)

        std::vector<OperationNode<Base>*> txArray(p1), tyArray(p1);
        for (size_t k = 0; k < p1; k++) {
            if (k == 0)
                txArray[k] = BaseAbstractAtomicFun<Base>::makeArray(*handler, taylor_x_local, order_up, k);
            else
                txArray[k] = BaseAbstractAtomicFun<Base>::makeSparseArray(*handler, taylor_x_local, order_up, k);
            tyArray[k] = BaseAbstractAtomicFun<Base>::makeZeroArray(*handler, m);
        }
        OperationNode<Base>* parArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, pars);

        std::vector<Argument<Base> > args(2 * p1 + 1);
        for (size_t k = 0; k < p1; k++) {
            args[0 * p1 + k] = *txArray[k];
            args[1 * p1 + k] = *tyArray[k];
        }
        args[2 * p1] = *parArray;

        OperationNode<Base>* atomicOp = handler->makeNode(CGOpCode::AtomicForward,{id_, order_low, order_up}, args);
        handler->registerAtomicFunction(*this);

        for (size_t k = 0; k < p1; k++) {
            for (size_t i = 0; i < m; i++) {
                size_t pos = i * p1 + k;
                if (vyLocal.size() == 0 || vyLocal[pos]) {
                    taylor_y[pos] = CGB(*handler->makeNode(CGOpCode::ArrayElement, {i}, {*tyArray[k], *atomicOp}));
                    if (valuesDefined) {
                        taylor_y[pos].setValue(tyb[pos]);
                    }
                } else {
                    CPPADCG_ASSERT_KNOWN(tyb.size() == 0 || IdenticalZero(tyb[pos]), "Invalid value")
                    taylor_y[pos] = Base(0.0); // not a variable (zero)
                }
            }
        }

        return true;
    }

    bool reverse(const CppAD::vector<CGB>& parameter_x,
                 const CppAD::vector<ad_type_enum>& type_x,
                 size_t order_up,
                 const CppAD::vector<CGB>& taylor_x,
                 const CppAD::vector<CGB>& taylor_y,
                 CppAD::vector<CGB>& partial_x,
                 const CppAD::vector<CGB>& partial_y) override {
        using CppAD::vector;

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() > size_dyn_ind(), "type_x must have the size of the sum of the number of independent variables and independent dynamic parameters")
        validateTypeX(type_x);

        bool allParameters = BaseAbstractAtomicFun<Base>::isParameters(taylor_x) &&
                             BaseAbstractAtomicFun<Base>::isParameters(taylor_y) &&
                             BaseAbstractAtomicFun<Base>::isParameters(partial_y);

        if (allParameters) {
            vector<Base> pxb;

            if (!evalReverseValues(parameter_x, type_x, order_up, taylor_x, taylor_y, pxb, partial_y))
                return false;

            CPPADCG_ASSERT_UNKNOWN(pxb.size() == partial_x.size())

            for (size_t i = 0; i < partial_x.size(); i++) {
                partial_x[i] = pxb[i];
            }
            return true;
        }

        /**
         * Use the Jacobian sparsity to determine which elements
         * will always be zero
         */

        size_t p1 = order_up + 1;
        // k == 0
        size_t m = taylor_y.size() / p1;
        size_t ndyn = taylor_x.size() / p1;
        size_t npar = this->size_dyn_ind();
        size_t n = ndyn - npar;

        vector<bool> select_y(m);
        for (size_t i = 0; i < m; i++) {
            select_y[i] = !partial_y[i * p1 + order_up].isIdenticalZero(); ////// TODO: verify index!
        }

        vector<bool> select_x(ndyn);
        for (size_t j = 0; j < n; j++) {
            select_x[j] = true;
        }
        for (size_t j = n; j < ndyn; j++) {
            select_x[j] = false;
        }

        vector<CGB> x(n);
        for (size_t j = 0; j < n; j++) {
            x[j] = taylor_x[j * p1];
        }

        sparse_rc<vector<size_t> > jac_pattern;
        bool good = this->jac_sparsity(parameter_x, type_x, false, select_x, select_y, jac_pattern);
        if (!good)
            return false;

        vector<bool> vxLocal(partial_x.size());

        for (size_t j = 0; j < n; j++) {
            vxLocal[j * p1 + order_up] = false;
        }

        size_t nnz = jac_pattern.nnz();
        for (size_t e = 0; e < nnz; e++) {
            size_t j = jac_pattern.col()[e];
            vxLocal[j * p1 + order_up] = true;
        }

        if (order_up >= 1) {
            /**
             * Use the Hessian sparsity to determine which elements
             * will always be zero
             */

            for (size_t j = 0; j < n; j++) {
                select_x[j] = !taylor_x[j * p1].isIdenticalZero();
            }

            for (size_t i = 0; i < m; i++) {
                select_y[i] = !partial_y[i * p1 + order_up].isIdenticalZero();
            }

            sparse_rc<vector<size_t> > pattern_hess;

            this->hes_sparsity(parameter_x, type_x, select_x, select_y, pattern_hess);

            nnz = pattern_hess.nnz();
            CppAD::vector<std::set<size_t>> pattern_hess_set;
            generateSparsitySet(pattern_hess, pattern_hess_set);
            for (size_t e = 0; e < nnz; e++) {
                size_t j = pattern_hess.row()[e];
                vxLocal[j * p1 + order_up - 1] = true;
            }

        }

        bool allZero = true;
        for (bool j : vxLocal) {
            if (j) {
                allZero = false;
                break;
            }
        }

        if (allZero) {
            for (size_t j = 0; j < partial_x.size(); j++) {
                partial_x[j] = Base(0.0);
            }
            return true;
        }

        bool valuesDefined = BaseAbstractAtomicFun<Base>::isValuesDefined(taylor_x) &&
                             BaseAbstractAtomicFun<Base>::isValuesDefined(taylor_y) &&
                             BaseAbstractAtomicFun<Base>::isValuesDefined(partial_y);

        vector<Base> pxb;
        if (valuesDefined) {
            if (!evalReverseValues(parameter_x, type_x, order_up, taylor_x, taylor_y, pxb, partial_y))
                return false;
        }

        CodeHandler<Base>* handler = findHandler(taylor_x);
        if (handler == nullptr) {
            handler = findHandler(taylor_y);
            if (handler == nullptr) {
                handler = findHandler(partial_y);
            }
        }
        CPPADCG_ASSERT_UNKNOWN(handler != nullptr)

        ArrayView<const CGB> taylor_x_local(taylor_x.data(), n * p1);
        ArrayView<const CGB> pars(parameter_x.data() + n, npar);

        std::vector<OperationNode<Base>*> txArray(p1), tyArray(p1), pxArray(p1), pyArray(p1);
        for (size_t k = 0; k <= order_up; k++) {
            if (k == 0)
                txArray[k] = BaseAbstractAtomicFun<Base>::makeArray(*handler, taylor_x_local, order_up, k);
            else
                txArray[k] = BaseAbstractAtomicFun<Base>::makeSparseArray(*handler, taylor_x_local, order_up, k);

            if (standAlone_) {
                tyArray[k] = BaseAbstractAtomicFun<Base>::makeEmptySparseArray(*handler, m);
            } else {
                tyArray[k] = BaseAbstractAtomicFun<Base>::makeSparseArray(*handler, taylor_y, order_up, k);
            }

            if (k == 0)
                pxArray[k] = BaseAbstractAtomicFun<Base>::makeZeroArray(*handler, n);
            else
                pxArray[k] = BaseAbstractAtomicFun<Base>::makeEmptySparseArray(*handler, n);

            if (k == 0)
                pyArray[k] = BaseAbstractAtomicFun<Base>::makeSparseArray(*handler, partial_y, order_up, k);
            else
                pyArray[k] = BaseAbstractAtomicFun<Base>::makeArray(*handler, partial_y, order_up, k);
        }

        OperationNode<Base>* parArray = BaseAbstractAtomicFun<Base>::makeArray(*handler, pars);

        std::vector<Argument<Base> > args(4 * p1 + 1);
        for (size_t k = 0; k <= order_up; k++) {
            args[0 * p1 + k] = *txArray[k];
            args[1 * p1 + k] = *tyArray[k];
            args[2 * p1 + k] = *pxArray[k];
            args[3 * p1 + k] = *pyArray[k];
        }
        args[4 * p1] = *parArray;

        OperationNode<Base>* atomicOp = handler->makeNode(CGOpCode::AtomicReverse,{id_, order_up}, args);
        handler->registerAtomicFunction(*this);

        for (size_t k = 0; k < p1; k++) {
            for (size_t j = 0; j < ndyn; j++) {
                size_t pos = j * p1 + k;
                if (vxLocal[pos]) {
                    partial_x[pos] = CGB(*handler->makeNode(CGOpCode::ArrayElement, {j}, {*pxArray[k], *atomicOp}));
                    if (valuesDefined) {
                        partial_x[pos].setValue(pxb[pos]);
                    }
                } else {
                    // CPPADCG_ASSERT_KNOWN(pxb.size() == 0 || IdenticalZero(pxb[j]), "Invalid value")
                    // pxb[j] might be non-zero but it is not required (it might have been used to determine other pxbs)
                    partial_x[pos] = Base(0); // not a variable (zero)
                }
            }
        }

        return true;
    }

    virtual sparse_rc<CppAD::vector<size_t>> jacobianSparsity(size_t m,
                                                              const CppAD::vector<CGB>& parameter_x,
                                                              const CppAD::vector<ad_type_enum>& type_x) {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() > size_dyn_ind(), "type_x must have the size of the sum of the number of independent variables and independent dynamic parameters")
        validateTypeX(type_x);

        size_t n = type_x.size();

        CppAD::vector<bool> select_x(n);
        CppAD::vector<bool> select_y(m);

        for (size_t j = 0; j < n; j++) {
            select_x[j]= true;
        }

        for (size_t i = 0; i < m; i++) {
            select_y[i]= true;
        }

        sparse_rc<CppAD::vector<size_t> > jac_pattern;

        bool dependency = false;

        bool good = this->jac_sparsity(parameter_x, type_x, dependency, select_x, select_y, jac_pattern);
        if (!good)
            throw CGException("Failed to compute Jacobian sparsity pattern for atomic function '", this->atomic_name(), "'");

        return jac_pattern;
    }

    virtual sparse_rc<CppAD::vector<size_t> > hessianSparsity(size_t m,
                                                              const CppAD::vector<CGB>& parameter_x,
                                                              const CppAD::vector<ad_type_enum>& type_x) {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        CPPADCG_ASSERT_KNOWN(type_x.size() > size_dyn_ind(), "type_x must have the size of the sum of the number of independent variables and independent dynamic parameters")
        validateTypeX(type_x);

        size_t n = type_x.size();

        CppAD::vector<bool> select_x(n);
        CppAD::vector<bool> select_y(m);

        for (size_t j = 0; j < n; j++) {
            select_x[j]= true;
        }

        for (size_t i = 0; i < m; i++) {
            select_y[i]= true;
        }

        sparse_rc<CppAD::vector<size_t> > pattern_hess;

        bool good = this->hes_sparsity(parameter_x, type_x, select_x, select_y, pattern_hess);
        if (!good)
            throw CGException("Failed to compute Hessian sparsity pattern for atomic function '",
                              this->atomic_name(), "'");

        return pattern_hess;
    }

    /**
     * Uses an internal counter to produce IDs for atomic functions.
     */
    static size_t createNewAtomicFunctionID() {
        CPPAD_ASSERT_FIRST_CALL_NOT_PARALLEL
        static size_t count = 0;
        count++;
        return count;
    }

protected:

    /**
     * Used to evaluate function values and forward mode function values and
     * derivatives.
     *
     * @param parameter_x The values, in afun(ax, ay), for arguments that are
     *                    parameters.
     * @param type_x
     * @param need_y
     * @param order_low Lowest order for this forward mode calculation.
     * @param order_up Highest order for this forward mode calculation.
     * @param taylor_x Taylor coefficients corresponding to \c x for this
     *                 calculation
     * @param taylor_y Taylor coefficient corresponding to \c y for this
     *                 calculation
     * @return true on success, false otherwise
     */
    virtual bool atomicForward(const CppAD::vector<Base>& parameter_x,
                               const CppAD::vector<CppAD::ad_type_enum>& type_x,
                               size_t need_y,
                               size_t order_low,
                               size_t order_up,
                               const CppAD::vector<Base>& taylor_x,
                               CppAD::vector<Base>& taylor_y) = 0;
    /**
     * Used to evaluate reverse mode function derivatives.
     *
     * @param parameter_x The values, in afun(ax, ay), for arguments that are
     *                    parameters.
     * @param type_x
     * @param order_up Highest order for this forward mode calculation.
     * @param taylor_x Taylor coefficients corresponding to \c x for this
     *                 calculation
     * @param taylor_y Taylor coefficient corresponding to \c y for this
     *                 calculation
     * @param partial_x Partials w.r.t. the \c x Taylor coefficients.
     * @param partial_y Partials w.r.t. the \c y Taylor coefficients
     * @return true on success, false otherwise
     */
    virtual bool atomicReverse(const CppAD::vector<Base>& parameter_x,
                               const CppAD::vector<CppAD::ad_type_enum>& type_x,
                               size_t order_up,
                               const CppAD::vector<Base>& taylor_x,
                               const CppAD::vector<Base>& taylor_y,
                               CppAD::vector<Base>& partial_x,
                               const CppAD::vector<Base>& partial_y) = 0;

    inline void validateTypeX(const CppAD::vector<ad_type_enum>& type_x) const {
        size_t ndyn = type_x.size();
        size_t npar = this->size_dyn_ind();

        CPPADCG_ASSERT_KNOWN(type_x.size() >= npar, "type_x must have at least the same size as type_x")

        for (size_t j = ndyn - npar; j < ndyn; ++j) {
            if (type_x[j] != ad_type_enum::constant_enum && type_x[j] != ad_type_enum::dynamic_enum) {
                throw CGException("Wrong variable type at ", j, ". Expected either constant_enum or dynamic_enum.");
            }
        }
    }

private:

    inline bool evalForwardValues(const CppAD::vector<CGB>& parameter_x,
                                  const CppAD::vector<CppAD::ad_type_enum>& type_x,
                                  size_t need_y,
                                  size_t order_low,
                                  size_t order_up,
                                  const CppAD::vector<CGB>& taylor_x,
                                  CppAD::vector<Base>& taylor_yb,
                                  size_t ty_size) {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        validateTypeX(type_x);

        CppAD::vector<Base> parb(parameter_x.size());
        CppAD::vector<Base> txb(taylor_x.size());
        taylor_yb.resize(ty_size);

        for (size_t i = 0; i < parb.size(); i++) {
            parb[i] = parameter_x[i].getValue();
        }
        for (size_t i = 0; i < txb.size(); i++) {
            txb[i] = taylor_x[i].getValue();
        }

        return atomicForward(parb, type_x, need_y, order_low, order_up, txb, taylor_yb);
    }

    inline bool evalReverseValues(const CppAD::vector<CGB>& parameter_x,
                                  const CppAD::vector<CppAD::ad_type_enum>& type_x,
                                  size_t order_up,
                                  const CppAD::vector<CGB>& taylor_x,
                                  const CppAD::vector<CGB>& taylor_y,
                                  CppAD::vector<Base>& pxb,
                                  const CppAD::vector<CGB>& py) {

        CPPADCG_ASSERT_KNOWN(type_x.size() == parameter_x.size(), "type_x must have the same size as parameter_x")
        validateTypeX(type_x);

        using CppAD::vector;

        vector<Base> parb(parameter_x.size());
        vector<Base> txb(taylor_x.size());
        vector<Base> tyb(taylor_y.size());
        pxb.resize(taylor_x.size());
        vector<Base> pyb(py.size());

        for (size_t i = 0; i < parb.size(); i++) {
            parb[i] = parameter_x[i].getValue();
        }
        for (size_t i = 0; i < txb.size(); i++) {
            txb[i] = taylor_x[i].getValue();
        }
        for (size_t i = 0; i < tyb.size(); i++) {
            tyb[i] = taylor_y[i].getValue();
        }
        for (size_t i = 0; i < pyb.size(); i++) {
            pyb[i] = py[i].getValue();
        }

        return atomicReverse(parb, type_x, order_up, txb, tyb, pxb, pyb);
    }

};

} // END cg namespace
} // END CppAD namespace

#endif
