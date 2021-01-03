#ifndef CPPAD_CG_ATOMIC_EXTERNAL_FUNCTION_INCLUDED
#define CPPAD_CG_ATOMIC_EXTERNAL_FUNCTION_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2018 Joao Leal
 *    Copyright (C) 2014 Ciengis
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

template<class Base>
class AtomicExternalFunctionWrapper : public ExternalFunctionWrapper<Base> {
private:
    atomic_three<Base>* atomic_;
public:

    inline explicit AtomicExternalFunctionWrapper(atomic_three<Base>& atomic) :
        atomic_(&atomic) {
    }

    inline virtual ~AtomicExternalFunctionWrapper() = default;

    bool forward(FunctorGenericModel<Base>& libModel,
                 int order_low,
                 int order_up,
                 const Array tx[],
                 const Array& params,
                 Array& ty) override {
        size_t m = ty.size;
        size_t n = tx[0].size;
        size_t np = params.size;
        size_t ndyn = n + np;

        prepareParameterX(libModel, n, params);
        prepareTypeX(libModel, n, np);
        auto need_y = size_t(ad_type_enum::variable_enum);

        convert(tx, libModel._tx, ndyn, n, order_up, order_up + 1);
        convertParams(params, libModel._tx, n, order_up);

        zero(libModel._ty, m * (order_up + 1));

        bool ret = atomic_->forward(libModel._par_x, libModel._type_x, need_y, order_low, order_up, libModel._tx, libModel._ty);

        convertAdd(libModel._ty, ty, m, order_up, order_up);

        return ret;
    }

    bool reverse(FunctorGenericModel<Base>& libModel,
                 int order_up,
                 const Array tx[],
                 const Array& params,
                 Array& px,
                 const Array py[]) override {
        size_t m = py[0].size;
        size_t n = tx[0].size;
        size_t np = params.size;
        size_t ndyn = n + np;

        prepareParameterX(libModel, n, params);
        prepareTypeX(libModel, n, np);

        convert(tx, libModel._tx, ndyn, n, order_up, order_up + 1);
        convertParams(params, libModel._tx, n, order_up);

        zero(libModel._ty, m * (order_up + 1));

        convert(py, libModel._py, m, m, order_up, order_up + 1);

        zero(libModel._px, ndyn * (order_up + 1));

#ifndef NDEBUG
        if (libModel._evalAtomicForwardOne4CppAD) {
            // only required in order to avoid an issue with a validation inside CppAD
            auto need_y = size_t(ad_type_enum::variable_enum);

            if (!atomic_->forward(libModel._par_x, libModel._type_x, need_y, order_up, order_up, libModel._tx, libModel._ty))
                return false;
        }
#endif
        bool ret = atomic_->reverse(libModel._par_x, libModel._type_x, order_up, libModel._tx, libModel._ty, libModel._px, libModel._py);

        convertAdd(libModel._px, px, n, order_up, 0); // k=0 for both order_up=0 and order_up=1

        return ret;
    }

private:

    inline void convert(const Array from[],
                        CppAD::vector<Base>& to,
                        size_t n,
                        size_t nMax,
                        size_t order_up,
                        size_t kmax) {
        assert(n >= nMax);

        size_t p1 = order_up + 1;
        to.resize(n * p1);

        for (size_t k = 0; k < kmax; k++) {
            Base* values = static_cast<Base*> (from[k].data);
            if (from[k].sparse) {
                if (order_up == 0) {
                    std::fill(to.data(), to.data() + nMax, Base(0));
                } else {
                    for (size_t j = 0; j < nMax; j++) {
                        to[j * p1 + k] = Base(0);
                    }
                }

                for (size_t e = 0; e < from[k].nnz; e++) {
                    size_t j = from[k].idx[e];
                    to[j * p1 + k] = values[e];
                }
            } else {
                for (size_t j = 0; j < nMax; j++) {
                    to[j * p1 + k] = values[j];
                }
            }
        }
    }

    inline void convertParams(const Array& from,
                              CppAD::vector<Base>& to,
                              size_t nMax,
                              size_t order_up) {

        size_t p1 = order_up + 1;

        Base* values = static_cast<Base*> (from.data);
        for (size_t j = 0; j < from.size; j++) {
            to[nMax * p1 + j * p1] = values[j];

            for (size_t k = 1; k <= order_up; k++) {
                to[nMax * p1 + j * p1 + k] = Base(0);
            }
        }
    }

    inline static void prepareParameterX(FunctorGenericModel<Base>& libModel,
                                         size_t n,
                                         const Array& params) {
        CPPADCG_ASSERT_UNKNOWN(params.sparse == 0)

        libModel._par_x.resize(n + params.size);
        auto* data = static_cast<Base*>(params.data);
        std::copy(data, data + params.size, libModel._par_x.data() + n);
    }

    inline static void prepareTypeX(FunctorGenericModel<Base>& libModel,
                                    size_t n,
                                    size_t np) {
        libModel._type_x.resize(n + np);
        std::fill(libModel._type_x.data(), libModel._type_x.data() + n, ad_type_enum::variable_enum);
        std::fill(libModel._type_x.data() + n, libModel._type_x.data() + n + np, ad_type_enum::dynamic_enum);
    }

    inline static void zero(CppAD::vector<Base>& to,
                            size_t size) {
        to.resize(size);
        std::fill(to.data(), to.data() + to.size(), Base(0));
    }

    inline static void convertDense(const Array& from,
                                    CppAD::vector<Base>& to) {
        CPPADCG_ASSERT_UNKNOWN(from.sparse == 0)
        Base* values = static_cast<Base*> (from.data);

        to.resize(from.size);

        std::copy(values, values + from.size, to.data());
    }

    inline static void convertAdd(const CppAD::vector<Base>& from,
                                  Array& to,
                                  size_t n,
                                  size_t p,
                                  size_t k) {
        CPPADCG_ASSERT_KNOWN(!to.sparse, "output must be a dense array")
        CPPADCG_ASSERT_KNOWN(to.size >= n, "invalid size")

        Base* values = static_cast<Base*> (to.data);

        if (p == 0) {
            std::copy(&from[0], &from[0] + n, values);
        } else {
            size_t p1 = p + 1;

            for (size_t j = 0; j < n; j++) {
                values[j] += from[j * p1 + k];
            }
        }

    }
};

} // END cg namespace
} // END CppAD namespace

#endif
