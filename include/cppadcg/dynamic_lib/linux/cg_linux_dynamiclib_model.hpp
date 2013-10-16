#ifndef CPPAD_CG_LINUX_DYNAMICLIB_MODEL_INCLUDED
#define CPPAD_CG_LINUX_DYNAMICLIB_MODEL_INCLUDED
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
#ifdef __linux__

#include <typeinfo>
#include <dlfcn.h>

namespace CppAD {

    /**
     * Useful class to call the compiled model in a dynamic library.
     * For the Linux Operating System only.
     * 
     * @author Joao Leal
     */
    template<class Base>
    class LinuxDynamicLibModel : public DynamicLibModel<Base> {
    protected:
        /// the model name
        const std::string _name;
        /// the dynamic library
        LinuxDynamicLib<Base>* _dynLib;
        size_t _m;
        size_t _n;
        std::vector<const Base*> _in;
        std::vector<const Base*> _inHess;
        std::vector<Base*> _out;
        CLangAtomicFun _atomicFuncArg;
        std::vector<atomic_base<Base>* > _atomic;
        size_t _missingAtomicFunctions;
        vector<Base> _tx, _ty, _px, _py;
        // original model function
        void (*_zero)(Base const*const*, Base * const*, CLangAtomicFun);
        // first order forward mode
        int (*_forwardOne)(Base const tx[], Base ty[], CLangAtomicFun);
        // first order reverse mode
        int (*_reverseOne)(Base const tx[], Base const ty[], Base px[], Base const py[], CLangAtomicFun);
        // second order reverse mode
        int (*_reverseTwo)(Base const tx[], Base const ty[], Base px[], Base const py[], CLangAtomicFun);
        // jacobian function in the dynamic library
        void (*_jacobian)(Base const*const*, Base * const*, CLangAtomicFun);
        // hessian function in the dynamic library
        void (*_hessian)(Base const*const*, Base * const*, CLangAtomicFun);
        //
        int (*_sparseForwardOne)(unsigned long, Base const *const *, Base * const *, CLangAtomicFun);
        //
        int (*_sparseReverseOne)(unsigned long, Base const *const *, Base * const *, CLangAtomicFun);
        //
        int (*_sparseReverseTwo)(unsigned long, Base const *const *, Base * const *, CLangAtomicFun);
        // sparse jacobian function in the dynamic library
        void (*_sparseJacobian)(Base const*const*, Base * const*, CLangAtomicFun);
        // sparse hessian function in the dynamic library
        void (*_sparseHessian)(Base const*const*, Base * const*, CLangAtomicFun);
        //
        void (*_forwardOneSparsity)(unsigned long, unsigned long const**, unsigned long*, CLangAtomicFun);
        //
        void (*_reverseOneSparsity)(unsigned long, unsigned long const**, unsigned long*, CLangAtomicFun);
        //
        void (*_reverseTwoSparsity)(unsigned long, unsigned long const**, unsigned long*, CLangAtomicFun);
        // jacobian sparsity function in the dynamic library
        void (*_jacobianSparsity)(unsigned long const** row,
                unsigned long const** col,
                unsigned long * nnz);
        // hessian sparsity function in the dynamic library
        void (*_hessianSparsity)(unsigned long const** row,
                unsigned long const** col,
                unsigned long * nnz);
        void (*_hessianSparsity2)(unsigned long i,
                unsigned long const** row,
                unsigned long const** col,
                unsigned long * nnz);
        void (*_atomicFunctions)(const char*** names,
                unsigned long * n);

    public:

        virtual const std::string& getName() const {
            return _name;
        }

        virtual bool addAtomicFunction(atomic_base<Base>& atomic) {
            const char** names;
            unsigned long n;
            (*_atomicFunctions)(&names, &n);

            assert(_atomic.size() == n);

            for (unsigned long i = 0; i < n; i++) {
                if (atomic.afun_name() == names[i]) {
                    if (_atomic[i] == NULL) {
                        _missingAtomicFunctions--;
                    }
                    _atomic[i] = &atomic;
                    return true;
                }
            }
            return false;
        }

        // Jacobian sparsity

        virtual std::vector<bool> JacobianSparsityBool() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            bool set_type = true;
            std::vector<bool> s;

            loadSparsity(set_type, s, _m, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > JacobianSparsitySet() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _m, _n, row, col, nnz);

            return s;
        }

        virtual void JacobianSparsity(std::vector<size_t>& rows, std::vector<size_t>& cols) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            rows.resize(nnz);
            cols.resize(nnz);

            std::copy(row, row + nnz, rows.begin());
            std::copy(col, col + nnz, cols.begin());
        }

        // Hessian sparsity 

        virtual std::vector<bool> HessianSparsityBool() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            bool set_type = true;
            std::vector<bool> s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > HessianSparsitySet() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual void HessianSparsity(std::vector<size_t>& rows, std::vector<size_t>& cols) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            rows.resize(nnz);
            cols.resize(nnz);

            std::copy(row, row + nnz, rows.begin());
            std::copy(col, col + nnz, cols.begin());
        }

        virtual std::vector<bool> HessianSparsityBool(size_t i) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity2 != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity2)(i, &row, &col, &nnz);

            bool set_type = true;
            std::vector<bool> s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > HessianSparsitySet(size_t i) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity2 != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity2)(i, &row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual void HessianSparsity(size_t i, std::vector<size_t>& rows, std::vector<size_t>& cols) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity2 != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity2)(i, &row, &col, &nnz);

            rows.resize(nnz);
            cols.resize(nnz);

            std::copy(row, row + nnz, rows.begin());
            std::copy(col, col + nnz, cols.begin());
        }

        /// number of independent variables

        virtual size_t Domain() const {
            return _n;
        }

        /// number of dependent variables

        virtual size_t Range() const {
            return _m;
        }

        /// calculate the dependent values (zero order)

        virtual void ForwardZero(const Base* x, size_t x_size, Base* dep, size_t dep_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(dep_size == _m, "Invalid dependent array size");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");


            _in[0] = x;
            _out[0] = dep;

            (*_zero)(&_in[0], &_out[0], _atomicFuncArg);
        }

        virtual void ForwardZero(const std::vector<const Base*> &x,
                                 Base* dep, size_t dep_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == x.size(), "The number of independent variable arrays is invalid");
            CPPADCG_ASSERT_KNOWN(dep_size == _m, "Invalid dependent array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            _out[0] = dep;

            (*_zero)(&x[0], &_out[0], _atomicFuncArg);
        }

        virtual void ForwardZero(const CppAD::vector<bool>& vx,
                                 CppAD::vector<bool>& vy,
                                 const CppAD::vector<Base> &tx,
                                 CppAD::vector<Base>& ty) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(tx.size() == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(ty.size() == _m, "Invalid dependent array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            _in[0] = &tx[0];
            _out[0] = &ty[0];

            (*_zero)(&_in[0], &_out[0], _atomicFuncArg);

            if (vx.size() > 0) {
                CPPADCG_ASSERT_KNOWN(vx.size() >= _n, "Invalid vx size");
                CPPADCG_ASSERT_KNOWN(vy.size() >= _m, "Invalid vy size");
                const std::vector<std::set<size_t> > jacSparsity = JacobianSparsitySet();
                for (size_t i = 0; i < _m; i++) {
                    std::set<size_t>::const_iterator it;
                    for (it = jacSparsity[i].begin(); it != jacSparsity[i].end(); ++it) {
                        size_t j = *it;
                        if (vx[j]) {
                            vy[i] = true;
                            break;
                        }
                    }
                }
            }
        }

        /// calculate entire Jacobian       

        virtual void Jacobian(const Base* x, size_t x_size,
                              Base* jac, size_t jac_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobian != NULL, "No Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(jac_size == _m * _n, "Invalid Jacobian array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");


            _in[0] = x;
            _out[0] = jac;

            (*_jacobian)(&_in[0], &_out[0], _atomicFuncArg);
        }

        /// calculate Hessian for one component of f

        virtual void Hessian(const Base* x, size_t x_size,
                             const Base* w, size_t w_size,
                             Base* hess) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(w_size == _m, "Invalid multiplier array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            _inHess[0] = x;
            _inHess[1] = w;
            _out[0] = hess;

            (*_hessian)(&_inHess[0], &_out[0], _atomicFuncArg);
        }

        virtual void ForwardOne(const Base tx[], size_t tx_size,
                                Base ty[], size_t ty_size) {
            const size_t k = 1;

            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_forwardOne != NULL, "No forward one function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(tx_size >= (k + 1) * _n, "Invalid tx size");
            CPPADCG_ASSERT_KNOWN(ty_size >= (k + 1) * _m, "Invalid ty size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            int ret = (*_forwardOne)(tx, ty, _atomicFuncArg);

            CPPADCG_ASSERT_KNOWN(ret == 0, "First-order forward mode failed."); // generic failure
        }

        virtual void ReverseOne(const Base tx[], size_t tx_size,
                                const Base ty[], size_t ty_size,
                                Base px[], size_t px_size,
                                const Base py[], size_t py_size) {
            const size_t k = 0;
            const size_t k1 = k + 1;

            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_reverseOne != NULL, "No reverse one function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(tx_size >= k1 * _n, "Invalid tx size");
            CPPADCG_ASSERT_KNOWN(ty_size >= k1 * _m, "Invalid ty size");
            CPPADCG_ASSERT_KNOWN(px_size >= k1 * _n, "Invalid px size");
            CPPADCG_ASSERT_KNOWN(py_size >= k1 * _m, "Invalid py size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            int ret = (*_reverseOne)(tx, ty, px, py, _atomicFuncArg);

            CPPADCG_ASSERT_KNOWN(ret == 0, "First-order reverse mode failed.");
        }

        virtual void ReverseTwo(const Base tx[], size_t tx_size,
                                const Base ty[], size_t ty_size,
                                Base px[], size_t px_size,
                                const Base py[], size_t py_size) {
            const size_t k = 1;
            const size_t k1 = k + 1;

            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_reverseTwo != NULL, "No sparse reverse two function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1");
            CPPADCG_ASSERT_KNOWN(tx_size >= k1 * _n, "Invalid tx size");
            CPPADCG_ASSERT_KNOWN(ty_size >= k1 * _m, "Invalid ty size");
            CPPADCG_ASSERT_KNOWN(px_size >= k1 * _n, "Invalid px size");
            CPPADCG_ASSERT_KNOWN(py_size >= k1 * _m, "Invalid py size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            int ret = (*_reverseTwo)(tx, ty, px, py, _atomicFuncArg);

            CPPADCG_ASSERT_KNOWN(ret != 1, "Second-order reverse mode failed: py[2*i] (i=0...m) must be zero.");
            CPPADCG_ASSERT_KNOWN(ret == 0, "Second-order reverse mode failed.");
        }

        /// calculate sparse Jacobians 

        virtual void SparseJacobian(const Base* x, size_t x_size,
                                    Base* jac, size_t jac_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* row;
            unsigned long const* col;
            unsigned long nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            CppAD::vector<Base> compressed(nnz);

            if (nnz > 0) {
                _in[0] = x;
                _out[0] = &compressed[0];

                (*_sparseJacobian)(&_in[0], &_out[0], _atomicFuncArg);
            }

            createDenseFromSparse(compressed,
                                  _m, _n,
                                  row, col,
                                  nnz,
                                  jac, jac_size);
        }

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac,
                                    std::vector<size_t>& row,
                                    std::vector<size_t>& col) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* drow;
            unsigned long const* dcol;
            unsigned long nnz;
            (*_jacobianSparsity)(&drow, &dcol, &nnz);

            jac.resize(nnz);
            row.resize(nnz);
            col.resize(nnz);

            if (nnz > 0) {
                _in[0] = &x[0];
                _out[0] = &jac[0];

                (*_sparseJacobian)(&_in[0], &_out[0], _atomicFuncArg);
                std::copy(drow, drow + nnz, row.begin());
                std::copy(dcol, dcol + nnz, col.begin());
            }
        }

        virtual void SparseJacobian(const Base* x, size_t x_size,
                                    Base* jac,
                                    size_t const** row,
                                    size_t const** col,
                                    size_t nnz) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* drow;
            unsigned long const* dcol;
            unsigned long K;
            (*_jacobianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Jacobian");
            *row = drow;
            *col = dcol;

            if (nnz > 0) {
                _in[0] = x;
                _out[0] = jac;

                (*_sparseJacobian)(&_in[0], &_out[0], _atomicFuncArg);
            }
        }

        virtual void SparseJacobian(const std::vector<const Base*>& x,
                                    Base* jac,
                                    size_t const** row,
                                    size_t const** col,
                                    size_t nnz) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == x.size(), "The number of independent variable arrays is invalid");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* drow;
            unsigned long const* dcol;
            unsigned long K;
            (*_jacobianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Jacobian");
            *row = drow;
            *col = dcol;

            if (nnz > 0) {
                _out[0] = jac;

                (*_sparseJacobian)(&x[0], &_out[0], _atomicFuncArg);
            }
        }

        /// calculate sparse Hessians 

        virtual void SparseHessian(const Base* x, size_t x_size,
                                   const Base* w, size_t w_size,
                                   Base* hess, size_t hess_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(w_size == _m, "Invalid multiplier array size");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* row, *col;
            unsigned long nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            CppAD::vector<Base> compressed(nnz);
            if (nnz > 0) {
                _inHess[0] = x;
                _inHess[1] = w;
                _out[0] = &compressed[0];

                (*_sparseHessian)(&_inHess[0], &_out[0], _atomicFuncArg);
            }

            createDenseFromSparse(compressed,
                                  _n, _n,
                                  row, col,
                                  nnz,
                                  hess, hess_size);
        }

        virtual void SparseHessian(const std::vector<Base> &x,
                                   const std::vector<Base> &w,
                                   std::vector<Base>& hess,
                                   std::vector<size_t>& row,
                                   std::vector<size_t>& col) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* drow, *dcol;
            unsigned long nnz;
            (*_hessianSparsity)(&drow, &dcol, &nnz);

            hess.resize(nnz);
            row.resize(nnz);
            col.resize(nnz);

            if (nnz > 0) {
                std::copy(drow, drow + nnz, row.begin());
                std::copy(dcol, dcol + nnz, col.begin());

                _inHess[0] = &x[0];
                _inHess[1] = &w[0];
                _out[0] = &hess[0];

                (*_sparseHessian)(&_inHess[0], &_out[0], _atomicFuncArg);
            }
        }

        virtual void SparseHessian(const Base* x, size_t x_size,
                                   const Base* w, size_t w_size,
                                   Base* hess,
                                   size_t const** row,
                                   size_t const** col,
                                   size_t nnz) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(w_size == _m, "Invalid multiplier array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* drow, *dcol;
            unsigned long K;
            (*_hessianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Hessian");
            *row = drow;
            *col = dcol;

            if (nnz > 0) {
                _inHess[0] = x;
                _inHess[1] = w;
                _out[0] = hess;

                (*_sparseHessian)(&_inHess[0], &_out[0], _atomicFuncArg);
            }
        }

        virtual void SparseHessian(const std::vector<const Base*>& x,
                                   const Base* w, size_t w_size,
                                   Base* hess,
                                   size_t const** row,
                                   size_t const** col,
                                   size_t nnz) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == x.size(), "The number of independent variable arrays is invalid");
            CPPADCG_ASSERT_KNOWN(w_size == _m, "Invalid multiplier array size");
            CPPADCG_ASSERT_KNOWN(_missingAtomicFunctions == 0, "Some atomic functions used by the compiled model have not been specified yet");

            unsigned long const* drow, *dcol;
            unsigned long K;
            (*_hessianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Hessian");
            *row = drow;
            *col = dcol;

            if (nnz > 0) {
                std::copy(x.begin(), x.end(), _inHess.begin());
                _inHess.back() = w;
                _out[0] = hess;

                (*_sparseHessian)(&_inHess[0], &_out[0], _atomicFuncArg);
            }
        }

        virtual ~LinuxDynamicLibModel() {
            if (_dynLib != NULL) {
                _dynLib->destroyed(this);
            }
        }

    protected:

        /**
         * Creates a new model 
         * 
         * @param name The model name
         */
        LinuxDynamicLibModel(LinuxDynamicLib<Base>* dynLib, const std::string& name) :
            _name(name),
            _dynLib(dynLib),
            _m(0),
            _n(0),
            _missingAtomicFunctions(0),
            _zero(NULL),
            _forwardOne(NULL),
            _reverseOne(NULL),
            _reverseTwo(NULL),
            _jacobian(NULL),
            _hessian(NULL),
            _sparseForwardOne(NULL),
            _sparseReverseOne(NULL),
            _sparseReverseTwo(NULL),
            _sparseJacobian(NULL),
            _sparseHessian(NULL),
            _forwardOneSparsity(NULL),
            _reverseOneSparsity(NULL),
            _reverseTwoSparsity(NULL),
            _jacobianSparsity(NULL),
            _hessianSparsity(NULL),
            _hessianSparsity2(NULL) {

            assert(_dynLib != NULL);

            // validate the dynamic library
            validate();

            // load functions from the dynamic library
            loadFunctions();
        }

        virtual void validate() {
            std::string error;

            /**
             * Check the data type
             */
            void (*infoFunc)(const char** baseName, unsigned long*, unsigned long*, unsigned int*, unsigned int*);
            *(void **) (&infoFunc) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_INFO, error);
            CPPADCG_ASSERT_KNOWN(error.empty(), error.c_str());

            // local
            const char* localBaseName = typeid (Base).name();
            std::string local = CLangCompileModelHelper<Base>::baseTypeName() + "  " + localBaseName;

            // from dynamic library
            const char* dynamicLibBaseName = NULL;
            unsigned int inSize = 0;
            unsigned int outSize = 0;
            (*infoFunc)(&dynamicLibBaseName, &_m, &_n, &inSize, &outSize);

            _in.resize(inSize);
            _inHess.resize(inSize + 1);
            _out.resize(outSize);

            CPPADCG_ASSERT_KNOWN(local == std::string(dynamicLibBaseName),
                                 (std::string("Invalid data type in dynamic library. Expected '") + local
                                 + "' but the library provided '" + dynamicLibBaseName + "'.").c_str());
            CPPADCG_ASSERT_KNOWN(inSize > 0,
                                 "Invalid dimension received from the dynamic library.");
            CPPADCG_ASSERT_KNOWN(outSize > 0,
                                 "Invalid dimension received from the dynamic library.");
        }

        virtual void loadFunctions() {
            *(void **) (&_zero) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_FORWAD_ZERO, false);
            *(void **) (&_forwardOne) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_FORWARD_ONE, false);
            *(void **) (&_reverseOne) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_REVERSE_ONE, false);
            *(void **) (&_reverseTwo) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_REVERSE_TWO, false);
            *(void **) (&_jacobian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN, false);
            *(void **) (&_hessian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_HESSIAN, false);
            *(void **) (&_sparseForwardOne) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_FORWARD_ONE, false);
            *(void **) (&_sparseReverseOne) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_REVERSE_ONE, false);
            *(void **) (&_sparseReverseTwo) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_REVERSE_TWO, false);
            *(void **) (&_sparseJacobian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_JACOBIAN, false);
            *(void **) (&_sparseHessian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_HESSIAN, false);
            *(void **) (&_forwardOneSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_FORWARD_ONE_SPARSITY, false);
            *(void **) (&_reverseOneSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_REVERSE_ONE_SPARSITY, false);
            *(void **) (&_reverseTwoSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_REVERSE_TWO_SPARSITY, false);
            *(void **) (&_jacobianSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN_SPARSITY, false);
            *(void **) (&_hessianSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_HESSIAN_SPARSITY, false);
            *(void **) (&_hessianSparsity2) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_HESSIAN_SPARSITY2, false);
            *(void **) (&_atomicFunctions) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_ATOMIC_FUNC_NAMES, true);

            CPPADCG_ASSERT_KNOWN((_sparseForwardOne == NULL) == (_forwardOneSparsity == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseForwardOne == NULL) == (_forwardOne == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseReverseOne == NULL) == (_reverseOneSparsity == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseReverseOne == NULL) == (_reverseOne == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseReverseTwo == NULL) == (_reverseTwoSparsity == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseReverseTwo == NULL) == (_reverseTwo == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseJacobian == NULL) || (_jacobianSparsity != NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseHessian == NULL) || (_hessianSparsity != NULL), "Missing functions in the dynamic library");

            /**
             * Prepare the atomic functions argument
             */
            const char** names;
            unsigned long n;
            (*_atomicFunctions)(&names, &n);
            _atomic.resize(n);

            _atomicFuncArg.libModel = this;
            _atomicFuncArg.forward = &atomicForward;
            _atomicFuncArg.reverse = &atomicReverse;

            _missingAtomicFunctions = n;
        }

        template <class VectorSet>
        inline void loadSparsity(bool set_type,
                                 VectorSet& s,
                                 unsigned long nrows, unsigned long ncols,
                                 unsigned long const* rows, unsigned long const* cols,
                                 unsigned long nnz) {
            s.resize(nrows * ncols, false);

            for (unsigned long i = 0; i < nnz; i++) {
                s[rows[i] * ncols + cols[i]] = true;
            }
        }

        template <class VectorSet>
        inline void loadSparsity(const std::set<size_t>& set_type,
                                 VectorSet& s,
                                 unsigned long nrows, unsigned long ncols,
                                 unsigned long const* rows, unsigned long const* cols,
                                 unsigned long nnz) {

            // dimension size of result vector
            s.resize(nrows);

            for (unsigned long i = 0; i < nnz; i++) {
                s[rows[i]].insert(cols[i]);
            }
        }

        inline void createDenseFromSparse(const CppAD::vector<Base>& compressed,
                                          unsigned long nrows, unsigned long ncols,
                                          unsigned long const* rows, unsigned long const* cols,
                                          unsigned long nnz,
                                          Base* mat, size_t mat_size) const {
            CPPADCG_ASSERT_KNOWN(mat_size == nrows * ncols, "Invalid matrix size");
            std::fill(mat, mat + mat_size, 0);

            for (size_t i = 0; i < nnz; i++) {
                mat[rows[i] * ncols + cols[i]] = compressed[i];
            }
        }

        void dynamicLibClosed() {
            _dynLib = NULL;
            _zero = NULL;
            _forwardOne = NULL;
            _reverseOne = NULL;
            _reverseTwo = NULL;
            _jacobian = NULL;
            _hessian = NULL;
            _sparseForwardOne = NULL;
            _sparseReverseOne = NULL;
            _sparseReverseTwo = NULL;
            _sparseJacobian = NULL;
            _sparseHessian = NULL;
            _forwardOneSparsity = NULL;
            _reverseOneSparsity = NULL;
            _reverseTwoSparsity = NULL;
            _jacobianSparsity = NULL;
            _hessianSparsity = NULL;
            _hessianSparsity2 = NULL;
        }

    private:

        static int atomicForward(void* libModelIn,
                                 int atomicIndex,
                                 int q,
                                 int p,
                                 const void* tx, unsigned long txSize,
                                 void* ty, unsigned long tySize) {
            LinuxDynamicLibModel<Base>* libModel = static_cast<LinuxDynamicLibModel<Base>*> (libModelIn);
            atomic_base<Base>* atomicBase = libModel->_atomic[atomicIndex];
            const Base* txb = static_cast<const Base*> (tx);
            Base* tyb = static_cast<Base*> (ty);

            vector<bool> vx, vy;
            libModel->_tx.resize(txSize);
            libModel->_ty.resize(tySize);
            std::copy(txb, txb + txSize, &libModel->_tx[0]);
            if (p > 0) {
                std::copy(tyb, tyb + tySize, &libModel->_ty[0]);
            }

            bool ret = atomicBase->forward(q, p, vx, vy, libModel->_tx, libModel->_ty);
            std::copy(&libModel->_ty[0], &libModel->_ty[0] + tySize, tyb);

            return ret;
        }

        static int atomicReverse(void* libModelIn,
                                 int atomicIndex,
                                 int p,
                                 const void* tx,
                                 const void* ty,
                                 void* px,
                                 const void* py,
                                 unsigned long xSize,
                                 unsigned long ySize) {
            LinuxDynamicLibModel<Base>* libModel = static_cast<LinuxDynamicLibModel<Base>*> (libModelIn);
            atomic_base<Base>* atomicBase = libModel->_atomic[atomicIndex];
            const Base* txb = static_cast<const Base*> (tx);
            const Base* tyb = static_cast<const Base*> (ty);
            Base* pxb = static_cast<Base*> (px);
            const Base* pyb = static_cast<const Base*> (py);

            libModel->_tx.resize(xSize);
            libModel->_ty.resize(ySize);
            libModel->_px.resize(xSize);
            libModel->_py.resize(ySize);

            std::copy(txb, txb + xSize, &libModel->_tx[0]);
            std::copy(tyb, tyb + ySize, &libModel->_ty[0]);
            if (p > 0) {
                std::copy(pxb, pxb + xSize, &libModel->_px[0]);
            }
            std::copy(pyb, pyb + ySize, &libModel->_py[0]);

#ifndef NDEBUG
            if (libModel->_evalAtomicForwardOne4CppAD) {
                // only required in order to avoid an issue with a validation inside CppAD 
                vector<bool> vx, vy;
                if (!atomicBase->forward(p, p, vx, vy, libModel->_tx, libModel->_ty))
                    return false;
            }
#endif

            bool ret = atomicBase->reverse(p, libModel->_tx, libModel->_ty, libModel->_px, libModel->_py);
            std::copy(&libModel->_px[0], &libModel->_px[0] + xSize, pxb);

            return ret;
        }

        LinuxDynamicLibModel(const LinuxDynamicLibModel&); // not implemented

        LinuxDynamicLibModel& operator=(const LinuxDynamicLibModel&); // not implemented

        friend class LinuxDynamicLib<Base>;
    };

}

#endif

#endif