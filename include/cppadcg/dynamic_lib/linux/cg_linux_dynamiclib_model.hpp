#ifndef CPPAD_CG_LINUX_DYNAMICLIB_MODEL_INCLUDED
#define	CPPAD_CG_LINUX_DYNAMICLIB_MODEL_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#ifdef __linux__

#include <typeinfo>
#include <dlfcn.h>

namespace CppAD {

    /**
     * Useful class to call the compiled model in a dynamic library.
     * For the Linux Operating System only.
     * 
     * \author Joao Leal
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
        void (*_zero)(double const*const*, double*const*);
        // jacobian function in the dynamic library
        void (*_jacobian)(double const*const*, double*const*);
        // hessian function in the dynamic library
        void (*_hessian)(double const*const*, double*const*);
        // sparse jacobian function in the dynamic library
        void (*_sparseJacobian)(double const*const*, double*const*);
        // sparse hessian function in the dynamic library
        void (*_sparseHessian)(double const*const*, double*const*);
        // jacobian sparsity function in the dynamic library
        void (*_jacobianSparsity)(unsigned long int const** row,
                unsigned long int const** col,
                unsigned long int * nnz);
        // hessian sparsity function in the dynamic library
        void (*_hessianSparsity)(unsigned long int const** row,
                unsigned long int const** col,
                unsigned long int * nnz);

    public:

        const std::string& getName() const {
            return _name;
        }

        // Jacobian sparsity

        virtual std::vector<bool> JacobianSparsityBool() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            bool set_type = true;
            std::vector<bool> s;

            loadSparsity(set_type, s, _m, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > JacobianSparsitySet() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _m, _n, row, col, nnz);

            return s;
        }

        virtual void JacobianSparsity(std::vector<size_t>& rows, std::vector<size_t>& cols) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
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

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            bool set_type = true;
            std::vector<bool> s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > HessianSparsitySet() {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual void HessianSparsity(std::vector<size_t>& rows, std::vector<size_t>& cols) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

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

        virtual std::vector<Base> ForwardZero(const std::vector<Base> &x) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x.size() == _n, "Invalid independent vector size");

            _in[0] = &x[0];
            std::vector<double> dep(_m);
            _out[0] = &dep[0];

            (*_zero)(&_in[0], &_out[0]);

            return dep;
        }

        virtual void ForwardZero(const std::vector<Base> &x, std::vector<Base>& dep) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x.size() == _n, "Invalid independent vector size");

            _in[0] = &x[0];
            dep.resize(_m);
            _out[0] = &dep[0];

            (*_zero)(&_in[0], &_out[0]);
        }

        virtual void ForwardZero(const Base* x, size_t x_size, Base* dep, size_t dep_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(dep_size == _m, "Invalid dependent array size");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");


            _in[0] = x;
            _out[0] = dep;

            (*_zero)(&_in[0], &_out[0]);
        }

        virtual void ForwardZero(const std::vector<const Base*> &x,
                                 Base* dep, size_t dep_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == x.size(), "The number of independent variable arrays is invalid");
            CPPADCG_ASSERT_KNOWN(dep_size == _m, "Invalid dependent array size");

            _out[0] = dep;

            (*_zero)(&x[0], &_out[0]);
        }

        /// calculate entire Jacobian       

        virtual std::vector<Base> Jacobian(const std::vector<Base> &x) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobian != NULL, "No Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x.size() == _n, "Invalid independent vector size");

            _in[0] = &x[0];
            std::vector<double> jac(_m * _n);
            _out[0] = &jac[0];

            (*_jacobian)(&_in[0], &_out[0]);

            return jac;
        }

        virtual void Jacobian(const std::vector<Base> &x, std::vector<Base>& jac) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobian != NULL, "No Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x.size() == _n, "Invalid independent vector size");

            _in[0] = &x[0];
            jac.resize(_m * _n);
            _out[0] = &jac[0];

            (*_jacobian)(&_in[0], &_out[0]);
        }

        virtual void Jacobian(const Base* x, size_t x_size,
                              Base* jac, size_t jac_size) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_jacobian != NULL, "No Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(jac_size == _m * _n, "Invalid Jacobian array size");


            _in[0] = x;
            _out[0] = jac;

            (*_jacobian)(&_in[0], &_out[0]);
        }

        /// calculate Hessian for one component of f

        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          const std::vector<Base> &w) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");

            _inHess[0] = &x[0];
            _inHess[1] = &w[0];
            std::vector<Base> hess(_n * _n);
            _out[0] = &hess[0];

            (*_hessian)(&_inHess[0], &_out[0]);

            return hess;
        }

        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          size_t i) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");

            std::vector<Base> w(_m);
            w[i] = 1.0;
            _inHess[0] = &x[0];
            _inHess[1] = &w[0];
            std::vector<Base> hess(_n * _n);
            _out[0] = &hess[0];

            (*_hessian)(&_inHess[0], &_out[0]);

            return hess;
        }

        virtual void Hessian(const std::vector<Base> &x,
                             const std::vector<Base> &w,
                             std::vector<Base>& hess) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");

            _inHess[0] = &x[0];
            _inHess[1] = &w[0];
            hess.resize(_n * _n);
            _out[0] = &hess[0];

            (*_hessian)(&_inHess[0], &_out[0]);
        }

        virtual void Hessian(const Base* x, size_t x_size,
                             const Base* w, size_t w_size,
                             Base* hess) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");
            CPPADCG_ASSERT_KNOWN(x_size == _n, "Invalid independent array size");
            CPPADCG_ASSERT_KNOWN(w_size == _m, "Invalid multiplier array size");

            _inHess[0] = x;
            _inHess[1] = w;
            _out[0] = hess;

            (*_hessian)(&_inHess[0], &_out[0]);
        }

        /// calculate sparse Jacobians 

        virtual std::vector<Base> SparseJacobian(const std::vector<Base> &x) {
            std::vector<Base> jac(_m * _n);
            SparseJacobian(x, jac);
            return jac;
        }

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");

            unsigned long int const* row;
            unsigned long int const* col;
            unsigned long int nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            std::vector<Base> compressed(nnz);

            _in[0] = &x[0];
            _out[0] = &compressed[0];

            (*_sparseJacobian)(&_in[0], &_out[0]);

            createDenseFromSparse(compressed,
                                  _m, _n,
                                  row, col,
                                  nnz,
                                  jac);
        }

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac,
                                    std::vector<size_t>& row,
                                    std::vector<size_t>& col) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");

            unsigned long int const* drow;
            unsigned long int const* dcol;
            unsigned long int nnz;
            (*_jacobianSparsity)(&drow, &dcol, &nnz);

            jac.resize(nnz);
            row.resize(nnz);
            col.resize(nnz);

            _in[0] = &x[0];
            _out[0] = &jac[0];

            (*_sparseJacobian)(&_in[0], &_out[0]);
            std::copy(drow, drow + nnz, row.begin());
            std::copy(dcol, dcol + nnz, col.begin());
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

            unsigned long int const* drow;
            unsigned long int const* dcol;
            unsigned long int K;
            (*_jacobianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Jacobian");
            *row = drow;
            *col = dcol;

            _in[0] = x;
            _out[0] = jac;

            (*_sparseJacobian)(&_in[0], &_out[0]);
        }

        virtual void SparseJacobian(const std::vector<const Base*>& x,
                                    Base* jac,
                                    size_t const** row,
                                    size_t const** col,
                                    size_t nnz) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse Jacobian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == x.size(), "The number of independent variable arrays is invalid");

            unsigned long int const* drow;
            unsigned long int const* dcol;
            unsigned long int K;
            (*_jacobianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Jacobian");
            *row = drow;
            *col = dcol;

            _out[0] = jac;

            (*_sparseJacobian)(&x[0], &_out[0]);
        }

        /// calculate sparse Hessians 

        virtual std::vector<Base> SparseHessian(const std::vector<Base> &x, const std::vector<Base> &w) {
            std::vector<Base> hess;
            SparseHessian(x, w, hess);
            return hess;
        }

        virtual void SparseHessian(const std::vector<Base> &x,
                                   const std::vector<Base> &w,
                                   std::vector<Base>& hess) {
            CPPADCG_ASSERT_KNOWN(_dynLib != NULL, "Dynamic library closed");
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");
            CPPADCG_ASSERT_KNOWN(_in.size() == 1, "The number of independent variable arrays is higher than 1,"
                                 " please use the variable size methods");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            std::vector<Base> compressed(nnz);

            _inHess[0] = &x[0];
            _inHess[1] = &w[0];
            _out[0] = &compressed[0];

            (*_sparseHessian)(&_inHess[0], &_out[0]);

            createDenseFromSparse(compressed,
                                  _n, _n,
                                  row, col,
                                  nnz,
                                  hess);
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

            unsigned long int const* drow, *dcol;
            unsigned long int nnz;
            (*_hessianSparsity)(&drow, &dcol, &nnz);

            hess.resize(nnz);
            row.resize(nnz);
            col.resize(nnz);

            _inHess[0] = &x[0];
            _inHess[1] = &w[0];
            _out[0] = &hess[0];

            (*_sparseHessian)(&_inHess[0], &_out[0]);
            std::copy(drow, drow + nnz, row.begin());
            std::copy(dcol, dcol + nnz, col.begin());
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

            unsigned long int const* drow, *dcol;
            unsigned long int K;
            (*_hessianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Hessian");
            *row = drow;
            *col = dcol;

            _inHess[0] = x;
            _inHess[1] = w;
            _out[0] = hess;

            (*_sparseHessian)(&_inHess[0], &_out[0]);
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

            unsigned long int const* drow, *dcol;
            unsigned long int K;
            (*_hessianSparsity)(&drow, &dcol, &K);
            CPPADCG_ASSERT_KNOWN(K == nnz, "Invalid number of non-zero elements in Hessian");
            *row = drow;
            *col = dcol;

            std::copy(x.begin(), x.end(), _inHess.begin());
            _inHess.back() = w;
            _out[0] = hess;

            (*_sparseHessian)(&_inHess[0], &_out[0]);
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
         * \param name The model name
         */
        LinuxDynamicLibModel(LinuxDynamicLib<Base>* dynLib, const std::string& name) :
            _name(name),
            _dynLib(dynLib),
            _m(0),
            _n(0),
            _zero(NULL),
            _jacobian(NULL),
            _hessian(NULL),
            _sparseJacobian(NULL),
            _sparseHessian(NULL),
            _jacobianSparsity(NULL),
            _hessianSparsity(NULL) {

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
            void (*infoFunc)(const char** baseName, unsigned long int*, unsigned long int*, unsigned int*, unsigned int*);
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
            *(void **) (&_jacobian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN, false);
            *(void **) (&_hessian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_HESSIAN, false);
            *(void **) (&_sparseJacobian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_JACOBIAN, false);
            *(void **) (&_sparseHessian) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_SPARSE_HESSIAN, false);
            *(void **) (&_jacobianSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_JACOBIAN_SPARSITY, false);
            *(void **) (&_hessianSparsity) = _dynLib->loadFunction(_name + "_" + CLangCompileModelHelper<Base>::FUNCTION_HESSIAN_SPARSITY, false);

            CPPADCG_ASSERT_KNOWN((_sparseJacobian == NULL) == (_jacobianSparsity == NULL), "Missing functions in the dynamic library");
            CPPADCG_ASSERT_KNOWN((_sparseHessian == NULL) == (_hessianSparsity == NULL), "Missing functions in the dynamic library");
        }

        template <class VectorSet>
        inline void loadSparsity(bool set_type,
                                 VectorSet& s,
                                 unsigned long int nrows, unsigned long int ncols,
                                 unsigned long int const* rows, unsigned long int const* cols,
                                 unsigned long int nnz) {
            s.resize(nrows * ncols, false);

            for (unsigned long int i = 0; i < nnz; i++) {
                s[rows[i] * ncols + cols[i]] = true;
            }
        }

        template <class VectorSet>
        inline void loadSparsity(const std::set<size_t>& set_type,
                                 VectorSet& s,
                                 unsigned long int nrows, unsigned long int ncols,
                                 unsigned long int const* rows, unsigned long int const* cols,
                                 unsigned long int nnz) {

            // dimension size of result vector
            s.resize(nrows);

            for (unsigned long int i = 0; i < nnz; i++) {
                s[rows[i]].insert(cols[i]);
            }
        }

        inline void createDenseFromSparse(const std::vector<Base>& compressed,
                                          unsigned long int nrows, unsigned long int ncols,
                                          unsigned long int const* rows, unsigned long int const* cols,
                                          unsigned long int nnz,
                                          std::vector<Base> mat) const {

            std::fill(mat.begin(), mat.end(), 0);
            mat.resize(nrows * ncols);

            for (size_t i = 0; i < nnz; i++) {
                mat[rows[i] * ncols + cols[i]] = compressed[i];
            }
        }

        void dynamicLibClosed() {
            _dynLib = NULL;
            _zero = NULL;
            _jacobian = NULL;
            _hessian = NULL;
            _sparseJacobian = NULL;
            _sparseHessian = NULL;
            _jacobianSparsity = NULL;
            _hessianSparsity = NULL;
        }

    private:
        LinuxDynamicLibModel(const LinuxDynamicLibModel&); // not implemented

        LinuxDynamicLibModel& operator=(const LinuxDynamicLibModel&); // not implemented

        friend class LinuxDynamicLib<Base>;
    };

}

#endif

#endif