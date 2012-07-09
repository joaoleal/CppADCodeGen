#ifndef CPPAD_CG_DYNAMICLIB_LINUX_INCLUDED
#define	CPPAD_CG_DYNAMICLIB_LINUX_INCLUDED
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
     * Useful class to call the compiled source code in a dynamic library.
     * For the Linux Operating System only.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class LinuxDynamicLib : public DynamicLib<Base> {
    protected:
        const std::string _dynLibName;
        /// the dynamic library handler
        void* _dynLibHandle;
        size_t _m;
        size_t _n;
        void (*_zero)(const double*, double*);
        // jacobian function in the dynamic library
        void (*_jacobian)(const double*, double*);
        // hessian function in the dynamic library
        void (*_hessian)(const double*, const double*, double*);
        // sparse jacobian function in the dynamic library
        void (*_sparseJacobian)(const double*, double*);
        // sparse hessian function in the dynamic library
        void (*_sparseHessian)(const double*, const double*, double*);
        // jacobian sparsity function in the dynamic library
        void (*_jacobianSparsity)(unsigned long int const** row,
                unsigned long int const** col,
                unsigned long int * nnz);
        // hessian sparsity function in the dynamic library
        void (*_hessianSparsity)(unsigned long int const** row,
                unsigned long int const** col,
                unsigned long int * nnz);

    public:

        LinuxDynamicLib(const std::string& dynLibName) :
            _dynLibName(dynLibName),
            _dynLibHandle(NULL),
            _m(0),
            _n(0),
            _zero(NULL),
            _jacobian(NULL),
            _hessian(NULL),
            _sparseJacobian(NULL),
            _sparseHessian(NULL),
            _jacobianSparsity(NULL),
            _hessianSparsity(NULL) {

            std::string path;
            if (dynLibName[0] == '/') {
                path = dynLibName; // absolute path
            } else {
                path = "./" + dynLibName; // relative path
            }

            // load the dynamic library
            _dynLibHandle = dlopen(path.c_str(), RTLD_NOW);
            CPPADCG_ASSERT_KNOWN(_dynLibHandle != NULL, ("Failed to dynamically load library '" + dynLibName + "'").c_str());

            // validate the dynamic library
            validate();

            // load functions from the dynamic library
            loadFunctions();
        }

        // Jacobian sparsity

        virtual std::vector<bool> JacobianSparsityBool() {
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            bool set_type;
            std::vector<bool> s;

            loadSparsity(set_type, s, _m, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > JacobianSparsitySet() {
            CPPADCG_ASSERT_KNOWN(_jacobianSparsity != NULL, "No Jacobian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _m, _n, row, col, nnz);

            return s;
        }

        // Hessian sparsity 

        virtual std::vector<bool> HessianSparsityBool() {
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            bool set_type;
            std::vector<bool> s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
        }

        virtual std::vector<std::set<size_t> > HessianSparsitySet() {
            CPPADCG_ASSERT_KNOWN(_hessianSparsity != NULL, "No Hessian sparsity function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            std::set<size_t> set_type;
            std::vector<std::set<size_t> > s;

            loadSparsity(set_type, s, _n, _n, row, col, nnz);

            return s;
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
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");

            std::vector<double> dep(_m);
            (*_zero)(&x[0], &dep[0]);
            return dep;
        }

        virtual void ForwardZero(const std::vector<Base> &x, std::vector<Base>& dep) {
            CPPADCG_ASSERT_KNOWN(_zero != NULL, "No zero order forward function defined in the dynamic library");

            dep.resize(_m);
            (*_zero)(&x[0], &dep[0]);
        }

        /// calculate entire Jacobian       

        virtual std::vector<Base> Jacobian(const std::vector<Base> &x) {
            CPPADCG_ASSERT_KNOWN(_jacobian != NULL, "No Jacobian function defined in the dynamic library");

            std::vector<double> jac(_m * _n);
            (*_jacobian)(&x[0], &jac[0]);
            return jac;
        }

        virtual void Jacobian(const std::vector<Base> &x, std::vector<Base>& jac) {
            CPPADCG_ASSERT_KNOWN(_jacobian != NULL, "No Jacobian function defined in the dynamic library");

            jac.resize(_m * _n);
            (*_jacobian)(&x[0], &jac[0]);
        }

        /// calculate Hessian for one component of f

        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          const std::vector<Base> &w) {
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");

            std::vector<Base> hess(_n * _n);
            (*_hessian)(&x[0], &w[0], &hess[0]);
            return hess;
        }

        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          size_t i) {
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");

            std::vector<Base> hess(_n * _n);
            std::vector<Base> w(_m);
            w[i] = 1.0;
            (*_hessian)(&x[0], &w[0], &hess[0]);
            return hess;
        }

        virtual void Hessian(const std::vector<Base> &x,
                             const std::vector<Base> &w,
                             std::vector<Base>& hess) {
            CPPADCG_ASSERT_KNOWN(_hessian != NULL, "No Hessian function defined in the dynamic library");

            hess.resize(_n * _n);
            (*_hessian)(&x[0], &w[0], &hess[0]);
        }

        /// calculate sparse Jacobians 

        virtual std::vector<Base> SparseJacobian(const std::vector<Base> &x) {
            std::vector<Base> jac(_m * _n);
            SparseJacobian(x, jac);
            return jac;
        }

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac) {
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse jacobian function defined in the dynamic library");

            unsigned long int const* row;
            unsigned long int const* col;
            unsigned long int nnz;
            (*_jacobianSparsity)(&row, &col, &nnz);

            std::vector<Base> compressed(nnz);

            (*_sparseJacobian)(&x[0], &compressed[0]);

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
            CPPADCG_ASSERT_KNOWN(_sparseJacobian != NULL, "No sparse Jacobian function defined in the dynamic library");

            unsigned long int const* drow;
            unsigned long int const* dcol;
            unsigned long int nnz;
            (*_jacobianSparsity)(&drow, &dcol, &nnz);

            jac.resize(nnz);
            row.resize(nnz);
            col.resize(nnz);

            (*_sparseJacobian)(&x[0], &jac[0]);
            std::copy(drow, drow + nnz, row.begin());
            std::copy(dcol, dcol + nnz, col.begin());
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
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");

            unsigned long int const* row, *col;
            unsigned long int nnz;
            (*_hessianSparsity)(&row, &col, &nnz);

            std::vector<Base> compressed(nnz);

            (*_sparseHessian)(&x[0], &w[0], &compressed[0]);

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
            CPPADCG_ASSERT_KNOWN(_sparseHessian != NULL, "No sparse Hessian function defined in the dynamic library");

            unsigned long int const* drow, *dcol;
            unsigned long int nnz;
            (*_hessianSparsity)(&drow, &dcol, &nnz);

            hess.resize(nnz);
            row.resize(nnz);
            col.resize(nnz);

            (*_sparseHessian)(&x[0], &w[0], &hess[0]);
            std::copy(drow, drow + nnz, row.begin());
            std::copy(dcol, dcol + nnz, col.begin());
        }

        virtual ~LinuxDynamicLib() {
            if (_dynLibHandle != NULL) {
                dlclose(_dynLibHandle);
                _dynLibHandle = NULL;
            }
        }

    protected:
        
        virtual void validate() {
            std::string error;
            
            /**
             * Check the version
             */
            unsigned long int (*versionFunc)();
            *(void **) (&versionFunc) = loadFunction(CLangCompileHelper<Base>::FUNCTION_VERSION, error);
            CPPADCG_ASSERT_KNOWN(error.empty(), error.c_str());
            
            unsigned long int version = (*versionFunc)();
            CPPADCG_ASSERT_KNOWN(CLangCompileHelper<Base>::API_VERSION == version,
                    "The version of the dynamic library is incompatible with the current version");
            
            /**
             * Check the data type
             */
            void (*infoFunc)(const char** baseName, unsigned long int*, unsigned long int*);
            *(void **) (&infoFunc) = loadFunction(CLangCompileHelper<Base>::FUNCTION_INFO, error);
            CPPADCG_ASSERT_KNOWN(error.empty(), error.c_str());

            // local
            const char* localBaseName = typeid (Base).name();
            std::string local = CLangCompileHelper<Base>::baseTypeName() + "  " + localBaseName;
            
            // from dynamic library
            const char* dynamicLibBaseName = NULL;
            (*infoFunc)(&dynamicLibBaseName, &_m, &_n);

            CPPADCG_ASSERT_KNOWN(local == std::string(dynamicLibBaseName),
                    (std::string("Invalid data type in dynamic library. Expected '") + local
                    + "' but the library provided '" + dynamicLibBaseName + "'.").c_str());
        }

        virtual void loadFunctions() {
            std::string error;
            
            *(void **) (&_zero) = loadFunction(CLangCompileHelper<Base>::FUNCTION_FORWAD_ZERO, error);
            *(void **) (&_jacobian) = loadFunction(CLangCompileHelper<Base>::FUNCTION_JACOBIAN, error);
            *(void **) (&_hessian) = loadFunction(CLangCompileHelper<Base>::FUNCTION_HESSIAN, error);
            *(void **) (&_sparseJacobian) = loadFunction(CLangCompileHelper<Base>::FUNCTION_SPARSE_JACOBIAN, error);
            *(void **) (&_sparseHessian) = loadFunction(CLangCompileHelper<Base>::FUNCTION_SPARSE_HESSIAN, error);
            *(void **) (&_jacobianSparsity) = loadFunction(CLangCompileHelper<Base>::FUNCTION_JACOBIAN_SPARSITY, error);
            *(void **) (&_hessianSparsity) = loadFunction(CLangCompileHelper<Base>::FUNCTION_HESSIAN_SPARSITY, error);

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

        inline void* loadFunction(const std::string& functionName, std::string& error) {
            void* functor = dlsym(_dynLibHandle, functionName.c_str());
            char *err;
            error.clear();
            if ((err = dlerror()) != NULL) {
                error += "Failed to load function '";
                error += functionName + ": " + err;
                return NULL;
            }
            return functor;
        }

    private:
        LinuxDynamicLib(const LinuxDynamicLib&); // not implemented

        LinuxDynamicLib& operator=(const LinuxDynamicLib&); // not implemented
    };

}

#endif

#endif