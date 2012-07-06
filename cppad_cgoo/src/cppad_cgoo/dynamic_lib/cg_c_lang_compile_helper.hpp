#ifndef CPPAD_CG_C_LANG_COMPILE_HELPER_INCLUDED
#define	CPPAD_CG_C_LANG_COMPILE_HELPER_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */
#include <assert.h>
#include <vector>

namespace CppAD {

    /**
     * Useful class for generating source code for the creation of a dynamic
     * library.
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CLangCompileHelper {
    public:
        static const std::string FUNCTION_FORWAD_ZERO;
        static const std::string FUNCTION_JACOBIAN;
        static const std::string FUNCTION_HESSIAN;
        static const std::string FUNCTION_SPARSE_JACOBIAN;
        static const std::string FUNCTION_SPARSE_HESSIAN;
        static const std::string FUNCTION_JACOBIAN_SPARSITY;
        static const std::string FUNCTION_HESSIAN_SPARSITY;
        static const std::string FUNCTION_INFO;
        static const std::string FUNCTION_VERSION;
        static const unsigned long int API_VERSION;
    protected:
        ADFun<CppAD::CG<Base> >* _fun;
        bool _zero;
        bool _jacobian;
        bool _hessian;
        bool _sparseJacobian;
        bool _sparseHessian;
        std::string _libraryName; // the path of the dynamic library to be created
        std::ostringstream _cache;
    public:

        CLangCompileHelper(ADFun<CppAD::CG<Base> >* fun) :
            _fun(fun),
            _zero(true),
            _jacobian(false),
            _hessian(false),
            _sparseJacobian(false),
            _sparseHessian(false),
            _libraryName("cppad_cg_model.so") {

            assert(_fun != NULL);
        }

        inline bool isCreateHessian() const {
            return _hessian;
        }

        inline void setCreateHessian(bool createFunction) {
            _hessian = createFunction;
        }

        inline bool isCreateJacobian() const {
            return _jacobian;
        }

        inline void setCreateJacobian(bool createFunction) {
            _jacobian = createFunction;
        }

        inline bool isCreateSparseHessian() const {
            return _sparseHessian;
        }

        inline void setCreateSparseHessian(bool createFunction) {
            _sparseHessian = createFunction;
        }

        inline bool isCreateSparseJacobian() const {
            return _sparseJacobian;
        }

        inline void setCreateSparseJacobian(bool createFunction) {
            _sparseJacobian = createFunction;
        }

        inline bool isCreateForwardZero() const {
            return _zero;
        }

        inline void setCreateForwardZero(bool createFunction) {
            _zero = createFunction;
        }

        DynamicLib<Base>* createDynamicLibrary(CLangCompiler<Base>& compiler);

        static inline const std::string baseTypeName();
        
    protected:

        virtual void generateVerionSource();
        
        virtual void generateInfoSource();

        virtual void generateZeroSource();

        virtual void generateJacobianSource();

        virtual void generateHessianSource();

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseHessianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparsitySource(const std::string& function, const std::vector<bool>& sparsityBool, size_t rowCount, size_t colCount);

        virtual DynamicLib<Base>* loadDynamicLibrary();

    private:
        CLangCompileHelper(const CLangCompileHelper&); // not implemented

        CLangCompileHelper& operator=(const CLangCompileHelper&); // not implemented
    };

}

#endif
