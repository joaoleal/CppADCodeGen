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
        static const std::string CONST;
    protected:
        ADFun<CppAD::CG<Base> >* _fun;
        const std::string _baseTypeName;
        bool _zero;
        bool _jacobian;
        bool _hessian;
        bool _sparseJacobian;
        bool _sparseHessian;
        std::vector<size_t> _custom_jac_row;
        std::vector<size_t> _custom_jac_col;
        std::vector<size_t> _custom_hess_row;
        std::vector<size_t> _custom_hess_col;
        std::string _libraryName; // the path of the dynamic library to be created
        std::ostringstream _cache;
        size_t _maxAssignPerFunc; // maximum number of assignments per function (~ lines)
    public:

        CLangCompileHelper(ADFun<CppAD::CG<Base> >* fun) :
            _fun(fun),
            _baseTypeName(CLangCompileHelper<Base>::baseTypeName()),
            _zero(true),
            _jacobian(false),
            _hessian(false),
            _sparseJacobian(false),
            _sparseHessian(false),
            _libraryName("cppad_cg_model.so"),
            _maxAssignPerFunc(20000) {

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

        inline void setCustomSparseJacobianElements(const std::vector<size_t>& row,
                                                    const std::vector<size_t>& col) {
            _custom_jac_row = row;
            _custom_jac_col = col;
        }

        inline void setCustomSparseHessianElements(const std::vector<size_t>& row,
                                                   const std::vector<size_t>& col) {
            _custom_hess_row = row;
            _custom_hess_col = col;
        }

        inline const std::string& getLibraryName() const {
            return _libraryName;
        }

        inline void setLibraryName(const std::string& libraryName) {
            _libraryName = libraryName;
        }

        inline size_t getMaxAssignmentsPerFunc() const {
            return _maxAssignPerFunc;
        }

        inline void setMaxAssignmentsPerFunc(size_t maxAssignPerFunc) {
            _maxAssignPerFunc = maxAssignPerFunc;
        }

        DynamicLib<Base>* createDynamicLibrary(CLangCompiler<Base>& compiler);

        static inline std::string baseTypeName();

    protected:

        virtual void generateVerionSource(std::map<std::string, std::string>& sources);

        virtual void generateInfoSource(std::map<std::string, std::string>& sources);

        virtual void generateZeroSource(std::map<std::string, std::string>& sources);

        virtual void generateJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateHessianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseJacobianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparseHessianSource(std::map<std::string, std::string>& sources);

        virtual void generateSparsitySource(const std::string& function,
                                            const std::vector<bool>& sparsity,
                                            size_t m, size_t n);

        virtual void generateSparsitySource(const std::string& function,
                                            const std::vector<size_t>& rows,
                                            const std::vector<size_t>& cols);

        virtual void generateSparsityIndexes(const std::vector<bool>& sparsity,
                                             size_t m, size_t n,
                                             std::vector<size_t>& rows,
                                             std::vector<size_t>& cols);

        virtual DynamicLib<Base>* loadDynamicLibrary();

    private:
        CLangCompileHelper(const CLangCompileHelper&); // not implemented

        CLangCompileHelper& operator=(const CLangCompileHelper&); // not implemented
    };

}

#endif
