#ifndef CPPAD_CG_DYNAMICLIB_MODEL_INCLUDED
#define CPPAD_CG_DYNAMICLIB_MODEL_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2012 Ciengis
 *
 *  CppADCodeGen is distributed under multiple licenses:
 *
 *   - Common Public License Version 1.0 (CPL1), and
 *   - GNU General Public License Version 2 (GPL2).
 *
 * CPL1 terms and conditions can be found in the file "epl-v10.txt", while
 * terms and conditions for the GPL2 can be found in the file "gpl2.txt".
 * ----------------------------------------------------------------------------
 * Author: Joao Leal
 */

namespace CppAD {

    /**
     * Abstract class used to call compiled code in a dynamic library
     * 
     * @author Joao Leal
     */
    template<class Base>
    class DynamicLibModel {
    protected:
        CGAtomicLibModel<Base>* _atomic;
    public:

        DynamicLibModel() :
            _atomic(NULL) {
        }

        /**
         * Provides the name for this model.
         * 
         * @return The model name
         */
        virtual const std::string& getName() const = 0;

        // Jacobian sparsity 
        virtual std::vector<std::set<size_t> > JacobianSparsitySet() = 0;
        virtual std::vector<bool> JacobianSparsityBool() = 0;
        virtual void JacobianSparsity(std::vector<size_t>& rows, std::vector<size_t>& cols) = 0;

        /**
         * Provides the sparsity of the sum of the hessian for each dependent 
         * variable.
         * 
         * @return The sparsity
         */
        virtual std::vector<std::set<size_t> > HessianSparsitySet() = 0;
        virtual std::vector<bool> HessianSparsityBool() = 0;
        virtual void HessianSparsity(std::vector<size_t>& rows, std::vector<size_t>& cols) = 0;

        /**
         * Provides the sparsity of the hessian for a dependent variable
         * 
         * @param i The index of the dependent variable
         * @return The sparsity
         */
        virtual std::vector<std::set<size_t> > HessianSparsitySet(size_t i) = 0;
        virtual std::vector<bool> HessianSparsityBool(size_t i) = 0;
        virtual void HessianSparsity(size_t i, std::vector<size_t>& rows, std::vector<size_t>& cols) = 0;

        /**
         * Provides the number of independent variables.
         * 
         * @return The number of independent variables
         */
        virtual size_t Domain() const = 0;

        /**
         * Provides the number of dependent variables.
         * 
         * @return The number of dependent variables.
         */
        virtual size_t Range() const = 0;

        /**
         * Defines an atomic function to be used by the compiled code.
         * It should match an atomic function name previously provided to
         * create the source.
         * 
         * @param atomic The atomic function. This object must only be deleted
         *               after the model.
         * @return true if the atomic function is required by the model, false
         *         if it will never be used.
         */
        virtual bool addAtomicFunction(atomic_base<Base>& atomic) = 0;

        /**
         * Evaluates the dependent model variables (zero-order).
         * This method considers that the dynamic library was compiled
         * using a single array for the independent variables (the default
         * behaviour).
         * 
         * @param x The independent variable vectior
         * @return The dependent variable vector
         */
        template<typename VectorBase>
        inline VectorBase ForwardZero(const VectorBase& x) {
            VectorBase dep(Range());
            this->ForwardZero(&x[0], x.size(), &dep[0], dep.size());
            return dep;
        }

        virtual void ForwardZero(const CppAD::vector<bool>& vx,
                                 CppAD::vector<bool>& vy,
                                 const CppAD::vector<Base> &tx,
                                 CppAD::vector<Base>& ty) = 0;

        /**
         * Evaluates the dependent model variables (zero-order).
         * This method considers that the dynamic library was compiled
         * using a single array for the independent variables (the default
         * behaviour).
         * 
         * @param x The independent variable vectior
         * @param dep The dependent variable vector
         */
        template<typename VectorBase>
        inline void ForwardZero(const VectorBase& x, VectorBase& dep) {
            dep.resize(Range());
            this->ForwardZero(&x[0], x.size(), &dep[0], dep.size());
        }

        virtual void ForwardZero(const Base* x, size_t x_size,
                                 Base* dep, size_t dep_size) = 0;

        /**
         * Determines the dependent variable values using a variable number of 
         * independent variable arrays.
         * This method can be useful if the dynamic library was compiled
         * considering that the independent variables are provided by several
         * arrays.
         * 
         * @param x Contains the several independent variable vectors
         * @param dep The values of the dependent variables
         * @param dep_size The number of dependent variables
         */
        virtual void ForwardZero(const std::vector<const Base*> &x,
                                 Base* dep, size_t dep_size) = 0;

        /// calculate dense Jacobian

        template<typename VectorBase>
        inline VectorBase Jacobian(const VectorBase& x) {
            VectorBase jac(Range() * Domain());
            Jacobian(&x[0], x.size(), &jac[0], jac.size());
            return jac;
        }

        template<typename VectorBase>
        inline void Jacobian(const VectorBase& x, VectorBase& jac) {
            jac.resize(Range() * Domain());
            Jacobian(&x[0], x.size(), &jac[0], jac.size());
        }

        virtual void Jacobian(const Base* x, size_t x_size,
                              Base* jac, size_t jac_size) = 0;

        /// calculate Hessian for one component of f

        template<typename VectorBase>
        inline VectorBase Hessian(const VectorBase& x,
                                  const VectorBase& w) {
            VectorBase hess(Domain() * Domain());
            this->Hessian(x, w, hess);
            return hess;
        }

        template<typename VectorBase>
        inline VectorBase Hessian(const VectorBase& x,
                                  size_t i) {
            CPPADCG_ASSERT_KNOWN(i < Range(), "Invalid equation index");

            VectorBase w(Range());
            w[i] = 1.0;
            VectorBase hess(Domain() * Domain());
            this->Hessian(x, w, hess);
            return hess;
        }

        template<typename VectorBase>
        inline void Hessian(const VectorBase& x,
                            const VectorBase& w,
                            VectorBase& hess) {
            this->Hessian(&x[0], x.size(), &w[0], w.size(), &hess[0]);
        }

        virtual void Hessian(const Base* x, size_t x_size,
                             const Base* w, size_t w_size,
                             Base* hess) = 0;

        /**
         * Computes results during a forward mode sweep. 
         * Computes the first-order Taylor coefficients for dependent variables
         * relative to a single independent variable.
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * @warning do not used it as a generic forward mode function!
         * 
         * @param tx The taylor coeficients of the independent variables 
         * @return The taylor coeficients of the dependent variables 
         */
        template<typename VectorBase>
        inline VectorBase ForwardOne(const VectorBase& tx) {
            size_t m = Range();
            const size_t k = 1;
            VectorBase ty((k + 1) * m);
            
            this->ForwardOne(tx, ty);

            VectorBase dy(m);
            for (size_t i = 0; i < m; i++) {
                dy[i] = ty[i * (k + 1) + k];
            }

            return dy;
        }

        /**
         * Computes results during a forward mode sweep. 
         * Computes the first-order Taylor coefficients for dependent variables
         * relative to a single independent variable.
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * @warning do not used it as a generic forward mode function!
         * 
         * @param tx The taylor coeficients of the independent variables 
         * @param ty The taylor coeficients of the dependent variables 
         */
        template<typename VectorBase>
        inline void ForwardOne(const VectorBase& tx,
                               VectorBase& ty) {
            this->ForwardOne(&tx[0], tx.size(), &ty[0], ty.size());
        }

        /**
         * Computes results during a forward mode sweep. 
         * Computes the first-order Taylor coefficients for dependent variables
         * relative to a single independent variable.
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * @warning do not used it as a generic forward mode function!
         * 
         * @param tx The taylor coeficients of the independent variables 
         * @param tx_size The size of tx
         * @param ty The taylor coeficients of the dependent variables 
         * @param ty_size The size of ty
         */
        virtual void ForwardOne(const Base tx[], size_t tx_size,
                                Base ty[], size_t ty_size) = 0;

        /**
         * Computes results during a reverse mode sweep for the evaluation of
         * the jacobian when the model is used through a user defined atomic AD
         * function.
         * @warning do not used it as a generic reverse mode function!
         * 
         * @param tx
         * @param ty
         * @param py
         * @return px
         */
        template<typename VectorBase>
        inline VectorBase ReverseOne(const VectorBase& tx,
                                     const VectorBase& ty,
                                     const VectorBase& py) {
            const size_t k = 0;
            VectorBase px((k + 1) * Domain());
            this->ReverseOne(tx, ty, px, py);
            return px;
        }

        /**
         * Computes results during a reverse mode sweep for the evaluation of
         * the jacobian when the model is used through a user defined atomic AD
         * function.
         * @warning do not used it as a generic reverse mode function!
         * 
         * @param tx
         * @param ty
         * @param px
         * @param py
         */
        template<typename VectorBase>
        inline void ReverseOne(const VectorBase& tx,
                               const VectorBase& ty,
                               VectorBase& px,
                               const VectorBase& py) {
            this->ReverseOne(&tx[0], tx.size(),
                             &ty[0], ty.size(),
                             &px[0], px.size(),
                             &py[0], py.size());
        }

        /**
         * Computes results during a reverse mode sweep for the evaluation of
         * the jacobian when the model is used through a user defined atomic AD
         * function.
         * @warning do not used it as a generic reverse mode function!
         * 
         * @param tx
         * @param ty
         * @param px
         * @param py
         */
        virtual void ReverseOne(const Base tx[], size_t tx_size,
                                const Base ty[], size_t ty_size,
                                Base px[], size_t px_size,
                                const Base py[], size_t py_size) = 0;

        /**
         * Computes second-order results during a reverse mode sweep (p = 2).
         * This method can be used during the evaluation of the hessian when
         * the model is used through a user defined atomic AD function.
         * @warning do not used it as a generic reverse mode function!
         * @warning only the values for px[j * (k+1)] are defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @param tx
         * @param ty
         * @param py
         * @return px
         */
        template<typename VectorBase>
        inline VectorBase ReverseTwo(const VectorBase& tx,
                                     const VectorBase& ty,
                                     const VectorBase& py) {
            const size_t k = 1;
            VectorBase px((k + 1) * Domain());
            this->ReverseTwo(tx, ty, px, py);
            return px;
        }

        /**
         * Computes second-order results during a reverse mode sweep (p = 2).
         * This method can be used during the evaluation of the hessian when
         * the model is used through a user defined atomic AD function.
         * @warning do not used it as a generic reverse mode function!
         * @warning only the values for px[j * (k+1)] are defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @param tx
         * @param ty
         * @param px
         * @param py
         */
        template<typename VectorBase>
        inline void ReverseTwo(const VectorBase& tx,
                               const VectorBase& ty,
                               VectorBase& px,
                               const VectorBase& py) {
            this->ReverseTwo(&tx[0], tx.size(),
                             &ty[0], ty.size(),
                             &px[0], px.size(),
                             &py[0], py.size());
        }

        /**
         * Computes second-order results during a reverse mode sweep (p = 2).
         * This method can be used during the evaluation of the hessian when
         * the model is used through a user defined atomic AD function.
         * @warning do not used it as a generic reverse mode function!
         * @warning only the values for px[j * (k+1)] are defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @param tx
         * @param ty
         * @param px
         * @param py
         */
        virtual void ReverseTwo(const Base tx[], size_t tx_size,
                                const Base ty[], size_t ty_size,
                                Base px[], size_t px_size,
                                const Base py[], size_t py_size) = 0;

        /// calculate sparse Jacobians 

        template<typename VectorBase>
        inline VectorBase SparseJacobian(const VectorBase& x) {
            VectorBase jac(Range() * Domain());
            SparseJacobian(x, jac);
            return jac;
        }

        template<typename VectorBase>
        inline void SparseJacobian(const VectorBase& x,
                                   VectorBase& jac) {
            jac.resize(Range() * Domain());
            SparseJacobian(&x[0], x.size(), &jac[0], jac.size());
        }

        virtual void SparseJacobian(const Base* x, size_t x_size,
                                    Base* jac, size_t jac_size) = 0;

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac,
                                    std::vector<size_t>& row,
                                    std::vector<size_t>& col) = 0;
        virtual void SparseJacobian(const Base* x, size_t x_size,
                                    Base* jac,
                                    size_t const** row,
                                    size_t const** col,
                                    size_t nnz) = 0;

        /**
         * Determines the sparse Jacobian using a variable number of independent 
         * variable arrays. This method can be useful if the dynamic library
         * was compiled considering that the independent variables are provided
         * by several arrays.
         * 
         * @param x Contains the several independent variable vectors
         * @param jac The values of the sparse Jacobian in the order provided by
         *            row and col
         * @param row The row indices of the Jacobian values
         * @param col The column indices of the Jacobian values
         * @param nnz The total number of non-zero elements
         */
        virtual void SparseJacobian(const std::vector<const Base*>& x,
                                    Base* jac,
                                    size_t const** row,
                                    size_t const** col,
                                    size_t nnz) = 0;

        /// calculate sparse Hessians 

        template<typename VectorBase>
        inline VectorBase SparseHessian(const VectorBase& x,
                                        const VectorBase& w) {
            VectorBase hess(Domain() * Domain());
            SparseHessian(x, w, hess);
            return hess;
        }

        template<typename VectorBase>
        inline void SparseHessian(const VectorBase& x,
                                  const VectorBase& w,
                                  VectorBase& hess) {
            hess.resize(Domain() * Domain());
            SparseHessian(&x[0], x.size(), &w[0], w.size(), &hess[0], hess.size());
        }

        virtual void SparseHessian(const Base* x, size_t x_size,
                                   const Base* w, size_t w_size,
                                   Base* hess, size_t hess_size) = 0;

        virtual void SparseHessian(const std::vector<Base> &x,
                                   const std::vector<Base> &w,
                                   std::vector<Base>& hess,
                                   std::vector<size_t>& row,
                                   std::vector<size_t>& col) = 0;
        virtual void SparseHessian(const Base* x, size_t x_size,
                                   const Base* w, size_t w_size,
                                   Base* hess,
                                   size_t const** row,
                                   size_t const** col,
                                   size_t nnz) = 0;

        /**
         * Determines the sparse Hessian using a variable number of independent 
         * variable arrays. This method can be useful if the dynamic library
         * was compiled considering that the independent variables are provided
         * by several arrays.
         * 
         * @param x Contains the several independent variable vectors
         * @param w The equation multipliers
         * @param w_size The number of equations
         * @param hess The values of the sparse hessian in the order provided by
         *             row and col
         * @param row The row indices of the hessian values
         * @param col The column indices of the hessian values
         * @param nnz The total number of non-zero elements
         */
        virtual void SparseHessian(const std::vector<const Base*>& x,
                                   const Base* w, size_t w_size,
                                   Base* hess,
                                   size_t const** row,
                                   size_t const** col,
                                   size_t nnz) = 0;

        /**
         * Provides a wrapper for this compiled model allowing it to be used as
         * an atomic function. The compiled model must not be deleted while
         * the atomic function is in used.
         * 
         * @return an atomic function wrapper for this model
         */
        virtual CGAtomicLibModel<Base>& asAtomic() {
            if (_atomic == NULL) {
                _atomic = new CGAtomicLibModel<Base>(*this);
            }
            return *_atomic;
        }

        inline virtual ~DynamicLibModel() {
            delete _atomic;
        };
    };

}

#endif


