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
    public:
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

        /// calculate the dependent values (zero order)
        virtual CppAD::vector<Base> ForwardZero(const CppAD::vector<Base> &x) = 0;
        virtual std::vector<Base> ForwardZero(const std::vector<Base> &x) = 0;

        virtual void ForwardZero(const CppAD::vector<bool>& vx,
                                 CppAD::vector<bool>& vy,
                                 const CppAD::vector<Base> &tx,
                                 CppAD::vector<Base>& ty) = 0;
        virtual void ForwardZero(const std::vector<Base> &x, std::vector<Base>& dep) = 0;
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

        /// calculate entire Jacobian
        virtual std::vector<Base> Jacobian(const std::vector<Base> &x) = 0;

        virtual void Jacobian(const std::vector<Base> &x, std::vector<Base>& jac) = 0;
        virtual void Jacobian(const Base* x, size_t x_size,
                              Base* jac, size_t jac_size) = 0;

        /// calculate Hessian for one component of f
        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          const std::vector<Base> &w) = 0;

        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          size_t i) = 0;

        virtual void Hessian(const std::vector<Base> &x,
                             const std::vector<Base> &w,
                             std::vector<Base>& hess) = 0;
        virtual void Hessian(const Base* x, size_t x_size,
                             const Base* w, size_t w_size,
                             Base* hess) = 0;

        /**
         * Computes results during a forward mode sweep. 
         * Computes the first-order Taylor coefficients for dependent variables
         * relative to a single independent variable.
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * 
         * @param tx
         * @param ty The values of the directional derivatives
         */
        virtual void SparseForwardOne(const CppAD::vector<Base>& tx,
                                      CppAD::vector<Base>& ty) = 0;

        /**
         * Computes results during a forward mode sweep. 
         * Computes the first-order Taylor coefficients for dependent variables
         * relative to a single independent variable.
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * 
         * @param tx
         * @return ty
         */
        virtual CppAD::vector<Base> SparseForwardOne(const CppAD::vector<Base>& tx) = 0;

        /**
         * Computes results during a reverse mode sweep. 
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * 
         * @param tx
         * @param ty
         * @param px
         * @param py
         */
        virtual void SparseReverseOne(const CppAD::vector<Base>& tx,
                                      const CppAD::vector<Base>& ty,
                                      CppAD::vector<Base>& px,
                                      const CppAD::vector<Base>& py) = 0;

        /**
         * Computes results during a reverse mode sweep. 
         * This method can be used during the evaluation of the jacobian when
         * the model is used through a user defined atomic AD function.
         * 
         * @param tx
         * @param ty
         * @param py
         * @return px
         */
        virtual CppAD::vector<Base> SparseReverseOne(const CppAD::vector<Base>& tx,
                                                     const CppAD::vector<Base>& ty,
                                                     const CppAD::vector<Base>& py) = 0;

        /**
         * Computes second-order results during a reverse mode sweep (p = 2).
         * This method can be used during the evaluation of the hessian when
         * the model is used through a user defined atomic AD function.
         * Warning: only the values for px[j * (k+1)] are defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @param tx
         * @param ty
         * @param px
         * @param py
         */
        virtual void SparseReverseTwo(const CppAD::vector<Base>& tx,
                                      const CppAD::vector<Base>& ty,
                                      CppAD::vector<Base>& px,
                                      const CppAD::vector<Base>& py) = 0;

        /**
         * Computes second-order results during a reverse mode sweep (p = 2).
         * This method can be used during the evaluation of the hessian when
         * the model is used through a user defined atomic AD function.
         * Warning: only the values for px[j * (k+1)] are defined, since
         *          px[j * (k+1) + 1] is not used during the hessian evaluation.
         * 
         * @param tx
         * @param ty
         * @param py
         * @return px
         */
        virtual CppAD::vector<Base> SparseReverseTwo(const CppAD::vector<Base>& tx,
                                                     const CppAD::vector<Base>& ty,
                                                     const CppAD::vector<Base>& py) = 0;

        /// calculate sparse Jacobians 
        virtual std::vector<Base> SparseJacobian(const std::vector<Base> &x) = 0;

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac) = 0;

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
         * @param jac The values of the sparse Jacobian in the order provided by row and col
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
        virtual std::vector<Base> SparseHessian(const std::vector<Base> &x,
                                                const std::vector<Base> &w) = 0;

        virtual void SparseHessian(const std::vector<Base> &x,
                                   const std::vector<Base> &w,
                                   std::vector<Base>& hess) = 0;

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
         * @param hess The values of the sparse hessian in the order provided by row and col
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

        inline virtual ~DynamicLibModel() {
        };
    };

}

#endif


