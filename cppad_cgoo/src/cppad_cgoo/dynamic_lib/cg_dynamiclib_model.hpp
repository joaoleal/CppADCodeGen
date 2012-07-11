#ifndef CPPAD_CG_DYNAMICLIB_MODEL_INCLUDED
#define	CPPAD_CG_DYNAMICLIB_MODEL_INCLUDED
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
     * Abstract class used to call compiled code in a dynamic library
     * 
     * \author Joao Leal
     */
    template<class Base>
    class DynamicLibModel {
    public:
        virtual const std::string& getName() const = 0;

        // Jacobian sparsity 
        virtual std::vector<std::set<size_t> > JacobianSparsitySet() = 0;
        virtual std::vector<bool> JacobianSparsityBool() = 0;

        // Hessian sparsity 
        virtual std::vector<std::set<size_t> > HessianSparsitySet() = 0;
        virtual std::vector<bool> HessianSparsityBool() = 0;

        /// number of independent variables
        virtual size_t Domain() const = 0;

        /// number of dependent variables
        virtual size_t Range() const = 0;

        /// calculate the dependent values (zero order)
        virtual std::vector<Base> ForwardZero(const std::vector<Base> &x) = 0;

        virtual void ForwardZero(const std::vector<Base> &x, std::vector<Base>& dep) = 0;

        /// calculate entire Jacobian
        virtual std::vector<Base> Jacobian(const std::vector<Base> &x) = 0;

        virtual void Jacobian(const std::vector<Base> &x, std::vector<Base>& jac) = 0;

        /// calculate Hessian for one component of f
        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          const std::vector<Base> &w) = 0;

        virtual std::vector<Base> Hessian(const std::vector<Base> &x,
                                          size_t i) = 0;

        virtual void Hessian(const std::vector<Base> &x,
                             const std::vector<Base> &w,
                             std::vector<Base>& hess) = 0;

        /// calculate sparse Jacobians 
        virtual std::vector<Base> SparseJacobian(const std::vector<Base> &x) = 0;

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac) = 0;

        virtual void SparseJacobian(const std::vector<Base> &x,
                                    std::vector<Base>& jac,
                                    std::vector<size_t>& row,
                                    std::vector<size_t>& col) = 0;

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
    };

}

#endif


