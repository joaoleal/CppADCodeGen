/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#ifndef CPPAD_CODEGEN_AD_FUN_CODE_GEN_INCLUDED
#define	CPPAD_CODEGEN_AD_FUN_CODE_GEN_INCLUDED

//#include <cppad/local/ad_fun.hpp>
#include <cppad_codegen/local/ad_code_gen_name_provider.hpp>

CPPAD_BEGIN_NAMESPACE

enum DiffMode {
    AUTO, FORWARD, REVERSE
};

static const std::set<size_t> CODEGEN_EMPTY_SET;

template <class Base>
class ADFunCodeGen : public ADFun<Base> {
    // ------------------------------------------------------------
private:
    /// number of taylor_ coefficient per variable (currently generated in source)
    size_t taylor_per_var_;
    DefaultCCodeGenNameProvider<Base> defaultCodeNameGen_;
    // Source code generation
    CodeGenNameProvider<Base>* nameGen_;
    // Indicates which taylor coefficients are always zero in the current 
    // source code generation
public:
    /// copy constructor

    ADFunCodeGen(const ADFunCodeGen<Base>& g) : ADFun<Base>(g) {
    }

    /// default constructor

    ADFunCodeGen() : ADFun<Base>() {
        nameGen_ = &defaultCodeNameGen_;
        taylor_per_var_ = 0;
    }

    /// sequence constructor

    template <typename ADvector>
    ADFunCodeGen(const ADvector &x, const ADvector &y) : ADFun<Base>(x, y) {
        nameGen_ = &defaultCodeNameGen_;
        taylor_per_var_ = 0;
    }

    // assignment operator
    // (see doxygen in fun_construct.hpp)

    void operator=(const ADFunCodeGen& f) {
        if (&f == this) {
            return;
        }

        ADFun<Base>::operator=(f);
        taylor_per_var_ = f.taylor_per_var_;
        nameGen_ = f.nameGen_;
    }

    /// destructor

    virtual ~ADFunCodeGen() {
    }

    CodeGenNameProvider<Base>* getCodeGenNameProvider() const {
        return nameGen_;
    }

    void ForwardCodeGen(size_t p, std::ostream& s_out);

    template <typename VectorBase>
    void ReverseCodeGen(size_t p, const VectorBase &w, std::ostream& s_out);

    template<class VectorBool>
    VectorBool SparseJacobianCodeGen(std::ostream& s_out, DiffMode& mode, bool compress = false, const std::set<size_t>& indepFilter = CODEGEN_EMPTY_SET);

    template<class VectorBool>
    void SparseJacobianCodeGen(std::ostream& s_out, const VectorBool& p, DiffMode& mode, bool compress = false, const std::set<size_t>& indepFilter = CODEGEN_EMPTY_SET);

private:

    template<class VectorBool>
    void SparseJacobianCaseCodeGen(std::ostream& s_out, DiffMode mode, const std::set<size_t>& set_type, const VectorBool& p, bool compress, const std::set<size_t>& indepFilter);

    template <class VectorBool>
    void SparseJacobianCaseCodeGen(std::ostream& s_out, DiffMode mode, bool set_type, const VectorBool& p, bool compress, const std::set<size_t>& indepFilter);
};

CPPAD_END_NAMESPACE

#endif	/* CPPAD_CODEGEN_AD_FUN_CODE_GEN_INCLUDED */

