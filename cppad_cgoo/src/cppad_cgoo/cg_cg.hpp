#ifndef CPPAD_ADCG_INCLUDED
#define	CPPAD_ADCG_INCLUDED
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2012 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#include <stddef.h>
#include <assert.h>
#include <string>

namespace CppAD {

    /**
     * 
     */
    template<class Base>
    class CG {
    private:
        // value for parameter
        Base* value_;
        // variable ID (zero means that it is either a temporary variable or a parameter)
        size_t id_;
        // the operations used to create this variable (temporary variables only)
        std::string operations_;
        // status of the operations
        OpContainement opTypes_;
        CodeHandler<Base>* handler_;

    public:
        // default constructor
        inline CG();
        // copy constructor
        inline CG(const CG<Base>& orig);
        //assignment operator
        inline CG& operator =(const CG<Base> &x);

        // construction and assignment from base type
        inline CG(const Base& orig);
        inline CG& operator=(const Base &b);

        //

        inline CodeHandler<Base>* getCodeHandler() const {
            return handler_;
        }

        //
        inline bool isVariable() const;
        inline bool isTemporaryVariable() const;
        inline bool isParameter() const;

        inline size_t getVariableID() const;
        inline const Base& getParameterValue() const throw (CGException);

        inline bool IdenticalZero() const throw (CGException);
        inline bool IdenticalOne() const throw (CGException);

        inline std::string operations() const throw (CGException);

        // computed assignment operators
        inline CG<Base>& operator+=(const CG<Base> &right);
        inline CG<Base>& operator-=(const CG<Base> &right);
        inline CG<Base>& operator*=(const CG<Base> &right);
        inline CG<Base>& operator/=(const CG<Base> &right);
        inline CG<Base>& operator+=(const Base &right);
        inline CG<Base>& operator-=(const Base &right);
        inline CG<Base>& operator*=(const Base &right);
        inline CG<Base>& operator/=(const Base &right);
        template< class T>
        inline CG<Base>& operator+=(const T &right);
        template<class T>
        inline CG<Base>& operator-=(const T &right);
        template<class T>
        inline CG<Base>& operator/=(const T &right);
        template<class T>
        inline CG<Base>& operator*=(const T &right);

        // unary operators
        inline CG<Base> operator+() const;
        inline CG<Base> operator-() const;

        // arithmetic binary operators
        inline CG<Base> operator+(const CG<Base> &right) const;
        inline CG<Base> operator-(const CG<Base> &right) const;
        inline CG<Base> operator*(const CG<Base> &right) const;
        inline CG<Base> operator/(const CG<Base> &right) const;

        //
        inline void makeParameter(const Base &b);

        inline void makeVariable(CodeHandler<Base>& handler);
        inline void makeVariable(CodeHandler<Base>& handler, size_t id);

        inline std::string createVariableName() const {
            assert(handler_ != NULL);
            return handler_->createVariableName(*this);
        }

        // destructor
        virtual ~CG();
    protected:

        inline void makeTemporaryVariable(CodeHandler<Base>& handler, const std::string& operations, OpContainement op);

        inline size_t createID() const {
            assert(handler_ != NULL);
            return handler_->createID();
        }

        inline void printOperationAssig(const std::string& var, const std::string& operations) const {
            assert(handler_ != NULL);
            return handler_->printOperationAssign(var, operations);
        }

    private:

        // comparison operators are not used to create code
        //        friend bool operator< <Base> (const CG<Base> &left, const CG<Base> &right);
        //        friend bool operator <= <Base> (const CG<Base> &left, const CG<Base> &right);
        //        friend bool operator> <Base> (const CG<Base> &left, const CG<Base> &right);
        //        friend bool operator >= <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator == <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator != <Base> (const CG<Base> &left, const CG<Base> &right);

        // order determining functions, see ordered.hpp
        friend bool GreaterThanZero <Base> (const CG<Base> &x);
        friend bool GreaterThanOrZero <Base> (const CG<Base> &x);
        friend bool LessThanZero <Base> (const CG<Base> &x);
        friend bool LessThanOrZero <Base> (const CG<Base> &x);
        friend bool abs_geq <Base> (const CG<Base>& x, const CG<Base>& y);

        // EqualOpSeq function
        friend bool EqualOpSeq <Base> (const CG<Base> &u, const CG<Base> &v);

        // NearEqual function
        friend bool NearEqual <Base> (const CG<Base> &x, const CG<Base> &y, const Base &r, const Base &a);

        friend bool NearEqual <Base> (const Base &x, const CG<Base> &y, const Base &r, const Base &a);

        friend bool NearEqual <Base> (const CG<Base> &x, const Base &y, const Base &r, const Base &a);

        // CondExp function
        friend CG<Base> CondExpOp <Base> (
                enum CompareOp cop,
                const CG<Base> &left,
                const CG<Base> &right,
                const CG<Base> &trueCase,
                const CG<Base> &falseCase
                );

        friend CG<Base> sign<Base>(const CG<Base> &x);

        // power function
        friend CG<Base> pow<Base>(const CG<Base> &x, const CG<Base> &y);

        friend CG<Base> abs<Base>(CG<Base> var);

        friend CG<Base> acos<Base>(CG<Base> var);

        friend CG<Base> asin<Base>(CG<Base> var);

        friend CG<Base> atan<Base>(CG<Base> var);

        friend CG<Base> cos<Base>(CG<Base> var);

        friend CG<Base> cosh<Base>(CG<Base> var);

        friend CG<Base> exp<Base>(CG<Base> var);

        friend CG<Base> log<Base>(CG<Base> var);

        friend CG<Base> sin<Base>(CG<Base> var);

        friend CG<Base> sinh<Base>(CG<Base> var);

        friend CG<Base> sqrt<Base>(CG<Base> var);

        friend CG<Base> tan<Base>(CG<Base> var);

        friend CG<Base> tanh<Base>(CG<Base> var);
    };

    template <class Base>
    int Integer(const CG<Base> &x) {
        if (x.isParameter()) {
            return Integer(x.getParameterValue());
        } else {
            assert(false);
        }
    }

    template<class Base>
    inline std::ostream& operator <<(
    /// stream to write to
    std::ostream& os,
    /// vector that is output
    const CppAD::CG<Base>& v) {
        if (v.isParameter()) {
            os << v.getParameterValue();
        } else if (v.isVariable()) {
            os << v.getCodeHandler()->createVariableName(v);
        } else if (v.isTemporaryVariable()) {
            os << v.operations();
        } else {
            assert(0);
        }
        return os;
    }

    template<class Base>
    inline std::ostringstream& operator <<(
    /// steam to write the vector to
    std::ostringstream& os,
    /// vector that is output
    const CppAD::CG<Base>& v) {
        if (v.isParameter()) {
            os << v.getParameterValue();
        } else if (v.isVariable()) {
            os << v.getCodeHandler()->createVariableName(v);
        } else if (v.isTemporaryVariable()) {
            os << v.operations();
        } else {
            assert(0);
        }
        return os;
    }

}

#endif	

