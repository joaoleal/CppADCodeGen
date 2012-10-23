#ifndef CPPAD_CG_CG_INCLUDED
#define	CPPAD_CG_CG_INCLUDED
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
     * A node in the operation graph
     * 
     * \author Joao Leal
     */
    template<class Base>
    class CG {
    private:
        // the source code handler (NULL for parameters)
        CodeHandler<Base>* handler_;
        // the source code that generated this variable (NULL for parameters)
        SourceCodeFragment<Base>* sourceCode_;
        // value (constant parameters only)
        Base * value_;

    public:
        // default constructor (creates a parameter with a zero value)
        inline CG();

        // copy constructor
        inline CG(const CG<Base>& orig);
        //assignment operator
        inline CG& operator =(const CG<Base> &rhs);

        // construction and assignment from base type
        inline CG(const Base& orig);
        inline CG& operator=(const Base &b);

        //

        inline CodeHandler<Base>* getCodeHandler() const {
            return handler_;
        }

        // variable classification methods
        inline bool isVariable() const;
        inline bool isParameter() const;

        inline const Base& getParameterValue() const throw (CGException);

        inline bool IdenticalZero() const throw (CGException);
        inline bool IdenticalOne() const throw (CGException);

        inline SourceCodeFragment<Base>* getSourceCodeFragment() const;

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

        // destructor
        virtual ~CG();
    protected:
        // creates a temporary variable
        inline CG(CodeHandler<Base>& handler, SourceCodeFragment<Base>* sourceCode);

        // 
        /**
         * creates a variable/parameter from an existing argument 
         * \param handler The code handler.
         * \param arg  An argument that may be a parameter or a variable. 
         *              (source code fragments are assumed to already be managed
         *              by the handler)
         */
        inline CG(CodeHandler<Base>& handler, const Argument<Base>& arg);

        //
        inline void makeParameter(const Base &b);
        inline void makeVariable(CodeHandler<Base>& handler, SourceCodeFragment<Base>* operation);

        inline Argument<Base> argument() const;
    private:

        friend class CodeHandler<Base>;

        template<class Base1, class Base2>
        friend class Evaluator;

        /**
         * arithmetic binary operators
         */
        friend CG<Base> CppAD::operator+ <Base>(const CG<Base> &left, const CG<Base> &right);
        friend CG<Base> CppAD::operator- <Base>(const CG<Base> &left, const CG<Base> &right);
        friend CG<Base> CppAD::operator* <Base>(const CG<Base> &left, const CG<Base> &right);
        friend CG<Base> CppAD::operator/ <Base>(const CG<Base> &left, const CG<Base> &right);

        /**
         * comparison operators are not used to create code
         */
        friend bool operator< <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator <= <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator> <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator >= <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator == <Base> (const CG<Base> &left, const CG<Base> &right);
        friend bool operator != <Base> (const CG<Base> &left, const CG<Base> &right);

        /**
         * order determining functions
         */
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
        friend CG<Base> CondExpOp <Base> (enum CompareOp cop,
                const CG<Base> &left,
                const CG<Base> &right,
                const CG<Base> &trueCase,
                const CG<Base> &falseCase);

        friend CG<Base> sign<Base>(const CG<Base> &x);

        /**
         * math functions
         */
        friend CG<Base> pow<Base>(const CG<Base> &x, const CG<Base> &y);

        friend CG<Base> abs<Base>(const CG<Base>& var);

        friend CG<Base> acos<Base>(const CG<Base>& var);

        friend CG<Base> asin<Base>(const CG<Base>& var);

        friend CG<Base> atan<Base>(const CG<Base>& var);

        friend CG<Base> cos<Base>(const CG<Base>& var);

        friend CG<Base> cosh<Base>(const CG<Base>& var);

        friend CG<Base> exp<Base>(const CG<Base>& var);

        friend CG<Base> log<Base>(const CG<Base>& var);

        friend CG<Base> sin<Base>(const CG<Base>& var);

        friend CG<Base> sinh<Base>(const CG<Base>& var);

        friend CG<Base> sqrt<Base>(const CG<Base>& var);

        friend CG<Base> tan<Base>(const CG<Base>& var);

        friend CG<Base> tanh<Base>(const CG<Base>& var);
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
    std::ostream& os, //< stream to write to
    const CppAD::CG<Base>& v//< vector that is output
    ) {
        if (v.isParameter()) {
            os << v.getParameterValue();
        } else {
            os << *v.getSourceCodeFragment();
        }
        return os;
    }

    template<class Base>
    inline std::ostringstream& operator <<(
    std::ostringstream& os, //< steam to write the vector to
    const CppAD::CG<Base>& v//< vector that is output
    ) {
        if (v.isParameter()) {
            os << v.getParameterValue();
        } else {
            os << *v.getSourceCodeFragment();
        }
        return os;
    }

}

#endif	

