#ifndef CPPAD_CG_DECLARE_CG_INCLUDED
#define CPPAD_CG_DECLARE_CG_INCLUDED
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

#include <cppad/local/define.hpp>
#include <cppad/local/cppad_assert.hpp>
#include <cppad/local/base_cond_exp.hpp>

// forward declarations
namespace CppAD {

    template<class Base>
    class vector;

    /***************************************************************************
     * Atomics
     **************************************************************************/
    template<class Base>
    class BaseAbstractAtomicFun;

    template<class Base>
    class CGAbstractAtomicFun;

    template<class Base>
    class CGAtomicFun;

    /***************************************************************************
     * Core
     **************************************************************************/
    template<class Base>
    class CodeHandler;

    template<class Base>
    class CG;

    template<class Base>
    class AD;

    template<class Base>
    struct OperationPathNode;

    template<class Base>
    class ScopePathElement;

    /**************************************************************************/
    class SparseForjacHessianWorkJac;
    class SparseForjacHessianWorkHes;
    class SparseForjacHessianWork;

    /***************************************************************************
     * Nodes
     **************************************************************************/
    template<class Base>
    class OperationNode;

    template<class Base>
    class IndexOperationNode;

    template<class Base>
    class IndexAssignOperationNode;

    template<class Base>
    class IndexDclrOperationNode;

    template<class Base>
    class LoopStartOperationNode;

    template<class Base>
    class LoopEndOperationNode;

    /***************************************************************************
     * Loops
     **************************************************************************/
    template<class Base>
    class EquationPattern;

    template<class Base>
    class DependentPatternMatcher;

    template<class Base>
    class Loop;

    template<class Base>
    class LoopFreeModel;

    template<class Base>
    class LoopModel;

    template<class Base>
    class IndexedDependentLoopInfo;

    class IndexPattern;

    /***************************************************************************
     * Languages
     **************************************************************************/

    template<class Base>
    class CLanguage;

    template<class Base>
    class VariableNameGenerator;

    template<class Base>
    class CLangDefaultVariableNameGenerator;

    template<class Base>
    class CLangCustomVariableNameGenerator;

    /***************************************************************************
     * Dynamic model compilation
     **************************************************************************/

    template<class Base>
    class CLangCompiler;

    template<class Base>
    class DynamicLibModel;

    template<class Base>
    class DynamicLib;

    template<class Base>
    class CGAtomicLibModel;

    template<class Base>
    class CLangCompileModelHelper;

    template<class Base>
    class CLangCompileDynamicHelper;

#ifdef __linux__
    template<class Base>
    class LinuxDynamicLibModel;

    template<class Base>
    class LinuxDynamicLib;
#endif

    /***************************************************************************
     * Index reduction classes
     **************************************************************************/
    template<class Base>
    class Enode;

    template<class Base>
    class Vnode;

    template<class Base, class Base2>
    class Evaluator;

    /***************************************************************************
     *  Utilities
     **************************************************************************/

    template<class Base>
    class SmartVectorPointer;

    template<class Key, class Value>
    class SmartMapValuePointer;

    /***************************************************************************
     * 
     **************************************************************************/
    // order determining functions, see ordered.hpp
    template<class Base>
    bool GreaterThanZero(const CG<Base> &x);

    template<class Base>
    bool GreaterThanOrZero(const CG<Base> &x);

    template<class Base>
    bool LessThanZero(const CG<Base> &x);

    template<class Base>
    bool LessThanOrZero(const CG<Base> &x);

    template<class Base>
    bool abs_geq(const CG<Base>& x, const CG<Base>& y);

    // The identical property functions, see identical.hpp
    template<class Base>
    inline bool IdenticalPar(const CG<Base>& x) throw (CGException);

    template<class Base>
    bool IdenticalZero(const CG<Base> &x) throw (CGException);

    template<class Base>
    bool IdenticalOne(const CG<Base> &x) throw (CGException);

    template<class Base>
    bool IdenticalEqualPar(const CG<Base> &x, const CG<Base> &y);

    // EqualOpSeq function
    template<class Base>
    bool EqualOpSeq(const CG<Base> &u, const CG<Base> &v);

    // NearEqual function
    template<class Base>
    bool NearEqual(const CG<Base> &x, const CG<Base> &y, const Base &r, const Base &a);

    template<class Base>
    bool NearEqual(const Base &x, const CG<Base> &y, const Base &r, const Base &a);

    template<class Base>
    bool NearEqual(const CG<Base> &x, const Base &y, const Base &r, const Base &a);

    template <class Base>
    inline bool isnan(const CG<Base> &s);

    // CondExp function
    //    template<class Base>
    //    AD<CG<Base> > CondExpOp(enum CompareOp cop, const AD<CG<Base> > &left, const AD<CG<Base> > &right, const AD<CG<Base> > &trueCase, const AD<CG<Base> > &falseCase);

    template<class Base>
    CG<Base> CondExpOp(enum CompareOp cop, const CG<Base> &left, const CG<Base> &right, const CG<Base> &trueCase, const CG<Base> &falseCase);

    /**
     * arithmetic
     */
    template<class Base>
    CG<Base> operator+(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    CG<Base> operator-(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    CG<Base> operator*(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    CG<Base> operator/(const CG<Base> &left, const CG<Base> &right);

    /**
     * comparisons
     */
    template<class Base>
    bool operator ==(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator !=(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator<(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator <=(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator>(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator >=(const CG<Base> &left, const CG<Base> &right);

    template<class Base>
    bool operator !=(const CG<Base> &left, double right);

    /**
     * Math functions
     */
    template<class Base>
    inline CG<Base> sign(const CG<Base> &x);

    // power function
    template<class Base>
    inline CG<Base> pow(const CG<Base> &x, const CG<Base> &y);

    template<class Base>
    inline CG<Base> abs(const CG<Base>& var);

    template<class Base>
    inline CG<Base> acos(const CG<Base>& var);

    template<class Base>
    inline CG<Base> asin(const CG<Base>& var);

    template<class Base>
    inline CG<Base> atan(const CG<Base>& var);

    template<class Base>
    inline CG<Base> cos(const CG<Base>& var);

    template<class Base>
    inline CG<Base> cosh(const CG<Base>& var);

    template<class Base>
    inline CG<Base> exp(const CG<Base>& var);

    template<class Base>
    inline CG<Base> log(const CG<Base>& var);

    template<class Base>
    inline CG<Base> sin(const CG<Base>& var);

    template<class Base>
    inline CG<Base> sinh(const CG<Base>& var);

    template<class Base>
    inline CG<Base> sqrt(const CG<Base>& var);

    template<class Base>
    inline CG<Base> tan(const CG<Base>& var);

    template<class Base>
    inline CG<Base> tanh(const CG<Base>& var);

    /***************************************************************************
     * Utility functions
     **************************************************************************/
    template<class VectorBool, class Base>
    inline VectorBool jacobianSparsity(ADFun<Base>& fun);

    template<class VectorBool, class VectorSize>
    inline void generateSparsityIndexes(const VectorBool& sparsity,
                                        size_t m,
                                        size_t n,
                                        VectorSize& row,
                                        VectorSize& col);

    template<class VectorSet, class VectorSize>
    inline void generateSparsityIndexes(const VectorSet& sparsity,
                                        VectorSize& row,
                                        VectorSize& col);

    /**
     * Index reduction functions
     */
    template<class Base>
    inline std::ostream& operator <<(std::ostream& os, const Enode<Base>& i);

    template<class Base>
    inline std::ostream& operator <<(std::ostream& os, const Vnode<Base>& j);

    /***************************************************************************
     * Enums
     **************************************************************************/

    /**
     * Automatic Differentiation modes used to determine the Jacobian
     */
    enum JacobianADMode {
        FORWARD, REVERSE, AUTOMATIC
    };

    /**
     * Index pattern types
     */
    enum IndexPatternType {
        LINEAR, // y = (x / dx) * dy + b
        SECTIONED, // several index patterns
        RANDOM,
        PLANE2D // y = f(x) + f(z)
    };
}

/**
 * loops namespace
 */
#include <cppadcg/cg_declare_cg_loops.hpp>

#endif

