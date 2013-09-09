#ifndef CPPAD_CG_BASE_ABSTRACT_ATOMIC_FUN_INCLUDED
#define CPPAD_CG_BASE_ABSTRACT_ATOMIC_FUN_INCLUDED
/* --------------------------------------------------------------------------
 *  CppADCodeGen: C++ Algorithmic Differentiation with Source Code Generation:
 *    Copyright (C) 2013 Ciengis
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
     * Contains some utility methods for atomic functions
     * 
     * @author Joao Leal
     */
    template <class Base>
    class BaseAbstractAtomicFun : public atomic_base<CppAD::CG<Base> > {
    public:
        typedef CppAD::CG<Base> CGB;
        typedef Argument<Base> Arg;
    protected:

        /**
         * Creates a new atomic function that is responsible for defining the
         * dependencies to calls of a user atomic function.
         * 
         * @param name The atomic function name.
         */
        BaseAbstractAtomicFun(const std::string& name) :
            atomic_base<CGB>(name) {
            CPPADCG_ASSERT_KNOWN(!name.empty(), "The atomic function name cannot be empty");
        }

    public:

        template <class ADVector>
        void operator()(const ADVector& ax, ADVector& ay, size_t id = 0) {
            this->atomic_base<CGB>::operator()(ax, ay, id);
        }

        virtual ~BaseAbstractAtomicFun() {
        }

    protected:

        static inline void appendAsArguments(typename std::vector<Arg>::iterator begin, const vector<CGB>& tx) {
            std::vector<Arg> arguments(tx.size());
            typename std::vector<Arg>::iterator it = begin;
            for (size_t i = 0; i < arguments.size(); i++, ++it) {
                if (tx[i].isParameter()) {
                    *it = Arg(tx[i].getValue());
                } else {
                    *it = Arg(*tx[i].getOperationNode());
                }
            }
        }

        static inline OperationNode<Base>* makeArray(CodeHandler<Base>& handler,
                                                     const vector<CGB>& tx) {
            if (false && tx.size() > 0) {
                /**
                 * Cannot reuse arrays yet! :(
                 * the usage of the new array elements would have to depend on
                 * the last usage of all the elements in tx
                 */
                OperationNode<Base>* op = tx[0].getOperationNode();
                if (op != NULL && op->getOperationType() == CGArrayElementOp) {
                    OperationNode<Base>* otherArray = op->getArguments()[0].getOperation();
                    bool reuseArray = true;
                    for (size_t i = 0; i < tx.size(); i++) {
                        op = tx[i].getOperationNode();
                        if (op == NULL ||
                                op->getOperationType() != CGArrayElementOp ||
                                op->getArguments()[0].getOperation() != otherArray ||
                                op->getInfo()[0] != i) {
                            reuseArray = false;
                            break;
                        }
                    }
                    if (reuseArray) {
                        return otherArray;
                    }
                }
            }
            std::vector<Arg> arrayArgs = asArguments(tx);
            std::vector<size_t> info; // empty
            OperationNode<Base>* array = new OperationNode<Base>(CGArrayCreationOp, info, arrayArgs);
            handler.manageOperationNode(array);
            return array;
        }

        static inline OperationNode<Base>* makeZeroArray(CodeHandler<Base>& handler,
                                                         const vector<CGB>& tx) {
            vector<CGB> tx2(tx.size());
            std::vector<Arg> arrayArgs = asArguments(tx2);
            std::vector<size_t> info; // empty
            OperationNode<Base>* array = new OperationNode<Base>(CGArrayCreationOp, info, arrayArgs);
            handler.manageOperationNode(array);
            return array;
        }

        static inline bool isParameters(const vector<CGB>& tx) {
            for (size_t i = 0; i < tx.size(); i++) {
                if (!tx[i].isParameter()) {
                    return false;
                }
            }
            return true;
        }

        static inline bool isValuesDefined(const vector<CGB>& tx) {
            for (size_t i = 0; i < tx.size(); i++) {
                if (!tx[i].isValueDefined()) {
                    return false;
                }
            }
            return true;
        }

    };

}

#endif